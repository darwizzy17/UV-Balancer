# uv_exposure.py
from dataclasses import dataclass
from typing import List, Optional, Tuple

MED_TABLE = {1:200.0,2:250.0,3:300.0,4:450.0,5:600.0,6:800.0}

def ery_irradiance_from_uvi(uvi: float) -> float:
    return 0.025 * uvi

def effective_protection_factor(spf: Optional[float], thickness_mg_cm2: Optional[float]) -> float:
    if not spf or spf <= 1 or not thickness_mg_cm2 or thickness_mg_cm2 <= 0:
        return 1.0
    return max(1.0, float(spf) ** (float(thickness_mg_cm2) / 2.0))

@dataclass
class UviPoint:
    t_sec: float
    uvi: float

@dataclass
class Inputs:
    uvi: Optional[float] = None
    uvi_series: Optional[List[UviPoint]] = None
    fitz: int = 3
    age: int = 30
    body_frac_exposed: float = 0.35
    spf: Optional[float] = None
    thickness_mg_cm2: float = 0.0
    covered_frac_of_exposed: float = 1.0
    vitaminD_fraction_of_MED: float = 0.25
    safety_factor: float = 0.8

@dataclass
class Outputs:
    time_vitaminD_on_exposed_sec: float
    effective_PF: float

def _time_to_dose_seconds(target_dose_jm2: float, attenuator: float, uvi: Optional[float], uvi_series: Optional[List[UviPoint]]) -> float:
    eps=1e-9
    if uvi_series and len(uvi_series)>=2:
        cum=0.0
        for i in range(1,len(uvi_series)):
            U0,U1= uvi_series[i-1].uvi, uvi_series[i].uvi
            t0,t1= uvi_series[i-1].t_sec, uvi_series[i].t_sec
            dt=max(0.0,t1-t0)
            Eavg=ery_irradiance_from_uvi(0.5*(U0+U1))/max(attenuator,eps)
            dD=Eavg*dt
            if cum+dD>=target_dose_jm2:
                remain=target_dose_jm2-cum
                frac=remain/max(dD,eps)
                return t0+frac*dt
            cum+=dD
        return float("inf")
    else:
        if uvi is None:
            raise ValueError("Need uvi or uvi_series")
        E=ery_irradiance_from_uvi(uvi)/max(attenuator,eps)
        return target_dose_jm2/max(E,eps)

def compute_uv_exposure(params: Inputs) -> Outputs:
    MED=MED_TABLE[params.fitz]
    D_vitD=params.vitaminD_fraction_of_MED*MED
    PF=effective_protection_factor(params.spf, params.thickness_mg_cm2)
    mix=(1.0-params.covered_frac_of_exposed)+(params.covered_frac_of_exposed/max(PF,1.0))
    atten=1.0/max(mix,1e-6)
    t_vitd=_time_to_dose_seconds(D_vitD,atten,params.uvi,params.uvi_series)
    return Outputs(time_vitaminD_on_exposed_sec=t_vitd,effective_PF=PF)

def skin_efficiency_factor(fitz:int)->float:
    return {1:1.0,2:1.0,3:0.8,4:0.6,5:0.4,6:0.25}.get(fitz,0.5)

def age_efficiency_factor(age:int)->float:
    return 0.7 if age>=65 else 1.0

def estimate_vitaminD_IU_for_session(params:Inputs, session_minutes:float)->Tuple[float,float,float]:
    if session_minutes<=0: return (0.0,0.0,0.0)
    U=params.uvi if params.uvi is not None else params.uvi_series[-1].uvi
    PF=effective_protection_factor(params.spf, params.thickness_mg_cm2)
    mix=(1.0-params.covered_frac_of_exposed)+(params.covered_frac_of_exposed/max(PF,1.0))
    atten=1.0/max(mix,1e-6)
    E=ery_irradiance_from_uvi(U)/atten
    t=session_minutes*60.0
    dose_exposed=E*t
    MED=MED_TABLE[params.fitz]
    fracMED_exposed=dose_exposed/MED
    whole_body_frac=fracMED_exposed*params.body_frac_exposed
    IU_MED_central, IU_MED_low, IU_MED_high=12000.0,10000.0,20000.0
    scale=skin_efficiency_factor(params.fitz)*age_efficiency_factor(params.age)
    return (IU_MED_central*whole_body_frac*scale,
            IU_MED_low*whole_body_frac*scale,
            IU_MED_high*whole_body_frac*scale)
