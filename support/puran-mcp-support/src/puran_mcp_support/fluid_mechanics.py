"""Source-bound API 520 relief sizing, API 521/IEC valve noise and liquid eductors.

SI inputs; no inferred relief load, fire exposure or correction factors. The
optional MIT ``fluids`` dependency is loaded only when a calculation is ready.
Upstream methods: https://fluids.readthedocs.io/fluids.safety_valve.html,
fluids.control_valve.html and fluids.jet_pump.html. Results are calculations,
not equipment approval or an API 521 relieving-scenario study.
"""
from __future__ import annotations
from importlib import import_module
from math import isfinite
from typing import Annotated, Literal
from pydantic import BaseModel, ConfigDict, Field
Positive = Annotated[float, Field(gt=0)]
Nonnegative = Annotated[float, Field(ge=0)]
Factor = Annotated[float, Field(gt=0, le=1)]
Text = Annotated[str, Field(min_length=1, pattern='\\S')]

class Input(BaseModel):
    model_config = ConfigDict(extra='forbid', allow_inf_nan=False)

def _library(module):
    try:
        return import_module('fluids.' + module)
    except ImportError as exc:
        raise ValueError('MECHANICAL_FLUIDS_UNAVAILABLE: install fluids (puran-mcp-support[fluids])') from exc

def _result(request, method, holds=(), **values):
    return dict(status='hold' if holds else 'calculated', method=method, inputs=request.model_dump(mode='json'), holds=list(holds), safety_factor=1.0, safety_factor_basis='No added multiplier; supplied duty and factors retained', **values)

class ReliefValveInput(Input):
    phase: Literal['gas', 'steam', 'liquid']
    material: Text | None = None
    set_pressure_pa_g: Positive | None = None
    relieving_case: Text | None = None
    case_basis: Text | None = None
    fire_case: bool | None = None
    fire_basis: Text | None = None
    back_pressure_pa_abs: Positive | None = None
    atmospheric_pressure_pa: Positive | None = None
    overpressure_fraction: Nonnegative | None = None
    mass_flow_kg_s: Positive | None = None
    temperature_k: Positive | None = None
    density_kg_m3: Positive | None = None
    compressibility: Positive | None = None
    molecular_weight_g_mol: Positive | None = None
    heat_capacity_ratio: Annotated[float, Field(gt=1)] | None = None
    discharge_coefficient: Factor | None = None
    backpressure_factor: Factor | None = None
    combination_factor: Factor | None = None
    viscosity_factor: Factor | None = None
    factors_basis: Text | None = None
    edition: Literal['7E', '10E'] | None = None

def relief_valve_size(request: ReliefValveInput | dict) -> dict:
    """API 520 gas/steam/liquid area and next API 526 effective orifice.

    API 521 scenario and fire-duty determination remain supplied engineering
    facts. Pressures: set gauge; back and atmosphere absolute, all Pa.
    Liquid viscosity/backpressure corrections are supplied explicitly.
    """
    r = ReliefValveInput.model_validate(request)
    needed = ['material', 'set_pressure_pa_g', 'relieving_case', 'case_basis', 'fire_case', 'back_pressure_pa_abs', 'atmospheric_pressure_pa', 'overpressure_fraction', 'mass_flow_kg_s', 'discharge_coefficient', 'backpressure_factor', 'combination_factor', 'factors_basis']
    if r.phase == 'gas':
        needed += ['temperature_k', 'compressibility', 'molecular_weight_g_mol', 'heat_capacity_ratio']
    elif r.phase == 'steam':
        needed += ['temperature_k', 'edition']
    else:
        needed += ['density_kg_m3', 'viscosity_factor', 'edition']
    holds = ['PSV_' + name.upper() + '_REQUIRED' for name in needed if getattr(r, name) is None]
    if (r.fire_case is True or (r.relieving_case and 'fire' in r.relieving_case.lower())) and (not r.fire_basis):
        holds.append('PSV_FIRE_BASIS_REQUIRED')
    method = 'API 520 Part I; API 521 scenario/load supplied; API 526 orifice rounding'
    if holds:
        return _result(r, method, holds)
    p1 = r.atmospheric_pressure_pa + r.set_pressure_pa_g * (1 + r.overpressure_fraction)
    if p1 <= r.back_pressure_pa_abs:
        raise ValueError('PSV_RELIEVING_PRESSURE_MUST_EXCEED_BACK_PRESSURE')
    f = _library('safety_valve')
    common = dict(m=r.mass_flow_kg_s, P1=p1, Kd=r.discharge_coefficient, Kc=r.combination_factor)
    if r.phase == 'gas':
        area = f.API520_A_g(**common, T=r.temperature_k, Z=r.compressibility, MW=r.molecular_weight_g_mol, k=r.heat_capacity_ratio, P2=r.back_pressure_pa_abs, Kb=r.backpressure_factor)
    elif r.phase == 'steam':
        area = f.API520_A_steam(**common, T=r.temperature_k, Kb=r.backpressure_factor, edition=r.edition)
    else:
        area = f.API520_A_l(**common, rho=r.density_kg_m3, P2=r.back_pressure_pa_abs, overpressure=r.overpressure_fraction, Kw=r.backpressure_factor, Kv=r.viscosity_factor, edition=r.edition)
    if not isfinite(area) or area <= 0:
        raise ValueError('PSV_NONPHYSICAL_AREA')
    selected = f.API520_round_size(area) if area <= max(f.API526_A) else None
    return _result(r, method, [] if selected else ['PSV_API526_RANGE_EXCEEDED'], material=r.material, relieving_pressure_pa_abs=p1, required_area_m2=area, selected_area_m2=selected, orifice_letter=f.API526_letters[f.API526_A.index(selected)] if selected else None)

class ReliefNoiseInput(Input):
    correlation: Literal['API521']
    mass_flow_kg_s: Positive
    inlet_pressure_pa_abs: Positive
    outlet_pressure_pa_abs: Positive
    sound_speed_m_s: Positive
    distance_m: Positive
    basis: Text

class ValveNoiseBase(Input):
    mass_flow_kg_s: Positive
    inlet_pressure_pa_abs: Positive
    outlet_pressure_pa_abs: Positive
    density_kg_m3: Positive
    kv_m3_h: Positive
    valve_diameter_m: Positive
    pipe_inside_diameter_m: Positive
    pipe_wall_m: Positive
    style_factor: Factor
    pipe_density_kg_m3: Positive
    pipe_sound_speed_m_s: Positive
    air_density_kg_m3: Positive
    air_sound_speed_m_s: Positive
    acoustic_efficiency: float
    basis: Text

class LiquidValveNoiseInput(ValveNoiseBase):
    correlation: Literal['IEC60534-8-4:2015']
    vapour_pressure_pa_abs: Nonnegative
    sound_speed_m_s: Positive
    recovery_factor: Factor
    cavitation_factor: Positive

class GasValveNoiseInput(ValveNoiseBase):
    correlation: Literal['IEC60534-8-3:2011']
    temperature_k: Positive
    outlet_temperature_k: Positive
    heat_capacity_ratio: Annotated[float, Field(gt=1)]
    molecular_weight_g_mol: Positive
    recovery_factor: Factor | None
    installed_recovery_factor: Factor | None
    piping_factor: Factor | None
    air_pressure_pa_abs: Positive
    strouhal_number: Positive
    beta: Positive
NoiseInput = Annotated[ReliefNoiseInput | LiquidValveNoiseInput | GasValveNoiseInput, Field(discriminator='correlation')]

def valve_noise(request: ReliefNoiseInput | LiquidValveNoiseInput | GasValveNoiseInput | dict) -> dict:
    """API 521 flare-tip noise at distance, or IEC downstream pipe valve noise.

    IEC returns A-weighted external sound pressure at its standard 1 m location;
    API 521 is unweighted flare-tip sound pressure, not an IEC valve prediction.
    """
    from pydantic import TypeAdapter
    r = TypeAdapter(NoiseInput).validate_python(request)
    if r.inlet_pressure_pa_abs <= r.outlet_pressure_pa_abs:
        raise ValueError('NOISE_INLET_PRESSURE_MUST_EXCEED_OUTLET')
    if isinstance(r, ReliefNoiseInput):
        level = _library('safety_valve').API521_noise(r.mass_flow_kg_s, r.inlet_pressure_pa_abs, r.outlet_pressure_pa_abs, r.sound_speed_m_s, r.distance_m)
        unit = 'dB (unweighted)'
    else:
        f = _library('control_valve')
        args = dict(m=r.mass_flow_kg_s, P1=r.inlet_pressure_pa_abs, P2=r.outlet_pressure_pa_abs, rho=r.density_kg_m3, Kv=r.kv_m3_h, d=r.valve_diameter_m, Di=r.pipe_inside_diameter_m, Fd=r.style_factor, t_pipe=r.pipe_wall_m, rho_pipe=r.pipe_density_kg_m3, c_pipe=r.pipe_sound_speed_m_s, rho_air=r.air_density_kg_m3, c_air=r.air_sound_speed_m_s, An=r.acoustic_efficiency)
        if isinstance(r, LiquidValveNoiseInput):
            if r.vapour_pressure_pa_abs >= r.outlet_pressure_pa_abs:
                raise ValueError('NOISE_FLASHING_LIQUID_OUTSIDE_MODEL')
            level = f.control_valve_noise_l_2015(**args, Psat=r.vapour_pressure_pa_abs, c=r.sound_speed_m_s, FL=r.recovery_factor, xFz=r.cavitation_factor)
        else:
            if r.recovery_factor is None and (r.installed_recovery_factor is None or r.piping_factor is None):
                raise ValueError('NOISE_RECOVERY_FACTORS_REQUIRED')
            level = f.control_valve_noise_g_2011(**args, T1=r.temperature_k, T2=r.outlet_temperature_k, gamma=r.heat_capacity_ratio, MW=r.molecular_weight_g_mol, FL=r.recovery_factor, FLP=r.installed_recovery_factor, FP=r.piping_factor, P_air=r.air_pressure_pa_abs, Stp=r.strouhal_number, beta=r.beta)
        unit = 'dB(A) at IEC 1 m location'
    if not isfinite(level):
        raise ValueError('NOISE_NONFINITE_RESULT')
    return _result(r, r.correlation, sound_pressure_level=level, unit=unit)

class EductorInput(Input):
    material: Text
    basis: Text
    motive_density_kg_m3: Positive
    suction_density_kg_m3: Positive
    nozzle_loss_coefficient: Nonnegative
    suction_loss_coefficient: Nonnegative
    mixing_loss_coefficient: Nonnegative
    diffuser_loss_coefficient: Nonnegative
    nozzle_retracted: bool
    diffuser_diameter_m: Positive
    nozzle_diameter_m: Positive | None = None
    mixing_diameter_m: Positive | None = None
    motive_flow_m3_s: Positive | None = None
    suction_flow_m3_s: Positive | None = None
    motive_pressure_pa_abs: Positive | None = None
    suction_pressure_pa_abs: Positive | None = None
    discharge_pressure_pa_abs: Positive | None = None

def eductor_size(request: EductorInput | dict) -> dict:
    """Single-phase liquid jet pump; exactly two unknowns, explicit diffuser.

    Uses fluids.jet_pump's energy/momentum model (ESDU), not a gas,
    flashing or cavitation-qualified ejector design.
    """
    r = EductorInput.model_validate(request)
    mapping = dict(d_nozzle='nozzle_diameter_m', d_mixing='mixing_diameter_m', Qp='motive_flow_m3_s', Qs='suction_flow_m3_s', P1='motive_pressure_pa_abs', P2='suction_pressure_pa_abs', P5='discharge_pressure_pa_abs')
    args = {key: getattr(r, value) for key, value in mapping.items()}
    if sum((value is None for value in args.values())) != 2:
        raise ValueError('EDUCTOR_EXACTLY_TWO_UNKNOWNS_REQUIRED')
    if r.nozzle_diameter_m is not None and r.mixing_diameter_m is not None and (r.nozzle_diameter_m >= r.mixing_diameter_m):
        raise ValueError('EDUCTOR_NOZZLE_MUST_BE_SMALLER_THAN_MIXER')
    result = _library('jet_pump').liquid_jet_pump(**args, rhop=r.motive_density_kg_m3, rhos=r.suction_density_kg_m3, Kp=r.nozzle_loss_coefficient, Ks=r.suction_loss_coefficient, Km=r.mixing_loss_coefficient, Kd=r.diffuser_loss_coefficient, d_diffuser=r.diffuser_diameter_m, nozzle_retracted=r.nozzle_retracted)
    if any((not isfinite(v) for v in result.values())) or not result['P1'] > result['P5'] > result['P2'] > 0 or (not 0 < result['d_nozzle'] < result['d_mixing'] <= result['d_diffuser']) or (not 0 < result['efficiency'] <= 1) or (min(result['Qp'], result['Qs']) <= 0):
        raise ValueError('EDUCTOR_NONPHYSICAL_SOLUTION')
    return _result(r, 'fluids.jet_pump.liquid_jet_pump; ESDU energy and momentum model', material=r.material, solution_si=result)
