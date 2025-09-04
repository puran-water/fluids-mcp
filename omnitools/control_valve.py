"""Unified control valve sizing for liquid and gas service."""

from typing import Optional, Literal
import inspect
from tools.liquid_control_valve import calculate_liquid_control_valve
from tools.gas_control_valve import calculate_gas_control_valve


def control_valve(
    phase: Literal["liquid", "gas"] = "liquid",
    
    # Flow parameters
    flow_rate: Optional[float] = None,
    flow_rate_gpm: Optional[float] = None,
    flow_rate_kg_s: Optional[float] = None,
    flow_rate_norm_m3_hr: Optional[float] = None,
    flow_rate_std_m3_hr: Optional[float] = None,
    
    # Pressure parameters
    inlet_pressure: Optional[float] = None,
    outlet_pressure: Optional[float] = None,
    inlet_pressure_psi: Optional[float] = None,
    outlet_pressure_psi: Optional[float] = None,
    
    # Fluid properties
    fluid_density: Optional[float] = None,
    fluid_density_lbft3: Optional[float] = None,
    fluid_viscosity: Optional[float] = None,
    fluid_viscosity_cp: Optional[float] = None,
    fluid_name: Optional[str] = None,
    temperature_c: Optional[float] = None,
    pressure_bar: Optional[float] = None,
    
    # Liquid-specific valve parameters
    fluid_saturation_pressure_pa: Optional[float] = None,
    fluid_saturation_pressure_psi: Optional[float] = None,
    fluid_critical_pressure_pa: Optional[float] = None,
    fluid_critical_pressure_psi: Optional[float] = None,
    valve_fl: Optional[float] = None,
    valve_fd: Optional[float] = None,
    
    # Gas-specific parameters
    inlet_temperature_c: Optional[float] = None,
    gas_mw: Optional[float] = None,
    gas_gamma: Optional[float] = None,
    gas_z_factor: Optional[float] = None,
    gas_viscosity: Optional[float] = None,
    valve_xt: Optional[float] = None,
    
    # Common valve parameters
    valve_type: Literal["globe", "butterfly", "ball", "gate", "other"] = "globe",
    size_units: Literal["m", "inch", "mm"] = "inch",
) -> str:
    """Unified control valve sizing calculations.
    
    Consolidates valve sizing for both liquid and gas service:
    - phase='liquid': Sizes valves for incompressible flow
    - phase='gas': Sizes valves for compressible flow
    
    Returns:
        JSON string with Cv/Kv values and recommended size
    """
    params = locals().copy()
    params.pop("phase")
    
    # Preflight checks for required parameters
    if phase == "gas" and inlet_temperature_c is None:
        return '{"error": "inlet_temperature_c is required for gas control valve calculations"}'
    
    if phase == "liquid":
        fn = calculate_liquid_control_valve
    elif phase == "gas":
        fn = calculate_gas_control_valve
    else:
        return f'{{"error": "Invalid phase: {phase}"}}'
    
    sig = inspect.signature(fn)
    allowed = set(sig.parameters.keys())
    forwarded = {k: v for k, v in params.items() if k in allowed and v is not None}
    
    return fn(**forwarded)