"""Unified parameter sweep for various calculations."""

from typing import Optional, Literal, Dict, Any
import inspect
from tools.pipe_pressure_drop import pipe_pressure_drop_sweep
from tools.gas_pipe_pressure_drop import gas_pipe_sweep
from tools.blower_compressor import blower_sweep


def parameter_sweep(
    calculation: Literal["pipe_liquid", "pipe_gas", "blower"],
    variable: str,
    start: float,
    stop: float,
    n: int,
    
    # Pipe flow parameters (liquid and gas)
    flow_rate: Optional[float] = None,
    flow_rate_gpm: Optional[float] = None,
    flow_rate_kg_s: Optional[float] = None,
    flow_rate_norm_m3_hr: Optional[float] = None,
    flow_rate_std_m3_hr: Optional[float] = None,
    pipe_diameter: Optional[float] = None,
    pipe_diameter_in: Optional[float] = None,
    pipe_length: Optional[float] = None,
    pipe_length_ft: Optional[float] = None,
    nominal_size_in: Optional[float] = None,
    schedule: str = "40",
    material: Optional[str] = None,
    pipe_roughness: Optional[float] = None,
    
    # Fluid properties
    fluid_density: Optional[float] = None,
    fluid_density_lbft3: Optional[float] = None,
    fluid_viscosity: Optional[float] = None,
    fluid_viscosity_cp: Optional[float] = None,
    fluid_name: Optional[str] = None,
    temperature_c: Optional[float] = None,
    pressure_bar: Optional[float] = None,
    
    # Gas-specific parameters
    inlet_pressure: Optional[float] = None,
    outlet_pressure: Optional[float] = None,
    gas_mw: Optional[float] = None,
    gas_gamma: Optional[float] = None,
    gas_z_factor: Optional[float] = None,
    gas_viscosity: Optional[float] = None,
    method: str = "Weymouth",
    
    # Blower-specific parameters
    inlet_pressure_psi: Optional[float] = None,
    outlet_pressure_psi: Optional[float] = None,
    inlet_temperature_c: Optional[float] = None,
    efficiency: Optional[float] = None,
    efficiency_type: str = "isentropic",
    
    # Compressor calculation flags
    allow_property_defaults: Optional[bool] = None,
    use_eos_calculations: Optional[bool] = None,
    solve_for: Optional[str] = None,
) -> str:
    """Unified parameter sweep calculations.
    
    Performs parameter sweeps for different calculation types:
    - calculation='pipe_liquid': Liquid pipe pressure drop sweep
    - calculation='pipe_gas': Gas pipe pressure drop sweep  
    - calculation='blower': Blower/compressor performance sweep
    
    Args:
        calculation: Type of calculation to sweep
        variable: Variable to sweep (flow_rate, pressure, etc.)
        start: Start value for sweep
        stop: End value for sweep
        n: Number of points in sweep
        
    Returns:
        JSON string with sweep results
    """
    params = locals().copy()
    params.pop("calculation")
    
    # Preflight checks and auto-defaults
    if calculation == "blower":
        if inlet_temperature_c is None:
            return '{"error": "inlet_temperature_c is required for blower sweeps"}'
        # Auto-set allow_property_defaults if no Z-factor or fluid provided
        if allow_property_defaults is None and gas_z_factor is None and fluid_name is None:
            params["allow_property_defaults"] = True
    
    if calculation == "pipe_liquid":
        fn = pipe_pressure_drop_sweep
    elif calculation == "pipe_gas":
        fn = gas_pipe_sweep
    elif calculation == "blower":
        fn = blower_sweep
    else:
        return f'{{"error": "Invalid calculation: {calculation}"}}'
    
    sig = inspect.signature(fn)
    allowed = set(sig.parameters.keys())
    forwarded = {k: v for k, v in params.items() if k in allowed and v is not None}
    
    return fn(**forwarded)