"""Unified pipe sizing optimization for liquid and gas flow."""

from typing import Optional, Literal
import inspect
from tools.pipe_sizing import select_liquid_pipe_size_for_dp, select_gas_pipe_size_for_dp


def pipe_sizing(
    phase: Literal["liquid", "gas"] = "liquid",
    
    # Flow parameters
    flow_rate: Optional[float] = None,
    flow_rate_gpm: Optional[float] = None,
    flow_rate_kg_s: Optional[float] = None,
    flow_rate_norm_m3_hr: Optional[float] = None,
    flow_rate_std_m3_hr: Optional[float] = None,
    
    # Pressure constraints
    allowable_dp_pa: Optional[float] = None,
    allowable_dp_psi: Optional[float] = None,
    inlet_pressure: Optional[float] = None,
    inlet_pressure_psi: Optional[float] = None,
    
    # Pipe specifications
    pipe_length: Optional[float] = None,
    pipe_length_ft: Optional[float] = None,
    schedule: str = "40",
    nps_candidates: Optional[list] = None,
    material: Optional[str] = None,
    pipe_roughness: Optional[float] = None,
    fittings: Optional[list] = None,
    
    # Fluid properties
    fluid_density: Optional[float] = None,
    fluid_density_lbft3: Optional[float] = None,
    fluid_viscosity: Optional[float] = None,
    fluid_viscosity_cp: Optional[float] = None,
    fluid_name: Optional[str] = None,
    temperature_c: Optional[float] = None,
    pressure_bar: Optional[float] = None,
    
    # Gas-specific parameters
    gas_mw: Optional[float] = None,
    gas_gamma: Optional[float] = None,
    gas_z_factor: Optional[float] = None,
    gas_viscosity: Optional[float] = None,
    method: str = "isothermal",
    
    # Constraints
    velocity_min_m_s: Optional[float] = None,
    velocity_max_m_s: Optional[float] = None,
    
    # Output options
    return_all_candidates: bool = False,
) -> str:
    """Unified pipe sizing optimization.
    
    Finds optimal pipe size based on pressure drop constraints:
    - phase='liquid': Uses Darcy-Weisbach for incompressible flow
    - phase='gas': Uses compressible flow with equivalent length iteration
    
    Returns:
        JSON string with recommended pipe size and evaluation results
    """
    params = locals().copy()
    params.pop("phase")
    
    if phase == "liquid":
        fn = select_liquid_pipe_size_for_dp
    elif phase == "gas":
        fn = select_gas_pipe_size_for_dp
    else:
        return f'{{"error": "Invalid phase: {phase}"}}'
    
    sig = inspect.signature(fn)
    allowed = set(sig.parameters.keys())
    forwarded = {k: v for k, v in params.items() if k in allowed and v is not None}
    
    return fn(**forwarded)