"""
Unified pipe flow calculations for both liquid and gas phases.
"""

from typing import Optional, Literal, Dict
import inspect
from tools.pipe_pressure_drop import calculate_pipe_pressure_drop
from tools.gas_pipe_pressure_drop import calculate_gas_pipe_pressure_drop


def pipe_flow(
    phase: Literal["liquid", "gas"] = "liquid",
    
    # Common parameters
    pipe_length: Optional[float] = None,
    pipe_diameter: Optional[float] = None,
    pipe_diameter_in: Optional[float] = None,
    pipe_length_ft: Optional[float] = None,
    nominal_size_in: Optional[float] = None,
    schedule: str = "40",
    material: Optional[str] = None,
    pipe_roughness: Optional[float] = None,
    fittings: Optional[list] = None,
    temperature_c: Optional[float] = None,
    fluid_name: Optional[str] = None,
    
    # Liquid-specific parameters
    flow_rate: Optional[float] = None,
    flow_rate_gpm: Optional[float] = None,
    fluid_density: Optional[float] = None,
    fluid_density_lbft3: Optional[float] = None,
    fluid_viscosity: Optional[float] = None,
    fluid_viscosity_cp: Optional[float] = None,
    pressure_drop: Optional[float] = None,
    pressure_drop_psi: Optional[float] = None,
    pressure_bar: Optional[float] = None,
    
    # Gas-specific parameters
    inlet_pressure: Optional[float] = None,
    outlet_pressure: Optional[float] = None,
    flow_rate_kg_s: Optional[float] = None,
    flow_rate_norm_m3_hr: Optional[float] = None,
    flow_rate_std_m3_hr: Optional[float] = None,
    gas_mw: Optional[float] = None,
    gas_gamma: Optional[float] = None,
    gas_z_factor: Optional[float] = None,
    gas_viscosity: Optional[float] = None,
    gas_composition_mol: Optional[Dict[str, float]] = None,
    method: Literal[
        "Weymouth", "Panhandle_A", "Panhandle_B", "IGT", "Oliphant",
        "Spitzglass_low", "Spitzglass_high", "isothermal_darcy"
    ] = "Weymouth",
) -> str:
    """
    Unified pipe flow calculations for liquid and gas phases.
    
    This omnitool consolidates pipe pressure drop calculations:
    - phase='liquid': Uses incompressible flow equations (Darcy-Weisbach)
    - phase='gas': Uses compressible flow equations (various methods)
    
    Args:
        phase: "liquid" or "gas" - determines calculation method
        
        Common parameters (both phases):
            pipe_length: Pipe length in meters
            pipe_length_ft: Pipe length in feet
            pipe_diameter: Inner diameter in meters
            pipe_diameter_in: Inner diameter in inches
            nominal_size_in: Nominal pipe size in inches
            schedule: Pipe schedule (default "40")
            material: Pipe material for roughness lookup
            pipe_roughness: Absolute roughness in meters
            fittings: List of fittings with type and quantity
            temperature_c: Temperature in Celsius
            fluid_name: Fluid name for property lookup
        
        Liquid-specific parameters:
            flow_rate: Flow rate in m³/s
            flow_rate_gpm: Flow rate in US GPM
            fluid_density: Density in kg/m³
            fluid_density_lbft3: Density in lb/ft³
            fluid_viscosity: Dynamic viscosity in Pa·s
            fluid_viscosity_cp: Viscosity in centipoise
            pressure_drop: Pressure drop in Pa (for solving)
            pressure_drop_psi: Pressure drop in psi
            pressure_bar: Pressure in bar (for properties)
        
        Gas-specific parameters:
            inlet_pressure: Inlet pressure in Pa (absolute)
            outlet_pressure: Outlet pressure in Pa (absolute)
            flow_rate_kg_s: Mass flow rate in kg/s
            flow_rate_norm_m3_hr: Normal volumetric flow in m³/hr
            flow_rate_std_m3_hr: Standard volumetric flow in m³/hr
            gas_mw: Molecular weight in kg/kmol
            gas_gamma: Specific heat ratio (Cp/Cv)
            gas_z_factor: Compressibility factor
            gas_viscosity: Dynamic viscosity in Pa·s
            gas_composition_mol: Composition as mol fractions
            method: Calculation method for gas flow
    
    Returns:
        JSON string with calculation results
    
    Examples:
        Liquid flow:
        >>> pipe_flow(phase="liquid", flow_rate_gpm=100, pipe_length_ft=100, 
        ...          nominal_size_in=2, fluid_name="Water", temperature_c=25)
        
        Gas flow:
        >>> pipe_flow(phase="gas", flow_rate_norm_m3_hr=150, inlet_pressure=101325,
        ...          pipe_length=25, nominal_size_in=3, fluid_name="Methane")
    """
    # Collect all local variables
    params = locals().copy()
    params.pop("phase")
    
    # Preflight checks for required parameters
    if phase == "gas" and temperature_c is None:
        return '{"error": "temperature_c is required for gas flow calculations"}'
    
    # Select the appropriate function based on phase
    if phase == "liquid":
        fn = calculate_pipe_pressure_drop
    elif phase == "gas":
        fn = calculate_gas_pipe_pressure_drop
    else:
        return f'{{"error": "Invalid phase: {phase}. Must be \'liquid\' or \'gas\'"}}'
    
    # Filter parameters to only those accepted by the target function
    sig = inspect.signature(fn)
    allowed = set(sig.parameters.keys())
    forwarded = {k: v for k, v in params.items() if k in allowed and v is not None}
    
    # Call the appropriate function
    return fn(**forwarded)