"""Unified machine requirements for pumps and compressors."""

from typing import Optional, Literal
import inspect
from tools.pump_requirements import calculate_pump_requirements
from tools.blower_compressor import calculate_blower_compressor_requirements
from tools.reynolds_number import calculate_reynolds_number
from tools.open_channel_flow_new import calculate_open_channel_flow


def machine_requirements(
    machine_type: Literal["pump", "compressor", "reynolds", "open_channel"] = "pump",
    
    # Common flow parameters
    flow_rate: Optional[float] = None,
    flow_rate_gpm: Optional[float] = None,
    flow_rate_kg_s: Optional[float] = None,
    flow_rate_norm_m3_hr: Optional[float] = None,
    flow_rate_std_m3_hr: Optional[float] = None,
    flow_rate_m3s: Optional[float] = None,
    
    # Pump-specific parameters
    pipe_diameter: Optional[float] = None,
    pipe_diameter_in: Optional[float] = None,
    suction_pipe_length: Optional[float] = None,
    suction_pipe_length_ft: Optional[float] = None,
    discharge_pipe_length: Optional[float] = None,
    discharge_pipe_length_ft: Optional[float] = None,
    static_suction_head: Optional[float] = None,
    static_suction_head_ft: Optional[float] = None,
    static_discharge_head: Optional[float] = None,
    static_discharge_head_ft: Optional[float] = None,
    suction_fittings: Optional[list] = None,
    discharge_fittings: Optional[list] = None,
    suction_nozzle_diameter: Optional[float] = None,
    discharge_nozzle_diameter: Optional[float] = None,
    suction_nozzle_diameter_in: Optional[float] = None,
    discharge_nozzle_diameter_in: Optional[float] = None,
    atmospheric_pressure: Optional[float] = 101325.0,
    atmospheric_pressure_psi: Optional[float] = None,
    suction_tank_pressure: Optional[float] = None,
    suction_tank_pressure_psi: Optional[float] = None,
    target_tdh: Optional[float] = None,
    target_tdh_ft: Optional[float] = None,
    
    # Compressor-specific parameters
    inlet_pressure: Optional[float] = None,
    outlet_pressure: Optional[float] = None,
    inlet_pressure_psi: Optional[float] = None,
    outlet_pressure_psi: Optional[float] = None,
    inlet_temperature_c: Optional[float] = None,
    efficiency: Optional[float] = None,
    efficiency_type: str = "isentropic",
    polytropic_n: Optional[float] = None,
    
    # Fluid properties
    fluid_density: Optional[float] = None,
    fluid_density_lbft3: Optional[float] = None,
    fluid_viscosity: Optional[float] = None,
    fluid_viscosity_cp: Optional[float] = None,
    fluid_vapor_pressure: Optional[float] = None,
    fluid_vapor_pressure_psi: Optional[float] = None,
    fluid_name: Optional[str] = None,
    temperature_c: Optional[float] = None,
    pressure_bar: Optional[float] = None,
    pipe_roughness: Optional[float] = None,
    nominal_size_in: Optional[float] = None,
    schedule: str = "40",
    material: Optional[str] = None,
    
    # Gas properties
    gas_mw: Optional[float] = None,
    gas_gamma: Optional[float] = None,
    gas_z_factor: Optional[float] = None,
    gas_viscosity: Optional[float] = None,
    gas_composition_mol: Optional[dict] = None,
    
    # Compressor calculation flags
    allow_property_defaults: Optional[bool] = None,
    use_eos_calculations: Optional[bool] = None,
    solve_for: Optional[str] = None,
    
    # Reynolds number parameters
    velocity: Optional[float] = None,
    velocity_ft_s: Optional[float] = None,
    characteristic_length: Optional[float] = None,
    characteristic_length_in: Optional[float] = None,
    fluid_kinematic_viscosity: Optional[float] = None,
    
    # Open channel parameters
    channel_type: Optional[str] = None,
    width_m: Optional[float] = None,
    diameter_m: Optional[float] = None,
    side_slope_z: Optional[float] = None,
    slope: Optional[float] = None,
    manning_n: Optional[float] = None,
    normal_depth_m: Optional[float] = None,
) -> str:
    """Unified machine requirements and hydraulic calculations.
    
    Consolidates pump, compressor, and hydraulic calculations:
    - machine_type='pump': Calculate pump TDH and NPSHA requirements
    - machine_type='compressor': Calculate compressor power and discharge temperature
    - machine_type='reynolds': Calculate Reynolds number
    - machine_type='open_channel': Calculate open channel flow
    
    Returns:
        JSON string with calculation results
    
    Examples:
        Pump requirements:
        >>> machine_requirements(machine_type="pump", flow_rate_gpm=100,
        ...                     static_discharge_head_ft=50, fluid_name="Water")
        
        Compressor requirements:
        >>> machine_requirements(machine_type="compressor", flow_rate_norm_m3_hr=150,
        ...                     inlet_pressure_psi=14.7, outlet_pressure_psi=30)
        
        Reynolds number:
        >>> machine_requirements(machine_type="reynolds", velocity=2.5,
        ...                     characteristic_length=0.05, fluid_name="Water")
    """
    params = locals().copy()
    params.pop("machine_type")
    
    if machine_type == "pump":
        fn = calculate_pump_requirements
    elif machine_type == "compressor":
        fn = calculate_blower_compressor_requirements
    elif machine_type == "reynolds":
        fn = calculate_reynolds_number
    elif machine_type == "open_channel":
        fn = calculate_open_channel_flow
    else:
        return f'{{"error": "Invalid machine_type: {machine_type}"}}'
    
    sig = inspect.signature(fn)
    allowed = set(sig.parameters.keys())
    forwarded = {k: v for k, v in params.items() if k in allowed and v is not None}
    
    return fn(**forwarded)