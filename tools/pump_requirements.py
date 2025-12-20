import json
import logging
import math
from typing import Dict, List, Optional, Union

# --- Third-Party Libraries ---
import fluids
import fluids.core
import fluids.friction
import fluids.piping
# Note: fluids.friction._roughness is imported dynamically within the function

# Import shared utilities
from utils.constants import (
    GPM_to_M3S, INCH_to_M, FT_to_M, PSI_to_PA, 
    LBFT3_to_KGM3, CENTIPOISE_to_PAS, DEFAULT_ROUGHNESS, DEG_C_to_K, G_GRAVITY
)
from utils.helpers import get_fitting_K
from utils.import_helpers import FLUIDPROP_AVAILABLE, FluidProperties, FLUID_SELECTION, COOLPROP_AVAILABLE, CP


# Configure logging
logger = logging.getLogger("fluids-mcp.calculate_pump_requirements")

# --- Main Function Definition Starts Below ---
def calculate_pump_requirements(
    # --- Core SI Inputs ---
    flow_rate: Optional[float] = None,          # Flow rate in m³/s
    pipe_diameter: Optional[float] = None,      # Pipe inner diameter in m
    suction_pipe_length: Optional[float] = None, # Suction pipe length in m
    discharge_pipe_length: Optional[float] = None, # Discharge pipe length in m
    static_suction_head: Optional[float] = None,  # Static suction head in m
    static_discharge_head: Optional[float] = None, # Static discharge head in m
    fluid_density: Optional[float] = None,      # Fluid density in kg/m³
    fluid_viscosity: Optional[float] = None,    # Fluid viscosity in Pa·s
    fluid_vapor_pressure: Optional[float] = None, # Fluid vapor pressure in Pa
    pipe_roughness: Optional[float] = None,     # Pipe roughness in m
    
    # --- Nozzle Diameters (for velocity head differential) ---
    suction_nozzle_diameter: Optional[float] = None,  # Pump suction nozzle diameter in m
    discharge_nozzle_diameter: Optional[float] = None, # Pump discharge nozzle diameter in m

    # --- Alternative Unit Inputs ---
    flow_rate_gpm: Optional[float] = None,        # Flow rate in US GPM
    pipe_diameter_in: Optional[float] = None,     # Pipe diameter in inches
    suction_pipe_length_ft: Optional[float] = None,# Suction pipe length in feet
    discharge_pipe_length_ft: Optional[float] = None, # Discharge pipe length in feet
    static_suction_head_ft: Optional[float] = None,# Static suction head in feet
    static_discharge_head_ft: Optional[float] = None, # Static discharge head in feet
    fluid_density_lbft3: Optional[float] = None,   # Fluid density in lb/ft³
    fluid_viscosity_cp: Optional[float] = None,    # Fluid viscosity in centipoise
    fluid_vapor_pressure_psi: Optional[float] = None, # Fluid vapor pressure in psi
    
    # --- Nozzle Diameters (Imperial) ---
    suction_nozzle_diameter_in: Optional[float] = None,  # Pump suction nozzle diameter in inches
    discharge_nozzle_diameter_in: Optional[float] = None, # Pump discharge nozzle diameter in inches

    # --- Property Lookup Inputs ---
    fluid_name: Optional[str] = None,            # e.g., "Water", "Air"
    temperature_c: Optional[float] = None,         # Temperature in Celsius
    pressure_bar: Optional[float] = None,          # Pressure in bar (optional for lookup)
    nominal_size_in: Optional[float] = None,       # Nominal pipe size in inches
    schedule: str = "40",                         # Pipe schedule
    material: Optional[str] = None,             # e.g., "Steel", "PVC"

    # --- Other Parameters ---
    atmospheric_pressure: float = 101325.0,       # Atmospheric pressure in Pa (can be overridden)
    atmospheric_pressure_psi: Optional[float] = None, # Allow overriding atm pressure in psi
    suction_tank_pressure: Optional[float] = None,    # Suction tank pressure in Pa (absolute, overrides atmospheric)
    suction_tank_pressure_psi: Optional[float] = None, # Suction tank pressure in psi (absolute)
    target_tdh: Optional[float] = None,            # Target TDH in m (for solving flow rate or diameter)
    target_tdh_ft: Optional[float] = None,         # Target TDH in ft (for solving flow rate or diameter) 
    suction_fittings: List[Dict[str, Union[str, int, float]]] = None,  # Suction line fittings
    discharge_fittings: List[Dict[str, Union[str, int, float]]] = None # Discharge line fittings
) -> str:
    """Calculate pump TDH and NPSHA requirements with flexible input options.

    Accepts inputs in SI units, common imperial units, or via property lookups.
    Priority: Explicit SI > Explicit Imperial > Lookup > Default.

    Args:
        flow_rate: Flow rate in m³/s
        pipe_diameter: Pipe inner diameter in m
        suction_pipe_length: Suction pipe length in m
        discharge_pipe_length: Discharge pipe length in m
        static_suction_head: Static suction head in m (relative to pump centerline, negative if below)
        static_discharge_head: Static discharge head in m (relative to pump centerline)
        fluid_density: Fluid density in kg/m³
        fluid_viscosity: Fluid viscosity in Pa·s
        fluid_vapor_pressure: Fluid vapor pressure in Pa
        pipe_roughness: Pipe roughness in m (default: lookup by material or 1.5e-5)
        flow_rate_gpm: Flow rate in US GPM
        pipe_diameter_in: Pipe diameter in inches
        suction_pipe_length_ft: Suction pipe length in feet
        discharge_pipe_length_ft: Discharge pipe length in feet
        static_suction_head_ft: Static suction head in feet
        static_discharge_head_ft: Static discharge head in feet
        fluid_density_lbft3: Fluid density in lb/ft³
        fluid_viscosity_cp: Fluid viscosity in centipoise
        fluid_vapor_pressure_psi: Fluid vapor pressure in psi
        fluid_name: Name of the fluid (e.g., "Water") for property lookup
        temperature_c: Temperature in Celsius (required with fluid_name)
        pressure_bar: Pressure in bar (optional for fluid lookup, defaults ~1 atm)
        nominal_size_in: Nominal pipe size in inches (for diameter lookup)
        schedule: Pipe schedule (default: "40", used with nominal_size_in)
        material: Pipe material name (e.g., "Steel", "PVC") for roughness lookup
        atmospheric_pressure: Atmospheric pressure in Pa (default: 101325.0)
        atmospheric_pressure_psi: Atmospheric pressure in psi (overrides Pa value if provided)
        suction_tank_pressure: Suction tank pressure in Pa (absolute, overrides atmospheric for closed systems)
        suction_tank_pressure_psi: Suction tank pressure in psi (absolute, overrides atmospheric for closed systems)
        target_tdh: Target TDH in m (for solving flow rate or diameter)
        target_tdh_ft: Target TDH in ft (for solving flow rate or diameter)
        suction_fittings: List of suction line fittings with type, quantity, and optionally K_value
        discharge_fittings: List of discharge line fittings with type, quantity, and optionally K_value
        suction_nozzle_diameter: Pump suction nozzle diameter in m (for velocity head calculation)
        discharge_nozzle_diameter: Pump discharge nozzle diameter in m (for velocity head calculation)
        suction_nozzle_diameter_in: Pump suction nozzle diameter in inches
        discharge_nozzle_diameter_in: Pump discharge nozzle diameter in inches

    Returns:
        JSON string with detailed pump sizing results including TDH, NPSHA, power estimate, and input resolution log.
    """
    results_log = [] # To track how values were obtained
    error_log = []
    fluid_info = {} # Store details from lookup
    pipe_info = {} # Store details from lookup

    try:
        # 1. Resolve Flow Rate (SI)
        local_flow_rate = None
        if flow_rate is not None:
            local_flow_rate = flow_rate
            results_log.append("Used provided SI flow_rate.")
        elif flow_rate_gpm is not None:
            local_flow_rate = flow_rate_gpm * GPM_to_M3S
            results_log.append(f"Converted flow_rate from {flow_rate_gpm} GPM.")
        else:
            error_log.append("Missing required input: flow_rate or flow_rate_gpm.")


        # 2. Resolve Pipe Diameter & Roughness (SI) - Reuse logic from calculate_pipe_pressure_drop
        local_pipe_diameter = None
        local_pipe_roughness = None
        # (Copying the exact resolution block from calculate_pipe_pressure_drop for consistency)
        if pipe_diameter is not None:
            local_pipe_diameter = pipe_diameter
            results_log.append("Used provided SI pipe_diameter.")
            if pipe_roughness is not None:
                local_pipe_roughness = pipe_roughness
                results_log.append("Used provided SI pipe_roughness.")
            elif material is not None:
                try:
                    from fluids.friction import _roughness
                    mat_lower = {k.lower(): k for k in _roughness.keys()}
                    if material.lower() in mat_lower:
                        actual_key = mat_lower[material.lower()]
                        local_pipe_roughness = _roughness[actual_key]
                        results_log.append(f"Looked up roughness for material '{actual_key}'.")
                    else:
                        local_pipe_roughness = 1.5e-5
                        results_log.append(f"Material '{material}' not found, using default roughness 1.5e-5 m.")
                        error_log.append(f"Warning: Material '{material}' not found for roughness lookup.")
                except Exception as lookup_e:
                    local_pipe_roughness = 1.5e-5
                    results_log.append(f"Roughness lookup failed ({lookup_e}), using default 1.5e-5 m.")
                    error_log.append(f"Warning: Roughness lookup failed: {lookup_e}")
            else:
                local_pipe_roughness = 1.5e-5
                results_log.append("Used default pipe_roughness 1.5e-5 m.")
        elif pipe_diameter_in is not None:
            local_pipe_diameter = pipe_diameter_in * INCH_to_M
            results_log.append(f"Converted pipe_diameter from {pipe_diameter_in} inches.")
            if pipe_roughness is not None:
                local_pipe_roughness = pipe_roughness
                results_log.append("Used provided SI pipe_roughness.")
            elif material is not None:
                try:
                    from fluids.friction import _roughness
                    mat_lower = {k.lower(): k for k in _roughness.keys()}
                    if material.lower() in mat_lower:
                        actual_key = mat_lower[material.lower()]
                        local_pipe_roughness = _roughness[actual_key]
                        results_log.append(f"Looked up roughness for material '{actual_key}'.")
                    else:
                        local_pipe_roughness = 1.5e-5
                        results_log.append(f"Material '{material}' not found, using default roughness 1.5e-5 m.")
                        error_log.append(f"Warning: Material '{material}' not found for roughness lookup.")
                except Exception as lookup_e:
                    local_pipe_roughness = 1.5e-5
                    results_log.append(f"Roughness lookup failed ({lookup_e}), using default 1.5e-5 m.")
                    error_log.append(f"Warning: Roughness lookup failed: {lookup_e}")
            else:
                local_pipe_roughness = 1.5e-5
                results_log.append("Used default pipe_roughness 1.5e-5 m.")
        elif nominal_size_in is not None:
            try:
                NPS, Di, Do, t = fluids.piping.nearest_pipe(NPS=nominal_size_in, schedule=schedule)
                local_pipe_diameter = Di
                pipe_info = {"NPS_in": NPS, "schedule": schedule, "Do_m": Do, "t_m": t, "Di_m": Di}
                results_log.append(f"Looked up pipe dimensions for NPS {nominal_size_in} Sch {schedule}.")
                if pipe_roughness is not None:
                    local_pipe_roughness = pipe_roughness
                    results_log.append("Used provided SI pipe_roughness.")
                elif material is not None:
                    try:
                        from fluids.friction import _roughness
                        mat_lower = {k.lower(): k for k in _roughness.keys()}
                        if material.lower() in mat_lower:
                            actual_key = mat_lower[material.lower()]
                            local_pipe_roughness = _roughness[actual_key]
                            results_log.append(f"Looked up roughness for material '{actual_key}'.")
                            pipe_info["roughness_m"] = local_pipe_roughness
                            pipe_info["roughness_source"] = "Material Lookup"
                        else:
                            local_pipe_roughness = 1.5e-5
                            results_log.append(f"Material '{material}' not found, using default roughness 1.5e-5 m.")
                            error_log.append(f"Warning: Material '{material}' not found for roughness lookup.")
                            pipe_info["roughness_m"] = local_pipe_roughness
                            pipe_info["roughness_source"] = "Default (Material Not Found)"
                    except Exception as lookup_e:
                        local_pipe_roughness = 1.5e-5
                        results_log.append(f"Roughness lookup failed ({lookup_e}), using default 1.5e-5 m.")
                        error_log.append(f"Warning: Roughness lookup failed: {lookup_e}")
                        pipe_info["roughness_m"] = local_pipe_roughness
                        pipe_info["roughness_source"] = "Default (Lookup Failed)"
                else:
                    local_pipe_roughness = 1.5e-5
                    results_log.append("Used default pipe_roughness 1.5e-5 m (material not specified).")
                    pipe_info["roughness_m"] = local_pipe_roughness
                    pipe_info["roughness_source"] = "Default (Material Not Specified)"
            except Exception as pipe_lookup_e:
                local_pipe_diameter = None
                local_pipe_roughness = None
                error_log.append(f"Failed to look up pipe NPS {nominal_size_in} Sch {schedule}: {pipe_lookup_e}")
        else:
            error_log.append("Missing required input: pipe_diameter, pipe_diameter_in, or nominal_size_in.")


        # 3. Resolve Pipe Lengths (SI)
        local_suction_pipe_length = None
        if suction_pipe_length is not None:
            local_suction_pipe_length = suction_pipe_length
            results_log.append("Used provided SI suction_pipe_length.")
        elif suction_pipe_length_ft is not None:
            local_suction_pipe_length = suction_pipe_length_ft * FT_to_M
            results_log.append(f"Converted suction_pipe_length from {suction_pipe_length_ft} ft.")
        else:
            error_log.append("Missing required input: suction_pipe_length or suction_pipe_length_ft.")


        local_discharge_pipe_length = None
        if discharge_pipe_length is not None:
            local_discharge_pipe_length = discharge_pipe_length
            results_log.append("Used provided SI discharge_pipe_length.")
        elif discharge_pipe_length_ft is not None:
            local_discharge_pipe_length = discharge_pipe_length_ft * FT_to_M
            results_log.append(f"Converted discharge_pipe_length from {discharge_pipe_length_ft} ft.")
        else:
            error_log.append("Missing required input: discharge_pipe_length or discharge_pipe_length_ft.")


        # 4. Resolve Static Heads (SI)
        local_static_suction_head = None
        if static_suction_head is not None:
            local_static_suction_head = static_suction_head
            results_log.append("Used provided SI static_suction_head.")
        elif static_suction_head_ft is not None:
            local_static_suction_head = static_suction_head_ft * FT_to_M
            results_log.append(f"Converted static_suction_head from {static_suction_head_ft} ft.")
        else:
            error_log.append("Missing required input: static_suction_head or static_suction_head_ft.")


        local_static_discharge_head = None
        if static_discharge_head is not None:
            local_static_discharge_head = static_discharge_head
            results_log.append("Used provided SI static_discharge_head.")
        elif static_discharge_head_ft is not None:
            local_static_discharge_head = static_discharge_head_ft * FT_to_M
            results_log.append(f"Converted static_discharge_head from {static_discharge_head_ft} ft.")
        else:
            error_log.append("Missing required input: static_discharge_head or static_discharge_head_ft.")


        # 5. Resolve Fluid Properties (SI) - Density, Viscosity, Vapor Pressure
        local_fluid_density = None
        local_fluid_viscosity = None
        local_fluid_vapor_pressure = None # Initialize


        # 5.1 Density and Viscosity (SI, Imperial, or Lookup)
        if fluid_density is not None and fluid_viscosity is not None:
            local_fluid_density = fluid_density
            local_fluid_viscosity = fluid_viscosity
            results_log.append("Used provided SI fluid density and viscosity.")
        elif fluid_density_lbft3 is not None and fluid_viscosity_cp is not None:
              local_fluid_density = fluid_density_lbft3 * LBFT3_to_KGM3
              local_fluid_viscosity = fluid_viscosity_cp * CENTIPOISE_to_PAS
              results_log.append(f"Converted density from {fluid_density_lbft3} lb/ft³ and viscosity from {fluid_viscosity_cp} cP.")
        elif fluid_name is not None and temperature_c is not None:
            if FLUIDPROP_AVAILABLE and FLUID_SELECTION is not None and FluidProperties is not None:
                try: # Fluid property lookup
                    try:
                        valid_fluids = [f[0] for f in FLUID_SELECTION if f is not None and hasattr(f, '__getitem__')]
                    except (TypeError, IndexError):
                        valid_fluids = []
                    actual_fluid_name = fluid_name
                    if not valid_fluids or fluid_name not in valid_fluids:
                        fluid_lower = fluid_name.lower()
                        match = next((f for f in valid_fluids if f.lower() == fluid_lower), None)
                        if match: actual_fluid_name = match
                        else: raise ValueError(f"Fluid '{fluid_name}' not found.")


                    p_bar = pressure_bar if pressure_bar is not None else 1.01325 # Default pressure
                    fluid_props = FluidProperties(coolprop_name=actual_fluid_name, T_in_deg_C=temperature_c, P_in_bar=p_bar)
                    local_fluid_density = float(fluid_props.rho[0])
                    local_fluid_viscosity = float(fluid_props.eta[0])
                    fluid_info = { # Store details
                        "name_used": actual_fluid_name,
                        "temperature_c": temperature_c,
                        "pressure_bar": p_bar,
                        "density_kgm3": local_fluid_density,
                        "viscosity_pas": local_fluid_viscosity,
                        "kinematic_viscosity_m2s": float(fluid_props.nu[0])
                    }
                    results_log.append(f"Looked up density/viscosity for {actual_fluid_name} at {temperature_c}°C, {p_bar} bar.")


                    # Attempt vapor pressure lookup ONLY if density/viscosity lookup succeeded
                    if COOLPROP_AVAILABLE:
                        try:
                            temp_k = temperature_c + DEG_C_to_K # Corrected constant usage
                            # Use saturation pressure at the given temperature
                            local_fluid_vapor_pressure = CP.PropsSI('P', 'T', temp_k, 'Q', 0, actual_fluid_name)
                            fluid_info["vapor_pressure_pa"] = local_fluid_vapor_pressure
                            fluid_info["vapor_pressure_source"] = "CoolProp Lookup"
                            results_log.append(f"Looked up vapor pressure ({local_fluid_vapor_pressure:.2f} Pa) for {actual_fluid_name} at {temperature_c}°C.")
                        except Exception as vapor_e:
                            error_log.append(f"Warning: Could not look up vapor pressure for {actual_fluid_name}: {vapor_e}")
                            fluid_info["vapor_pressure_source"] = "Lookup Failed"
                    else:
                        error_log.append("Warning: Vapor pressure lookup skipped: CoolProp package not available.")


                except Exception as fluid_lookup_e:
                    local_fluid_density = None
                    local_fluid_viscosity = None
                    error_log.append(f"Failed to look up base fluid properties for {fluid_name} at {temperature_c}°C: {fluid_lookup_e}")
            else:
                error_log.append("Fluid property lookup skipped: fluidprop package not available.")
        else:
            # This error is raised only if density/viscosity cannot be determined
            error_log.append("Missing inputs for fluid density/viscosity: Need (SI pair) OR (Imperial pair) OR (fluid_name AND temperature_c).")


        # 5.2 Vapor Pressure (SI or Imperial, only if not found via lookup)
        if local_fluid_vapor_pressure is None:
            if fluid_vapor_pressure is not None:
                local_fluid_vapor_pressure = fluid_vapor_pressure
                results_log.append("Used provided SI fluid_vapor_pressure.")
                if "vapor_pressure_source" not in fluid_info: fluid_info["vapor_pressure_source"] = "Provided SI"
            elif fluid_vapor_pressure_psi is not None:
                local_fluid_vapor_pressure = fluid_vapor_pressure_psi * PSI_to_PA
                results_log.append(f"Converted fluid_vapor_pressure from {fluid_vapor_pressure_psi} psi.")
                if "vapor_pressure_source" not in fluid_info: fluid_info["vapor_pressure_source"] = "Provided Imperial"
            else:
                # This error is critical only if density/viscosity *were* found previously
                if local_fluid_density is not None:
                    error_log.append("Missing required input for NPSHa: fluid_vapor_pressure or fluid_vapor_pressure_psi (could not be looked up).")
                # If density/visc also failed, the previous error message is sufficient


        # 6. Resolve Suction Side Pressure (Atmospheric or Tank Pressure)
        local_suction_side_pressure = atmospheric_pressure # Default SI
        suction_pressure_source = "Atmospheric (default)"
        
        # Check for tank pressure first (overrides atmospheric)
        if suction_tank_pressure is not None:
            local_suction_side_pressure = suction_tank_pressure
            suction_pressure_source = "Tank pressure (SI)"
            results_log.append(f"Used provided suction_tank_pressure ({suction_tank_pressure} Pa).")
        elif suction_tank_pressure_psi is not None:
            local_suction_side_pressure = suction_tank_pressure_psi * PSI_to_PA
            suction_pressure_source = "Tank pressure (converted from psi)"
            results_log.append(f"Used provided suction_tank_pressure_psi ({suction_tank_pressure_psi} psi).")
        elif atmospheric_pressure_psi is not None:
            local_suction_side_pressure = atmospheric_pressure_psi * PSI_to_PA
            suction_pressure_source = "Atmospheric (converted from psi)"
            results_log.append(f"Used provided atmospheric_pressure_psi ({atmospheric_pressure_psi} psi).")
        else:
            suction_pressure_source = "Atmospheric (default SI)"
            results_log.append(f"Used default atmospheric_pressure ({atmospheric_pressure} Pa).")

        # 7. Resolve Target TDH (for solving scenarios)
        local_target_tdh = None
        if target_tdh is not None:
            local_target_tdh = target_tdh
            results_log.append("Used provided SI target_tdh.")
        elif target_tdh_ft is not None:
            local_target_tdh = target_tdh_ft * FT_to_M
            results_log.append(f"Converted target_tdh from {target_tdh_ft} ft.")

        # 8. Determine solve_for variable (basic capability)
        specified_vars = {
            'flow_rate': local_flow_rate,
            'pipe_diameter': local_pipe_diameter
        }
        unknown_vars = [k for k, v in specified_vars.items() if v is None]
        
        solve_for = None
        if local_target_tdh is not None and len(unknown_vars) == 1:
            solve_for = unknown_vars[0]
            results_log.append(f"Target TDH specified - solving for: {solve_for}")
        elif len(unknown_vars) > 1 and local_target_tdh is not None:
            error_log.append(f"Cannot solve: too many unknowns with target TDH. Specify flow_rate OR pipe_diameter when using target_tdh.")

        # --- Check for any errors before proceeding ---
        missing_critical = any(e.startswith("Missing required input") for e in error_log)
        if missing_critical:
            return json.dumps({"errors": error_log, "log": results_log})


        # Check if any critical *local* variables failed resolution
        critical_vars = [
            local_flow_rate, local_pipe_diameter, local_suction_pipe_length,
            local_discharge_pipe_length, local_static_suction_head, local_static_discharge_head,
            local_fluid_density, local_fluid_viscosity, local_fluid_vapor_pressure,
            local_pipe_roughness, local_suction_side_pressure
        ]
        if any(v is None for v in critical_vars):
            error_log.append("Critical input parameter resolution failed after attempting conversions/lookups.")
            # Add specific missing items for better debug
            missing_items = [
                name for name, val in zip([
                    "flow_rate", "pipe_diameter", "suction_length", "discharge_length",
                    "static_suction_head", "static_discharge_head", "density", "viscosity",
                    "vapor_pressure", "roughness", "suction_side_pressure"
                ], critical_vars) if val is None
            ]
            error_log.append(f"Missing resolved values for: {', '.join(missing_items)}")
            return json.dumps({"errors": error_log, "log": results_log})


        # --- Resolve Nozzle Diameters for Velocity Head Differential ---
        local_suction_nozzle_diameter = None
        local_discharge_nozzle_diameter = None
        
        if suction_nozzle_diameter is not None:
            local_suction_nozzle_diameter = suction_nozzle_diameter
            results_log.append("Used provided SI suction_nozzle_diameter.")
        elif suction_nozzle_diameter_in is not None:
            local_suction_nozzle_diameter = suction_nozzle_diameter_in * INCH_to_M
            results_log.append(f"Converted suction_nozzle_diameter from {suction_nozzle_diameter_in} inches.")
        else:
            # Default: assume nozzle diameter equals pipe diameter
            local_suction_nozzle_diameter = local_pipe_diameter
            results_log.append("Using pipe diameter as suction nozzle diameter (default assumption).")
            
        if discharge_nozzle_diameter is not None:
            local_discharge_nozzle_diameter = discharge_nozzle_diameter
            results_log.append("Used provided SI discharge_nozzle_diameter.")
        elif discharge_nozzle_diameter_in is not None:
            local_discharge_nozzle_diameter = discharge_nozzle_diameter_in * INCH_to_M
            results_log.append(f"Converted discharge_nozzle_diameter from {discharge_nozzle_diameter_in} inches.")
        else:
            # Default: assume nozzle diameter equals pipe diameter
            local_discharge_nozzle_diameter = local_pipe_diameter
            results_log.append("Using pipe diameter as discharge nozzle diameter (default assumption).")

        # --- Core Calculation Logic ---
        suction_fittings_list = suction_fittings if suction_fittings else []
        discharge_fittings_list = discharge_fittings if discharge_fittings else []

        # --- Calculate Pipe and Nozzle Properties ---
        pipe_area = math.pi * (local_pipe_diameter / 2) ** 2
        if pipe_area == 0:
            return json.dumps({"error": "Pipe diameter resolved to zero.", "log": results_log})
        pipe_velocity = local_flow_rate / pipe_area
        
        suction_nozzle_area = math.pi * (local_suction_nozzle_diameter / 2) ** 2
        discharge_nozzle_area = math.pi * (local_discharge_nozzle_diameter / 2) ** 2
        
        if suction_nozzle_area == 0 or discharge_nozzle_area == 0:
            return json.dumps({"error": "Nozzle diameter resolved to zero.", "log": results_log})
            
        suction_nozzle_velocity = local_flow_rate / suction_nozzle_area
        discharge_nozzle_velocity = local_flow_rate / discharge_nozzle_area

        Re = fluids.core.Reynolds(V=pipe_velocity, D=local_pipe_diameter, rho=local_fluid_density, mu=local_fluid_viscosity)
        fd = fluids.friction.friction_factor(Re=Re, eD=local_pipe_roughness / local_pipe_diameter)


        # Calculate Suction Line Losses using helper
        K_suction_pipe = fd * local_suction_pipe_length / local_pipe_diameter
        K_suction_fittings = 0
        suction_fitting_details = []
        for fitting in suction_fittings_list:
            fitting_type = fitting.get("type", "").lower()
            quantity = fitting.get("quantity", 1)
            K_value_provided = fitting.get("K_value")
            K_s_fitting = 0
            source = "Error"
            try:
                if K_value_provided is not None:
                    K_s_fitting = float(K_value_provided)
                    source = "Provided"
                else:
                    K_s_fitting = get_fitting_K(fitting_type, local_pipe_diameter, Re, local_flow_rate)
                    source = "Calculated"
                K_suction_fittings += K_s_fitting * quantity
                suction_fitting_details.append({"type": fitting_type, "quantity": quantity, "K_source": source, "K_individual": round(K_s_fitting,4), "K_total": round(K_s_fitting * quantity, 4)})
            except ValueError:
                error_log.append(f"Invalid K_value '{K_value_provided}' for suction fitting type '{fitting_type}'. Skipping fitting.")
                suction_fitting_details.append({"type": fitting_type, "quantity": quantity, "K_source": "Error", "error": "Invalid K_value provided"})


        # Calculate Discharge Line Losses using helper
        K_discharge_pipe = fd * local_discharge_pipe_length / local_pipe_diameter
        K_discharge_fittings = 0
        discharge_fitting_details = []
        for fitting in discharge_fittings_list:
            fitting_type = fitting.get("type", "").lower()
            quantity = fitting.get("quantity", 1)
            K_value_provided = fitting.get("K_value")
            K_d_fitting = 0
            source = "Error"
            try:
                if K_value_provided is not None:
                    K_d_fitting = float(K_value_provided)
                    source = "Provided"
                else:
                    K_d_fitting = get_fitting_K(fitting_type, local_pipe_diameter, Re, local_flow_rate)
                    source = "Calculated"
                K_discharge_fittings += K_d_fitting * quantity
                discharge_fitting_details.append({"type": fitting_type, "quantity": quantity, "K_source": source, "K_individual": round(K_d_fitting,4), "K_total": round(K_d_fitting * quantity, 4)})
            except ValueError:
                error_log.append(f"Invalid K_value '{K_value_provided}' for discharge fitting type '{fitting_type}'. Skipping fitting.")
                discharge_fitting_details.append({"type": fitting_type, "quantity": quantity, "K_source": "Error", "error": "Invalid K_value provided"})


        K_suction_total = K_suction_pipe + K_suction_fittings
        K_discharge_total = K_discharge_pipe + K_discharge_fittings


        suction_head_loss = fluids.head_from_K(K=K_suction_total, V=pipe_velocity)
        discharge_head_loss = fluids.head_from_K(K=K_discharge_total, V=pipe_velocity)
        
        # Calculate velocity heads at pump nozzles
        suction_velocity_head = suction_nozzle_velocity ** 2 / (2 * G_GRAVITY)
        discharge_velocity_head = discharge_nozzle_velocity ** 2 / (2 * G_GRAVITY)
        velocity_head_differential = discharge_velocity_head - suction_velocity_head

        # TDH includes static head, friction losses, AND velocity head differential
        tdh = ((local_static_discharge_head - local_static_suction_head) + 
               discharge_head_loss + suction_head_loss + 
               velocity_head_differential)
        
        results_log.append(f"TDH components: Static={local_static_discharge_head - local_static_suction_head:.3f}m, "
                          f"Friction={discharge_head_loss + suction_head_loss:.3f}m, "
                          f"Velocity={velocity_head_differential:.3f}m")


        # NPSHA = Suction side pressure - Vapor pressure - Suction friction losses - Suction velocity head
        suction_pressure_head = local_suction_side_pressure / (local_fluid_density * G_GRAVITY)
        vapor_pressure_head = local_fluid_vapor_pressure / (local_fluid_density * G_GRAVITY)
        
        npsha = (suction_pressure_head +                     # Pressure head on liquid surface (atmospheric or tank)
                 local_static_suction_head -                 # Elevation relative to pump suction
                 vapor_pressure_head -                       # Vapor pressure head
                 suction_head_loss -                         # Friction losses in suction
                 suction_velocity_head)                      # Velocity head at suction nozzle
        
        results_log.append(f"NPSHA calculation: {suction_pressure_source}")
        results_log.append(f"NPSHA components: Suction_pressure={suction_pressure_head:.3f}m, "
                          f"Static={local_static_suction_head:.3f}m, "
                          f"Vapor={-vapor_pressure_head:.3f}m, "
                          f"Friction={-suction_head_loss:.3f}m, "
                          f"Velocity={-suction_velocity_head:.3f}m")


        # Approximate power - using a default efficiency
        pump_eff = 0.70 # Assumed efficiency - could be an input parameter
        brake_power_watts = (local_fluid_density * G_GRAVITY * local_flow_rate * tdh) / pump_eff if pump_eff > 0 else 0
        results_log.append(f"Power calculation based on assumed efficiency: {pump_eff*100:.1f}%")


        # --- Consolidate Results ---
        final_result = {
            "inputs_resolved": results_log,
            "warnings": error_log if error_log else None,
            "pipe_details": pipe_info if pipe_info else None,
            "fluid_details": fluid_info if fluid_info else None,
            # Core results with unit conversions
            "flow_rate_m3s": round(local_flow_rate, 6),
            "flow_rate_gpm": round(local_flow_rate / GPM_to_M3S, 2),
            "pipe_velocity_m_s": round(pipe_velocity, 3),
            "suction_nozzle_velocity_m_s": round(suction_nozzle_velocity, 3),
            "discharge_nozzle_velocity_m_s": round(discharge_nozzle_velocity, 3),
            "suction_nozzle_diameter_m": round(local_suction_nozzle_diameter, 4),
            "discharge_nozzle_diameter_m": round(local_discharge_nozzle_diameter, 4),
            "reynolds_number": round(Re, 2),
            "friction_factor": round(fd, 6),
            "suction_velocity_head_m": round(suction_velocity_head, 4),
            "discharge_velocity_head_m": round(discharge_velocity_head, 4),
            "velocity_head_differential_m": round(velocity_head_differential, 4),
            "suction_head_loss_m": round(suction_head_loss, 3),
            "suction_head_loss_ft": round(suction_head_loss / FT_to_M, 3),
            "suction_fitting_details": suction_fitting_details,
            "discharge_head_loss_m": round(discharge_head_loss, 3),
            "discharge_head_loss_ft": round(discharge_head_loss / FT_to_M, 3),
            "discharge_fitting_details": discharge_fitting_details,
            "suction_velocity_head_ft": round(suction_velocity_head / FT_to_M, 4),
            "discharge_velocity_head_ft": round(discharge_velocity_head / FT_to_M, 4),
            "velocity_head_differential_ft": round(velocity_head_differential / FT_to_M, 4),
            "total_dynamic_head_m": round(tdh, 3),
            "total_dynamic_head_ft": round(tdh / FT_to_M, 3),
            "npsha_m": round(npsha, 3),
            "npsha_ft": round(npsha / FT_to_M, 3),
            "estimated_power_kw": round(brake_power_watts / 1000, 2),
            "estimated_power_hp": round(brake_power_watts / 745.7, 2),
            "assumed_pump_efficiency": pump_eff
        }


        # Remove null keys if empty
        if not final_result["warnings"]: del final_result["warnings"]
        if not final_result["pipe_details"]: del final_result["pipe_details"]
        if not final_result["fluid_details"]: del final_result["fluid_details"]


        return json.dumps(final_result)
    except Exception as e:
        logger.error(f"Error in calculate_pump_requirements: {e}", exc_info=True)
        return json.dumps({"error": f"Calculation error: {str(e)}", "log": results_log, "errors_occurred": error_log})
