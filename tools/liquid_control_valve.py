import json
import logging
import math
from typing import Optional, Literal # Added Literal

# --- Third-Party Libraries ---
import fluids
import fluids.control_valve # Explicitly import submodule used

# Import shared utilities
from utils.constants import (
    GPM_to_M3S, INCH_to_M, FT_to_M, PSI_to_PA, 
    LBFT3_to_KGM3, CENTIPOISE_to_PAS, DEFAULT_ROUGHNESS, DEG_C_to_K
)
from utils.helpers import get_fitting_K
from utils.import_helpers import FLUIDPROP_AVAILABLE, FluidProperties, FLUID_SELECTION, COOLPROP_AVAILABLE, CP
from utils.fluid_aliases import map_fluid_name

# Configure logging
logger = logging.getLogger("fluids-mcp.calculate_liquid_control_valve")

def calculate_liquid_control_valve(
    # --- Core SI Inputs ---
    flow_rate: Optional[float] = None,             # Flow rate in m³/s
    inlet_pressure: Optional[float] = None,        # Inlet pressure in Pa
    outlet_pressure: Optional[float] = None,       # Outlet pressure in Pa
    fluid_density: Optional[float] = None,         # Fluid density in kg/m³
    fluid_viscosity: Optional[float] = None,       # Fluid viscosity in Pa·s
    # Optional SI critical props (unlikely user input, but allow override)
    fluid_saturation_pressure_pa: Optional[float] = None,
    fluid_critical_pressure_pa: Optional[float] = None,

    # --- Alternative Unit Inputs ---
    flow_rate_gpm: Optional[float] = None,           # Flow rate in US GPM
    inlet_pressure_psi: Optional[float] = None,    # Inlet pressure in psi
    outlet_pressure_psi: Optional[float] = None,   # Outlet pressure in psi
    fluid_density_lbft3: Optional[float] = None,   # Fluid density in lb/ft³
    fluid_viscosity_cp: Optional[float] = None,    # Fluid viscosity in centipoise
    fluid_saturation_pressure_psi: Optional[float] = None,
    fluid_critical_pressure_psi: Optional[float] = None,

    # --- Property Lookup Inputs ---
    fluid_name: Optional[str] = None,               # e.g., "Water"
    temperature_c: Optional[float] = None,          # Temperature in Celsius
    pressure_bar: Optional[float] = None,           # Pressure in bar (context for lookup)

    # --- Valve & Sizing Parameters ---
    valve_type: Literal["globe", "butterfly", "ball", "gate", "other"] = "globe",  # Valve type for estimating factors
    # Provide known factors if 'other' or to override defaults
    valve_fl: Optional[float] = None,                 # Liquid pressure recovery factor
    valve_fd: Optional[float] = None,                 # Valve style modifier (often near 1 for line size valves)
    # Output formatting
    size_units: Literal["m", "inch", "mm"] = "inch" # Default to inches for valve size output
) -> str:
    """Calculate liquid control valve coefficient (Cv/Kv) and recommended size with flexible input options.

    Uses an iterative approach to determine a suitable standard valve size.
    Accepts inputs in SI units, common imperial units, or via property lookups.
    Priority: Explicit SI > Explicit Imperial > Lookup > Default.

    Args:
        # --- SI Units ---
        flow_rate: Optional - Flow rate in m³/s
        inlet_pressure: Optional - Inlet pressure in Pa
        outlet_pressure: Optional - Outlet pressure in Pa
        fluid_density: Optional - Fluid density in kg/m³
        fluid_viscosity: Optional - Fluid viscosity in Pa·s
        fluid_saturation_pressure_pa: Optional - Fluid saturation pressure in Pa
        fluid_critical_pressure_pa: Optional - Fluid critical pressure in Pa
        # --- Imperial Units ---
        flow_rate_gpm: Optional - Flow rate in US GPM
        inlet_pressure_psi: Optional - Inlet pressure in psi
        outlet_pressure_psi: Optional - Outlet pressure in psi
        fluid_density_lbft3: Optional - Fluid density in lb/ft³
        fluid_viscosity_cp: Optional - Fluid viscosity in centipoise
        fluid_saturation_pressure_psi: Optional - Fluid saturation pressure in psi
        fluid_critical_pressure_psi: Optional - Fluid critical pressure in psi
        # --- Property Lookups ---
        fluid_name: Optional - Name of the fluid (e.g., "Water") for property lookup
        temperature_c: Optional - Temperature in Celsius (required with fluid_name for lookup)
        pressure_bar: Optional - Pressure in bar (context for fluid property lookup, defaults ~1 atm)
        # --- Valve & Sizing ---
        valve_type: Type of valve ("globe", "butterfly", "ball", "gate", "other") for estimating factors (default: globe)
        valve_fl: Optional - Override calculated/default Liquid Pressure Recovery Factor (FL)
        valve_fd: Optional - Override calculated/default Valve Style Modifier (Fd)
        size_units: Units for recommended valve size output ("m", "inch", "mm", default: "inch")

    Returns:
        JSON string with valve sizing results (Cv, Kv), recommended size, and input resolution log.
    """
    results_log = []
    error_log = []
    fluid_info = {} # Store details from lookup

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

        # 2. Resolve Pressures (SI)
        local_inlet_pressure = None
        if inlet_pressure is not None:
            local_inlet_pressure = inlet_pressure
            results_log.append("Used provided SI inlet_pressure.")
        elif inlet_pressure_psi is not None:
            local_inlet_pressure = inlet_pressure_psi * PSI_to_PA
            results_log.append(f"Converted inlet_pressure from {inlet_pressure_psi} psi.")
        else:
            error_log.append("Missing required input: inlet_pressure or inlet_pressure_psi.")

        local_outlet_pressure = None
        if outlet_pressure is not None:
            local_outlet_pressure = outlet_pressure
            results_log.append("Used provided SI outlet_pressure.")
        elif outlet_pressure_psi is not None:
            local_outlet_pressure = outlet_pressure_psi * PSI_to_PA
            results_log.append(f"Converted outlet_pressure from {outlet_pressure_psi} psi.")
        else:
            error_log.append("Missing required input: outlet_pressure or outlet_pressure_psi.")

        # 3. Resolve Fluid Properties (Density, Viscosity, Psat, Pc)
        local_fluid_density = None
        local_fluid_viscosity = None
        local_Psat = None
        local_Pc = None
        prop_lookup_success = False

        # 3.1 Density and Viscosity
        if fluid_density is not None and fluid_viscosity is not None:
            local_fluid_density = fluid_density
            local_fluid_viscosity = fluid_viscosity
            results_log.append("Used provided SI fluid density and viscosity.")
            prop_lookup_success = True # Assume if density/visc provided, user might provide Psat/Pc too
        elif fluid_density_lbft3 is not None and fluid_viscosity_cp is not None:
              local_fluid_density = fluid_density_lbft3 * LBFT3_to_KGM3
              local_fluid_viscosity = fluid_viscosity_cp * CENTIPOISE_to_PAS
              results_log.append(f"Converted density from {fluid_density_lbft3} lb/ft³ and viscosity from {fluid_viscosity_cp} cP.")
              prop_lookup_success = True
        elif fluid_name is not None and temperature_c is not None:
            if FLUIDPROP_AVAILABLE and FLUID_SELECTION is not None and FluidProperties is not None:
                try: # Fluid property lookup
                    # First try mapping common aliases
                    mapped_fluid_name = map_fluid_name(fluid_name)

                    try:
                        valid_fluids = [f[0] for f in FLUID_SELECTION if f is not None and hasattr(f, '__getitem__')]
                    except (TypeError, IndexError):
                        valid_fluids = []
                    actual_fluid_name = mapped_fluid_name
                    
                    # Handle incompressible fluids (glycols)
                    if mapped_fluid_name.startswith('INCOMP::'):
                        # For incompressible fluids, use directly with CoolProp
                        actual_fluid_name = mapped_fluid_name
                        # Note: FluidProperties may not handle INCOMP:: fluids well
                        # This is a known limitation
                    elif not valid_fluids or mapped_fluid_name not in valid_fluids:
                        fluid_lower = mapped_fluid_name.lower()
                        match = next((f for f in valid_fluids if f.lower() == fluid_lower), None)
                        if match: actual_fluid_name = match
                        else: raise ValueError(f"Fluid '{fluid_name}' (mapped to '{mapped_fluid_name}') not found.")


                    p_bar = pressure_bar if pressure_bar is not None else 1.01325 # Default pressure
                    fluid_props = FluidProperties(coolprop_name=actual_fluid_name, T_in_deg_C=temperature_c, P_in_bar=p_bar)
                    local_fluid_density = float(fluid_props.rho[0])
                    local_fluid_viscosity = float(fluid_props.eta[0])
                    fluid_info = { # Store details
                        "name_used": actual_fluid_name, "temperature_c": temperature_c, "pressure_bar": p_bar,
                        "density_kgm3": local_fluid_density, "viscosity_pas": local_fluid_viscosity
                    }
                    results_log.append(f"Looked up density/viscosity for {actual_fluid_name} at {temperature_c}°C, {p_bar} bar.")
                    prop_lookup_success = True


                    # Attempt Psat/Pc lookup ONLY if base props succeeded & CoolProp available
                    if COOLPROP_AVAILABLE:
                        try:
                            temp_k = temperature_c + DEG_C_to_K
                            # Use saturation pressure at the given temperature
                            local_Psat = CP.PropsSI('P', 'T', temp_k, 'Q', 0, actual_fluid_name)
                            fluid_info["saturation_pressure_pa"] = local_Psat
                            results_log.append(f"Looked up saturation pressure ({local_Psat:.2f} Pa).")
                        except Exception as psat_e:
                            error_log.append(f"Warning: Could not look up Psat for {actual_fluid_name}: {psat_e}")
                        try:
                            # Critical pressure might depend less on T/P context, use generic lookup
                            local_Pc = CP.PropsSI('PCRIT', actual_fluid_name)
                            fluid_info["critical_pressure_pa"] = local_Pc
                            results_log.append(f"Looked up critical pressure ({local_Pc:.2f} Pa).")
                        except Exception as pc_e:
                                error_log.append(f"Warning: Could not look up Pc for {actual_fluid_name}: {pc_e}")

                except Exception as fluid_lookup_e:
                    local_fluid_density = None
                    local_fluid_viscosity = None
                    error_log.append(f"Failed base fluid property lookup for {fluid_name} at {temperature_c}°C: {fluid_lookup_e}")
                    prop_lookup_success = False
            else:
                error_log.append("Fluid property lookup skipped: fluidprop package not available.")
        else:
            # This error is raised only if density/viscosity cannot be determined
            error_log.append("Missing inputs for fluid density/viscosity: Need (SI pair) OR (Imperial pair) OR (fluid_name AND temperature_c).")

        # 3.2 Saturation Pressure (if not found via lookup)
        if local_Psat is None:
            if fluid_saturation_pressure_pa is not None:
                local_Psat = fluid_saturation_pressure_pa
                results_log.append("Used provided SI fluid_saturation_pressure_pa.")
            elif fluid_saturation_pressure_psi is not None:
                local_Psat = fluid_saturation_pressure_psi * PSI_to_PA
                results_log.append(f"Converted fluid_saturation_pressure from {fluid_saturation_pressure_psi} psi.")
            else:
                # Use water defaults only for water; other fluids get warning about missing Psat
                if prop_lookup_success:
                    is_water = fluid_name is not None and fluid_name.lower() in ['water', 'h2o']
                    if is_water:
                        local_Psat = 3170.0  # Default for water at 25C
                        results_log.append(f"Using default saturation pressure: {local_Psat} Pa (Water @ 25C).")
                    else:
                        # For non-water fluids, use a conservative estimate but warn user
                        local_Psat = 1000.0  # Low conservative estimate
                        results_log.append(f"Using conservative Psat estimate: {local_Psat} Pa for non-water fluid.")
                        error_log.append(f"Warning: Saturation pressure unknown for '{fluid_name}'. Cavitation check may be inaccurate. Provide fluid_saturation_pressure for accurate results.")
                # else: error already logged for density/viscosity

        # 3.3 Critical Pressure (if not found via lookup)
        if local_Pc is None:
            if fluid_critical_pressure_pa is not None:
                local_Pc = fluid_critical_pressure_pa
                results_log.append("Used provided SI fluid_critical_pressure_pa.")
            elif fluid_critical_pressure_psi is not None:
                local_Pc = fluid_critical_pressure_psi * PSI_to_PA
                results_log.append(f"Converted fluid_critical_pressure from {fluid_critical_pressure_psi} psi.")
            else:
                # Use water defaults only for water; other fluids get warning about missing Pc
                if prop_lookup_success:
                    is_water = fluid_name is not None and fluid_name.lower() in ['water', 'h2o']
                    if is_water:
                        local_Pc = 22064000.0  # Default for water
                        results_log.append(f"Using default critical pressure: {local_Pc} Pa (Water).")
                    else:
                        # For non-water fluids, use a conservative estimate but warn user
                        local_Pc = 5000000.0  # Conservative estimate (~50 bar)
                        results_log.append(f"Using conservative Pc estimate: {local_Pc} Pa for non-water fluid.")
                        error_log.append(f"Warning: Critical pressure unknown for '{fluid_name}'. Cavitation check may be inaccurate. Provide fluid_critical_pressure for accurate results.")
                # else: error already logged for density/viscosity

        # 4. Resolve Valve Factors (FL, Fd)
        local_FL = None
        local_Fd = None
        FL_values = {"globe": 0.9, "butterfly": 0.7, "ball": 0.6, "gate": 0.9, "other": 0.9} # Default for 'other'
        Fd_values = {"globe": 0.9, "butterfly": 0.7, "ball": 0.98, "gate": 0.9, "other": 1.0} # Default for 'other' is 1

        if valve_fl is not None:
            local_FL = valve_fl
            results_log.append("Used provided valve_fl.")
        else:
            local_FL = FL_values.get(valve_type, 0.9) # Default to globe/other if type invalid
            results_log.append(f"Used default FL={local_FL} for valve type '{valve_type}'.")

        if valve_fd is not None:
            local_Fd = valve_fd
            results_log.append("Used provided valve_fd.")
        else:
            local_Fd = Fd_values.get(valve_type, 1.0) # Default to 1.0 for 'other'
            results_log.append(f"Used default Fd={local_Fd} for valve type '{valve_type}'.")

        # --- Check for any errors before proceeding ---
        missing_critical = any(e.startswith("Missing required input") for e in error_log)
        if missing_critical:
            return json.dumps({"errors": error_log, "log": results_log})

        # Check resolved critical local variables
        critical_vars_valve = [
            local_flow_rate, local_inlet_pressure, local_outlet_pressure,
            local_fluid_density, local_fluid_viscosity, local_Psat, local_Pc,
            local_FL, local_Fd
        ]
        if any(v is None for v in critical_vars_valve):
            error_log.append("Critical input parameter resolution failed for valve calculation.")
            missing_items = [
                name for name, val in zip([
                    "flow_rate", "inlet_P", "outlet_P", "density", "viscosity",
                    "Psat", "Pc", "FL", "Fd"
                ], critical_vars_valve) if val is None
            ]
            error_log.append(f"Missing resolved values for: {', '.join(missing_items)}")
            return json.dumps({"errors": error_log, "log": results_log})

        # --- Check Pressure Drop ---
        delta_P = local_inlet_pressure - local_outlet_pressure
        if delta_P <= 0:
            error_log.append("Inlet pressure must be greater than outlet pressure.")
            # Return error immediately as calculation cannot proceed
            return json.dumps({"errors": error_log, "log": results_log})

        # --- Iterative Sizing Logic ---
        standard_sizes_inch = [
            0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 24 # Add more if needed
        ]
        current_size_inch = 2.0 # Initial guess
        previous_size_inch = 0.0
        iteration = 0
        max_iterations = 10
        calculated_size_inch = 0.0 # Initialize
        Kv = 0.0 # Initialize

        while abs(current_size_inch - previous_size_inch) > 1e-6 and iteration < max_iterations:
            previous_size_inch = current_size_inch
            current_size_m = current_size_inch * INCH_to_M

            # Assume line size valve for sizing calculation (d=D1=D2)
            d = D1 = D2 = current_size_m

            try:
                Kv = fluids.control_valve.size_control_valve_l(
                    rho=local_fluid_density,
                    Psat=local_Psat,
                    Pc=local_Pc,
                    mu=local_fluid_viscosity,
                    P1=local_inlet_pressure,
                    P2=local_outlet_pressure,
                    Q=local_flow_rate,
                    D1=D1, # Pipe dia = valve dia
                    D2=D2, # Pipe dia = valve dia
                    d=d,   # Valve bore = valve dia
                    FL=local_FL,
                    Fd=local_Fd
                )

                # Estimate required size based on Kv
                # Using the formula: Kv ≈ 100 * (DN/25)^2 => DN_mm = 25 * sqrt(Kv/100)
                # DN_inches = DN_mm / 25.4
                if Kv > 0:
                    calculated_size_mm = 25.0 * math.sqrt(Kv / 100.0)
                    calculated_size_inch = calculated_size_mm / 25.4
                else:
                    calculated_size_inch = 0 # Handle negative Kv case

                # Find nearest *larger or equal* standard size
                selected_size_inch = next((size for size in standard_sizes_inch if size >= calculated_size_inch), standard_sizes_inch[-1])

                results_log.append(f"Size iteration {iteration+1}: Assumed NPS={current_size_inch}\", Kv={Kv:.3f}, Required NPS~{calculated_size_inch:.3f}\", Next NPS={selected_size_inch}\"")
                current_size_inch = selected_size_inch # Update for next iteration

            except Exception as valve_calc_e:
                error_log.append(f"Error during valve Kv calculation (iteration {iteration+1}): {valve_calc_e}")
                # Stop iteration if calculation fails
                iteration = max_iterations # Force exit
                break # Exit loop

            iteration += 1

        converged = abs(current_size_inch - previous_size_inch) <= 1e-6 and iteration <= max_iterations

        # --- Final Calculation with Converged Size ---
        final_size_inch = current_size_inch
        final_size_m = final_size_inch * INCH_to_M
        d = D1 = D2 = final_size_m

        try:
            # Recalculate Kv/Cv with the final determined standard size
            final_Kv = fluids.control_valve.size_control_valve_l(
                rho=local_fluid_density, Psat=local_Psat, Pc=local_Pc, mu=local_fluid_viscosity,
                P1=local_inlet_pressure, P2=local_outlet_pressure, Q=local_flow_rate,
                D1=D1, D2=D2, d=d, FL=local_FL, Fd=local_Fd
            )
            final_Cv = final_Kv / 0.865
        except Exception as final_calc_e:
            error_log.append(f"Error during final valve coefficient calculation: {final_calc_e}")
            final_Kv = None
            final_Cv = None

        # --- Format Output Size ---
        if size_units == "inch":
            size_value = final_size_inch
            size_description = f"{final_size_inch} inch"
        elif size_units == "mm":
            size_value = final_size_inch * 25.4
            size_description = f"{round(size_value)} mm (NPS {final_size_inch})"
        else: # meters
            size_value = final_size_m
            size_description = f"{round(size_value, 4)} m (NPS {final_size_inch})"

        # --- Consolidate Results ---
        final_result = {
            "inputs_resolved": results_log,
            "warnings": error_log if error_log else None,
            "fluid_details": fluid_info if fluid_info else None,
            # Core Results
            "flow_rate_m3s": round(local_flow_rate, 6),
            "flow_rate_gpm": round(local_flow_rate / GPM_to_M3S, 2),
            "inlet_pressure_pa": round(local_inlet_pressure, 2),
            "outlet_pressure_pa": round(local_outlet_pressure, 2),
            "pressure_drop_pa": round(delta_P, 2),
            "pressure_drop_psi": round(delta_P / PSI_to_PA, 3),
            "recommended_valve_size": size_description,
            "recommended_size_value": round(size_value, 4),
            "size_units": size_units,
            "valve_type": valve_type,
            "cv_us": round(final_Cv, 2) if final_Cv is not None else "Error",
            "kv_metric": round(final_Kv, 2) if final_Kv is not None else "Error",
            "FL_used": local_FL,
            "Fd_used": local_Fd,
            "sizing_iterations": iteration,
            "sizing_converged": converged
        }

        # Remove null keys if empty
        if not final_result["warnings"]: del final_result["warnings"]
        if not final_result["fluid_details"]: del final_result["fluid_details"]

        return json.dumps(final_result)

    except Exception as e:
        logger.error(f"Error in calculate_liquid_control_valve: {e}", exc_info=True)
        return json.dumps({"error": f"Calculation error: {str(e)}", "log": results_log, "errors_occurred": error_log})

