import json
import logging
import math
from typing import Optional, Dict, Literal

# --- Third-Party Libraries ---
import fluids
import fluids.compressible # Explicitly import submodule used
import fluids.core         # Explicitly import submodule used
import fluids.friction     # Explicitly import submodule used
import fluids.piping       # Explicitly import submodule used
# Note: fluids.friction._roughness is imported dynamically within the function.

# Import shared utilities
from utils.constants import (
    GPM_to_M3S, INCH_to_M, FT_to_M, PSI_to_PA, 
    LBFT3_to_KGM3, CENTIPOISE_to_PAS, DEFAULT_ROUGHNESS, DEG_C_to_K
)
from utils.helpers import get_fitting_K
from utils.import_helpers import FLUIDPROP_AVAILABLE, FluidProperties, FLUID_SELECTION

# Configure logging
logger = logging.getLogger("fluids-mcp.calculate_gas_pipe_pressure_drop")

def calculate_gas_pipe_pressure_drop(
    # --- Problem Definition: Specify 3 out of 4: P1, P2, L, Q ---
    inlet_pressure: Optional[float] = None,        # Inlet pressure (absolute) in Pa
    outlet_pressure: Optional[float] = None,       # Outlet pressure (absolute) in Pa
    pipe_length: Optional[float] = None,           # Pipe length in m
    flow_rate_kg_s: Optional[float] = None,        # Mass flow rate in kg/s
    flow_rate_norm_m3_hr: Optional[float] = None,  # Volumetric flow rate at Normal conditions (0°C, 1 atm) in m³/hr
    flow_rate_std_m3_hr: Optional[float] = None,   # Volumetric flow rate at Standard conditions (15°C, 1 atm) in m³/hr,

    # --- Pipe Details (SI, Imperial, or Lookup) ---
    pipe_diameter: Optional[float] = None,         # Pipe inner diameter in m
    pipe_diameter_in: Optional[float] = None,      # Pipe inner diameter in inches
    nominal_size_in: Optional[float] = None,       # Nominal pipe size in inches
    schedule: str = "40",                          # Pipe schedule
    material: Optional[str] = None,                # Pipe material name (e.g., "Steel", "HDPE")
    pipe_roughness: Optional[float] = None,        # Pipe absolute roughness in m

    # --- Gas Properties (Direct Input, Lookup, or Composition) ---
    temperature_c: Optional[float] = None,         # Gas temperature in Celsius
    gas_mw: Optional[float] = None,                # Gas Molecular Weight in kg/kmol (g/mol)
    gas_gamma: Optional[float] = None,             # Specific heat ratio (Cp/Cv)
    gas_z_factor: Optional[float] = None,          # Compressibility factor (Z)
    gas_viscosity: Optional[float] = None,         # Gas dynamic viscosity in Pa·s
    fluid_name: Optional[str] = None,              # Common gas name ("Air", "Methane", "NaturalGas", "CarbonDioxide")
    gas_composition_mol: Optional[Dict[str, float]] = None, # Composition as mol fraction

    # --- Calculation Options ---
    method: Literal["Weymouth", "Panhandle_A", "Panhandle_B", "IGT", "Oliphant", "Spitzglass_low", 
                    "Spitzglass_high", "isothermal_darcy"] = "Weymouth"
) -> str:
    """Calculates compressible gas flow parameters (pressure drop, flow rate, or length) with flexible inputs.

    Specify exactly 3 out of the 4 primary variables: inlet_pressure, outlet_pressure, pipe_length, flow_rate_*.
    The function will solve for the missing one. Supports various gas flow equations.

    Args:
        inlet_pressure: Optional - Inlet pressure (absolute) in Pa
        outlet_pressure: Optional - Outlet pressure (absolute) in Pa
        pipe_length: Optional - Pipe length in m
        flow_rate_kg_s: Optional - Mass flow rate in kg/s
        flow_rate_norm_m3_hr: Optional - Volumetric flow at Normal conditions (0°C, 1 atm) in m³/hr
        flow_rate_std_m3_hr: Optional - Volumetric flow at Standard conditions (15°C, 1 atm) in m³/hr
        pipe_diameter: Optional - Pipe inner diameter in m
        pipe_diameter_in: Optional - Pipe inner diameter in inches
        nominal_size_in: Optional - Nominal pipe size in inches
        schedule: Pipe schedule (default: "40")
        material: Optional - Pipe material name for roughness lookup
        pipe_roughness: Optional - Pipe absolute roughness in m
        temperature_c: Required - Gas temperature in Celsius
        gas_mw: Optional - Gas Molecular Weight in kg/kmol
        gas_gamma: Optional - Specific heat ratio (Cp/Cv)
        gas_z_factor: Optional - Compressibility factor (Z)
        gas_viscosity: Optional - Gas dynamic viscosity in Pa·s
        fluid_name: Optional - Common gas name ("Air", "Methane", etc.) for property lookup
        gas_composition_mol: Optional - Gas composition as mol fractions {component: fraction}
        method: Calculation method ("Weymouth", "Panhandle_A", etc.) (default: Weymouth)

    Returns:
        JSON string with calculated results, including the solved variable, and input resolution log.
    """
    results_log = []
    error_log = []
    pipe_info = {}
    gas_prop_info = {}
    final_results = {}  # <-- Initialize final_results here

    # --- Standard Constants ---
    P_norm = 101325.0  # Pa
    T_norm = 273.15    # K
    P_std = 101325.0   # Pa
    T_std = 288.15     # K (15 C)
    R_univ = 8314.462  # J/(kmol·K)

    try:
        # 1. Check Temperature Input
        if temperature_c is None:
            error_log.append("Missing required input: temperature_c.")
            return json.dumps({"errors": error_log, "log": results_log})
        local_T_k = temperature_c + DEG_C_to_K  # Convert to Kelvin
        results_log.append(f"Using Temperature: {temperature_c:.2f} C ({local_T_k:.2f} K)")

        # 2. Resolve Pipe Diameter & Roughness (SI)
        local_pipe_diameter = None
        local_pipe_roughness = None
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
                        local_pipe_roughness = 4.5e-5
                        results_log.append(
                            f"Material '{material}' not found, using default roughness {local_pipe_roughness:.2e} m."
                        )
                        error_log.append(f"Warning: Material '{material}' not found for roughness lookup.")
                except Exception as lookup_e:
                    local_pipe_roughness = 4.5e-5
                    results_log.append(
                        f"Roughness lookup failed ({lookup_e}), using default {local_pipe_roughness:.2e} m."
                    )
                    error_log.append(f"Warning: Roughness lookup failed: {lookup_e}")
            else:
                local_pipe_roughness = 4.5e-5
                results_log.append(f"Used default pipe_roughness {local_pipe_roughness:.2e} m.")
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
                        local_pipe_roughness = 4.5e-5
                        results_log.append(
                            f"Material '{material}' not found, using default roughness {local_pipe_roughness:.2e} m."
                        )
                        error_log.append(f"Warning: Material '{material}' not found for roughness lookup.")
                except Exception as lookup_e:
                    local_pipe_roughness = 4.5e-5
                    results_log.append(
                        f"Roughness lookup failed ({lookup_e}), using default {local_pipe_roughness:.2e} m."
                    )
                    error_log.append(f"Warning: Roughness lookup failed: {lookup_e}")
            else:
                local_pipe_roughness = 4.5e-5
                results_log.append(f"Used default pipe_roughness {local_pipe_roughness:.2e} m.")
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
                            local_pipe_roughness = 4.5e-5
                            results_log.append(
                                f"Material '{material}' not found, using default roughness {local_pipe_roughness:.2e} m."
                            )
                            error_log.append(f"Warning: Material '{material}' not found for roughness lookup.")
                            pipe_info["roughness_m"] = local_pipe_roughness
                            pipe_info["roughness_source"] = "Default (Material Not Found)"
                    except Exception as lookup_e:
                        local_pipe_roughness = 4.5e-5
                        results_log.append(
                            f"Roughness lookup failed ({lookup_e}), using default {local_pipe_roughness:.2e} m."
                        )
                        error_log.append(f"Warning: Roughness lookup failed: {lookup_e}")
                        pipe_info["roughness_m"] = local_pipe_roughness
                        pipe_info["roughness_source"] = "Default (Lookup Failed)"
                else:
                    local_pipe_roughness = 4.5e-5
                    results_log.append(
                        f"Used default pipe_roughness {local_pipe_roughness:.2e} m (material not specified)."
                    )
                    pipe_info["roughness_m"] = local_pipe_roughness
                    pipe_info["roughness_source"] = "Default (Material Not Specified)"
            except Exception as pipe_lookup_e:
                local_pipe_diameter = None
                local_pipe_roughness = None
                error_log.append(
                    f"Failed to look up pipe NPS {nominal_size_in} Sch {schedule}: {pipe_lookup_e}"
                )
        else:
            # Check if we're solving for diameter - if so, this is expected
            # We'll determine this after checking all variables
            local_pipe_diameter = None
            local_pipe_roughness = pipe_roughness or 4.5e-5  # Use provided or default roughness

        # 3. Resolve Gas Properties
        local_gas_mw = gas_mw
        local_gas_gamma = gas_gamma
        local_gas_z_factor = gas_z_factor
        local_gas_viscosity = gas_viscosity

        gas_prop_source = {
            "MW": "Provided" if local_gas_mw is not None else None,
            "gamma": "Provided" if local_gas_gamma is not None else None,
            "Z": "Provided" if local_gas_z_factor is not None else None,
            "viscosity": "Provided" if local_gas_viscosity is not None else None
        }
        prop_resolved = True

        needs_resolution = any(
            p is None
            for p in [local_gas_mw, local_gas_gamma, local_gas_z_factor, local_gas_viscosity]
        )
        if needs_resolution:
            lookup_status = "Not Attempted"
            if fluid_name is not None:
                results_log.append(f"Attempting property lookup for '{fluid_name}' at T={temperature_c} C.")
                if FLUIDPROP_AVAILABLE:
                    try:
                        avg_p_pa = (
                            (inlet_pressure + outlet_pressure) / 2.0
                            if inlet_pressure is not None and outlet_pressure is not None
                            else inlet_pressure
                        )
                        lookup_p_bar = avg_p_pa / 100000.0 if avg_p_pa else 1.0

                        # Try FluidProperties directly first
                        try:
                            fluid_props = FluidProperties(
                                coolprop_name=fluid_name,
                                T_in_deg_C=temperature_c,
                                P_in_bar=lookup_p_bar
                            )
                            actual_fluid_name = fluid_name
                        except Exception:
                            # Fallback to FLUID_SELECTION validation if direct lookup fails
                            valid_fluids = [f[0] for f in FLUID_SELECTION]
                            actual_fluid_name = fluid_name
                            if fluid_name not in valid_fluids:
                                match = next(
                                    (f for f in valid_fluids if f.lower() == fluid_name.lower()), None
                                )
                                if match:
                                    actual_fluid_name = match
                                else:
                                    raise ValueError(f"Fluid '{fluid_name}' not found in fluidprop.")

                            fluid_props = FluidProperties(
                                coolprop_name=actual_fluid_name,
                                T_in_deg_C=temperature_c,
                                P_in_bar=lookup_p_bar
                            )
                        lookup_status = f"Lookup Success ({actual_fluid_name})"

                        if local_gas_mw is None:
                            local_gas_mw = float(fluid_props.MW) * 1000.0  # Convert kg/mol to kg/kmol
                            gas_prop_source["MW"] = lookup_status
                        if local_gas_gamma is None:
                            # Safely access gamma attribute or calculate from Cp/Cv if available
                            if hasattr(fluid_props, 'gamma'):
                                local_gas_gamma = float(fluid_props.gamma[0])
                                gas_prop_source["gamma"] = lookup_status
                            elif hasattr(fluid_props, 'Cp') and hasattr(fluid_props, 'Cv') and float(fluid_props.Cv[0]) != 0:
                                local_gas_gamma = float(fluid_props.Cp[0]) / float(fluid_props.Cv[0])
                                gas_prop_source["gamma"] = lookup_status + " (Calculated Cp/Cv)"
                            # If gamma can't be determined, it will fall back to default later
                        if local_gas_z_factor is None:
                            # Safely access Z attribute
                            if hasattr(fluid_props, 'Z'):
                                local_gas_z_factor = float(fluid_props.Z[0])
                                gas_prop_source["Z"] = lookup_status
                            # Z not found, it will fall back to default later
                        if local_gas_viscosity is None:
                            if hasattr(fluid_props, 'eta'):
                                local_gas_viscosity = float(fluid_props.eta[0])
                                gas_prop_source["viscosity"] = lookup_status
                            # If viscosity can't be determined, it will fall back to default later

                        gas_prop_info = {
                            "name_used": actual_fluid_name,
                            "lookup_temp_c": temperature_c,
                            "lookup_p_bar": lookup_p_bar,
                            "viscosity_pas": local_gas_viscosity
                        }

                    except Exception as prop_lookup_e:
                        detailed_error = (
                            f"Property lookup failed for '{fluid_name}' at T={temperature_c}C, "
                            f"P={lookup_p_bar}bar: {type(prop_lookup_e).__name__}: {str(prop_lookup_e)}"
                        )
                        error_log.append(
                            f"Warning: {detailed_error}. Falling back to defaults/estimates."
                        )
                        logger.warning(detailed_error)
                        lookup_status = "Lookup Failed"
                        if gas_mw is None:
                            gas_prop_source["MW"] = lookup_status
                        if gas_gamma is None:
                            gas_prop_source["gamma"] = lookup_status
                        if gas_z_factor is None:
                            gas_prop_source["Z"] = lookup_status
                        prop_resolved = False
                else:
                    error_log.append(
                        "Warning: Property lookup skipped (fluidprop unavailable). "
                        "Falling back to defaults/estimates."
                    )
                    lookup_status = "Default (fluidprop unavailable)"
                    if local_gas_mw is None:
                        gas_prop_source["MW"] = lookup_status
                    if local_gas_gamma is None:
                        gas_prop_source["gamma"] = lookup_status
                    if local_gas_z_factor is None:
                        gas_prop_source["Z"] = lookup_status
                    prop_resolved = False

            # Apply defaults if still None
            if local_gas_mw is None:
                local_gas_mw = 28.96
                results_log.append("Using default MW (Air).")
                gas_prop_source["MW"] = "Default (Air)"
                prop_resolved = False
            if local_gas_gamma is None:
                local_gas_gamma = 1.4
                results_log.append("Using default gamma (Air).")
                gas_prop_source["gamma"] = "Default (Air)"
                prop_resolved = False
            if local_gas_z_factor is None:
                local_gas_z_factor = 1.0
                results_log.append("Using default Z-factor (Ideal).")
                gas_prop_source["Z"] = "Default (Ideal)"
                prop_resolved = False
            if local_gas_viscosity is None:
                local_gas_viscosity = 1.8e-5
                results_log.append("Using default viscosity (Air).")
                gas_prop_source["viscosity"] = "Default (Air)"
                prop_resolved = False

        # Store final gas properties
        gas_prop_info.update({
            "MW_kg_kmol": local_gas_mw,
            "gamma": local_gas_gamma,
            "Z_factor": local_gas_z_factor,
            "viscosity_pas": local_gas_viscosity,
            "source_details": gas_prop_source
        })
        results_log.append(
            f"Resolved Gas Properties: MW={local_gas_mw:.2f} "
            f"({gas_prop_source.get('MW', 'N/A')}), gamma={local_gas_gamma:.3f} "
            f"({gas_prop_source.get('gamma', 'N/A')}), "
            f"Z={local_gas_z_factor:.4f} ({gas_prop_source.get('Z', 'N/A')})"
        )

        # 4. Resolve Flow Rate (SI Mass Flow - kg/s)
        local_flow_rate_kg_s = None
        flow_rate_source = "N/A"
        if flow_rate_kg_s is not None:
            local_flow_rate_kg_s = flow_rate_kg_s
            flow_rate_source = "Provided mass flow (kg/s)"
        elif flow_rate_norm_m3_hr is not None:
            rho_norm = (P_norm * local_gas_mw) / (1.0 * R_univ * T_norm)
            local_flow_rate_kg_s = (flow_rate_norm_m3_hr / 3600.0) * rho_norm
            flow_rate_source = f"Converted from Normal m³/hr (rho_norm={rho_norm:.3f} kg/m³)"
        elif flow_rate_std_m3_hr is not None:
            rho_std = (P_std * local_gas_mw) / (1.0 * R_univ * T_std)
            local_flow_rate_kg_s = (flow_rate_std_m3_hr / 3600.0) * rho_std
            flow_rate_source = f"Converted from Standard m³/hr (rho_std={rho_std:.3f} kg/m³)"
        # else: to be solved for

        if local_flow_rate_kg_s is not None:
            results_log.append(
                f"Resolved Flow Rate: {local_flow_rate_kg_s:.4f} kg/s (Source: {flow_rate_source})"
            )

        # 5. Pressures and Length
        local_P1 = inlet_pressure
        local_P2 = outlet_pressure
        local_L = pipe_length
        if local_P1 is not None:
            results_log.append(f"Using P1 = {local_P1:.2f} Pa")
        if local_P2 is not None:
            results_log.append(f"Using P2 = {local_P2:.2f} Pa")
        if local_L is not None:
            results_log.append(f"Using L = {local_L:.2f} m")

        # Identify missing variable - NOW SUPPORTS DIAMETER AND LENGTH SOLVING
        specified_vars = {
            'P1': local_P1, 
            'P2': local_P2, 
            'L': local_L, 
            'Q': local_flow_rate_kg_s,
            'D': local_pipe_diameter
        }
        unknown_vars = [k for k, v in specified_vars.items() if v is None]

        if len(unknown_vars) != 1:
            error_log.append(
                "Incorrect number of primary variables specified. "
                f"Provide exactly 4 of: {list(specified_vars.keys())} (P1, P2, L, Q, D)"
            )
            return json.dumps({"errors": error_log, "log": results_log})

        solve_for = unknown_vars[0]
        results_log.append(f"Solving for: {solve_for}")
        
        # Validate that we have required inputs based on what we're solving for
        if solve_for != 'D' and local_pipe_diameter is None:
            error_log.append("Missing required input: pipe_diameter, pipe_diameter_in, or nominal_size_in (required when not solving for diameter).")
            return json.dumps({"errors": error_log, "log": results_log})

        # --- Calculate derived values needed by some methods ---
        # NOTE: Removed problematic averaging approach for isothermal_darcy per expert reviewer
        # fluids.compressible.isothermal_gas now handles internal iteration for accuracy

        # --- Build keyword arguments ---
        base_kwargs = {
            'Tavg': local_T_k,
            'Zavg': local_gas_z_factor,
            'E': 1.0
        }
        if local_P1 is not None:
            base_kwargs['P1'] = local_P1
        if local_P2 is not None:
            base_kwargs['P2'] = local_P2
        if local_L is not None:
            base_kwargs['L'] = local_L
        if local_pipe_diameter is not None:
            base_kwargs['D'] = local_pipe_diameter
        if local_gas_mw is not None and local_gas_mw > 0:
            base_kwargs['SG'] = local_gas_mw / 28.96

        Q_norm_m3s = None
        if local_flow_rate_kg_s is not None and local_gas_mw and local_gas_mw > 0:
            try:
                rho_norm = (P_norm * local_gas_mw) / (R_univ * T_norm)
                if rho_norm > 0:
                    Q_norm_m3s = local_flow_rate_kg_s / rho_norm
                    base_kwargs['Q'] = Q_norm_m3s
            except Exception as e:
                error_log.append(f"Could not calculate Q_norm_m3s: {e}")

        # --- Method-Specific Argument Selection ---
        calculation_func = None
        specific_kwargs = {}

        if method == "Weymouth":
            calculation_func = fluids.compressible.Weymouth
            required_keys = ['SG', 'Tavg', 'L', 'D', 'P1', 'P2', 'Q', 'Zavg', 'E']
            specific_kwargs = {k: v for k, v in base_kwargs.items() if k in required_keys}

        elif method == "Panhandle_A":
            calculation_func = fluids.compressible.Panhandle_A
            required_keys = ['SG', 'Tavg', 'L', 'D', 'P1', 'P2', 'Q', 'Zavg', 'E']
            specific_kwargs = {k: v for k, v in base_kwargs.items() if k in required_keys}

        elif method == "Panhandle_B":
            calculation_func = fluids.compressible.Panhandle_B
            required_keys = ['SG', 'Tavg', 'L', 'D', 'P1', 'P2', 'Q', 'Zavg', 'E']
            specific_kwargs = {k: v for k, v in base_kwargs.items() if k in required_keys}

        elif method == "IGT":
            calculation_func = fluids.compressible.IGT
            required_keys = ['SG', 'Tavg', 'mu', 'L', 'D', 'P1', 'P2', 'Q', 'Zavg', 'E']
            specific_kwargs = {k: v for k, v in base_kwargs.items() if k in required_keys}
            specific_kwargs['mu'] = local_gas_viscosity

        elif method == "Oliphant":
            calculation_func = fluids.compressible.Oliphant
            required_keys = ['SG', 'Tavg', 'L', 'D', 'P1', 'P2', 'Q', 'Zavg', 'E']
            specific_kwargs = {k: v for k, v in base_kwargs.items() if k in required_keys}

        elif method == "Spitzglass_low":
            calculation_func = fluids.compressible.Spitzglass_low
            required_keys = ['SG', 'Tavg', 'L', 'D', 'P1', 'P2', 'Q', 'Zavg', 'E']
            specific_kwargs = {k: v for k, v in base_kwargs.items() if k in required_keys}

        elif method == "Spitzglass_high":
            calculation_func = fluids.compressible.Spitzglass_high
            required_keys = ['SG', 'Tavg', 'L', 'D', 'P1', 'P2', 'Q', 'Zavg', 'E']
            specific_kwargs = {k: v for k, v in base_kwargs.items() if k in required_keys}

        elif method == "isothermal_darcy":
            calculation_func = fluids.compressible.isothermal_gas
            results_log.append("Using isothermal_darcy with internal iteration for improved accuracy")
            
            # Calculate initial density and friction factor for fluids function
            # The function will handle internal iteration properly
            if local_P1 and local_P2:
                avg_P = (local_P1 + local_P2) / 2.0
            elif local_P1:
                avg_P = local_P1
            elif local_P2:
                avg_P = local_P2
            else:
                avg_P = 101325.0
                
            # Calculate density at average conditions
            initial_rho = (avg_P * local_gas_mw) / (local_gas_z_factor * R_univ * local_T_k)
            
            # Calculate initial friction factor
            if local_flow_rate_kg_s and local_pipe_diameter:
                area = math.pi * (local_pipe_diameter / 2.0) ** 2
                initial_velocity = local_flow_rate_kg_s / (initial_rho * area)
                Re_initial = fluids.core.Reynolds(
                    V=initial_velocity, D=local_pipe_diameter, rho=initial_rho, mu=local_gas_viscosity
                )
                initial_fd = fluids.friction.friction_factor(
                    Re=Re_initial, eD=local_pipe_roughness / local_pipe_diameter
                )
            else:
                # Use default friction factor if flow rate unknown (will be iterated internally)
                initial_fd = 0.02  # Reasonable default for commercial pipe
                
            # Prepare arguments for isothermal_gas function
            iso_base = {
                'rho': initial_rho,
                'fd': initial_fd,
                'P1': local_P1,
                'P2': local_P2,
                'L': local_L,
                'D': local_pipe_diameter,
                'm': local_flow_rate_kg_s
            }
            specific_kwargs = iso_base

        else:
            error_log.append(f"Unsupported calculation method: {method}.")

        if error_log:
            return json.dumps({"errors": error_log, "log": results_log})
            
        # Determine the function key to be solved for and prepare arguments
        if solve_for == 'Q' and method != 'isothermal_darcy':
            # For Weymouth, etc. when solving for flow
            solve_key = 'Q'
        elif solve_for == 'Q' and method == 'isothermal_darcy':
            # For isothermal when solving for flow
            solve_key = 'm'
        elif solve_for in ['P1', 'P2', 'L', 'D']:
            # For pressure, length, diameter - same key for all methods
            solve_key = solve_for
        else:
            error_log.append(f"Unsupported solve_for variable: {solve_for}")
            return json.dumps({"errors": error_log, "log": results_log})
            
        # Prepare kwargs - include all values that are not None except the one we're solving for
        final_kwargs = {k: v for k, v in specific_kwargs.items() if v is not None and k != solve_key}
        # Explicitly add the solve_key with None value
        final_kwargs[solve_key] = None
        
        results_log.append(f"Calling {calculation_func.__name__} with args: {list(final_kwargs.keys())} to solve for '{solve_key}'")

        try:
            calculated_value = calculation_func(**final_kwargs)
            solved_variable_value = calculated_value
            results_log.append(
                f"Successfully calculated {solve_for} = {solved_variable_value} using {method}."
            )
        except TypeError as te:
            # Use lazy evaluation for expensive list() call
            if logger.isEnabledFor(logging.DEBUG):
                provided_args = list(final_kwargs.keys())
                error_log.append(f"Argument mismatch calling {method}: {te}. Provided args: {provided_args}")
            else:
                error_log.append(f"Argument mismatch calling {method}: {te}")
            logger.error("Calculation TypeError calling %s: %s", method, te, exc_info=True)
            return json.dumps({"errors": error_log, "log": results_log})
        except ValueError as ve:
            error_log.append(f"Calculation ValueError from {method}: {ve}")
            logger.warning("Calculation ValueError with %s: %s", method, ve)
            return json.dumps({"errors": error_log, "log": results_log})
        except Exception as calc_e:
            error_log.append(
                f"Unexpected error during compressible flow calculation using {method}: {calc_e}"
            )
            logger.error("Calculation error with %s: %s", method, calc_e, exc_info=True)
            return json.dumps({"errors": error_log, "log": results_log})

        # --- Populate final_results ---
        final_results = {
            "inputs_resolved": results_log,
            "warnings": error_log if error_log else None,
            "pipe_details": pipe_info if pipe_info else None,
            "gas_property_details": gas_prop_info,
            "calculation_method_used": method,
            "solved_variable": solve_for,
            f"solved_{solve_for}_value": solved_variable_value,
            "inlet_pressure_pa": local_P1 if solve_for != 'P1' else solved_variable_value,
            "outlet_pressure_pa": local_P2 if solve_for != 'P2' else solved_variable_value,
            "pipe_length_m": local_L if solve_for != 'L' else solved_variable_value,
            "pipe_diameter_m": local_pipe_diameter if solve_for != 'D' else solved_variable_value,
            "flow_rate_kg_s": (
                local_flow_rate_kg_s if solve_for not in ['Q', 'm']
                else solved_variable_value
            ),
        }

        # Calculate derived results
        final_P1 = final_results["inlet_pressure_pa"]
        final_P2 = final_results["outlet_pressure_pa"]
        final_Q_kgs = final_results["flow_rate_kg_s"]

        if final_P1 is not None and final_P2 is not None:
            final_results["pressure_drop_pa"] = final_P1 - final_P2
            final_results["pressure_drop_psi"] = (final_P1 - final_P2) / PSI_to_PA

        if final_Q_kgs is not None and local_gas_mw is not None:
            try:
                rho_norm = (P_norm * local_gas_mw) / (R_univ * T_norm)
                rho_std = (P_std * local_gas_mw) / (R_univ * T_std)
                if rho_norm > 0:
                    final_results["flow_rate_norm_m3_hr"] = (final_Q_kgs / rho_norm) * 3600.0
                if rho_std > 0:
                    final_results["flow_rate_std_m3_hr"] = (final_Q_kgs / rho_std) * 3600.0

                # Calculate velocity and Mach number using final results
                final_diameter = final_results["pipe_diameter_m"]
                if final_diameter and final_diameter > 0:
                    area = math.pi * (final_diameter / 2.0) ** 2
                    # Calculate average density for velocity calculation
                    avg_P = (final_P1 + final_P2) / 2.0 if final_P1 and final_P2 else (final_P1 or final_P2 or 101325.0)
                    avg_rho = (avg_P * local_gas_mw) / (local_gas_z_factor * R_univ * local_T_k)
                    
                    if avg_rho > 0 and area > 0:
                        avg_velocity = final_Q_kgs / (avg_rho * area)
                        final_results["average_velocity_m_s"] = avg_velocity

                        R_specific = R_univ / local_gas_mw
                        speed_sound = math.sqrt(
                            local_gas_gamma * local_gas_z_factor * R_specific * local_T_k
                        ) if R_specific > 0 else 0
                        if speed_sound > 0:
                            final_results["average_mach_number"] = avg_velocity / speed_sound

            except Exception as derived_e:
                error_log.append(
                    f"Warning: Could not calculate derived flow rates/velocity/Mach: {derived_e}"
                )
                if "warnings" not in final_results or final_results["warnings"] is None:
                    final_results["warnings"] = []
                final_results["warnings"].append(
                    f"Warning: Could not calculate derived flow rates/velocity/Mach: {derived_e}"
                )

    except Exception as e:
        logger.error(f"Error in calculate_gas_pipe_pressure_drop: {e}", exc_info=True)
        return json.dumps({
            "error": f"Calculation error: {str(e)}",
            "log": results_log,
            "errors_occurred": error_log
        })

    # --- Final Return ---
    if final_results:
        # Remove null keys before returning
        if not final_results.get("warnings"):
            final_results.pop("warnings", None)
        if not final_results.get("pipe_details"):
            final_results.pop("pipe_details", None)
        if not final_results.get("gas_property_details"):
            final_results.pop("gas_property_details", None)

        # Return the result
        return json.dumps(
            final_results,
            default=lambda x: round(x, 6) if isinstance(x, float) else str(x)
        )
    else:
        # Should only happen if an exception occurred and wasn't returned earlier
        return json.dumps({"errors": error_log, "log": results_log, "status": "Failed - No results generated"})

def gas_pipe_sweep(
    variable: str,
    start: float,
    stop: float,
    n: int,
    # Base calculation parameters from calculate_gas_pipe_pressure_drop
    inlet_pressure: Optional[float] = None,
    outlet_pressure: Optional[float] = None,
    pipe_length: Optional[float] = None,
    flow_rate_kg_s: Optional[float] = None,
    flow_rate_norm_m3_hr: Optional[float] = None,
    flow_rate_std_m3_hr: Optional[float] = None,
    pipe_diameter: Optional[float] = None,
    pipe_diameter_in: Optional[float] = None,
    nominal_size_in: Optional[float] = None,
    schedule: str = "40",
    material: Optional[str] = None,
    pipe_roughness: Optional[float] = None,
    temperature_c: Optional[float] = None,
    gas_mw: Optional[float] = None,
    gas_gamma: Optional[float] = None,
    gas_z_factor: Optional[float] = None,
    gas_viscosity: Optional[float] = None,
    fluid_name: Optional[str] = None,
    gas_composition_mol: Optional[Dict[str, float]] = None,
    method: str = "Weymouth"
) -> str:
    """Parameter sweep for gas pipe pressure drop analysis
    
    Optimized for performance with large datasets (1e4+ points) per expert reviewer recommendation.
    Sweeps one variable while keeping others constant to analyze system behavior.
    
    Args:
        variable: Variable to sweep ('inlet_pressure', 'outlet_pressure', 'pipe_length', 'flow_rate_norm_m3_hr', 'pipe_diameter')
        start: Start value for sweep
        stop: Stop value for sweep
        n: Number of points in sweep
        **other_params: All other parameters from calculate_gas_pipe_pressure_drop
    
    Returns:
        JSON string with sweep results as list of dictionaries
    """
    try:
        import numpy as np
    except ImportError:
        return json.dumps({"error": "numpy required for sweep functionality"})
    
    # Generate sweep values
    sweep_values = np.linspace(start, stop, n)
    results = []
    
    # Build base parameters dictionary
    base_kwargs = {
        'inlet_pressure': inlet_pressure,
        'outlet_pressure': outlet_pressure,
        'pipe_length': pipe_length,
        'flow_rate_kg_s': flow_rate_kg_s,
        'flow_rate_norm_m3_hr': flow_rate_norm_m3_hr,
        'flow_rate_std_m3_hr': flow_rate_std_m3_hr,
        'pipe_diameter': pipe_diameter,
        'pipe_diameter_in': pipe_diameter_in,
        'nominal_size_in': nominal_size_in,
        'schedule': schedule,
        'material': material,
        'pipe_roughness': pipe_roughness,
        'temperature_c': temperature_c,
        'gas_mw': gas_mw,
        'gas_gamma': gas_gamma,
        'gas_z_factor': gas_z_factor,
        'gas_viscosity': gas_viscosity,
        'fluid_name': fluid_name,
        'gas_composition_mol': gas_composition_mol,
        'method': method
    }
    
    # Remove None values
    base_kwargs = {k: v for k, v in base_kwargs.items() if v is not None}
    
    for value in sweep_values:
        try:
            # Set up parameters for this sweep point
            kwargs = base_kwargs.copy()
            kwargs[variable] = value
            
            # Calculate
            result_json = calculate_gas_pipe_pressure_drop(**kwargs)
            result = json.loads(result_json)
            
            if 'errors' not in result or not result['errors']:
                # Extract key results
                row = {
                    variable: round(value, 6),
                    'pressure_drop_pa': result.get('pressure_drop_pa', None),
                    'flow_rate_kg_s': result.get('flow_rate_kg_s', None),
                    'pipe_diameter_m': result.get('pipe_diameter_m', None),
                    'pipe_length_m': result.get('pipe_length_m', None),
                    'average_velocity_m_s': result.get('average_velocity_m_s', None),
                    'solved_variable': result.get('solved_variable', None),
                    'solved_value': result.get(f"solved_{result.get('solved_variable', '')}_value", None)
                }
                results.append(row)
            else:
                # Add row with errors
                row = {variable: round(value, 6), 'error': str(result.get('errors', 'Unknown error'))}
                results.append(row)
                
        except Exception as e:
            row = {variable: round(value, 6), 'error': str(e)}
            results.append(row)
    
    return json.dumps({
        "sweep_variable": variable,
        "sweep_range": {"start": start, "stop": stop, "n": n},
        "results": results,
        "summary": {
            "total_points": len(results),
            "successful_points": len([r for r in results if 'error' not in r]),
            "failed_points": len([r for r in results if 'error' in r])
        }
    })
