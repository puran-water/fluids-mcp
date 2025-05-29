import json
import logging
from typing import Optional, Dict, Literal

# --- Third-Party Libraries ---
try:
    import fluids
    import fluids.constants # Needed for MWs lookup if available
    FLUIDS_AVAILABLE = True
except ImportError:
    fluids = None # Define fluids as None if not available
    FLUIDS_AVAILABLE = False
    # logging.warning("Fluids library not found. Composition MW lookup may be limited.")

# Import shared utilities
from utils.constants import (
    GPM_to_M3S, INCH_to_M, FT_to_M, PSI_to_PA, 
    LBFT3_to_KGM3, CENTIPOISE_to_PAS, DEFAULT_ROUGHNESS, DEG_C_to_K
)
from utils.helpers import get_fitting_K
from utils.import_helpers import FLUIDPROP_AVAILABLE, FluidProperties, FLUID_SELECTION

# Configure logging
logger = logging.getLogger("fluids-mcp.calculate_blower_compressor_requirements")

def calculate_blower_compressor_requirements(
    # --- Flow Rate ---
    flow_rate_kg_s: Optional[float] = None,        # Mass flow rate in kg/s
    flow_rate_norm_m3_hr: Optional[float] = None,  # Volumetric flow rate at Normal conditions (0°C, 1 atm) in m³/hr
    flow_rate_std_m3_hr: Optional[float] = None,   # Volumetric flow rate at Standard conditions (15°C, 1 atm) in m³/hr

    # --- Pressures (Absolute) ---
    inlet_pressure: Optional[float] = None,        # Inlet pressure (absolute) in Pa
    outlet_pressure: Optional[float] = None,       # Required outlet pressure (absolute) in Pa
    inlet_pressure_psi: Optional[float] = None,    # Inlet pressure (absolute) in psi
    outlet_pressure_psi: Optional[float] = None,   # Required outlet pressure (absolute) in psi

    # --- Temperature ---
    inlet_temperature_c: Optional[float] = None,   # Inlet temperature in Celsius

    # --- Gas Properties (Direct Input, Lookup, or Composition) ---
    gas_mw: Optional[float] = None,                # Gas Molecular Weight in kg/kmol (g/mol)
    gas_gamma: Optional[float] = None,             # Specific heat ratio (Cp/Cv)
    gas_z_factor: Optional[float] = None,          # Compressibility factor (Z) - use average if possible, else inlet
    fluid_name: Optional[str] = None,              # Common gas name ("Air", "Methane", "Biogas") for property lookup
    gas_composition_mol: Optional[Dict[str, float]] = None, # Mol fraction composition {component: fraction}

    # --- Efficiency ---
    efficiency_type: Literal['isentropic', 'polytropic'] = 'isentropic', # Type of efficiency provided
    efficiency: Optional[float] = None,              # Machine efficiency (decimal, e.g., 0.7 for 70%)

    # --- Polytropic Calculation Specific (Optional) ---
    polytropic_n: Optional[float] = None             # Polytropic exponent (if efficiency_type is 'polytropic' and efficiency is based on 'n')

) -> str:
    """Estimates discharge temperature and power requirements for a gas compressor or blower.

    Calculates based on isentropic or polytropic compression relationships.
    Accepts flexible inputs for flow, pressure, temperature, and gas properties.

    Args:
        # --- Flow ---
        flow_rate_kg_s: Optional - Mass flow rate in kg/s
        flow_rate_norm_m3_hr: Optional - Volumetric flow at Normal conditions (0°C, 1 atm) in m³/hr
        flow_rate_std_m3_hr: Optional - Volumetric flow at Standard conditions (15°C, 1 atm) in m³/hr
        # --- Pressures ---
        inlet_pressure: Optional - Inlet pressure (absolute) in Pa
        outlet_pressure: Optional - Outlet pressure (absolute) in Pa
        inlet_pressure_psi: Optional - Inlet pressure (absolute) in psi
        outlet_pressure_psi: Optional - Outlet pressure (absolute) in psi
        # --- Temperature ---
        inlet_temperature_c: Required - Inlet temperature in Celsius
        # --- Gas Properties ---
        gas_mw: Optional - Gas Molecular Weight in kg/kmol (g/mol)
        gas_gamma: Optional - Specific heat ratio (Cp/Cv)
        gas_z_factor: Optional - Compressibility factor (Z). If not provided, lookup/default=1.0 used. Average Z recommended if known.
        fluid_name: Optional - Common gas name ("Air", "Methane", "Biogas", "NaturalGas", etc.)
        gas_composition_mol: Optional - Gas composition as mol fractions {component: fraction}
        # --- Efficiency ---
        efficiency_type: Type of efficiency ('isentropic' or 'polytropic', default: 'isentropic')
        efficiency: Optional - Machine efficiency (0 to 1). Default depends on type (e.g., 0.7 isentropic, 0.75 polytropic).
        polytropic_n: Optional - Polytropic exponent 'n'. Used if efficiency_type='polytropic'. If not given, it's estimated from gamma and efficiency.

    Returns:
        JSON string with calculated discharge temperature, power requirements, and input resolution log.
    """
    results_log = []
    error_log = []
    gas_prop_info = {}

    # --- Standard Constants ---
    P_norm = 101325.0 # Pa
    T_norm = 273.15 # K
    P_std = 101325.0 # Pa
    T_std = 288.15 # K (15 C)
    R_univ = 8314.462 # J/(kmol·K)

    try:
        # 1. Resolve Inlet Temperature (K)
        if inlet_temperature_c is None:
            error_log.append("Missing required input: inlet_temperature_c.")
            # Cannot proceed without temperature
            return json.dumps({"errors": error_log, "log": results_log})
        local_T1_k = inlet_temperature_c + DEG_C_to_K # Convert to Kelvin
        results_log.append(f"Using Inlet Temperature: {inlet_temperature_c:.2f} C ({local_T1_k:.2f} K)")

        # 2. Resolve Pressures (Pa Absolute)
        local_P1 = None
        if inlet_pressure is not None:
            local_P1 = inlet_pressure
            results_log.append("Used provided SI inlet_pressure.")
        elif inlet_pressure_psi is not None:
            local_P1 = inlet_pressure_psi * PSI_to_PA
            results_log.append(f"Converted inlet_pressure from {inlet_pressure_psi} psi.")
        else:
            error_log.append("Missing required input: inlet_pressure or inlet_pressure_psi.")

        local_P2 = None
        if outlet_pressure is not None:
            local_P2 = outlet_pressure
            results_log.append("Used provided SI outlet_pressure.")
        elif outlet_pressure_psi is not None:
            local_P2 = outlet_pressure_psi * PSI_to_PA
            results_log.append(f"Converted outlet_pressure from {outlet_pressure_psi} psi.")
        else:
            error_log.append("Missing required input: outlet_pressure or outlet_pressure_psi.")

        # Check P2 > P1 after resolution attempt
        if local_P1 is not None and local_P2 is not None and local_P2 <= local_P1:
            error_log.append("Outlet pressure must be greater than inlet pressure for compression.")
            # Return now if pressures are invalid
            if error_log: return json.dumps({"errors": error_log, "log": results_log})

        # 3. Resolve Gas Properties (MW, gamma, Z) - Reusing logic structure
        local_gas_mw = gas_mw
        local_gas_gamma = gas_gamma
        # Use provided Z if available (assumed average Z), else lookup at inlet, else default 1.0
        local_gas_z_factor = gas_z_factor
        # Initialize local_gas_viscosity with None
        local_gas_viscosity = None
        gas_prop_source = "Provided"
        prop_resolved = True # Flag if base properties are found

        if any(p is None for p in [local_gas_mw, local_gas_gamma]): # Z has fallback default
            if fluid_name is not None:
                gas_prop_source = f"Lookup ({fluid_name} @ {inlet_temperature_c} C)" # Use inlet conditions for lookup
                results_log.append(f"Attempting property lookup for '{fluid_name}'.")
                if FLUIDPROP_AVAILABLE:
                    try:
                        # Use average pressure for lookup if P1/P2 are known, else use 1 atm default
                        avg_p_pa = 101325.0
                        if inlet_pressure is not None and outlet_pressure is not None:
                            avg_p_pa = (inlet_pressure + outlet_pressure) / 2.0
                        elif inlet_pressure is not None:
                            avg_p_pa = inlet_pressure # Approximate with inlet if outlet unknown
                        # Convert average pressure from Pa to bar for lookup
                        lookup_p_bar = avg_p_pa / 100000.0

                        # Try FluidProperties directly first
                        try:
                            fluid_props = FluidProperties(coolprop_name=fluid_name, T_in_deg_C=inlet_temperature_c, P_in_bar=lookup_p_bar)
                            actual_fluid_name = fluid_name
                        except Exception:
                            # Fallback to FLUID_SELECTION validation if direct lookup fails
                            valid_fluids = [f[0] for f in FLUID_SELECTION]
                            actual_fluid_name = fluid_name
                            if fluid_name not in valid_fluids:
                                match = next((f for f in valid_fluids if f.lower() == fluid_name.lower()), None)
                                if match: actual_fluid_name = match
                                else: raise ValueError(f"Fluid '{fluid_name}' not found in fluidprop.")

                            fluid_props = FluidProperties(coolprop_name=actual_fluid_name, T_in_deg_C=inlet_temperature_c, P_in_bar=lookup_p_bar)
                        if local_gas_mw is None: local_gas_mw = float(fluid_props.MW) * 1000.0  # Convert kg/mol to kg/kmol
                        if local_gas_gamma is None:
                            # Safely access gamma attribute or calculate from Cp/Cv if available
                            if hasattr(fluid_props, 'gamma'):
                                local_gas_gamma = float(fluid_props.gamma[0])
                            elif hasattr(fluid_props, 'Cp') and hasattr(fluid_props, 'Cv') and float(fluid_props.Cv[0]) != 0:
                                local_gas_gamma = float(fluid_props.Cp[0]) / float(fluid_props.Cv[0])
                        if local_gas_z_factor is None:
                            # Safely access Z attribute
                            if hasattr(fluid_props, 'Z'):
                                local_gas_z_factor = float(fluid_props.Z[0])
                        gas_prop_info = {"name_used": actual_fluid_name, "lookup_temp_c": inlet_temperature_c, "lookup_p_bar": lookup_p_bar}
                        prop_resolved = True

                    except Exception as prop_lookup_e:
                        error_log.append(f"Warning: Property lookup failed for '{fluid_name}': {prop_lookup_e}. Falling back to defaults.")
                        gas_prop_source = "Lookup Failed - Using Defaults"
                        prop_resolved = False # Mark as failed
                else:
                    error_log.append("Warning: Property lookup skipped (fluidprop unavailable). Falling back to defaults.")
                    gas_prop_source = "Defaults (fluidprop unavailable)"
                    prop_resolved = False

            elif gas_composition_mol is not None:
                gas_prop_source = f"Composition Calculation"
                results_log.append("Attempting property calculation from composition.")
                # --- Placeholder for actual mixture calculation ---
                # This needs a proper implementation using mixing rules or a library
                # Example: Use fluids.constants.MWs for MW, estimate others
                temp_mw = None
                try:
                    # Check if fluids.constants.MWs is available (might not be directly accessible depending on install/version)
                    if hasattr(fluids, 'constants') and hasattr(fluids.constants, 'MWs'):
                        temp_mw = sum(fluids.constants.MWs[comp] * frac for comp, frac in gas_composition_mol.items() if comp in fluids.constants.MWs)
                    else: # Fallback if constants not easily accessible
                        # Very rough MW estimation based on common biogas components
                        mw_map = {'CH4': 16040, 'CO2': 44010, 'N2': 28010, 'H2S': 34080, 'O2': 32000, 'H2': 2016}  # kg/kmol
                        temp_mw = sum(mw_map.get(comp, 28960) * frac for comp, frac in gas_composition_mol.items()) # Default to air MW if unknown

                    if temp_mw > 0 and local_gas_mw is None:
                        local_gas_mw = temp_mw
                        prop_resolved = True # At least MW resolved
                    elif local_gas_mw is None:
                        error_log.append("Warning: Could not calculate MW from composition.")
                        prop_resolved = False
                except Exception as comp_e:
                    error_log.append(f"Error during composition MW calculation: {comp_e}")
                    prop_resolved = False

                # Estimate others (highly approximate - replace with better methods)
                if local_gas_gamma is None: local_gas_gamma = 1.3 # Typical biogas/natural gas range
                if local_gas_z_factor is None: local_gas_z_factor = 0.98 # Estimate
                # Viscosity would also need mixing rules - skipping estimate here, will default later if needed
                gas_prop_info = {"composition_mol": gas_composition_mol}
                error_log.append("Warning: Using estimated gamma & Z for mixture based on composition.")
                # --- End Placeholder ---

            # Apply defaults if properties still missing or lookup/composition failed
            if not prop_resolved or local_gas_mw is None:
                local_gas_mw = 28.96; results_log.append("Using default MW (Air)."); gas_prop_source += " + MW Default"
            if not prop_resolved or local_gas_gamma is None:
                local_gas_gamma = 1.4; results_log.append("Using default gamma (Air)."); gas_prop_source += " + Gamma Default"

            # Handle Z factor default separately if it wasn't set by lookup/input
            if local_gas_z_factor is None:
                local_gas_z_factor = 1.0; results_log.append("Using default Z-factor (Ideal)."); gas_prop_source += " + Z Default"

        gas_prop_info.update({
            "MW_kg_kmol": local_gas_mw, "gamma": local_gas_gamma,
            "Z_factor_used": local_gas_z_factor, # Clarify Z used (might be inlet Z or average estimate)
            "viscosity": local_gas_viscosity, # Placeholder for viscosity if needed
            "source": gas_prop_source
        })
        results_log.append(f"Resolved Gas Properties: MW={local_gas_mw:.2f}, gamma={local_gas_gamma:.3f}, Z={local_gas_z_factor:.4f} (Source: {gas_prop_source})")

        # 4. Resolve Flow Rate (SI Mass Flow - kg/s)
        local_flow_rate_kg_s = None
        flow_rate_source = "N/A"
        if flow_rate_kg_s is not None:
            local_flow_rate_kg_s = flow_rate_kg_s
            flow_rate_source = "Provided mass flow (kg/s)"
        elif flow_rate_norm_m3_hr is not None:
            if local_gas_mw is None: error_log.append("Cannot convert Norm m3/hr without MW."); # Continue check below
            else:
                rho_norm = (P_norm * local_gas_mw) / (1.0 * R_univ * T_norm)
                local_flow_rate_kg_s = (flow_rate_norm_m3_hr / 3600.0) * rho_norm
                flow_rate_source = f"Converted from Normal m³/hr ({rho_norm=:.3f} kg/m³)"
        elif flow_rate_std_m3_hr is not None:
            if local_gas_mw is None: error_log.append("Cannot convert Std m3/hr without MW."); # Continue check below
            else:
                rho_std = (P_std * local_gas_mw) / (1.0 * R_univ * T_std)
                local_flow_rate_kg_s = (flow_rate_std_m3_hr / 3600.0) * rho_std
                flow_rate_source = f"Converted from Standard m³/hr ({rho_std=:.3f} kg/m³)"
        else:
            error_log.append("Missing required input: flow_rate_kg_s or flow_rate_norm_m3_hr or flow_rate_std_m3_hr.")

        if local_flow_rate_kg_s is not None:
            results_log.append(f"Resolved Flow Rate: {local_flow_rate_kg_s:.4f} kg/s (Source: {flow_rate_source})")
        else:
            # If flow rate resolution failed (e.g., needed MW but MW failed) add specific error
            if "Missing required input: flow_rate" not in str(error_log): # Avoid duplicate
                error_log.append("Could not resolve mass flow rate (e.g., missing MW for conversion).")

        # 5. Resolve Efficiency
        local_efficiency = efficiency
        local_polytropic_n = polytropic_n
        if local_efficiency is None:
            if efficiency_type == 'isentropic':
                local_efficiency = 0.70 # Default isentropic efficiency
            else: # polytropic
                local_efficiency = 0.75 # Default polytropic efficiency
            results_log.append(f"Using default {efficiency_type} efficiency: {local_efficiency:.2f}")
        else:
            results_log.append(f"Using provided {efficiency_type} efficiency: {local_efficiency:.2f}")

        if not (0 < local_efficiency <= 1.0):
            error_log.append("Efficiency must be between 0 (exclusive) and 1 (inclusive).")

        # --- Check for critical errors before calculation ---
        critical_inputs = [
            local_T1_k, local_P1, local_P2, local_gas_mw, local_gas_gamma,
            local_gas_z_factor, local_flow_rate_kg_s, local_efficiency
        ]
        # Check if any critical input values are still None OR if previous critical errors occurred
        if any(v is None for v in critical_inputs) or any(e.startswith("Missing required input") or "pressure must be greater" in e for e in error_log):
            if any(v is None for v in critical_inputs):
                missing_items = [name for name, val in zip([
                    "Inlet Temp", "Inlet Pressure", "Outlet Pressure", "MW", "gamma",
                    "Z Factor", "Mass Flow", "Efficiency"
                    ], critical_inputs) if val is None]
                error_log.append(f"Critical input parameter resolution failed for: {', '.join(missing_items)}")
            return json.dumps({"errors": error_log, "log": results_log})

        # --- Core Calculation Logic ---
        P_ratio = local_P2 / local_P1
        k = local_gas_gamma
        Z_avg = local_gas_z_factor # Using inlet Z as average Z estimate - improve if possible
        R_specific = R_univ / local_gas_mw # J/(kg·K)

        T2_s = 0.0 # Isentropic discharge temp (K)
        W_s = 0.0  # Isentropic work (J/kg)
        T2_p = 0.0 # Polytropic discharge temp (K)
        W_p = 0.0  # Polytropic work (J/kg)
        n_poly_calc = None # Calculated polytropic n

        # Isentropic Calculations
        k_exp = (k - 1.0) / k
        if k_exp == 0 or P_ratio < 0: # Avoid math errors
            error_log.append("Invalid inputs for isentropic calculation (k=1 or negative P_ratio).")
            return json.dumps({"errors": error_log, "log": results_log})

        T2_s = local_T1_k * (P_ratio)**k_exp
        # Isentropic Work: Ws = Z_avg*R_spec*T1*k/(k-1) * [ (P2/P1)^((k-1)/k) - 1 ]
        W_s = (Z_avg * R_specific * local_T1_k * k / (k - 1.0)) * ( (P_ratio)**k_exp - 1.0 )

        # Actual Discharge Temp and Work based on efficiency type
        W_actual = 0.0
        T2_actual_k = 0.0

        if efficiency_type == 'isentropic':
            if local_efficiency == 0: error_log.append("Isentropic efficiency cannot be zero."); return json.dumps({"errors": error_log, "log": results_log})
            W_actual = W_s / local_efficiency
            # W_actual = Cp * (T2_actual - T1) = (k*R_spec / (k-1)) * (T2_actual - T1)
            # T2_actual = T1 + W_actual * (k-1) / (k*R_spec)
            # Check k*R_specific is not zero
            if abs(k * R_specific) < 1e-12:
                error_log.append("Invalid gas properties (k or R_specific is zero) for temperature calculation.")
                return json.dumps({"errors": error_log, "log": results_log})
            T2_actual_k = local_T1_k + (W_actual * (k - 1.0)) / (k * R_specific)
            results_log.append("Calculated using isentropic efficiency.")

        elif efficiency_type == 'polytropic':
            # Determine polytropic exponent n
            n = local_polytropic_n # Use provided n if available
            if n is None:
                # Estimate n from polytropic efficiency
                eta_p = local_efficiency
                k = local_gas_gamma
                # eta_p = [(k-1)/k]/[(n-1)/n]
                # eta_p = (k-1)/k * n/(n-1) => (n-1)/n = (k-1)/(k*eta_p) => 1 - 1/n = (k-1)/(k*eta_p) => 1/n = 1 - (k-1)/(k*eta_p)
                # 1/n = (k*eta_p - k + 1) / (k*eta_p) => n = k*eta_p / (k*eta_p - k + 1)
                n_num = k * eta_p
                n_den = k * eta_p - k + 1.0
                if abs(n_den) < 1e-9:
                    error_log.append("Cannot calculate polytropic exponent 'n' from efficiency (denominator near zero). Check inputs.")
                    return json.dumps({"errors": error_log, "log": results_log})
                n = n_num / n_den
                n_poly_calc = n # Store calculated n
                results_log.append(f"Calculated polytropic exponent n = {n:.4f} from efficiency.")

            if n <= 1.0: # Check if calculated n is physically plausible
                error_log.append(f"Calculated/Provided polytropic exponent n={n:.4f} is not valid (must be > 1). Check inputs/efficiency.")
                return json.dumps({"errors": error_log, "log": results_log})

            # Polytropic Discharge Temp: T2p/T1 = (P2/P1)^((n-1)/n)
            n_exp = (n - 1.0) / n
            if n_exp == 0 or P_ratio < 0: # Avoid math errors
                error_log.append("Invalid inputs for polytropic calculation (n=1 or negative P_ratio).")
                return json.dumps({"errors": error_log, "log": results_log})

            T2_p = local_T1_k * (P_ratio)**n_exp
            # Polytropic Work: Wp = Z_avg*R_spec*T1*n/(n-1) * [ (P2/P1)^((n-1)/n) - 1 ]
            W_p = (Z_avg * R_specific * local_T1_k * n / (n - 1.0)) * ( (P_ratio)**n_exp - 1.0 )
            # Assuming the definition where W_p IS the actual work based on the path defined by 'n'
            W_actual = W_p
            T2_actual_k = T2_p
            results_log.append(f"Calculated using polytropic relations with n={n:.4f}.")
        else:
            # Should not happen due to Literal check
            error_log.append(f"Unknown efficiency type: {efficiency_type}")
            return json.dumps({"errors": error_log, "log": results_log})

        # Total Power (Watts) = mass_flow_rate (kg/s) * work (J/kg)
        Power_watts = local_flow_rate_kg_s * W_actual

        # --- Consolidate Results ---
        final_result = {
            "inputs_resolved": results_log,
            "warnings": error_log if error_log else None,
            "gas_property_details": gas_prop_info,
            # Input Conditions Summary
            "inlet_pressure_pa": round(local_P1, 2),
            "outlet_pressure_pa": round(local_P2, 2),
            "pressure_ratio": round(P_ratio, 4),
            "inlet_temperature_c": round(inlet_temperature_c, 2),
            "inlet_temperature_k": round(local_T1_k, 2),
            "flow_rate_kg_s": round(local_flow_rate_kg_s, 4),
            # Efficiency Details
            "efficiency_type_used": efficiency_type,
            "efficiency_used": round(local_efficiency, 4),
            "polytropic_n_calculated": round(n_poly_calc, 4) if n_poly_calc is not None else None,
            "polytropic_n_input": local_polytropic_n,
            # Calculated Results
            "isentropic_discharge_temp_c": round(T2_s - DEG_C_to_K, 2),
            "isentropic_discharge_temp_k": round(T2_s, 2),
            "isentropic_work_kj_kg": round(W_s / 1000.0, 3),
            "actual_discharge_temp_c": round(T2_actual_k - DEG_C_to_K, 2),
            "actual_discharge_temp_k": round(T2_actual_k, 2),
            "actual_specific_work_kj_kg": round(W_actual / 1000.0, 3),
            "required_power_kw": round(Power_watts / 1000.0, 3),
            "required_power_hp": round(Power_watts / 745.7, 3)
        }

        # Add derived values if possible
        if local_flow_rate_kg_s is not None and local_gas_mw is not None:
            try:
                # Convert back to volumetric flows
                rho_norm = (P_norm * local_gas_mw) / (1.0 * R_univ * T_norm)
                rho_std = (P_std * local_gas_mw) / (1.0 * R_univ * T_std)
                final_result["flow_rate_norm_m3_hr"] = (local_flow_rate_kg_s / rho_norm) * 3600.0
                final_result["flow_rate_std_m3_hr"] = (local_flow_rate_kg_s / rho_std) * 3600.0
            except Exception as derived_e:
                error_log.append(f"Warning: Could not calculate derived flow rates: {derived_e}")

        # Remove null keys
        if not final_result.get("warnings"): final_result.pop("warnings", None)
        if final_result.get("polytropic_n_calculated") is None: final_result.pop("polytropic_n_calculated", None)
        if final_result.get("polytropic_n_input") is None: final_result.pop("polytropic_n_input", None)

        return json.dumps(final_result, default=lambda x: round(x, 6) if isinstance(x, float) else x) # Round floats for cleaner output

    except Exception as e:
        logger.error(f"Error in calculate_blower_compressor_requirements: {e}", exc_info=True)
        # Add context to the error message
        return json.dumps({
            "error": f"Calculation error: {str(e)}",
            "log": results_log,
            "errors_occurred": error_log
        })
