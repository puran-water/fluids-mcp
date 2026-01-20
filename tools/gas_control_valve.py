import json
import logging
import math
from typing import Optional, Dict, Literal

# --- Third-Party Libraries ---
try:
    import fluids
    import fluids.control_valve # Explicitly import submodule used
    import fluids.constants     # Potentially used for MW lookup
    FLUIDS_AVAILABLE = True
except ImportError:
    fluids = None # Define fluids as None if not available
    FLUIDS_AVAILABLE = False
    # logging.warning("Fluids library not found. Some functionality might be limited.")

# Import shared utilities
from utils.constants import (
    GPM_to_M3S, INCH_to_M, FT_to_M, PSI_to_PA, 
    LBFT3_to_KGM3, CENTIPOISE_to_PAS, DEFAULT_ROUGHNESS, DEG_C_to_K
)
from utils.helpers import get_fitting_K
from utils.import_helpers import FLUIDPROP_AVAILABLE, FluidProperties, FLUID_SELECTION
from utils.fluid_aliases import map_fluid_name

# Configure logging
logger = logging.getLogger("fluids-mcp.gas_control_valve")

def calculate_gas_control_valve(
    # --- Flow Rate ---
    flow_rate_kg_s: Optional[float] = None,        # Mass flow rate in kg/s
    flow_rate_norm_m3_hr: Optional[float] = None,  # Volumetric flow rate at Normal conditions (0°C, 1 atm) in m³/hr
    flow_rate_std_m3_hr: Optional[float] = None,   # Volumetric flow rate at Standard conditions (15°C, 1 atm) in m³/hr,

    # --- Pressures (Absolute) ---
    inlet_pressure: Optional[float] = None,        # Inlet pressure (absolute) in Pa
    outlet_pressure: Optional[float] = None,       # Outlet pressure (absolute) in Pa
    inlet_pressure_psi: Optional[float] = None,    # Inlet pressure (absolute) in psi
    outlet_pressure_psi: Optional[float] = None,   # Outlet pressure (absolute) in psi

    # --- Temperature ---
    inlet_temperature_c: Optional[float] = None,   # Inlet temperature in Celsius

    # --- Gas Properties (Direct Input, Lookup, or Composition) ---
    gas_mw: Optional[float] = None,                # Gas Molecular Weight in kg/kmol (g/mol)
    gas_gamma: Optional[float] = None,             # Specific heat ratio (Cp/Cv)
    gas_z_factor: Optional[float] = None,          # Compressibility factor (Z)
    gas_viscosity: Optional[float] = None,         # Gas dynamic viscosity in Pa·s
    fluid_name: Optional[str] = None,              # Common gas name ("Air", "Methane", "Biogas") for property lookup
    gas_composition_mol: Optional[Dict[str, float]] = None, # Mol fraction composition

    # --- Valve Parameters ---
    valve_type: Literal["globe", "butterfly", "ball", "gate", "other"] = "globe", # Valve type
    valve_xt: Optional[float] = None,              # Pressure drop ratio factor (overrides default)
    valve_fl: Optional[float] = None,              # Liquid pressure recovery factor (usually 1 for gas)
    
    # --- Output Formatting ---
    size_units: Literal["m", "inch", "mm"] = "inch"
) -> str:
    """Calculate gas control valve coefficient (Cv/Kv) and estimate size with flexible inputs."""
    results_log = []
    error_log = []
    gas_prop_info = {}

    # --- Standard Constants ---
    P_norm = 101325.0 # Pa
    T_norm = 273.15   # K
    P_std = 101325.0  # Pa
    T_std = 288.15    # K (15 C)
    R_univ = 8314.462 # J/(kmol·K)

    try:
        # 1. Resolve Inlet Temperature (K)
        if inlet_temperature_c is None:
            error_log.append("Missing required input: inlet_temperature_c.")
            return json.dumps({"errors": error_log, "log": results_log})
        local_T1_k = inlet_temperature_c + DEG_C_to_K
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

        # Check validity: P1 must be > P2 for a control valve
        if local_P1 is not None and local_P2 is not None:
            if local_P2 >= local_P1:
                error_log.append("Inlet pressure must be greater than outlet pressure for a control valve.")
                if error_log: 
                    return json.dumps({"errors": error_log, "log": results_log})
            delta_P = local_P1 - local_P2
        else:
            delta_P = None

        # 3. Resolve Gas Properties
        local_gas_mw = gas_mw
        local_gas_gamma = gas_gamma
        local_gas_z_factor = gas_z_factor
        local_gas_viscosity = gas_viscosity
        # FIX: Initialize as a dictionary, not a string
        gas_prop_source = {}  
        prop_resolved = True

        # If any property is missing, attempt lookup or composition approach
        if any(x is None for x in [local_gas_mw, local_gas_gamma, local_gas_z_factor, local_gas_viscosity]):
            # Attempt fluid-name-based lookup if fluidprop is available
            if fluid_name is not None:
                gas_prop_source["source_info"] = f"Lookup ({fluid_name} @ {inlet_temperature_c} C)"
                results_log.append(f"Attempting property lookup for '{fluid_name}'.")
                if FLUIDPROP_AVAILABLE and FLUID_SELECTION is not None and FluidProperties is not None:
                    try:
                        lookup_p_bar = (local_P1 / 100000.0) if local_P1 else 1.01325

                        # Try FluidProperties directly first
                        # Map fluid name through aliasing system first
                        actual_fluid_name = map_fluid_name(fluid_name)

                        try:
                            fluid_props = FluidProperties(
                                coolprop_name=actual_fluid_name,
                                T_in_deg_C=inlet_temperature_c,
                                P_in_bar=lookup_p_bar
                            )
                        except Exception:
                            # If mapped name fails, try FLUID_SELECTION validation
                            try:
                                valid_fluids = [f[0] for f in FLUID_SELECTION if f is not None and hasattr(f, '__getitem__')]
                            except (TypeError, IndexError):
                                valid_fluids = []
                            if not valid_fluids or actual_fluid_name not in valid_fluids:
                                match = next((f for f in valid_fluids if f.lower() == actual_fluid_name.lower()), None)
                                if match:
                                    actual_fluid_name = match
                                else:
                                    raise ValueError(f"Fluid '{fluid_name}' (mapped to '{actual_fluid_name}') not found in fluidprop.")

                            fluid_props = FluidProperties(
                                coolprop_name=actual_fluid_name,
                                T_in_deg_C=inlet_temperature_c,
                                P_in_bar=lookup_p_bar
                            )
                        if local_gas_mw is None: 
                            local_gas_mw = float(fluid_props.MW)  # Already in kg/kmol (numerically equal to g/mol)
                            gas_prop_source["MW"] = "Lookup"
                        if local_gas_gamma is None:
                            # Safely access gamma attribute or calculate from Cp/Cv if available
                            if hasattr(fluid_props, 'gamma'):
                                local_gas_gamma = float(fluid_props.gamma[0])
                                gas_prop_source["gamma"] = "Lookup"
                            elif hasattr(fluid_props, 'Cp') and hasattr(fluid_props, 'Cv') and float(fluid_props.Cv[0]) != 0:
                                local_gas_gamma = float(fluid_props.Cp[0]) / float(fluid_props.Cv[0])
                                gas_prop_source["gamma"] = "Lookup (Calculated Cp/Cv)"
                        if local_gas_z_factor is None:
                            # Safely access Z attribute
                            if hasattr(fluid_props, 'Z'):
                                local_gas_z_factor = float(fluid_props.Z[0])
                                gas_prop_source["Z"] = "Lookup"
                        if local_gas_viscosity is None: 
                            local_gas_viscosity = float(fluid_props.eta[0])
                            gas_prop_source["viscosity"] = "Lookup"

                        gas_prop_info.update({
                            "name_used": actual_fluid_name,
                            "lookup_temp_c": inlet_temperature_c,
                            "lookup_p_bar": lookup_p_bar,
                            "viscosity_pas": local_gas_viscosity
                        })
                        prop_resolved = True
                    except Exception as prop_lookup_e:
                        error_log.append(f"Warning: Property lookup failed for '{fluid_name}': {prop_lookup_e}. Falling back to defaults.")
                        gas_prop_source["source_info"] = "Lookup Failed - Using Defaults"
                        prop_resolved = False
                else:
                    error_log.append("Warning: Property lookup skipped (fluidprop unavailable). Falling back to defaults.")
                    gas_prop_source["source_info"] = "Defaults (fluidprop unavailable)"
                    prop_resolved = False

            # If composition is provided, do a rudimentary mixture approach
            elif gas_composition_mol is not None:
                gas_prop_source["source_info"] = "Composition Calculation"
                results_log.append("Attempting property calculation from composition.")
                try:
                    # Simplified approach using built-in or custom partial MW map
                    temp_mw = None
                    if hasattr(fluids, 'constants') and hasattr(fluids.constants, 'MWs'):
                        temp_mw = sum(
                            fluids.constants.MWs[comp] * frac
                            for comp, frac in gas_composition_mol.items()
                            if comp in fluids.constants.MWs
                        )
                    else:
                        mw_map = {'CH4': 16.04, 'CO2': 44.01, 'N2': 28.01, 'H2S': 34.08, 'O2': 32.00, 'H2': 2.016}  # kg/kmol
                        temp_mw = sum(mw_map.get(comp, 28.96) * frac for comp, frac in gas_composition_mol.items())

                    if temp_mw and temp_mw > 0 and local_gas_mw is None:
                        local_gas_mw = temp_mw
                        gas_prop_source["MW"] = "Composition"
                        prop_resolved = True
                    elif local_gas_mw is None:
                        error_log.append("Warning: Could not calculate MW from composition.")
                        prop_resolved = False
                except Exception as comp_e:
                    error_log.append(f"Error during composition MW calculation: {comp_e}")
                    prop_resolved = False

                # Estimate gamma, Z, viscosity if still None
                if local_gas_gamma is None: 
                    local_gas_gamma = 1.3
                    gas_prop_source["gamma"] = "Estimated"
                if local_gas_z_factor is None:
                    local_gas_z_factor = 0.98
                    gas_prop_source["Z"] = "Estimated"
                if local_gas_viscosity is None:
                    local_gas_viscosity = 1.1e-5
                    gas_prop_source["viscosity"] = "Estimated"
                gas_prop_info["composition_mol"] = gas_composition_mol
                error_log.append("Warning: Using estimated gamma, Z, viscosity for mixture.")

            # If still missing viscosity, apply a default
            if local_gas_viscosity is None:
                local_gas_viscosity = 1.8e-5
                gas_prop_source["viscosity"] = "Default (Air)"

            # Apply MW, gamma, Z defaults if needed
            if not prop_resolved or local_gas_mw is None:
                local_gas_mw = 28.96
                results_log.append("Using default MW (Air).")
                gas_prop_source["MW"] = "Default (Air)"
            if not prop_resolved or local_gas_gamma is None:
                local_gas_gamma = 1.4
                results_log.append("Using default gamma (Air).")
                gas_prop_source["gamma"] = "Default (Air)"
            if local_gas_z_factor is None:
                local_gas_z_factor = 1.0
                results_log.append("Using default Z-factor (Ideal).")
                gas_prop_source["Z"] = "Default (Ideal)"

        else:
            # All properties were provided explicitly
            gas_prop_source["MW"] = "Provided"
            gas_prop_source["gamma"] = "Provided"
            gas_prop_source["Z"] = "Provided"
            gas_prop_source["viscosity"] = "Provided"
            gas_prop_source["source_info"] = "All properties provided by user"

        # Store final gas properties
        gas_prop_info.update({
            "MW_kg_kmol": local_gas_mw,
            "gamma": local_gas_gamma,
            "Z_factor": local_gas_z_factor,
            "viscosity_pas": local_gas_viscosity,
            "source_details": gas_prop_source
        })
        results_log.append(
            f"Resolved Gas Properties: MW={local_gas_mw:.2f}, "
            f"gamma={local_gas_gamma:.3f}, Z={local_gas_z_factor:.4f}, "
            f"mu={local_gas_viscosity:.2e} Pa.s (Source: {gas_prop_source})"
        )

        # 4. Resolve Flow Rate (SI Mass Flow - kg/s)
        local_flow_rate_kg_s = None
        flow_rate_source = "N/A"
        if flow_rate_kg_s is not None:
            local_flow_rate_kg_s = flow_rate_kg_s
            flow_rate_source = "Provided mass flow (kg/s)"
        elif flow_rate_norm_m3_hr is not None:
            if local_gas_mw is None:
                error_log.append("Cannot convert Norm m3/hr without MW.")
            else:
                rho_norm = (P_norm * local_gas_mw) / (R_univ * T_norm)  # Z=1 at normal
                local_flow_rate_kg_s = (flow_rate_norm_m3_hr / 3600.0) * rho_norm
                flow_rate_source = f"Converted from Normal m³/hr (rho_norm={rho_norm:.3f} kg/m³)"
        elif flow_rate_std_m3_hr is not None:
            if local_gas_mw is None:
                error_log.append("Cannot convert Std m3/hr without MW.")
            else:
                rho_std = (P_std * local_gas_mw) / (R_univ * T_std)  # Z=1 at standard
                local_flow_rate_kg_s = (flow_rate_std_m3_hr / 3600.0) * rho_std
                flow_rate_source = f"Converted from Standard m³/hr (rho_std={rho_std:.3f} kg/m³)"
        else:
            error_log.append("Missing required input: flow_rate_kg_s or flow_rate_norm_m3_hr or flow_rate_std_m3_hr.")

        if local_flow_rate_kg_s is not None:
            results_log.append(
                f"Resolved Flow Rate: {local_flow_rate_kg_s:.4f} kg/s (Source: {flow_rate_source})"
            )
        else:
            error_log.append("Could not resolve mass flow rate.")

        # --- We will skip pipe-diameter-related checks for control valves;
        #     not all code references are needed here. This block is in the snippet,
        #     but for brevity, we omit advanced friction/pipe checks. ---

        # 5. Resolve Valve Factors (xt, FL, Fk)
        local_xt = valve_xt
        local_FL = valve_fl
        local_Fk = None

        xt_defaults = {"globe": 0.75, "butterfly": 0.35, "ball": 0.60, "gate": 0.80, "other": 0.70}
        if local_xt is None:
            local_xt = xt_defaults.get(valve_type, 0.70)
            results_log.append(f"Used default xt={local_xt} for valve type '{valve_type}'.")
        else:
            results_log.append(f"Used provided valve_xt={local_xt}.")

        if local_FL is None:
            local_FL = 1.0
            results_log.append("Used default FL=1.0 for gas service.")
        else:
            results_log.append(f"Used provided valve_fl={local_FL}.")

        if local_gas_gamma is not None and local_gas_gamma > 1.0:
            local_Fk = local_gas_gamma / (local_gas_gamma - 1.0)
            results_log.append(f"Calculated Fk={local_Fk:.4f} from gamma={local_gas_gamma:.3f}.")
        else:
            error_log.append("Cannot calculate Fk: Invalid or missing gamma.")

        # Verify critical inputs for calculation
        critical_inputs_valve = [
            local_flow_rate_kg_s, local_T1_k, local_gas_mw, local_gas_z_factor,
            local_gas_gamma, local_P1, local_P2, local_xt, local_Fk, local_FL
        ]
        if any(v is None for v in critical_inputs_valve) or any(
            e.startswith("Missing required input") or "pressure must be greater" in e for e in error_log
        ):
            if any(v is None for v in critical_inputs_valve):
                missing_items = [
                    name for name, val in zip(
                        ["Mass Flow", "Inlet Temp", "MW", "Z Factor", "gamma",
                         "Inlet Pressure", "Outlet Pressure", "xt", "Fk", "FL"],
                        critical_inputs_valve
                    ) if val is None
                ]
                error_log.append(f"Critical input parameter resolution failed for: {', '.join(missing_items)}")
            return json.dumps({"errors": error_log, "log": results_log})

        # --- Valve sizing calculation with fluids.control_valve.size_control_valve_g ---
        try:
            # Convert mass flow to normal m³/hr for size_control_valve_g
            rho_norm = (P_norm * local_gas_mw) / (R_univ * T_norm)  # Z=1 at normal conditions
            if rho_norm <= 0:
                raise ValueError("Calculated normal density is non-positive.")
            local_Q_norm_m3hr = (local_flow_rate_kg_s / rho_norm) * 3600.0
            results_log.append(f"Converted mass flow {local_flow_rate_kg_s:.4f} kg/s to normal flow {local_Q_norm_m3hr:.2f} m³/hr")

            calculated_Kv = fluids.control_valve.size_control_valve_g(
                Q=local_Q_norm_m3hr / 3600.0,  # Convert m³/hr to m³/s for fluids library
                T=local_T1_k,
                MW=local_gas_mw,  # MW already in kg/kmol (fluids library expects kg/kmol)
                Z=local_gas_z_factor,
                gamma=local_gas_gamma,
                P1=local_P1,
                P2=local_P2,
                mu=local_gas_viscosity,
                xT=local_xt,  # Make sure to use 'xT' not 'xt' for fluids library
                FL=local_FL
            )
            calculated_Cv = calculated_Kv * 1.156  # Convert Kv to Cv (Cv = Kv × 1.156)
        except Exception as calc_e:
            error_log.append(f"Error during gas valve sizing calculation: {calc_e}")
            logger.error(f"Calculation error in size_control_valve_g: {calc_e}", exc_info=True)
            return json.dumps({"errors": error_log, "log": results_log})

        # --- Estimate Valve Size (rough) ---
        cv_factor = 10.0 if valve_type in ["globe", "gate", "other"] else 15.0
        estimated_nps_inch = math.sqrt(calculated_Cv / cv_factor) if calculated_Cv > 0 else 0
        standard_sizes_inch = [0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 24]
        selected_size_inch = next((sz for sz in standard_sizes_inch if sz >= estimated_nps_inch), standard_sizes_inch[-1])

        size_value = None
        size_description = "N/A"
        try:
            if size_units == "inch":
                size_value = selected_size_inch
                size_description = f"{selected_size_inch} inch"
            elif size_units == "mm":
                size_value = selected_size_inch * 25.4
                size_description = f"{round(size_value)} mm (NPS {selected_size_inch})"
            else:  # meters
                size_value = selected_size_inch * INCH_to_M
                size_description = f"{round(size_value, 4)} m (NPS {selected_size_inch})"
            results_log.append(
                f"Estimated required NPS based on Cv: ~{estimated_nps_inch:.2f}\". "
                f"Suggesting standard size: {size_description}"
            )
        except Exception as size_fmt_e:
            error_log.append(f"Could not format estimated size: {size_fmt_e}")

        # --- Consolidate Final Results ---
        final_result = {
            "inputs_resolved": results_log,
            "warnings": error_log if error_log else None,
            "gas_property_details": gas_prop_info,
            "flow_rate_kg_s": round(local_flow_rate_kg_s, 5),
            "inlet_pressure_pa": round(local_P1, 2) if local_P1 is not None else None,
            "outlet_pressure_pa": round(local_P2, 2) if local_P2 is not None else None,
            "inlet_temperature_c": round(inlet_temperature_c, 2),
            "pressure_drop_pa": round(delta_P, 2) if delta_P is not None else None,
            # Valve Factors
            "valve_type": valve_type,
            "xt_used": round(local_xt, 3),
            "FL_used": round(local_FL, 3),
            "Fk_used": round(local_Fk, 4) if local_Fk is not None else None,
            # Sizing Results
            "required_cv_us": round(calculated_Cv, 2),
            "required_kv_metric": round(calculated_Kv, 2),
            "estimated_valve_size": size_description,
            "estimated_size_value": size_value,
            "estimated_size_units": size_units,
            "note_on_size": (
                "Estimated size is approximate. Select valve from manufacturer data "
                "using the calculated Cv/Kv."
            )
        }

        # Remove null/empty fields
        if not final_result.get("warnings"):
            final_result.pop("warnings", None)
        if final_result.get("pressure_drop_pa") is None:
            final_result.pop("pressure_drop_pa", None)

        return json.dumps(final_result)

    except Exception as e:
        logger.error(f"Error in calculate_gas_control_valve: {e}", exc_info=True)
        return json.dumps({
            "error": f"Calculation error: {str(e)}",
            "log": results_log,
            "errors_occurred": error_log
        })
