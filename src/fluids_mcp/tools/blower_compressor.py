import json
import logging
import math
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

# Import thermo library for advanced EOS calculations
try:
    from thermo import Chemical
    from scipy.integrate import quad
    import numpy as np
    THERMO_AVAILABLE = True
except ImportError:
    Chemical = None
    quad = None
    np = None
    THERMO_AVAILABLE = False

# Import shared utilities
from fluids_mcp.errors import CalculationError
from fluids_mcp.utils.constants import (
    GPM_to_M3S, INCH_to_M, FT_to_M, PSI_to_PA,
    LBFT3_to_KGM3, CENTIPOISE_to_PAS, DEFAULT_ROUGHNESS, DEG_C_to_K
)
from fluids_mcp.utils.helpers import get_fitting_K

# Configure logging
logger = logging.getLogger("fluids-mcp.calculate_blower_compressor_requirements")


def _structured_error(
    code: str,
    message: str,
    log_messages: list,
    details: Optional[dict] = None,
) -> str:
    error = CalculationError(code, message, details=details)
    return json.dumps({
        "error_code": error.code,
        "error": error.message,
        "errors": [error.to_public_dict()],
        "log": log_messages,
    }, allow_nan=False)


def _finite_number(value) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _component_mw(component: str) -> float:
    mw_map = {
        "CH4": 16.04,
        "METHANE": 16.04,
        "CO2": 44.01,
        "CARBONDIOXIDE": 44.01,
        "CARBON DIOXIDE": 44.01,
        "N2": 28.01,
        "NITROGEN": 28.01,
        "H2S": 34.08,
        "O2": 32.00,
        "OXYGEN": 32.00,
        "H2": 2.016,
        "HYDROGEN": 2.016,
    }
    key = component.strip().upper()
    return mw_map.get(key, 28.96)


def _composition_mw(gas_composition_mol: Dict[str, float]) -> float:
    total_fraction = 0.0
    weighted_mw = 0.0
    for component, fraction in gas_composition_mol.items():
        if not _finite_number(fraction) or fraction < 0:
            raise ValueError(f"Invalid mole fraction for {component}: {fraction!r}")
        total_fraction += float(fraction)
        weighted_mw += _component_mw(str(component)) * float(fraction)
    if total_fraction <= 0:
        raise ValueError("Composition mole fractions must sum to a positive value.")
    return weighted_mw / total_fraction


def calculate_integral_Z(chemical, P1, P2, T1):
    """Calculate pressure-averaged Z factor using Peng-Robinson EOS"""
    if not THERMO_AVAILABLE or chemical is None:
        return 1.0  # Fallback to ideal gas

    try:
        def Z_at_pressure(P):
            temp_chemical = Chemical(chemical.ID, T=T1, P=P)
            return temp_chemical.Z

        Z_avg, _ = quad(Z_at_pressure, P1, P2)
        Z_avg /= (P2 - P1)
        return Z_avg
    except Exception:
        return 1.0  # Fallback to ideal gas

def calculate_integral_Cp(chemical, P1, P2, T1):
    """Calculate pressure-averaged Cp using EOS"""
    if not THERMO_AVAILABLE or chemical is None:
        # Fallback to ideal gas Cp from gamma
        return None

    try:
        def Cp_at_pressure(P):
            temp_chemical = Chemical(chemical.ID, T=T1, P=P)
            return temp_chemical.Cp

        Cp_avg, _ = quad(Cp_at_pressure, P1, P2)
        Cp_avg /= (P2 - P1)
        return Cp_avg
    except Exception:
        return None

def get_chemical_properties(fluid_name, T, P):
    """Get thermodynamic properties using thermo Chemical class"""
    if not THERMO_AVAILABLE:
        return None

    # Map common fluid names to thermo chemical names
    chemical_name_map = {
        'Air': 'nitrogen',  # Use nitrogen as approximation for air
        'air': 'nitrogen',
        'Methane': 'methane',
        'methane': 'methane',
        'CH4': 'methane',
        'NaturalGas': 'methane',  # Approximate with methane
        'Biogas': 'methane',  # Approximate with methane
        'CarbonDioxide': 'carbon dioxide',
        'CO2': 'carbon dioxide',
        'Nitrogen': 'nitrogen',
        'N2': 'nitrogen',
        'Oxygen': 'oxygen',
        'O2': 'oxygen',
        'Hydrogen': 'hydrogen',
        'H2': 'hydrogen',
        'Water': 'water',
        'H2O': 'water'
    }

    try:
        # Map the fluid name to thermo chemical name
        thermo_name = chemical_name_map.get(fluid_name, fluid_name.lower())
        chemical = Chemical(thermo_name, T=T, P=P)

        # Calculate gamma more carefully
        gamma_val = None
        cv_val = None

        try:
            if hasattr(chemical, 'Cvg'):  # Gas-phase Cv
                cv_val = chemical.Cvg
            elif hasattr(chemical, 'Cv'):
                cv_val = chemical.Cv

            if cv_val is not None and cv_val > 0 and hasattr(chemical, 'Cp'):
                gamma_val = chemical.Cp / cv_val
        except Exception as e:
            logger.debug("Gamma calculation failed: %s", e)
            gamma_val = None

        return {
            'MW': chemical.MW,  # Already in g/mol which numerically equals kg/kmol
            'gamma': gamma_val,
            'Z': chemical.Z,
            'Cp': chemical.Cp if hasattr(chemical, 'Cp') else None,
            'Cv': cv_val,
            'chemical_obj': chemical
        }
    except Exception as e:
        logger.debug(f"get_chemical_properties failed for {fluid_name}: {e}")
        return None

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
    inlet_pressure_barg: Optional[float] = None,   # Inlet gauge pressure in bar(g)
    outlet_pressure_barg: Optional[float] = None,  # Required outlet gauge pressure in bar(g)
    atmospheric_pressure: Optional[float] = None, # Absolute reference pressure for bar(g) conversion, Pa

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
    compressor_type: Literal['turbo', 'roots_pd'] = 'turbo',

    # --- Polytropic Calculation Specific (Optional) ---
    polytropic_n: Optional[float] = None,             # Polytropic exponent (if efficiency_type is 'polytropic' and efficiency is based on 'n')

    # --- Default Handling ---
    allow_property_defaults: bool = False,            # Allow using default gas properties (MW=28.96, gamma=1.4, Z=1.0) which can cause >10% error for non-air gases

    # --- Advanced EOS Options ---
    use_eos_calculations: bool = True,                # Use equation of state for accurate Z-factor and Cp calculations (requires thermo library)

) -> str:
    """Estimate discharge temperature and shaft power for a compressor/blower.

    ``compressor_type='turbo'`` uses isentropic or polytropic compression.
    ``compressor_type='roots_pd'`` uses the positive-displacement inlet-volume
    basis ``shaft power = Q1 * (P2-P1) / efficiency``. The latter follows the
    inlet actual-volume and shaft-power reporting basis of ISO 1217 for Roots/
    lobe blower duty; vendor guaranteed specific power remains the final
    selection basis.

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
        inlet_pressure_barg: Optional - Inlet pressure in bar(g)
        outlet_pressure_barg: Optional - Outlet pressure in bar(g)
        atmospheric_pressure: Absolute local atmospheric pressure in Pa used to convert bar(g);
            omitted values resolve to 101325 Pa with a warning
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
        compressor_type: ``turbo`` or external-compression ``roots_pd``.
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
        if compressor_type not in ("turbo", "roots_pd"):
            return _structured_error(
                "INVALID_COMPRESSOR_TYPE",
                "compressor_type must be one of: turbo, roots_pd.",
                results_log,
                {"compressor_type": compressor_type},
            )

        # 1. Resolve Inlet Temperature (K)
        if inlet_temperature_c is None:
            error_log.append("Missing required input: inlet_temperature_c.")
            # Cannot proceed without temperature
            return json.dumps({"errors": error_log, "log": results_log})
        local_T1_k = inlet_temperature_c + DEG_C_to_K # Convert to Kelvin
        results_log.append(f"Using Inlet Temperature: {inlet_temperature_c:.2f} C ({local_T1_k:.2f} K)")

        # 2. Resolve Pressures (Pa Absolute)
        pressure_inputs = {
            "inlet": (inlet_pressure, inlet_pressure_psi, inlet_pressure_barg),
            "outlet": (outlet_pressure, outlet_pressure_psi, outlet_pressure_barg),
        }
        for pressure_name, candidates in pressure_inputs.items():
            if sum(value is not None for value in candidates) > 1:
                return _structured_error(
                    "AMBIGUOUS_PRESSURE_INPUT",
                    f"Provide only one {pressure_name} pressure basis: absolute Pa, absolute psi, or bar(g).",
                    results_log,
                    {"pressure": pressure_name},
                )

        local_atmospheric_pressure = 101325.0 if atmospheric_pressure is None else atmospheric_pressure
        if not _finite_number(local_atmospheric_pressure) or local_atmospheric_pressure <= 0:
            return _structured_error(
                "INVALID_PRESSURE_INPUT",
                "atmospheric_pressure must be a positive finite absolute pressure in Pa.",
                results_log,
                {"atmospheric_pressure": atmospheric_pressure},
            )
        if atmospheric_pressure is None and (
            inlet_pressure_barg is not None or outlet_pressure_barg is not None
        ):
            error_log.append(
                "Warning: atmospheric_pressure was omitted; used 101325 Pa to convert bar(g). "
                "Supply site atmospheric pressure for site-specific blower duties."
            )

        pressure_basis = {}
        local_P1 = None
        if inlet_pressure is not None:
            local_P1 = inlet_pressure
            results_log.append("Used provided SI inlet_pressure.")
            pressure_basis["inlet"] = "absolute_pa"
        elif inlet_pressure_psi is not None:
            local_P1 = inlet_pressure_psi * PSI_to_PA
            results_log.append(f"Converted inlet_pressure from {inlet_pressure_psi} psi.")
            pressure_basis["inlet"] = "absolute_psi"
        elif inlet_pressure_barg is not None:
            local_P1 = local_atmospheric_pressure + inlet_pressure_barg * 100000.0
            results_log.append(
                f"Converted inlet_pressure from {inlet_pressure_barg} bar(g) using "
                f"atmospheric_pressure={local_atmospheric_pressure:.2f} Pa."
            )
            pressure_basis["inlet"] = "gauge_bar"
        else:
            error_log.append(
                "Missing required input: inlet_pressure, inlet_pressure_psi, or inlet_pressure_barg."
            )

        local_P2 = None
        if outlet_pressure is not None:
            local_P2 = outlet_pressure
            results_log.append("Used provided SI outlet_pressure.")
            pressure_basis["outlet"] = "absolute_pa"
        elif outlet_pressure_psi is not None:
            local_P2 = outlet_pressure_psi * PSI_to_PA
            results_log.append(f"Converted outlet_pressure from {outlet_pressure_psi} psi.")
            pressure_basis["outlet"] = "absolute_psi"
        elif outlet_pressure_barg is not None:
            local_P2 = local_atmospheric_pressure + outlet_pressure_barg * 100000.0
            results_log.append(
                f"Converted outlet_pressure from {outlet_pressure_barg} bar(g) using "
                f"atmospheric_pressure={local_atmospheric_pressure:.2f} Pa."
            )
            pressure_basis["outlet"] = "gauge_bar"
        else:
            error_log.append(
                "Missing required input: outlet_pressure, outlet_pressure_psi, or outlet_pressure_barg."
            )

        for pressure_name, pressure_value in (("inlet", local_P1), ("outlet", local_P2)):
            if pressure_value is None:
                continue
            if not _finite_number(pressure_value) or pressure_value <= 0:
                return _structured_error(
                    "INVALID_PRESSURE_INPUT",
                    f"Resolved {pressure_name} absolute pressure must be positive and finite.",
                    results_log,
                    {f"{pressure_name}_pressure_pa": pressure_value},
                )
            if pressure_basis[pressure_name].startswith("absolute") and pressure_value < 0.5 * 101325.0:
                return _structured_error(
                    "AMBIGUOUS_PRESSURE_BASIS",
                    f"{pressure_name}_pressure resolved below 0.5 atm absolute and looks like a gauge "
                    f"value. Absolute pressure fields require Pa/psi absolute; use {pressure_name}_pressure_barg "
                    "for blower pressures quoted in bar(g).",
                    results_log,
                    {
                        f"{pressure_name}_pressure_pa": pressure_value,
                        "minimum_plausibility_threshold_pa": 0.5 * 101325.0,
                    },
                )

        # Check P2 > P1 after resolution attempt
        if local_P1 is not None and local_P2 is not None and local_P2 <= local_P1:
            error_log.append("Outlet pressure must be greater than inlet pressure for compression.")
            # Return now if pressures are invalid
            if error_log: return json.dumps({"errors": error_log, "log": results_log})

        # 3. Resolve Gas Properties (MW, gamma, Z) - Enhanced with EOS calculations
        local_gas_mw = gas_mw
        local_gas_gamma = gas_gamma
        local_gas_z_factor = gas_z_factor
        local_gas_viscosity = None
        local_chemical = None  # Store Chemical object for EOS calculations
        integral_z_factor = None
        integral_cp = None
        gas_prop_source = "Provided"
        prop_resolved = all(p is not None for p in [local_gas_mw, local_gas_gamma])

        # Try EOS-based calculations first if enabled and thermo available
        if use_eos_calculations and THERMO_AVAILABLE and fluid_name is not None:
            results_log.append(f"Attempting EOS-based property calculation for '{fluid_name}'.")
            try:
                # Use average pressure for initial property lookup
                avg_p_pa = 101325.0
                if local_P1 is not None and local_P2 is not None:
                    avg_p_pa = (local_P1 + local_P2) / 2.0
                elif local_P1 is not None:
                    avg_p_pa = local_P1

                eos_props = get_chemical_properties(fluid_name, local_T1_k, avg_p_pa)
                if eos_props:
                    local_chemical = eos_props['chemical_obj']
                    if local_gas_mw is None:
                        local_gas_mw = eos_props['MW']
                    if local_gas_gamma is None and eos_props['gamma'] is not None:
                        local_gas_gamma = eos_props['gamma']
                    if local_gas_z_factor is None:
                        local_gas_z_factor = eos_props['Z']

                    gas_prop_source = f"EOS ({fluid_name} @ {inlet_temperature_c}C, {avg_p_pa/100000:.1f}bar)"
                    results_log.append(f"Successfully obtained properties using EOS calculations.")
                    prop_resolved = True
                else:
                    results_log.append("EOS calculation failed, falling back to FluidProp lookup.")
            except Exception as eos_e:
                results_log.append(f"EOS calculation failed ({eos_e}), falling back to FluidProp lookup.")

        # If EOS failed or disabled, try centralized property resolver (CoolProp → fluidprop → thermo)
        if any(p is None for p in [local_gas_mw, local_gas_gamma]) and fluid_name is not None and not prop_resolved:
            gas_prop_source = f"Lookup ({fluid_name} @ {inlet_temperature_c} C)"
            results_log.append(f"Attempting property lookup for '{fluid_name}'.")
            try:
                from fluids_mcp.utils.resolve_properties import resolve_gas_properties
                avg_p_pa = 101325.0
                if local_P1 is not None and local_P2 is not None:
                    avg_p_pa = (local_P1 + local_P2) / 2.0
                elif local_P1 is not None:
                    avg_p_pa = local_P1
                lookup_p_bar = avg_p_pa / 100000.0

                gas_props = resolve_gas_properties(fluid_name, inlet_temperature_c, lookup_p_bar)
                if gas_props is not None:
                    if local_gas_mw is None and gas_props.mw is not None:
                        local_gas_mw = gas_props.mw
                    if local_gas_gamma is None and gas_props.gamma is not None:
                        local_gas_gamma = gas_props.gamma
                    if local_gas_z_factor is None and gas_props.z_factor is not None:
                        local_gas_z_factor = gas_props.z_factor
                    gas_prop_info = {"source": gas_props.source, "lookup_temp_c": inlet_temperature_c, "lookup_p_bar": lookup_p_bar}
                    gas_prop_source = f"Resolver ({gas_props.source})"
                    prop_resolved = True
                else:
                    error_log.append(f"Warning: Property lookup failed for '{fluid_name}'. Falling back to defaults.")
                    gas_prop_source = "Lookup Failed - Using Defaults"
                    prop_resolved = False
            except Exception as prop_lookup_e:
                error_log.append(f"Warning: Property lookup failed for '{fluid_name}': {prop_lookup_e}. Falling back to defaults.")
                gas_prop_source = "Lookup Failed - Using Defaults"
                prop_resolved = False

        if any(p is None for p in [local_gas_mw, local_gas_gamma, local_gas_z_factor]) and gas_composition_mol is not None:
            gas_prop_source = "Composition Calculation"
            results_log.append("Attempting property calculation from composition.")
            try:
                if local_gas_mw is None:
                    local_gas_mw = _composition_mw(gas_composition_mol)
                if local_gas_gamma is None:
                    local_gas_gamma = 1.3
                if local_gas_z_factor is None:
                    local_gas_z_factor = 0.98
                gas_prop_info = {"composition_mol": gas_composition_mol}
                error_log.append("Warning: Using estimated gamma & Z for mixture based on composition.")
                prop_resolved = local_gas_mw is not None and local_gas_gamma is not None
            except Exception as comp_e:
                error_log.append(f"Error during composition property calculation: {comp_e}")
                prop_resolved = False

        if local_gas_mw is None:
            if allow_property_defaults:
                local_gas_mw = 28.96; results_log.append("WARNING: Using default MW (Air) - may cause >10% error for non-air gases."); gas_prop_source += " + MW Default"
            else:
                error_log.append("Missing gas molecular weight. Provide gas_mw, fluid_name, gas_composition_mol, or set allow_property_defaults=True (may cause >10% error for non-air gases).")
        if local_gas_gamma is None:
            if allow_property_defaults:
                local_gas_gamma = 1.4; results_log.append("WARNING: Using default gamma (Air) - may cause >10% error for non-air gases."); gas_prop_source += " + Gamma Default"
            else:
                error_log.append("Missing gas specific heat ratio. Provide gas_gamma, fluid_name, gas_composition_mol, or set allow_property_defaults=True (may cause >10% error for non-air gases).")

        if local_gas_z_factor is None:
            if allow_property_defaults or prop_resolved:
                local_gas_z_factor = 1.0
                if prop_resolved:
                    results_log.append("INFO: Z-factor not available from property lookup, using ideal gas assumption (Z=1.0)")
                    gas_prop_source += " + Z Ideal"
                else:
                    results_log.append("WARNING: Using default Z-factor (Ideal) - may cause >5% error at high pressures.")
                    gas_prop_source += " + Z Default"
            else:
                error_log.append("Missing gas compressibility factor. Provide gas_z_factor, fluid_name, gas_composition_mol, or set allow_property_defaults=True (may cause errors at high pressures).")

        gas_prop_info.update({
            "MW_kg_kmol": local_gas_mw, "gamma": local_gas_gamma,
            "Z_factor_used": local_gas_z_factor, # Clarify Z used (might be inlet Z or average estimate)
            "viscosity": local_gas_viscosity, # Placeholder for viscosity if needed
            "source": gas_prop_source
        })
        results_log.append(f"Resolved Gas Properties: MW={local_gas_mw or 'None'}, gamma={local_gas_gamma or 'None'}, Z={local_gas_z_factor or 'None'} (Source: {gas_prop_source})")

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
            if local_flow_rate_kg_s <= 0:
                return _structured_error(
                    "INVALID_FLOW_RATE",
                    "Compressor/blower flow rate must be positive.",
                    results_log,
                    {
                        "flow_rate_kg_s": local_flow_rate_kg_s,
                        "flow_rate_source": flow_rate_source,
                    },
                )
        else:
            # If flow rate resolution failed (e.g., needed MW but MW failed) add specific error
            if "Missing required input: flow_rate" not in str(error_log): # Avoid duplicate
                error_log.append("Could not resolve mass flow rate (e.g., missing MW for conversion).")

        # 5. Resolve Efficiency
        local_efficiency = efficiency
        local_polytropic_n = polytropic_n
        if local_efficiency is None:
            if compressor_type == 'roots_pd':
                local_efficiency = 0.70
            elif efficiency_type == 'isentropic':
                local_efficiency = 0.70 # Default isentropic efficiency
            else: # polytropic
                local_efficiency = 0.75 # Default polytropic efficiency
            efficiency_label = "overall shaft" if compressor_type == 'roots_pd' else efficiency_type
            results_log.append(f"Using default {efficiency_label} efficiency: {local_efficiency:.2f}")
        else:
            efficiency_label = "overall shaft" if compressor_type == 'roots_pd' else efficiency_type
            results_log.append(f"Using provided {efficiency_label} efficiency: {local_efficiency:.2f}")

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

        # Calculate integral-averaged Z and Cp if EOS available
        if use_eos_calculations and local_chemical is not None and local_P1 is not None and local_P2 is not None:
            try:
                integral_z_factor = calculate_integral_Z(local_chemical, local_P1, local_P2, local_T1_k)
                integral_cp = calculate_integral_Cp(local_chemical, local_P1, local_P2, local_T1_k)
                results_log.append(f"Calculated integral-averaged Z = {integral_z_factor:.4f}")
                if integral_cp is not None:
                    results_log.append(f"Calculated integral-averaged Cp = {integral_cp:.1f} J/(mol·K)")
            except Exception as integral_e:
                results_log.append(f"Integral averaging failed ({integral_e}), using point properties.")

        # Use integral Z if available, otherwise use point Z
        Z_avg = integral_z_factor if integral_z_factor is not None else local_gas_z_factor
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
        inlet_actual_flow_m3_s = None
        power_method = None

        if compressor_type == 'roots_pd':
            rho_inlet = local_P1 * local_gas_mw / (
                Z_avg * R_univ * local_T1_k
            )
            if rho_inlet <= 0:
                return _structured_error(
                    "INVALID_GAS_DENSITY",
                    "Resolved inlet gas density must be positive for Roots/PD power calculation.",
                    results_log,
                    {"rho_inlet_kg_m3": rho_inlet},
                )
            inlet_actual_flow_m3_s = local_flow_rate_kg_s / rho_inlet
            ideal_pd_power_w = inlet_actual_flow_m3_s * (local_P2 - local_P1)
            W_actual = ideal_pd_power_w / (local_efficiency * local_flow_rate_kg_s)
            cp_ideal = k * R_specific / (k - 1.0)
            T2_actual_k = local_T1_k + W_actual / cp_ideal
            power_method = "ISO 1217 Roots/PD inlet-volume basis: Q1*dP/overall_efficiency"
            results_log.append(
                "Calculated Roots/PD shaft power from actual inlet volume flow and differential pressure "
                "per the ISO 1217 reporting basis."
            )
        elif efficiency_type == 'isentropic':
            if local_efficiency == 0: error_log.append("Isentropic efficiency cannot be zero."); return json.dumps({"errors": error_log, "log": results_log})
            W_actual = W_s / local_efficiency
            # W_actual = Cp * (T2_actual - T1) = (k*R_spec / (k-1)) * (T2_actual - T1)
            # T2_actual = T1 + W_actual * (k-1) / (k*R_spec)
            # Check k*R_specific is not zero
            if abs(k * R_specific) < 1e-12:
                error_log.append("Invalid gas properties (k or R_specific is zero) for temperature calculation.")
                return json.dumps({"errors": error_log, "log": results_log})
            T2_actual_k = local_T1_k + (W_actual * (k - 1.0)) / (k * R_specific)
            power_method = "isentropic compression divided by isentropic efficiency"
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
            power_method = "polytropic compression path"
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
            "atmospheric_pressure_pa": round(local_atmospheric_pressure, 2),
            "pressure_input_basis": pressure_basis,
            "compressor_type": compressor_type,
            "power_method": power_method,
            "differential_pressure_pa": round(local_P2 - local_P1, 2),
            "pressure_ratio": round(P_ratio, 4),
            "inlet_temperature_c": round(inlet_temperature_c, 2),
            "inlet_temperature_k": round(local_T1_k, 2),
            "flow_rate_kg_s": round(local_flow_rate_kg_s, 4),
            # Efficiency Details
            "efficiency_type_used": "overall_shaft" if compressor_type == "roots_pd" else efficiency_type,
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
            "required_power_hp": round(Power_watts / 745.7, 3),
            "shaft_power_kw": round(Power_watts / 1000.0, 3),
            "shaft_power_hp": round(Power_watts / 745.7, 3),
        }
        if inlet_actual_flow_m3_s is not None:
            final_result["inlet_actual_flow_m3_s"] = inlet_actual_flow_m3_s

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

        return json.dumps(
            final_result,
            default=lambda x: round(x, 6) if isinstance(x, float) else x,
            allow_nan=False,
        ) # Round floats for cleaner output

    except Exception as e:
        logger.error(f"Error in calculate_blower_compressor_requirements: {e}", exc_info=True)
        return _structured_error(
            "BLOWER_COMPRESSOR_CALCULATION_FAILED",
            f"Calculation error: {str(e)}",
            results_log,
            {"errors_occurred": error_log},
        )

def blower_sweep(
    variable: str,
    start: float,
    stop: float,
    n: int,
    # Base calculation parameters from calculate_blower_compressor_requirements
    flow_rate_kg_s: Optional[float] = None,
    flow_rate_norm_m3_hr: Optional[float] = None,
    flow_rate_std_m3_hr: Optional[float] = None,
    inlet_pressure: Optional[float] = None,
    outlet_pressure: Optional[float] = None,
    inlet_pressure_psi: Optional[float] = None,
    outlet_pressure_psi: Optional[float] = None,
    inlet_temperature_c: Optional[float] = None,
    gas_mw: Optional[float] = None,
    gas_gamma: Optional[float] = None,
    gas_z_factor: Optional[float] = None,
    fluid_name: Optional[str] = None,
    gas_composition_mol: Optional[Dict[str, float]] = None,
    efficiency_type: Literal['isentropic', 'polytropic'] = 'isentropic',
    efficiency: Optional[float] = None,
    polytropic_n: Optional[float] = None,
    allow_property_defaults: bool = False,
    use_eos_calculations: bool = True,
) -> str:
    """Parameter sweep for blower compressor analysis

    Optimized for performance with large datasets (1e4+ points) per expert reviewer recommendation.
    Sweeps one variable while keeping others constant to analyze system behavior.

    Args:
        variable: Variable to sweep ('pressure_ratio', 'flow_rate_kg_s', 'efficiency', 'inlet_temperature_c')
        start: Start value
        stop: Stop value
        n: Number of points
        **other_params: All other parameters from calculate_blower_compressor_requirements

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
        'flow_rate_kg_s': flow_rate_kg_s,
        'flow_rate_norm_m3_hr': flow_rate_norm_m3_hr,
        'flow_rate_std_m3_hr': flow_rate_std_m3_hr,
        'inlet_pressure': inlet_pressure,
        'outlet_pressure': outlet_pressure,
        'inlet_pressure_psi': inlet_pressure_psi,
        'outlet_pressure_psi': outlet_pressure_psi,
        'inlet_temperature_c': inlet_temperature_c,
        'gas_mw': gas_mw,
        'gas_gamma': gas_gamma,
        'gas_z_factor': gas_z_factor,
        'fluid_name': fluid_name,
        'gas_composition_mol': gas_composition_mol,
        'efficiency_type': efficiency_type,
        'efficiency': efficiency,
        'polytropic_n': polytropic_n,
        'allow_property_defaults': allow_property_defaults,
        'use_eos_calculations': use_eos_calculations,
    }

    # Remove None values
    base_kwargs = {k: v for k, v in base_kwargs.items() if v is not None}

    for value in sweep_values:
        try:
            # Set up parameters for this sweep point
            kwargs = base_kwargs.copy()

            if variable == 'pressure_ratio':
                if 'inlet_pressure' in kwargs:
                    kwargs['outlet_pressure'] = kwargs['inlet_pressure'] * value
                else:
                    kwargs['inlet_pressure'] = 101325.0
                    kwargs['outlet_pressure'] = 101325.0 * value
            elif variable == 'flow_rate_kg_s':
                kwargs['flow_rate_kg_s'] = value
            elif variable == 'efficiency':
                kwargs['efficiency'] = value
            elif variable == 'inlet_temperature_c':
                kwargs['inlet_temperature_c'] = value
            else:
                # For any other variable, set it directly
                kwargs[variable] = value

            # Calculate
            result_json = calculate_blower_compressor_requirements(**kwargs)
            result = json.loads(result_json)

            if 'errors' not in result or not result['errors']:
                # Extract key results
                row = {
                    variable: round(value, 6),
                    'power_kw': result.get('required_power_kw', None),
                    'discharge_temp_c': result.get('actual_discharge_temp_c', None),
                    'pressure_ratio': result.get('pressure_ratio', None),
                    'efficiency_used': result.get('efficiency_used', None),
                    'z_factor': result.get('gas_property_details', {}).get('Z_factor_used', None),
                    'specific_work_kj_kg': result.get('actual_specific_work_kj_kg', None)
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
