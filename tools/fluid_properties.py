"""
Fluid properties tools.

This module provides tools to retrieve fluid properties and list available fluids.
Delegates to the centralized resolver (utils/resolve_properties.py) for core property
resolution, supplemented by CachedFluidProperties for extended thermal properties.
"""

import json
import logging
from utils.import_helpers import COOLPROP_AVAILABLE, get_coolprop_fluids_list
from utils.fluid_aliases import map_fluid_name

# Configure logging
logger = logging.getLogger("fluids-mcp.fluid_properties")

def get_fluid_properties(
    fluid_name: str,           # Name of the fluid (e.g., "Water", "Air", "Nitrogen")
    temperature_c: float,      # Temperature in degrees Celsius
    pressure_bar: float = 1.0, # Pressure in bar (default: 1.0 bar)
) -> str:
    """Retrieve thermodynamic properties of a fluid at specified temperature and pressure.

    Args:
        fluid_name: Name of the fluid from the available list: Water (121), Air (2),
                   Nitrogen (63), Methane (43), Ammonia (3), CarbonDioxide (6),
                   Ethane (21), Ethanol (22), Hydrogen (30), Oxygen (69), etc.
                   For a complete list, use "list_available_fluids" tool.
        temperature_c: Temperature in degrees Celsius
        pressure_bar: Pressure in bar (default: 1.0 bar)

    Returns:
        Detailed fluid properties including density, viscosity, thermal properties, etc.
    """
    try:
        # First map the fluid name through aliasing system
        mapped_fluid_name = map_fluid_name(fluid_name)

        # Validate fluid name against known list
        valid_fluids = get_coolprop_fluids_list()
        if valid_fluids:
            if not mapped_fluid_name.startswith('INCOMP::') and mapped_fluid_name not in valid_fluids:
                fluid_lower = mapped_fluid_name.lower()
                match = next((f for f in valid_fluids if f.lower() == fluid_lower), None)
                if match:
                    mapped_fluid_name = match
                else:
                    return json.dumps({
                        "error": f"Fluid '{fluid_name}' (mapped to '{mapped_fluid_name}') not found",
                        "available_fluids": valid_fluids[:10],
                        "note": "Use 'list_available_fluids' tool for a complete list of valid fluids."
                    })

        # Use centralized resolver for core properties (CoolProp -> fluidprop -> thermo)
        from utils.resolve_properties import resolve_liquid_properties, resolve_gas_properties

        liquid_props = resolve_liquid_properties(mapped_fluid_name, temperature_c, pressure_bar)
        gas_props = resolve_gas_properties(mapped_fluid_name, temperature_c, pressure_bar)

        rho = getattr(liquid_props, 'density', None) if liquid_props else None
        eta = getattr(liquid_props, 'viscosity', None) if liquid_props else None
        nu = getattr(liquid_props, 'kinematic_viscosity', None) if liquid_props else None
        mw = getattr(gas_props, 'mw', None) if gas_props else None
        cp = getattr(gas_props, 'cp', None) if gas_props else None
        cv = getattr(gas_props, 'cv', None) if gas_props else None

        # Extended properties via CachedFluidProperties (thermal conductivity, Prandtl, etc.)
        k = None      # thermal conductivity
        alpha = None   # thermal expansion coefficient
        kappa = None   # thermal diffusivity
        comp = None    # isothermal compressibility
        pr = None      # Prandtl number

        if COOLPROP_AVAILABLE:
            try:
                from utils.property_cache import CachedFluidProperties
                fluid = CachedFluidProperties(
                    coolprop_name=mapped_fluid_name,
                    T_in_deg_C=temperature_c,
                    P_in_bar=pressure_bar
                )
                def _to_float(x):
                    try:
                        return float(x)
                    except Exception:
                        return None

                k = _to_float(fluid.lambda_[0])
                alpha = _to_float(fluid.alpha[0])
                # Fill in any gaps from the resolver with CoolProp-direct values
                if rho is None:
                    rho = _to_float(fluid.rho[0])
                if eta is None:
                    eta = _to_float(fluid.eta[0])
                if nu is None:
                    nu = _to_float(fluid.nu[0])
                if cp is None:
                    cp = _to_float(fluid.Cp[0])
                if cv is None:
                    cv = _to_float(fluid.Cv[0])
                if mw is None:
                    mw = _to_float(fluid.MW)
            except Exception as cp_err:
                logger.warning(f"CachedFluidProperties failed for {mapped_fluid_name}: {cp_err}")

        # Format the output
        result = {
            "fluid_name": mapped_fluid_name,
            "temperature_c": temperature_c,
            "pressure_bar": pressure_bar,
            # Physical properties
            "density_kg_m3": rho,
            "dynamic_viscosity_pa_s": eta,
            "kinematic_viscosity_m2_s": nu,
            # Thermal properties
            "thermal_conductivity_w_m_k": k,
            "thermal_expansion_coefficient_1_k": alpha,
            "thermal_diffusivity_m2_s": kappa,
            "specific_heat_cp_j_kg_k": cp,
            "specific_heat_cv_j_kg_k": cv,
            # Other properties
            "isothermal_compressibility_1_pa": comp,
            "prandtl_number": pr,
            "molecular_weight_kg_kmol": mw
        }

        return json.dumps(result)
    except Exception as e:
        logger.error(f"Error in get_fluid_properties: {e}", exc_info=True)
        return json.dumps({"error": f"Error retrieving fluid properties: {str(e)}"})

def list_available_fluids() -> str:
    """List all available fluids supported by the system.

    Returns:
        JSON string containing all available fluid names and their indices.
    """
    # Get fluids directly from CoolProp
    if COOLPROP_AVAILABLE:
        try:
            # Get the complete list from CoolProp
            fluids_list = get_coolprop_fluids_list()

            if fluids_list:
                # Format the fluids list with indices
                fluids_dict = {i: name for i, name in enumerate(fluids_list)}

                result = {
                    "available_fluids": fluids_dict,
                    "total_count": len(fluids_dict),
                    "source": "CoolProp direct lookup",
                    "usage_example": "Use the fluid name (e.g., 'Water') when calling get_fluid_properties"
                }

                return json.dumps(result)
        except Exception as e:
            logger.error(f"Error getting fluid list from CoolProp: {e}", exc_info=True)

    # If CoolProp unavailable
    return json.dumps({
        "error": "Fluid property lookup is not available. CoolProp is not properly configured."
    })
