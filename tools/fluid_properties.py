"""
Fluid properties tools.

This module provides tools to retrieve fluid properties and list available fluids.
"""

import json
import logging
from utils.import_helpers import (
    FLUIDPROP_AVAILABLE, FluidProperties, FLUID_SELECTION, 
    COOLPROP_AVAILABLE, CP, get_coolprop_fluids_list
)
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
    if not FLUIDPROP_AVAILABLE or FLUID_SELECTION is None:
        return json.dumps({
            "error": "Fluid property lookup is not available. The fluidprop package is not installed or not properly initialized."
        })
        
    try:
        # First map the fluid name through aliasing system
        mapped_fluid_name = map_fluid_name(fluid_name)
        
        # Then try getting complete list from CoolProp
        valid_fluids = get_coolprop_fluids_list()
        
        # If CoolProp list is not available, fall back to FLUID_SELECTION
        if not valid_fluids and FLUID_SELECTION:
            valid_fluids = [f[0] for f in FLUID_SELECTION]
            
        if not valid_fluids:
            return json.dumps({
                "error": "Could not retrieve fluid list from either CoolProp or FluidProp"
            })
            
        # Check if INCOMP:: prefix fluids (glycols) should bypass validation
        if mapped_fluid_name.startswith('INCOMP::'):
            # Skip validation for incompressible fluids
            pass
        elif mapped_fluid_name not in valid_fluids:
            # Try case-insensitive match
            fluid_lower = mapped_fluid_name.lower()
            match = next((f for f in valid_fluids if f.lower() == fluid_lower), None)
            if match:
                mapped_fluid_name = match
            else:
                # Return available fluids if not found
                return json.dumps({
                    "error": f"Fluid '{fluid_name}' (mapped to '{mapped_fluid_name}') not found",
                    "available_fluids": valid_fluids[:10],  # Show first 10 fluids
                    "note": "Use 'list_available_fluids' tool for a complete list of valid fluids."
                })
        
        # Get fluid properties
        fluid = FluidProperties(
            coolprop_name=mapped_fluid_name,
            T_in_deg_C=temperature_c,
            P_in_bar=pressure_bar
        )
        
        # Extract values
        def _to_float(x):
            try:
                return float(x)
            except Exception:
                return None

        def _is_bad(x):
            try:
                import math
                return x is None or (isinstance(x, float) and math.isnan(x))
            except Exception:
                return False

        rho = _to_float(fluid.rho[0])
        eta = _to_float(fluid.eta[0])
        nu = _to_float(fluid.nu[0])
        k = _to_float(fluid.lambda_[0])
        alpha = _to_float(fluid.alpha[0])
        kappa = _to_float(fluid.kappa[0])
        cp = _to_float(fluid.Cp[0])
        cv = _to_float(fluid.Cv[0])
        comp = _to_float(fluid.comp[0])
        pr = _to_float(fluid.Pr[0])
        mw = _to_float(fluid.MW)

        # Fallback for NaN/missing values using thermo.Chemical (better SO2 support)
        need_fallback = any(_is_bad(v) for v in [eta, cp, cv, mw])
        if need_fallback:
            try:
                from thermo.chemical import Chemical
                # Attempt with multiple identifiers if needed
                cand = [mapped_fluid_name]
                try:
                    import re as _re
                    spaced = _re.sub(r'(?<!^)(?=[A-Z])', ' ', mapped_fluid_name)
                    cand.extend([spaced, spaced.lower(), mapped_fluid_name.lower()])
                except Exception:
                    pass
                # Common formula fallbacks
                formula_map = {
                    'sulfurdioxide': 'SO2', 'sulfur dioxide': 'SO2', 'sulphur dioxide': 'SO2',
                    'carbondioxide': 'CO2', 'carbon dioxide': 'CO2',
                    'hydrogensulfide': 'H2S', 'hydrogen sulfide': 'H2S',
                }
                cand.extend([formula_map.get(c, c) for c in list(cand)])

                chem = None
                last_err = None
                for ident in cand:
                    try:
                        chem = Chemical(ident, T=temperature_c + 273.15, P=pressure_bar*1e5)
                        if chem and chem.P is not None:
                            break
                    except Exception as ee:
                        chem = None
                        last_err = ee

                if chem is not None:
                    # MW g/mol (== kg/kmol numerically)
                    if _is_bad(mw) and getattr(chem, 'MW', None):
                        mw = float(chem.MW)
                    # eta
                    if _is_bad(eta) and getattr(chem, 'mug', None):
                        eta = float(chem.mug)
                    # cp, cv
                    if _is_bad(cp) and getattr(chem, 'Cpg', None):
                        cp = float(chem.Cpg)
                    if _is_bad(cv) and getattr(chem, 'Cvg', None):
                        cv = float(chem.Cvg)
                    # If kinematic viscosity missing but rho and mu known
                    if _is_bad(nu) and rho and eta:
                        try:
                            nu = eta / rho
                        except Exception:
                            pass
                else:
                    logger.warning(f"thermo.Chemical fallback failed for {mapped_fluid_name}: {last_err}")
            except Exception as fb_e:
                logger.warning(f"Fallback to thermo.Chemical failed: {fb_e}")

        # Format the output
        result = {
            "fluid_name": fluid.coolprop_name,
            "formula": fluid.formula,
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
            "molecular_weight_kg_mol": mw
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
    # Method 1: Get fluids directly from CoolProp
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
            # Fall through to next method if this fails
    
    # Method 2: Fall back to FLUID_SELECTION from FluidProp
    if FLUIDPROP_AVAILABLE and FLUID_SELECTION:
        try:
            # Format the fluids list with indices
            fluids_dict = {i: name for i, (name, _) in enumerate(FLUID_SELECTION)}
            
            result = {
                "available_fluids": fluids_dict,
                "total_count": len(fluids_dict),
                "source": "FluidProp FLUID_SELECTION",
                "usage_example": "Use the fluid name (e.g., 'Water') or its index when calling get_fluid_properties",
                "note": "This is a limited list. Full CoolProp fluid list may be available with an update."
            }
            
            return json.dumps(result)
        except Exception as e:
            logger.error(f"Error using FLUID_SELECTION: {e}", exc_info=True)
    
    # If both methods fail
    return json.dumps({
        "error": "Fluid property lookup is not available. Neither CoolProp nor FluidProp is properly configured."
    })
