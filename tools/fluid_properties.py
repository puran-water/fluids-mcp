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
    if not FLUIDPROP_AVAILABLE:
        return json.dumps({
            "error": "Fluid property lookup is not available. The fluidprop package is not installed."
        })
        
    try:        
        # First try getting complete list from CoolProp
        valid_fluids = get_coolprop_fluids_list()
        
        # If CoolProp list is not available, fall back to FLUID_SELECTION
        if not valid_fluids and FLUID_SELECTION:
            valid_fluids = [f[0] for f in FLUID_SELECTION]
            
        if not valid_fluids:
            return json.dumps({
                "error": "Could not retrieve fluid list from either CoolProp or FluidProp"
            })
            
        if fluid_name not in valid_fluids:
            # Try case-insensitive match
            fluid_lower = fluid_name.lower()
            match = next((f for f in valid_fluids if f.lower() == fluid_lower), None)
            if match:
                fluid_name = match
            else:
                # Return available fluids if not found
                return json.dumps({
                    "error": f"Fluid '{fluid_name}' not found",
                    "available_fluids": valid_fluids[:10],  # Show first 10 fluids
                    "note": "Use 'list_available_fluids' tool for a complete list of valid fluids."
                })
        
        # Get fluid properties
        fluid = FluidProperties(
            coolprop_name=fluid_name,
            T_in_deg_C=temperature_c,
            P_in_bar=pressure_bar
        )
        
        # Format the output
        result = {
            "fluid_name": fluid.coolprop_name,
            "formula": fluid.formula,
            "temperature_c": temperature_c,
            "pressure_bar": pressure_bar,
            # Physical properties
            "density_kg_m3": float(fluid.rho[0]),
            "dynamic_viscosity_pa_s": float(fluid.eta[0]),
            "kinematic_viscosity_m2_s": float(fluid.nu[0]),
            # Thermal properties
            "thermal_conductivity_w_m_k": float(fluid.lambda_[0]),
            "thermal_expansion_coefficient_1_k": float(fluid.alpha[0]),
            "thermal_diffusivity_m2_s": float(fluid.kappa[0]),
            "specific_heat_cp_j_kg_k": float(fluid.Cp[0]),
            "specific_heat_cv_j_kg_k": float(fluid.Cv[0]),
            # Other properties
            "isothermal_compressibility_1_pa": float(fluid.comp[0]),
            "prandtl_number": float(fluid.Pr[0]),
            "molecular_weight_kg_mol": float(fluid.MW)
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
