"""
Import helpers for optional dependencies.

This module provides functions to gracefully handle optional dependencies
and provide fallback values when packages are not available.
"""

import logging

logger = logging.getLogger("fluids-mcp.imports")

# FluidProperties availability check
FLUIDPROP_AVAILABLE = False
FLUID_SELECTION = None
FluidProperties = None

try:
    from fluidprop import FluidProperties as _FluidProperties, FLUID_SELECTION as _FLUID_SELECTION
    FLUIDPROP_AVAILABLE = True
    FLUID_SELECTION = _FLUID_SELECTION
    FluidProperties = _FluidProperties
    logger.info("FluidProperties package successfully imported")
except ImportError:
    logger.info("fluidprop module not available. CoolProp will be used as primary property backend.")

# CoolProp availability check
COOLPROP_AVAILABLE = False
CP = None

try:
    import CoolProp.CoolProp as _CP
    COOLPROP_AVAILABLE = True
    CP = _CP
    logger.info("CoolProp package successfully imported")
except ImportError:
    logger.warning("CoolProp module not available. Vapor pressure lookups will be disabled.")

def get_coolprop_fluids_list():
    """Get the complete list of fluids available in CoolProp.
    
    Returns:
        List of fluid names or empty list if CoolProp not available
    """
    if not COOLPROP_AVAILABLE:
        return []
        
    try:
        # Get the full list of fluids from CoolProp
        fluids_string = CP.get_global_param_string("FluidsList")
        return fluids_string.split(',')
    except Exception as e:
        logger.error("Error getting fluid list from CoolProp: %s", e)
        return []

