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
    logger.warning("fluidprop module not available. Fluid property lookups will be disabled.")

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

def get_fluid_properties(fluid_name, temperature_c, pressure_bar=1.01325):
    """Attempt to get fluid properties with graceful fallback to defaults.
    
    Args:
        fluid_name: Name of fluid
        temperature_c: Temperature in Celsius
        pressure_bar: Pressure in bar
        
    Returns:
        Dictionary of fluid properties or None if lookup failed
    """
    if not FLUIDPROP_AVAILABLE:
        logger.warning("Cannot lookup fluid properties: fluidprop package not available")
        return None
        
    try:
        # Use cached version for better performance
        from .property_cache import CachedFluidProperties
        fluid_props = CachedFluidProperties(
            coolprop_name=fluid_name,
            T_in_deg_C=temperature_c,
            P_in_bar=pressure_bar
        )
        
        return {
            "name": fluid_name,
            "temperature_c": temperature_c,
            "pressure_bar": pressure_bar,
            "density_kgm3": float(fluid_props.rho[0]),
            "viscosity_pas": float(fluid_props.eta[0]),
            "kinematic_viscosity_m2s": float(fluid_props.nu[0])
        }
    except Exception as e:
        logger.error("Error getting fluid properties for %s: %s", fluid_name, e)
        return None

def get_saturation_pressure(fluid_name, temperature_c):
    """Attempt to get saturation pressure with graceful fallback.
    
    Args:
        fluid_name: Name of fluid
        temperature_c: Temperature in Celsius
        
    Returns:
        Saturation pressure in Pa or None if lookup failed
    """
    if not COOLPROP_AVAILABLE:
        logger.warning("Cannot lookup saturation pressure: CoolProp package not available")
        return None
        
    try:
        # Use cached version for better performance
        from .property_cache import cached_saturation_pressure
        return cached_saturation_pressure(fluid_name, temperature_c)
    except Exception as e:
        logger.error("Error getting saturation pressure for %s: %s", fluid_name, e)
        return None

def get_critical_pressure(fluid_name):
    """Attempt to get critical pressure with graceful fallback.
    
    Args:
        fluid_name: Name of fluid
        
    Returns:
        Critical pressure in Pa or None if lookup failed
    """
    if not COOLPROP_AVAILABLE:
        logger.warning("Cannot lookup critical pressure: CoolProp package not available")
        return None
        
    try:
        # Use cached version for better performance
        from .property_cache import cached_critical_pressure
        return cached_critical_pressure(fluid_name)
    except Exception as e:
        logger.error("Error getting critical pressure for %s: %s", fluid_name, e)
        return None
