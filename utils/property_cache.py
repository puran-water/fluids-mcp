"""
Performance optimization module for caching expensive CoolProp property calculations.

This module provides cached versions of CoolProp calls to significantly improve
performance for repeated property lookups with the same conditions.
"""

import logging
import math
from functools import lru_cache
from typing import Optional, Tuple, Any, Dict
import time

from .json_helpers import is_valid_number

logger = logging.getLogger("fluids-mcp.property_cache")


def validate_coolprop_result(value: float, prop_name: str, fluid_name: str,
                              T_K: float = None, P_Pa: float = None) -> float:
    """
    Validate a CoolProp property result for physically meaningful values.

    CoolProp returns _HUGE (approximately 1e308) for invalid states,
    inf for some error conditions, and nan for failed calculations.

    Args:
        value: The property value from CoolProp
        prop_name: Name of the property (for error messages)
        fluid_name: Name of the fluid (for error messages)
        T_K: Temperature in Kelvin (optional, for error messages)
        P_Pa: Pressure in Pascals (optional, for error messages)

    Returns:
        The validated value

    Raises:
        ValueError: If the value is invalid (inf, nan, or _HUGE)
    """
    if not is_valid_number(value):
        conditions = []
        if T_K is not None:
            conditions.append(f"T={T_K-273.15:.1f}°C")
        if P_Pa is not None:
            conditions.append(f"P={P_Pa/1e5:.2f} bar")
        cond_str = f" at {', '.join(conditions)}" if conditions else ""

        raise ValueError(
            f"CoolProp returned invalid value for {prop_name} of {fluid_name}{cond_str}. "
            f"Value={value}. This may indicate the fluid is outside its valid state range "
            f"(e.g., supercritical, two-phase, or beyond property limits)."
        )
    return value


def validate_properties_dict(properties: Dict[str, float], fluid_name: str,
                             T_K: float = None, P_Pa: float = None) -> Dict[str, float]:
    """
    Validate a dictionary of fluid properties.

    Args:
        properties: Dictionary of property name -> value
        fluid_name: Name of the fluid (for error messages)
        T_K: Temperature in Kelvin (optional)
        P_Pa: Pressure in Pascals (optional)

    Returns:
        The validated properties dictionary

    Raises:
        ValueError: If any property value is invalid
    """
    invalid_props = []
    for name, value in properties.items():
        if value is not None and not is_valid_number(value):
            invalid_props.append(f"{name}={value}")

    if invalid_props:
        conditions = []
        if T_K is not None:
            conditions.append(f"T={T_K-273.15:.1f}°C")
        if P_Pa is not None:
            conditions.append(f"P={P_Pa/1e5:.2f} bar")
        cond_str = f" at {', '.join(conditions)}" if conditions else ""

        raise ValueError(
            f"Invalid CoolProp properties for {fluid_name}{cond_str}: {', '.join(invalid_props)}. "
            f"The fluid may be outside its valid state range."
        )
    return properties

# Cache statistics for monitoring performance
_cache_stats = {
    "hits": 0,
    "misses": 0,
    "total_calls": 0,
    "total_time_saved": 0.0
}

def get_cache_stats() -> Dict[str, Any]:
    """Get cache performance statistics."""
    total = _cache_stats["hits"] + _cache_stats["misses"]
    hit_rate = (_cache_stats["hits"] / total * 100) if total > 0 else 0
    
    return {
        "hit_rate_percent": hit_rate,
        "total_calls": total,
        "cache_hits": _cache_stats["hits"],
        "cache_misses": _cache_stats["misses"],
        "estimated_time_saved_seconds": _cache_stats["total_time_saved"]
    }

def clear_cache_stats():
    """Clear cache statistics."""
    global _cache_stats
    _cache_stats = {"hits": 0, "misses": 0, "total_calls": 0, "total_time_saved": 0.0}

@lru_cache(maxsize=1000)
def cached_props_si(prop: str, input1_name: str, input1_val: float, 
                   input2_name: str, input2_val: float, fluid_name: str) -> float:
    """
    Cached version of CoolProp.PropsSI for individual property lookups.
    
    Args:
        prop: Property to calculate (e.g., 'P', 'PCRIT', 'TCRIT')
        input1_name: First input property name (e.g., 'T')
        input1_val: First input property value
        input2_name: Second input property name (e.g., 'Q')
        input2_val: Second input property value
        fluid_name: CoolProp fluid name
        
    Returns:
        Property value from CoolProp
    """
    global _cache_stats
    
    # Check if this is a cache hit (for statistics)
    cache_info = cached_props_si.cache_info()
    prev_hits = cache_info.hits
    
    start_time = time.time()
    
    try:
        from CoolProp.CoolProp import PropsSI
        result = PropsSI(prop, input1_name, input1_val, input2_name, input2_val, fluid_name)

        # Validate the result for _HUGE, inf, nan values
        T_K = input1_val if input1_name == 'T' else (input2_val if input2_name == 'T' else None)
        P_Pa = input1_val if input1_name == 'P' else (input2_val if input2_name == 'P' else None)
        result = validate_coolprop_result(result, prop, fluid_name, T_K, P_Pa)

        # Update statistics
        elapsed_time = time.time() - start_time
        new_cache_info = cached_props_si.cache_info()

        if new_cache_info.hits > prev_hits:
            # This was a cache hit
            _cache_stats["hits"] += 1
            _cache_stats["total_time_saved"] += elapsed_time * 0.9  # Estimate 90% time saved
        else:
            # This was a cache miss
            _cache_stats["misses"] += 1

        return result
        
    except Exception as e:
        logger.warning("CoolProp calculation failed: %s", e)
        raise

@lru_cache(maxsize=500)
def cached_saturation_pressure(fluid_name: str, temperature_c: float) -> float:
    """
    Cached saturation pressure lookup.
    
    Args:
        fluid_name: CoolProp fluid name
        temperature_c: Temperature in Celsius
        
    Returns:
        Saturation pressure in Pa
    """
    temperature_k = temperature_c + 273.15
    return cached_props_si('P', 'T', temperature_k, 'Q', 0, fluid_name)

@lru_cache(maxsize=100)
def cached_critical_pressure(fluid_name: str) -> float:
    """
    Cached critical pressure lookup.
    
    Args:
        fluid_name: CoolProp fluid name
        
    Returns:
        Critical pressure in Pa
    """
    return cached_props_si('PCRIT', '', 0, '', 0, fluid_name)

@lru_cache(maxsize=100)
def cached_critical_temperature(fluid_name: str) -> float:
    """
    Cached critical temperature lookup.
    
    Args:
        fluid_name: CoolProp fluid name
        
    Returns:
        Critical temperature in K
    """
    return cached_props_si('TCRIT', '', 0, '', 0, fluid_name)

@lru_cache(maxsize=1000)
def cached_fluid_properties(fluid_name: str, temperature_c: float, pressure_bar: float) -> Tuple[float, ...]:
    """
    Cached comprehensive fluid properties lookup.
    
    This replaces the FluidProperties class instantiation with a cached version
    that returns the most commonly used properties as a tuple.
    
    Args:
        fluid_name: CoolProp fluid name  
        temperature_c: Temperature in Celsius
        pressure_bar: Pressure in bar
        
    Returns:
        Tuple of (density_kg_m3, viscosity_pas, thermal_conductivity_w_mk, 
                 specific_heat_cp_j_kgk, molecular_weight_kg_kmol, 
                 kinematic_viscosity_m2s, thermal_diffusivity_m2s)
    """
    global _cache_stats
    
    # Check if this is a cache hit (for statistics)
    cache_info = cached_fluid_properties.cache_info()
    prev_hits = cache_info.hits
    
    start_time = time.time()
    
    try:
        from .import_helpers import FluidProperties

        if FluidProperties is None:
            raise ImportError("FluidProperties is not available (fluidprop package not installed)")

        # Create FluidProperties object (this is expensive)
        fluid_props = FluidProperties(
            coolprop_name=fluid_name,
            T_in_deg_C=temperature_c,
            P_in_bar=pressure_bar
        )
        
        # Extract commonly used properties
        result = (
            float(fluid_props.rho[0]),      # density kg/m³
            float(fluid_props.eta[0]),      # dynamic viscosity Pa·s
            float(fluid_props.lambda_[0]),  # thermal conductivity W/(m·K)
            float(fluid_props.Cp[0]),       # specific heat J/(kg·K)
            float(fluid_props.MW),          # molecular weight kg/kmol
            float(fluid_props.nu[0]),       # kinematic viscosity m²/s
            float(fluid_props.alpha[0])     # thermal diffusivity m²/s
        )
        
        # Update statistics
        elapsed_time = time.time() - start_time
        new_cache_info = cached_fluid_properties.cache_info()
        
        if new_cache_info.hits > prev_hits:
            # This was a cache hit
            _cache_stats["hits"] += 1
            _cache_stats["total_time_saved"] += elapsed_time * 0.95  # Estimate 95% time saved
        else:
            # This was a cache miss  
            _cache_stats["misses"] += 1
            
        return result
        
    except Exception as e:
        logger.warning("FluidProperties calculation failed: %s", e)
        raise

class CachedFluidProperties:
    """
    Cached replacement for FluidProperties class.

    This provides the same interface as FluidProperties but uses CoolProp
    directly for cross-platform compatibility (avoids fluidprop package issues).
    """

    def __init__(self, coolprop_name: str, T_in_deg_C: float, P_in_bar: float):
        self.coolprop_name = coolprop_name
        self.T_in_deg_C = T_in_deg_C
        self.P_in_bar = P_in_bar

        # Use CoolProp directly via cached_props_si
        try:
            T_K = T_in_deg_C + 273.15
            P_Pa = P_in_bar * 1e5

            # Get properties using cached PropsSI calls (validation happens in cached_props_si)
            rho = cached_props_si('D', 'T', T_K, 'P', P_Pa, coolprop_name)  # density kg/m³
            eta = cached_props_si('V', 'T', T_K, 'P', P_Pa, coolprop_name)  # viscosity Pa·s
            lambda_ = cached_props_si('L', 'T', T_K, 'P', P_Pa, coolprop_name)  # conductivity W/(m·K)
            Cp = cached_props_si('C', 'T', T_K, 'P', P_Pa, coolprop_name)  # specific heat J/(kg·K)
            MW = cached_props_si('M', 'T', T_K, 'P', P_Pa, coolprop_name) * 1000  # mol weight kg/kmol

            # Additional validation: ensure physically reasonable values
            if rho <= 0:
                raise ValueError(f"Invalid density {rho} kg/m³ for {coolprop_name}")
            if eta <= 0:
                raise ValueError(f"Invalid viscosity {eta} Pa·s for {coolprop_name}")
            if Cp <= 0:
                raise ValueError(f"Invalid specific heat {Cp} J/(kg·K) for {coolprop_name}")

            # Derived properties
            nu = eta / rho  # kinematic viscosity m²/s
            alpha = lambda_ / (rho * Cp) if lambda_ > 0 else 1e-7  # thermal diffusivity m²/s

            # Store as lists to match FluidProperties interface
            self.rho = [rho]
            self.eta = [eta]
            self.lambda_ = [lambda_]
            self.Cp = [Cp]
            self.MW = MW
            self.nu = [nu]
            self.alpha = [alpha]

            # Derived properties for gas calculations
            R_specific = 8314.46 / MW if MW > 0 else 287.0  # J/(kg·K)
            Cv = Cp - R_specific  # Cv = Cp - R for ideal gas approximation
            self.Cv = [Cv]
            self.gamma = [Cp / Cv] if Cv > 0 else [1.4]
            self.Z = [1.0]  # Compressibility factor (approximation)

        except Exception as e:
            logger.error("Failed to get CoolProp properties for %s: %s", coolprop_name, e)
            raise

def get_cache_info() -> Dict[str, Any]:
    """Get detailed cache information for all cached functions."""
    return {
        "props_si_cache": {
            "hits": cached_props_si.cache_info().hits,
            "misses": cached_props_si.cache_info().misses,
            "maxsize": cached_props_si.cache_info().maxsize,
            "currsize": cached_props_si.cache_info().currsize
        },
        "saturation_pressure_cache": {
            "hits": cached_saturation_pressure.cache_info().hits,
            "misses": cached_saturation_pressure.cache_info().misses,
            "maxsize": cached_saturation_pressure.cache_info().maxsize,
            "currsize": cached_saturation_pressure.cache_info().currsize
        },
        "critical_pressure_cache": {
            "hits": cached_critical_pressure.cache_info().hits,
            "misses": cached_critical_pressure.cache_info().misses,
            "maxsize": cached_critical_pressure.cache_info().maxsize,
            "currsize": cached_critical_pressure.cache_info().currsize
        },
        "fluid_properties_cache": {
            "hits": cached_fluid_properties.cache_info().hits,
            "misses": cached_fluid_properties.cache_info().misses,
            "maxsize": cached_fluid_properties.cache_info().maxsize,
            "currsize": cached_fluid_properties.cache_info().currsize
        },
        "overall_stats": get_cache_stats()
    }

def clear_all_caches():
    """Clear all property caches."""
    cached_props_si.cache_clear()
    cached_saturation_pressure.cache_clear()
    cached_critical_pressure.cache_clear()
    cached_fluid_properties.cache_clear()
    clear_cache_stats()
    logger.info("All property caches cleared")