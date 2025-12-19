"""
Performance optimization module for caching expensive CoolProp property calculations.

This module provides cached versions of CoolProp calls to significantly improve
performance for repeated property lookups with the same conditions.
"""

import logging
from functools import lru_cache
from typing import Optional, Tuple, Any, Dict
import time

logger = logging.getLogger("fluids-mcp.property_cache")

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
        import CoolProp as CP
        result = CP.PropsSI(prop, input1_name, input1_val, input2_name, input2_val, fluid_name)
        
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
    
    This provides the same interface as FluidProperties but uses cached lookups
    for significant performance improvement.
    """
    
    def __init__(self, coolprop_name: str, T_in_deg_C: float, P_in_bar: float):
        self.coolprop_name = coolprop_name
        self.T_in_deg_C = T_in_deg_C
        self.P_in_bar = P_in_bar
        
        # Get cached properties
        try:
            props = cached_fluid_properties(coolprop_name, T_in_deg_C, P_in_bar)
            
            # Unpack properties (maintaining FluidProperties interface)
            self.rho = [props[0]]       # density kg/m³
            self.eta = [props[1]]       # dynamic viscosity Pa·s  
            self.lambda_ = [props[2]]   # thermal conductivity W/(m·K)
            self.Cp = [props[3]]        # specific heat J/(kg·K)
            self.MW = props[4]          # molecular weight kg/kmol
            self.nu = [props[5]]        # kinematic viscosity m²/s
            self.alpha = [props[6]]     # thermal diffusivity m²/s
            
            # Derived properties
            # Calculate specific gas constant from molecular weight
            R_specific = 8314.46 / self.MW if self.MW > 0 else 287.0  # J/(kg·K)
            self.Cv = [self.Cp[0] - R_specific]  # Cv = Cp - R for ideal gas approximation
            self.gamma = [self.Cp[0] / self.Cv[0]] if self.Cv[0] > 0 else [1.4]
            self.Z = [1.0]  # Compressibility factor (approximation)
            
        except Exception as e:
            logger.error("Failed to get cached properties for %s: %s", coolprop_name, e)
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