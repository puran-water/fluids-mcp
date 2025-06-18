"""
Vectorized calculation utilities for high-performance batch processing.

This module provides numpy-based vectorized versions of common fluid calculations
for significant performance improvements when processing multiple data points.
"""

import numpy as np
import logging
from typing import Union, List, Tuple, Optional, Dict, Any
from functools import lru_cache

logger = logging.getLogger("fluids-mcp.vectorized_calcs")

def vectorized_pipe_pressure_drop(
    flow_rates: np.ndarray,
    pipe_diameters: np.ndarray,
    pipe_lengths: np.ndarray,
    fluid_densities: np.ndarray,
    fluid_viscosities: np.ndarray,
    pipe_roughnesses: np.ndarray
) -> Dict[str, np.ndarray]:
    """
    Vectorized pressure drop calculation for multiple pipe segments.
    
    Args:
        flow_rates: Array of flow rates in m³/s
        pipe_diameters: Array of pipe diameters in m
        pipe_lengths: Array of pipe lengths in m  
        fluid_densities: Array of fluid densities in kg/m³
        fluid_viscosities: Array of fluid viscosities in Pa·s
        pipe_roughnesses: Array of pipe roughnesses in m
        
    Returns:
        Dictionary containing vectorized results:
        - velocities: Flow velocities in m/s
        - reynolds_numbers: Reynolds numbers
        - friction_factors: Darcy friction factors
        - pressure_drops: Pressure drops in Pa
    """
    
    # Ensure all inputs are numpy arrays of same shape
    arrays = [flow_rates, pipe_diameters, pipe_lengths, 
              fluid_densities, fluid_viscosities, pipe_roughnesses]
    arrays = [np.asarray(arr) for arr in arrays]
    
    # Check shapes are compatible
    shapes = [arr.shape for arr in arrays]
    if not all(shape == shapes[0] for shape in shapes):
        raise ValueError("All input arrays must have the same shape")
    
    # Vectorized calculations
    areas = np.pi * (pipe_diameters / 2.0) ** 2
    velocities = flow_rates / areas
    
    # Reynolds numbers
    reynolds_numbers = fluid_densities * velocities * pipe_diameters / fluid_viscosities
    
    # Relative roughness
    relative_roughness = pipe_roughnesses / pipe_diameters
    
    # Friction factors using Colebrook-White approximation (vectorized)
    # Using explicit approximation for vectorization
    friction_factors = vectorized_friction_factor(reynolds_numbers, relative_roughness)
    
    # Pressure drops
    pressure_drops = friction_factors * (pipe_lengths / pipe_diameters) * \
                    (fluid_densities * velocities**2 / 2.0)
    
    return {
        "velocities": velocities,
        "reynolds_numbers": reynolds_numbers, 
        "friction_factors": friction_factors,
        "pressure_drops": pressure_drops,
        "areas": areas,
        "relative_roughness": relative_roughness
    }

def vectorized_friction_factor(reynolds: np.ndarray, relative_roughness: np.ndarray) -> np.ndarray:
    """
    Vectorized friction factor calculation using Swamee-Jain approximation.
    
    This is more accurate than Blasius and works for both smooth and rough pipes.
    Valid for: 5000 < Re < 10^8 and 10^-6 < ε/D < 10^-2
    
    Args:
        reynolds: Array of Reynolds numbers
        relative_roughness: Array of relative roughness (ε/D)
        
    Returns:
        Array of Darcy friction factors
    """
    
    # Laminar flow
    laminar_mask = reynolds < 2300
    f_laminar = 64.0 / reynolds
    
    # Turbulent flow - Swamee-Jain approximation
    turbulent_mask = reynolds >= 2300
    
    # Avoid division by zero and log(0)
    eps_d_safe = np.maximum(relative_roughness, 1e-12)
    re_safe = np.maximum(reynolds, 1.0)
    
    # Swamee-Jain equation: f = 0.25 / [log10(ε/D/3.7 + 5.74/Re^0.9)]^2
    log_term = np.log10(eps_d_safe / 3.7 + 5.74 / (re_safe ** 0.9))
    f_turbulent = 0.25 / (log_term ** 2)
    
    # Combine laminar and turbulent
    friction_factors = np.where(laminar_mask, f_laminar, f_turbulent)
    
    return friction_factors

def vectorized_pump_curves(
    flow_rates: np.ndarray,
    pump_coefficients: Tuple[float, float, float],
    pump_speed_rpm: float = 1750,
    pump_diameter_in: float = 10
) -> Dict[str, np.ndarray]:
    """
    Vectorized pump performance curve calculations.
    
    Uses typical centrifugal pump affinity laws and performance curves.
    
    Args:
        flow_rates: Array of flow rates in m³/s
        pump_coefficients: Tuple of (a, b, c) for head curve H = a - b*Q - c*Q²
        pump_speed_rpm: Pump speed in RPM
        pump_diameter_in: Impeller diameter in inches
        
    Returns:
        Dictionary with pump performance arrays:
        - heads: Total head in m
        - powers: Shaft power in kW  
        - efficiencies: Pump efficiency (0-1)
        - npsh_required: Required NPSH in m
    """
    
    flow_rates = np.asarray(flow_rates)
    a, b, c = pump_coefficients
    
    # Convert flow rates to GPM for typical pump curve equations
    flow_gpm = flow_rates * 15850.3  # m³/s to GPM
    
    # Head curve: H = a - b*Q - c*Q²
    heads_ft = a - b * flow_gpm - c * flow_gpm**2
    heads_m = heads_ft * 0.3048  # ft to m
    
    # Power curve (typical quadratic): P = d + e*Q + f*Q²
    # Using typical coefficients based on pump size
    d = pump_diameter_in * 0.5  # Base power
    e = pump_diameter_in * 0.001
    f = pump_diameter_in * 0.000001
    
    power_hp = d + e * flow_gpm + f * flow_gpm**2
    power_kw = power_hp * 0.7457  # HP to kW
    
    # Efficiency curve (typical parabolic): η = g - h*(Q - Q_bep)²
    q_bep = np.sqrt(a / (2 * c)) if c > 0 else flow_gpm.mean()  # Best efficiency point
    eta_max = 0.85  # Typical max efficiency
    efficiency_drop_factor = 0.000005
    
    efficiencies = eta_max - efficiency_drop_factor * (flow_gpm - q_bep)**2
    efficiencies = np.clip(efficiencies, 0.1, 0.9)  # Reasonable bounds
    
    # NPSH Required (typical curve): NPSH_R = j + k*Q²
    j = 2.0  # Base NPSH in ft
    k = 0.00001  # Flow dependency
    npsh_required_ft = j + k * flow_gpm**2
    npsh_required_m = npsh_required_ft * 0.3048
    
    return {
        "heads": heads_m,
        "powers": power_kw,
        "efficiencies": efficiencies,
        "npsh_required": npsh_required_m,
        "flow_gpm": flow_gpm
    }

def vectorized_gas_compression(
    inlet_pressures: np.ndarray,
    outlet_pressures: np.ndarray, 
    inlet_temperatures: np.ndarray,
    gas_molecular_weights: np.ndarray,
    specific_heat_ratios: np.ndarray,
    compression_efficiencies: np.ndarray = None
) -> Dict[str, np.ndarray]:
    """
    Vectorized gas compression calculations for multiple operating points.
    
    Args:
        inlet_pressures: Array of inlet pressures in Pa
        outlet_pressures: Array of outlet pressures in Pa
        inlet_temperatures: Array of inlet temperatures in K
        gas_molecular_weights: Array of molecular weights in kg/kmol
        specific_heat_ratios: Array of specific heat ratios (γ = Cp/Cv)
        compression_efficiencies: Array of isentropic efficiencies (default 0.8)
        
    Returns:
        Dictionary with compression results:
        - pressure_ratios: Compression ratios
        - outlet_temperatures_isentropic: Isentropic outlet temperatures in K
        - outlet_temperatures_actual: Actual outlet temperatures in K  
        - work_isentropic: Isentropic work in J/kg
        - work_actual: Actual work in J/kg
        - power_per_mass_flow: Power per unit mass flow in W/(kg/s)
    """
    
    # Convert to numpy arrays
    arrays = [inlet_pressures, outlet_pressures, inlet_temperatures, 
              gas_molecular_weights, specific_heat_ratios]
    arrays = [np.asarray(arr) for arr in arrays]
    
    if compression_efficiencies is None:
        compression_efficiencies = np.full_like(inlet_pressures, 0.8)
    else:
        compression_efficiencies = np.asarray(compression_efficiencies)
    
    # Universal gas constant
    R_universal = 8314.462  # J/(kmol·K)
    
    # Specific gas constants
    R_specific = R_universal / gas_molecular_weights
    
    # Pressure ratios
    pressure_ratios = outlet_pressures / inlet_pressures
    
    # Isentropic exponent for compression
    k = specific_heat_ratios
    k_exp = (k - 1.0) / k
    
    # Isentropic outlet temperatures
    outlet_temperatures_isentropic = inlet_temperatures * (pressure_ratios ** k_exp)
    
    # Isentropic work per unit mass
    work_isentropic = R_specific * inlet_temperatures * (k / (k - 1.0)) * \
                     (pressure_ratios ** k_exp - 1.0)
    
    # Actual work accounting for efficiency
    work_actual = work_isentropic / compression_efficiencies
    
    # Actual outlet temperatures
    outlet_temperatures_actual = inlet_temperatures + \
                               (outlet_temperatures_isentropic - inlet_temperatures) / compression_efficiencies
    
    # Power per unit mass flow
    power_per_mass_flow = work_actual
    
    return {
        "pressure_ratios": pressure_ratios,
        "outlet_temperatures_isentropic": outlet_temperatures_isentropic,
        "outlet_temperatures_actual": outlet_temperatures_actual,
        "work_isentropic": work_isentropic,
        "work_actual": work_actual, 
        "power_per_mass_flow": power_per_mass_flow,
        "specific_gas_constants": R_specific
    }

def vectorized_reynolds_numbers(
    velocities: np.ndarray,
    characteristic_lengths: np.ndarray,
    fluid_densities: np.ndarray,
    fluid_viscosities: np.ndarray
) -> np.ndarray:
    """
    Vectorized Reynolds number calculation.
    
    Args:
        velocities: Array of velocities in m/s
        characteristic_lengths: Array of characteristic lengths (e.g., diameter) in m
        fluid_densities: Array of fluid densities in kg/m³
        fluid_viscosities: Array of dynamic viscosities in Pa·s
        
    Returns:
        Array of Reynolds numbers (dimensionless)
    """
    
    arrays = [velocities, characteristic_lengths, fluid_densities, fluid_viscosities]
    arrays = [np.asarray(arr) for arr in arrays]
    
    reynolds_numbers = fluid_densities * velocities * characteristic_lengths / fluid_viscosities
    
    return reynolds_numbers

@lru_cache(maxsize=100)
def get_batch_property_template(fluid_name: str, n_points: int) -> Dict[str, Any]:
    """
    Get a template for batch property calculations to optimize memory allocation.
    
    Args:
        fluid_name: Fluid name for CoolProp
        n_points: Number of calculation points
        
    Returns:
        Dictionary with pre-allocated arrays for property storage
    """
    
    return {
        "fluid_name": fluid_name,
        "n_points": n_points,
        "densities": np.zeros(n_points),
        "viscosities": np.zeros(n_points),
        "thermal_conductivities": np.zeros(n_points),
        "specific_heats": np.zeros(n_points),
        "temperatures": np.zeros(n_points),
        "pressures": np.zeros(n_points)
    }

def batch_fluid_properties(
    fluid_name: str,
    temperatures_c: np.ndarray,
    pressures_bar: np.ndarray
) -> Dict[str, np.ndarray]:
    """
    Batch fluid property calculations using vectorized operations where possible.
    
    For now, this is a framework for future vectorized CoolProp implementation.
    CoolProp doesn't natively support vectorized calls, but this provides the
    interface for when that becomes available or we implement workarounds.
    
    Args:
        fluid_name: CoolProp fluid name
        temperatures_c: Array of temperatures in Celsius
        pressures_bar: Array of pressures in bar
        
    Returns:
        Dictionary with property arrays
    """
    
    temperatures_c = np.asarray(temperatures_c)
    pressures_bar = np.asarray(pressures_bar)
    
    # For now, we use the cached property system
    # Future: implement true vectorization when CoolProp supports it
    n_points = len(temperatures_c)
    
    densities = np.zeros(n_points)
    viscosities = np.zeros(n_points) 
    thermal_conductivities = np.zeros(n_points)
    specific_heats = np.zeros(n_points)
    
    # Use cached properties for each point (still faster than uncached)
    from .property_cache import cached_fluid_properties
    
    for i, (temp, press) in enumerate(zip(temperatures_c, pressures_bar)):
        try:
            props = cached_fluid_properties(fluid_name, temp, press)
            densities[i] = props[0]  # density
            viscosities[i] = props[1]  # viscosity
            thermal_conductivities[i] = props[2]  # thermal conductivity
            specific_heats[i] = props[3]  # specific heat
        except Exception as e:
            logger.warning("Failed to get properties for %s at T=%f, P=%f: %s", 
                         fluid_name, temp, press, e)
            # Use reasonable defaults
            densities[i] = 1000.0  # kg/m³
            viscosities[i] = 0.001  # Pa·s
            thermal_conductivities[i] = 0.6  # W/(m·K)
            specific_heats[i] = 4180.0  # J/(kg·K)
    
    return {
        "temperatures_c": temperatures_c,
        "pressures_bar": pressures_bar,
        "densities": densities,
        "viscosities": viscosities,
        "thermal_conductivities": thermal_conductivities,
        "specific_heats": specific_heats,
        "fluid_name": fluid_name
    }

def demo_vectorized_performance():
    """
    Demonstration function showing performance benefits of vectorized calculations.
    """
    import time
    
    # Generate test data
    n_points = 1000
    flow_rates = np.linspace(0.001, 0.1, n_points)  # m³/s
    pipe_diameters = np.full(n_points, 0.1)  # 100mm pipe
    pipe_lengths = np.full(n_points, 100.0)  # 100m length
    fluid_densities = np.full(n_points, 1000.0)  # Water density
    fluid_viscosities = np.full(n_points, 0.001)  # Water viscosity
    pipe_roughnesses = np.full(n_points, 4.5e-5)  # Steel roughness
    
    # Time vectorized calculation
    start_time = time.time()
    results = vectorized_pipe_pressure_drop(
        flow_rates, pipe_diameters, pipe_lengths,
        fluid_densities, fluid_viscosities, pipe_roughnesses
    )
    vectorized_time = time.time() - start_time
    
    logger.info("Vectorized calculation of %d points completed in %.4f seconds", 
                n_points, vectorized_time)
    logger.info("Average pressure drop: %.2f Pa", np.mean(results["pressure_drops"]))
    logger.info("Performance: %.1f calculations per second", n_points / vectorized_time)
    
    return results, vectorized_time