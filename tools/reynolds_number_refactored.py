"""
Refactored Reynolds number calculation using the new InputResolver system.

This demonstrates how the new architecture eliminates code duplication
and provides consistent input handling across all tools.
"""

import json
import logging
from typing import Optional

# Import the new InputResolver
from utils.input_resolver import InputResolver

# Configure logging
logger = logging.getLogger("fluids-mcp.reynolds_number_refactored")


def calculate_reynolds_number_refactored(
    # Flow inputs (any one required)
    flow_rate_m3_s: Optional[float] = None,
    flow_rate_gpm: Optional[float] = None,
    velocity_m_s: Optional[float] = None,
    
    # Pipe dimension inputs (any one required)
    pipe_diameter_m: Optional[float] = None,
    pipe_diameter_in: Optional[float] = None,
    
    # Fluid property inputs
    fluid_density_kg_m3: Optional[float] = None,
    fluid_density_lbft3: Optional[float] = None,
    fluid_viscosity_pas: Optional[float] = None,
    fluid_viscosity_cp: Optional[float] = None,
    
    # Property lookup option
    fluid_name: Optional[str] = None,
    temperature_c: Optional[float] = None,
    pressure_bar: Optional[float] = None,
    
) -> str:
    """
    Calculate Reynolds number using refactored InputResolver system.
    
    This version demonstrates the new architecture with:
    - Centralized input validation and unit conversion
    - Consistent error handling and logging
    - Elimination of code duplication
    """
    
    # Initialize the input resolver
    resolver = InputResolver("reynolds_number_refactored")
    
    try:
        # Resolve flow rate (only if needed for velocity calculation)
        flow_rate = None
        if velocity_m_s is None:
            flow_rate = resolver.resolve_flow_rate(
                m3_s=flow_rate_m3_s,
                gpm=flow_rate_gpm
            )
        
        # Resolve pipe diameter  
        pipe_diameter = resolver.resolve_dimension(
            "Pipe diameter",
            m=pipe_diameter_m,
            inch=pipe_diameter_in
        )
        
        # Resolve fluid density
        fluid_density = None
        if fluid_density_kg_m3 is not None:
            fluid_density = fluid_density_kg_m3
            resolver.results_log.append("Fluid density: SI (kg/m³)")
        elif fluid_density_lbft3 is not None:
            from utils.constants import LBFT3_to_KGM3
            fluid_density = fluid_density_lbft3 * LBFT3_to_KGM3
            resolver.results_log.append(f"Fluid density: Imperial (lb/ft³) -> {fluid_density:.3f} kg/m³")
        elif fluid_name and temperature_c is not None:
            # Try property lookup
            try:
                from utils.import_helpers import FLUIDPROP_AVAILABLE, FluidProperties
                if FLUIDPROP_AVAILABLE:
                    lookup_p_bar = pressure_bar or 1.0
                    fluid_props = FluidProperties(
                        coolprop_name=fluid_name,
                        T_in_deg_C=temperature_c,
                        P_in_bar=lookup_p_bar
                    )
                    fluid_density = float(fluid_props.rho[0])
                    resolver.results_log.append(f"Fluid density: Lookup ({fluid_name}) -> {fluid_density:.3f} kg/m³")
                else:
                    resolver.error_log.append("FluidProp unavailable for density lookup")
            except Exception as e:
                resolver.error_log.append(f"Fluid density lookup failed: {e}")
        
        if fluid_density is None:
            resolver.error_log.append("Missing fluid density")
        
        # Resolve fluid viscosity
        fluid_viscosity = None
        if fluid_viscosity_pas is not None:
            fluid_viscosity = fluid_viscosity_pas
            resolver.results_log.append("Fluid viscosity: SI (Pa·s)")
        elif fluid_viscosity_cp is not None:
            from utils.constants import CENTIPOISE_to_PAS
            fluid_viscosity = fluid_viscosity_cp * CENTIPOISE_to_PAS
            resolver.results_log.append(f"Fluid viscosity: Imperial (cP) -> {fluid_viscosity:.6f} Pa·s")
        elif fluid_name and temperature_c is not None:
            # Try property lookup
            try:
                from utils.import_helpers import FLUIDPROP_AVAILABLE, FluidProperties
                if FLUIDPROP_AVAILABLE:
                    lookup_p_bar = pressure_bar or 1.0
                    fluid_props = FluidProperties(
                        coolprop_name=fluid_name,
                        T_in_deg_C=temperature_c,
                        P_in_bar=lookup_p_bar
                    )
                    if hasattr(fluid_props, 'eta'):
                        fluid_viscosity = float(fluid_props.eta[0])
                        resolver.results_log.append(f"Fluid viscosity: Lookup ({fluid_name}) -> {fluid_viscosity:.6f} Pa·s")
                    else:
                        resolver.error_log.append("Viscosity not available from property lookup")
                else:
                    resolver.error_log.append("FluidProp unavailable for viscosity lookup")
            except Exception as e:
                resolver.error_log.append(f"Fluid viscosity lookup failed: {e}")
        
        if fluid_viscosity is None:
            resolver.error_log.append("Missing fluid viscosity")
        
        # Check for early errors
        logs = resolver.get_logs()
        if logs["errors"]:
            return json.dumps({
                "errors": logs["errors"],
                "log": logs["log"]
            })
        
        # Calculate velocity if not provided
        if velocity_m_s is None and flow_rate is not None and pipe_diameter is not None:
            import math
            area = math.pi * (pipe_diameter / 2.0) ** 2
            velocity_m_s = flow_rate / area
            resolver.results_log.append(f"Calculated velocity: {velocity_m_s:.3f} m/s")
        elif velocity_m_s is None:
            resolver.error_log.append("Missing velocity (provide directly or via flow_rate + pipe_diameter)")
            return json.dumps({
                "errors": resolver.error_log,
                "log": resolver.results_log
            })
        
        # Calculate Reynolds number
        reynolds_number = (fluid_density * velocity_m_s * pipe_diameter) / fluid_viscosity
        
        # Determine flow regime
        if reynolds_number < 2300:
            flow_regime = "Laminar"
        elif reynolds_number < 4000:
            flow_regime = "Transitional"
        else:
            flow_regime = "Turbulent"
        
        resolver.results_log.append(f"Reynolds number calculated: {reynolds_number:.1f}")
        resolver.results_log.append(f"Flow regime: {flow_regime}")
        
        # Return results
        final_logs = resolver.get_logs()
        return json.dumps({
            "reynolds_number": reynolds_number,
            "flow_regime": flow_regime,
            "inputs_used": {
                "velocity_m_s": velocity_m_s,
                "pipe_diameter_m": pipe_diameter,
                "fluid_density_kg_m3": fluid_density,
                "fluid_viscosity_pas": fluid_viscosity
            },
            "log": final_logs["log"],
            "errors": final_logs["errors"]
        })
        
    except Exception as e:
        resolver.error_log.append(f"Calculation error: {str(e)}")
        final_logs = resolver.get_logs()
        return json.dumps({
            "errors": final_logs["errors"],
            "log": final_logs["log"]
        })