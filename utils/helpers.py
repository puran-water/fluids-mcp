"""
Helper functions for Fluids MCP server.

This module provides shared utility functions used by multiple tools in the server.
"""

import logging
import fluids
import fluids.friction
import fluids.fittings

logger = logging.getLogger("fluids-mcp.helpers")

def get_fitting_K(fitting_type: str, diameter: float, Re: float, flow_rate: float) -> float:
    """Maps simplified fitting names to fluids library calls.
    
    Args:
        fitting_type: Simplified name of the fitting type (e.g., "90_elbow", "gate_valve")
        diameter: Pipe inner diameter in meters
        Re: Reynolds number
        flow_rate: Flow rate in m³/s
    
    Returns:
        K-value for the fitting
    """
    fitting_type = fitting_type.lower()
    # Mapping fitting types to fluids library functions
    if fitting_type in ["90_elbow", "90_lr_elbow", "elbow_90"]:
        return fluids.fittings.bend_rounded(Di=diameter, angle=90, Re=Re, method='Crane')
    elif fitting_type in ["45_elbow", "elbow_45"]:
        return fluids.fittings.bend_rounded(Di=diameter, angle=45, Re=Re, method='Crane')
    elif fitting_type in ["gate_valve", "gate_valve_open"]:
        return fluids.fittings.K_gate_valve_Crane(D1=diameter, D2=diameter, angle=0)
    elif fitting_type == "gate_valve_half_open":
        return fluids.fittings.K_gate_valve_Crane(D1=diameter, D2=diameter, angle=45)
    elif fitting_type == "gate_valve_quarter_open":
        return fluids.fittings.K_gate_valve_Crane(D1=diameter, D2=diameter, angle=68)
    elif fitting_type == "globe_valve":
        return fluids.fittings.K_globe_valve_Crane(D1=diameter, D2=diameter)
    elif fitting_type == "check_valve_swing":
        return fluids.fittings.K_swing_check_valve_Crane(D=diameter)
    elif fitting_type == "check_valve_lift":
        return fluids.fittings.K_lift_check_valve_Crane(D=diameter)
    elif fitting_type == "ball_valve":
        return fluids.fittings.K_ball_valve_Crane(D1=diameter, D2=diameter, angle=0)
    elif fitting_type == "ball_valve_half_open":
        return fluids.fittings.K_ball_valve_Crane(D1=diameter, D2=diameter, angle=30)
    elif fitting_type == "butterfly_valve":
        return fluids.fittings.K_butterfly_valve_Crane(D=diameter)
    elif fitting_type in ["tee_run", "tee_run_through"]: 
        # Assuming equal split for Crane K factors
        return fluids.fittings.K_run_diverging_Crane(
            D_run=diameter, D_branch=diameter, Q_run=flow_rate*0.5, Q_branch=flow_rate*0.5
        )
    elif fitting_type in ["tee_branch", "tee_branch_flow"]:
        return fluids.fittings.K_branch_diverging_Crane(
            D_run=diameter, D_branch=diameter, Q_run=flow_rate*0.5, Q_branch=flow_rate*0.5
        )
    elif fitting_type == "entrance_sharp":
        return fluids.fittings.entrance_sharp()
    elif fitting_type == "entrance_rounded":
        return fluids.fittings.entrance_rounded(Di=diameter, rc=0.1*diameter)
    elif fitting_type == "entrance_beveled":
        return fluids.fittings.entrance_beveled(Di=diameter, l=0.1*diameter, angle=45)
    elif fitting_type == "exit_normal":
        return fluids.fittings.exit_normal()
    elif fitting_type == "expansion":
        # Assuming 50% expansion by default
        D_larger = diameter * 1.5
        return fluids.fittings.contraction_sharp(Di1=diameter, Di2=D_larger)
    elif fitting_type == "contraction":
        # Assuming 50% contraction by default
        D_smaller = diameter * 0.5
        return fluids.fittings.contraction_sharp(Di1=diameter, Di2=D_smaller)
    # Add more mappings here as needed
    else:
        logger.warning(f"K value calculation not implemented for fitting type: '{fitting_type}'. Assuming K=0.")
        return 0.0

def get_pipe_roughness(material: str = None, pipe_roughness: float = None) -> tuple:
    """Get pipe roughness based on material or default.
    
    Args:
        material: Optional pipe material name
        pipe_roughness: Optional explicit roughness value
        
    Returns:
        Tuple of (roughness value, source description)
    """
    if pipe_roughness is not None:
        return pipe_roughness, "Provided"
        
    if material is not None:
        try:
            from fluids.friction import _roughness
            mat_lower = {k.lower(): k for k in _roughness.keys()}
            if material.lower() in mat_lower:
                actual_key = mat_lower[material.lower()]
                return _roughness[actual_key], f"Material Lookup '{actual_key}'"
            else:
                from utils.constants import DEFAULT_ROUGHNESS
                logger.warning(f"Material '{material}' not found, using default roughness {DEFAULT_ROUGHNESS}")
                return DEFAULT_ROUGHNESS, f"Default (Material '{material}' Not Found)"
        except Exception as e:
            from utils.constants import DEFAULT_ROUGHNESS
            logger.warning(f"Roughness lookup failed: {e}")
            return DEFAULT_ROUGHNESS, "Default (Lookup Failed)"
            
    # Default case
    from utils.constants import DEFAULT_ROUGHNESS
    return DEFAULT_ROUGHNESS, "Default"

def check_fluidprop_available():
    """Check if fluidprop package is available.
    
    Returns:
        Tuple of (available, module if available or None)
    """
    try:
        from fluidprop import FluidProperties, FLUID_SELECTION
        return True, (FluidProperties, FLUID_SELECTION)
    except ImportError:
        logger.warning("fluidprop module not available. Fluid property lookups will be disabled.")
        return False, None

def check_coolprop_available():
    """Check if CoolProp package is available.
    
    Returns:
        Tuple of (available, module if available or None)
    """
    try:
        import CoolProp.CoolProp as CP
        return True, CP
    except ImportError:
        logger.warning("CoolProp module not available. Vapor pressure lookups will be disabled.")
        return False, None
