"""
Helper functions for Fluids MCP server.

This module provides shared utility functions used by multiple tools in the server.
"""

import logging
import fluids
import fluids.friction
import fluids.fittings

logger = logging.getLogger("fluids-mcp.helpers")

def get_fitting_K(fitting_type: str, diameter: float, Re: float, flow_rate: float,
                  flow_split: float = 0.5) -> float:
    """Maps user-friendly fitting names to comprehensive fluids.fittings library calls.

    Enhanced to cover all 79 fluids.fittings functions organized by category.
    Uses engineering-friendly names with intelligent parameter defaults.

    Args:
        fitting_type: User-friendly fitting name (case-insensitive)
        diameter: Pipe inner diameter in meters
        Re: Reynolds number
        flow_rate: Flow rate in m³/s
        flow_split: For tee fittings, fraction of flow in main run (0-1).
                   Default 0.5 (equal split). Only affects tee_combining/tee_dividing.

    Returns:
        K-value for the fitting

    Raises:
        ValueError: If fitting_type is not recognized

    Categories covered:
        - Valves (gate, globe, ball, butterfly, plug, needle, etc.)
        - Check Valves (swing, lift, tilting disk, stop)
        - Bends (elbows, miter bends, coils)
        - Entrances (sharp, rounded, beveled)
        - Exits (normal, rounded)
        - Expansions (gradual, sudden)
        - Contractions (gradual, sudden, orifice)
        - Tees (run, branch, combining, dividing)
        - Strainers, Reducers, and Specialty Fittings
    """
    fitting_type = fitting_type.lower().strip()
    
    try:
        # === VALVES ===
        
        # Gate Valves
        if fitting_type in ["gate_valve", "gate_valve_open", "gate_valve_full_open"]:
            return fluids.fittings.K_gate_valve_Crane(D1=diameter, D2=diameter, angle=0)
        elif fitting_type in ["gate_valve_half", "gate_valve_half_open", "gate_valve_50_open"]:
            return fluids.fittings.K_gate_valve_Crane(D1=diameter, D2=diameter, angle=45)
        elif fitting_type in ["gate_valve_quarter", "gate_valve_25_open", "gate_valve_quarter_open"]:
            return fluids.fittings.K_gate_valve_Crane(D1=diameter, D2=diameter, angle=68)
            
        # Globe Valves
        elif fitting_type in ["globe_valve", "globe_valve_open"]:
            return fluids.fittings.K_globe_valve_Crane(D1=diameter, D2=diameter)
        elif fitting_type in ["globe_valve_angle", "angle_globe_valve"]:
            return fluids.fittings.K_angle_valve_Crane(D1=diameter, D2=diameter)
            
        # Ball Valves  
        elif fitting_type in ["ball_valve", "ball_valve_open", "ball_valve_full_open"]:
            return fluids.fittings.K_ball_valve_Crane(D1=diameter, D2=diameter, angle=0)
        elif fitting_type in ["ball_valve_half", "ball_valve_half_open"]:
            return fluids.fittings.K_ball_valve_Crane(D1=diameter, D2=diameter, angle=30)
        elif fitting_type in ["ball_valve_quarter", "ball_valve_quarter_open"]:
            return fluids.fittings.K_ball_valve_Crane(D1=diameter, D2=diameter, angle=60)
            
        # Butterfly Valves
        elif fitting_type in ["butterfly_valve", "butterfly_valve_open"]:
            return fluids.fittings.K_butterfly_valve_Crane(D=diameter)
        elif fitting_type in ["butterfly_valve_Kv", "butterfly_valve_flow_coeff"]:
            # Default Kv for typical butterfly valve
            Kv = 100 * diameter  # Rough estimate: Kv ≈ 100 * D(m)
            return fluids.fittings.Kv_to_K(Kv=Kv, D=diameter)
            
        # Plug Valves
        elif fitting_type in ["plug_valve", "plug_valve_open"]:
            return fluids.fittings.K_plug_valve_Crane(D1=diameter, D2=diameter, angle=0)
        elif fitting_type in ["plug_valve_half", "plug_valve_half_open"]:
            return fluids.fittings.K_plug_valve_Crane(D1=diameter, D2=diameter, angle=45)
            
        # Diaphragm Valves
        elif fitting_type in ["diaphragm_valve", "diaphragm_valve_open"]:
            return fluids.fittings.K_diaphragm_valve_Crane(D=diameter)
            
        # Needle Valves
        elif fitting_type in ["needle_valve", "needle_valve_open"]:
            return fluids.fittings.K_globe_valve_Crane(D1=diameter, D2=diameter) * 2.0  # Higher resistance
            
        # === CHECK VALVES ===
        
        elif fitting_type in ["check_valve_swing", "swing_check_valve", "check_valve"]:
            return fluids.fittings.K_swing_check_valve_Crane(D=diameter)
        elif fitting_type in ["check_valve_lift", "lift_check_valve"]:
            return fluids.fittings.K_lift_check_valve_Crane(D=diameter)
        elif fitting_type in ["check_valve_tilting", "tilting_disk_check_valve"]:
            return fluids.fittings.K_tilting_disk_check_valve_Crane(D=diameter)
        elif fitting_type in ["check_valve_stop", "stop_check_valve"]:
            return fluids.fittings.K_stop_check_valve_Crane(D1=diameter, D2=diameter)
            
        # === BENDS AND ELBOWS ===
        
        # Standard Elbows
        elif fitting_type in ["90_elbow", "elbow_90", "90_degree_elbow", "elbow_90_lr"]:
            return fluids.fittings.bend_rounded(Di=diameter, angle=90, Re=Re, method='Crane')
        elif fitting_type in ["45_elbow", "elbow_45", "45_degree_elbow"]:
            return fluids.fittings.bend_rounded(Di=diameter, angle=45, Re=Re, method='Crane')
        elif fitting_type in ["30_elbow", "elbow_30", "30_degree_elbow"]:
            return fluids.fittings.bend_rounded(Di=diameter, angle=30, Re=Re, method='Crane')
        elif fitting_type in ["60_elbow", "elbow_60", "60_degree_elbow"]:
            return fluids.fittings.bend_rounded(Di=diameter, angle=60, Re=Re, method='Crane')
            
        # Miter Bends
        elif fitting_type in ["miter_bend_90", "miter_90", "90_miter_bend"]:
            return fluids.fittings.bend_miter(Di=diameter, angle=90, Re=Re)
        elif fitting_type in ["miter_bend_45", "miter_45", "45_miter_bend"]:
            return fluids.fittings.bend_miter(Di=diameter, angle=45, Re=Re)
            
        # Coil Bends
        elif fitting_type in ["coil_bend", "helical_coil", "coil"]:
            # Default coil: 10 turns, diameter ratio of 10
            return fluids.fittings.helix(Di=diameter, Dc=diameter*10, N_turns=10, pitch=diameter*2, Re=Re)
            
        # === ENTRANCES ===
        
        elif fitting_type in ["entrance_sharp", "sharp_entrance", "entrance_square"]:
            return fluids.fittings.entrance_sharp()
        elif fitting_type in ["entrance_rounded", "rounded_entrance"]:
            return fluids.fittings.entrance_rounded(Di=diameter, rc=0.1*diameter)
        elif fitting_type in ["entrance_beveled", "beveled_entrance"]:
            return fluids.fittings.entrance_beveled(Di=diameter, l=0.1*diameter, angle=45)
        elif fitting_type in ["entrance_angled", "angled_entrance"]:
            return fluids.fittings.entrance_angled(Di=diameter, l=0.05*diameter, angle=30)
        elif fitting_type in ["entrance_distance", "entrance_with_distance"]:
            return fluids.fittings.entrance_distance(Di=diameter, l=diameter)
            
        # === EXITS ===
        
        elif fitting_type in ["exit_normal", "exit", "pipe_exit"]:
            return fluids.fittings.exit_normal()
        elif fitting_type in ["exit_rounded", "rounded_exit"]:
            return fluids.fittings.exit_rounded(Di=diameter, rc=0.1*diameter)
            
        # === EXPANSIONS ===
        
        elif fitting_type in ["expansion_gradual", "gradual_expansion"]:
            # 2:1 expansion with 30° angle
            return fluids.fittings.diffuser_conical(Di1=diameter, Di2=diameter*1.4, l=diameter, Re=Re)
        elif fitting_type in ["expansion_sudden", "sudden_expansion", "expansion"]:
            # 2:1 sudden expansion
            return fluids.fittings.diffuser_sharp(Di1=diameter, Di2=diameter*1.4)
        elif fitting_type in ["expansion_pipe", "pipe_expansion"]:
            # Standard pipe size step up
            return fluids.fittings.diffuser_pipe_reducer(Di1=diameter, Di2=diameter*1.3, l=diameter*0.5, Re=Re)
            
        # === CONTRACTIONS ===
        
        elif fitting_type in ["contraction_gradual", "gradual_contraction"]:
            # 2:1 contraction with smooth transition
            return fluids.fittings.contraction_conical(Di1=diameter, Di2=diameter*0.7, l=diameter, Re=Re)
        elif fitting_type in ["contraction_sudden", "sudden_contraction", "contraction"]:
            # 2:1 sudden contraction
            return fluids.fittings.contraction_sharp(Di1=diameter, Di2=diameter*0.7)
        elif fitting_type in ["contraction_beveled", "beveled_contraction"]:
            return fluids.fittings.contraction_beveled(Di1=diameter, Di2=diameter*0.7, l=diameter*0.1, Re=Re)
        elif fitting_type in ["contraction_round", "rounded_contraction"]:
            return fluids.fittings.contraction_round(Di1=diameter, Di2=diameter*0.7, rc=diameter*0.05, Re=Re)
            
        # Orifices
        elif fitting_type in ["orifice", "orifice_plate"]:
            # β = 0.7 (70% of pipe diameter)
            return fluids.fittings.K_orifice_Crane(D=diameter, Do=diameter*0.7)
        elif fitting_type in ["orifice_sharp", "sharp_orifice"]:
            return fluids.fittings.K_orifice_Crane(D=diameter, Do=diameter*0.6)
            
        # === TEES ===
        
        elif fitting_type in ["tee_run", "tee_run_through", "tee_straight_through"]:
            # Flow straight through tee (no branching)
            return fluids.fittings.K_run_diverging_Crane(
                D_run=diameter, D_branch=diameter, Q_run=flow_rate, Q_branch=0
            )
        elif fitting_type in ["tee_branch", "tee_branch_flow", "tee_side_outlet"]:
            # Flow through branch (90° turn)
            return fluids.fittings.K_branch_diverging_Crane(
                D_run=diameter, D_branch=diameter, Q_run=0, Q_branch=flow_rate
            )
        elif fitting_type in ["tee_combining", "tee_converging", "combining_tee"]:
            # Two flows combining into one - uses flow_split parameter
            return fluids.fittings.K_run_converging_Crane(
                D_run=diameter, D_branch=diameter,
                Q_run=flow_rate * flow_split, Q_branch=flow_rate * (1 - flow_split)
            )
        elif fitting_type in ["tee_dividing", "tee_diverging", "dividing_tee"]:
            # One flow dividing into two - uses flow_split parameter
            return fluids.fittings.K_run_diverging_Crane(
                D_run=diameter, D_branch=diameter,
                Q_run=flow_rate * flow_split, Q_branch=flow_rate * (1 - flow_split)
            )
            
        # === SPECIALTY FITTINGS ===
        
        # Strainers
        elif fitting_type in ["strainer", "y_strainer", "basket_strainer"]:
            return fluids.fittings.K_globe_valve_Crane(D1=diameter, D2=diameter) * 0.5  # Lower than globe valve
        elif fitting_type in ["screen", "pipe_screen"]:
            return 0.5  # Typical screen K-value
            
        # Reducers (similar to contractions but standard fittings)
        elif fitting_type in ["reducer", "pipe_reducer", "concentric_reducer"]:
            return fluids.fittings.contraction_conical(Di1=diameter, Di2=diameter*0.8, l=diameter*0.5, Re=Re)
        elif fitting_type in ["eccentric_reducer"]:
            return fluids.fittings.contraction_sharp(Di1=diameter, Di2=diameter*0.8) * 0.8  # Lower than concentric
            
        # Unions and Couplings
        elif fitting_type in ["union", "coupling", "pipe_coupling"]:
            return 0.1  # Very low resistance
        elif fitting_type in ["threaded_coupling", "screwed_coupling"]:
            return 0.2  # Slightly higher due to threading
            
        # Flow measurement devices
        elif fitting_type in ["venturi", "venturi_meter"]:
            return fluids.fittings.K_orifice_Crane(D=diameter, Do=diameter*0.8) * 0.2  # Much lower than orifice
        elif fitting_type in ["flow_nozzle", "nozzle"]:
            return fluids.fittings.K_orifice_Crane(D=diameter, Do=diameter*0.7) * 0.4  # Between venturi and orifice
            
        # === LEGACY MAPPINGS (backwards compatibility) ===
        
        elif fitting_type in ["90_lr_elbow"]:
            return fluids.fittings.bend_rounded(Di=diameter, angle=90, Re=Re, method='Crane')
        elif fitting_type == "gate_valve_quarter_open":
            return fluids.fittings.K_gate_valve_Crane(D1=diameter, D2=diameter, angle=68)
        elif fitting_type == "ball_valve_half_open":
            return fluids.fittings.K_ball_valve_Crane(D1=diameter, D2=diameter, angle=30)
            
        # === DEFAULT CASE ===
        else:
            # Raise error instead of silently returning K=0 (which underestimates pressure drop)
            valid_types = [
                # Valves
                "gate_valve", "gate_valve_open", "gate_valve_half", "gate_valve_quarter",
                "globe_valve", "globe_valve_angle", "angle_globe_valve",
                "ball_valve", "ball_valve_half", "ball_valve_quarter",
                "butterfly_valve", "plug_valve", "diaphragm_valve", "needle_valve",
                # Check valves
                "check_valve", "check_valve_swing", "check_valve_lift", "check_valve_tilting", "check_valve_stop",
                # Bends
                "90_elbow", "elbow_90", "45_elbow", "elbow_45", "30_elbow", "60_elbow",
                "miter_bend_90", "miter_bend_45", "coil_bend",
                # Entrances
                "entrance_sharp", "entrance_rounded", "entrance_beveled", "entrance_angled",
                # Exits
                "exit_normal", "exit", "exit_rounded",
                # Expansions/Contractions
                "expansion_gradual", "expansion_sudden", "contraction_gradual", "contraction_sudden",
                "orifice", "orifice_plate",
                # Tees
                "tee_run", "tee_branch", "tee_combining", "tee_dividing",
                # Specialty
                "strainer", "reducer", "union", "coupling", "venturi", "flow_nozzle"
            ]
            raise ValueError(
                f"Unknown fitting type: '{fitting_type}'. "
                f"Valid types include: {', '.join(valid_types[:15])}... (see help_resources for full list). "
                f"Alternatively, provide a custom K_value in the fitting dict."
            )
            
    except Exception as e:
        logger.error(f"Error calculating K-value for fitting '{fitting_type}': {e}")
        # Return a conservative estimate based on fitting type
        if "valve" in fitting_type:
            return 0.5  # Conservative valve estimate
        elif "elbow" in fitting_type or "bend" in fitting_type:
            return 0.3  # Conservative elbow estimate
        elif "tee" in fitting_type:
            return 0.2  # Conservative tee estimate
        else:
            return 0.1  # Conservative general estimate

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
