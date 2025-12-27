"""Help resources omnitool - lists available fittings, fluids, materials, and methods."""

import json
import logging
from typing import Literal, Optional

logger = logging.getLogger("fluids-mcp.help_resources")


# Comprehensive list of available fitting types grouped by category
FITTING_TYPES = {
    "valves": {
        "gate_valves": [
            {"name": "gate_valve", "aliases": ["gate_valve_open"], "description": "Fully open gate valve"},
            {"name": "gate_valve_half", "aliases": ["gate_valve_half_open"], "description": "Half open gate valve"},
            {"name": "gate_valve_quarter", "aliases": ["gate_valve_quarter_open"], "description": "Quarter open gate valve"},
        ],
        "globe_valves": [
            {"name": "globe_valve", "aliases": ["globe_valve_open"], "description": "Standard globe valve"},
            {"name": "globe_valve_angle", "aliases": ["angle_globe_valve"], "description": "Angle globe valve"},
        ],
        "ball_valves": [
            {"name": "ball_valve", "aliases": ["ball_valve_open"], "description": "Fully open ball valve"},
            {"name": "ball_valve_half", "aliases": ["ball_valve_half_open"], "description": "Half open ball valve"},
            {"name": "ball_valve_quarter", "aliases": ["ball_valve_quarter_open"], "description": "Quarter open ball valve"},
        ],
        "butterfly_valves": [
            {"name": "butterfly_valve", "aliases": ["butterfly_valve_open"], "description": "Butterfly valve"},
        ],
        "other_valves": [
            {"name": "plug_valve", "aliases": ["plug_valve_open"], "description": "Plug valve"},
            {"name": "plug_valve_half", "aliases": ["plug_valve_half_open"], "description": "Half open plug valve"},
            {"name": "diaphragm_valve", "aliases": ["diaphragm_valve_open"], "description": "Diaphragm valve"},
            {"name": "needle_valve", "aliases": ["needle_valve_open"], "description": "Needle valve"},
        ],
    },
    "check_valves": [
        {"name": "check_valve_swing", "aliases": ["swing_check_valve", "check_valve"], "description": "Swing check valve"},
        {"name": "check_valve_lift", "aliases": ["lift_check_valve"], "description": "Lift check valve"},
        {"name": "check_valve_tilting", "aliases": ["tilting_disk_check_valve"], "description": "Tilting disk check valve"},
        {"name": "check_valve_stop", "aliases": ["stop_check_valve"], "description": "Stop check valve"},
    ],
    "elbows_and_bends": [
        {"name": "90_elbow", "aliases": ["elbow_90", "90_degree_elbow", "elbow_90_lr"], "description": "90-degree elbow"},
        {"name": "45_elbow", "aliases": ["elbow_45", "45_degree_elbow"], "description": "45-degree elbow"},
        {"name": "30_elbow", "aliases": ["elbow_30", "30_degree_elbow"], "description": "30-degree elbow"},
        {"name": "60_elbow", "aliases": ["elbow_60", "60_degree_elbow"], "description": "60-degree elbow"},
        {"name": "miter_bend_90", "aliases": ["miter_90", "90_miter_bend"], "description": "90-degree miter bend"},
        {"name": "miter_bend_45", "aliases": ["miter_45", "45_miter_bend"], "description": "45-degree miter bend"},
        {"name": "coil_bend", "aliases": ["helical_coil", "coil"], "description": "Helical coil (10 turns, D_coil=10*D_pipe)"},
    ],
    "entrances": [
        {"name": "entrance_sharp", "aliases": ["sharp_entrance", "entrance_square"], "description": "Sharp-edged entrance (K~0.5)"},
        {"name": "entrance_rounded", "aliases": ["rounded_entrance"], "description": "Rounded entrance (rc=0.1D)"},
        {"name": "entrance_beveled", "aliases": ["beveled_entrance"], "description": "Beveled entrance (45 deg, l=0.1D)"},
        {"name": "entrance_angled", "aliases": ["angled_entrance"], "description": "Angled entrance"},
        {"name": "entrance_distance", "aliases": ["entrance_with_distance"], "description": "Entrance with distance from wall"},
    ],
    "exits": [
        {"name": "exit_normal", "aliases": ["exit", "pipe_exit"], "description": "Normal pipe exit (K=1.0)"},
        {"name": "exit_rounded", "aliases": ["rounded_exit"], "description": "Rounded pipe exit"},
    ],
    "expansions": [
        {"name": "expansion_gradual", "aliases": ["gradual_expansion"], "description": "Gradual expansion (2:1 area, 30 deg)"},
        {"name": "expansion_sudden", "aliases": ["sudden_expansion", "expansion"], "description": "Sudden expansion (2:1 area)"},
        {"name": "expansion_pipe", "aliases": ["pipe_expansion"], "description": "Standard pipe size step-up"},
    ],
    "contractions": [
        {"name": "contraction_gradual", "aliases": ["gradual_contraction"], "description": "Gradual contraction (0.7:1 dia)"},
        {"name": "contraction_sudden", "aliases": ["sudden_contraction", "contraction"], "description": "Sudden contraction"},
        {"name": "contraction_beveled", "aliases": ["beveled_contraction"], "description": "Beveled contraction"},
        {"name": "contraction_round", "aliases": ["rounded_contraction"], "description": "Rounded contraction"},
    ],
    "orifices": [
        {"name": "orifice", "aliases": ["orifice_plate"], "description": "Orifice plate (beta=0.7)"},
        {"name": "orifice_sharp", "aliases": ["sharp_orifice"], "description": "Sharp orifice (beta=0.6)"},
    ],
    "tees": [
        {"name": "tee_run", "aliases": ["tee_run_through", "tee_straight_through"], "description": "Flow straight through tee"},
        {"name": "tee_branch", "aliases": ["tee_branch_flow", "tee_side_outlet"], "description": "Flow through branch (90 deg turn)"},
        {"name": "tee_combining", "aliases": ["tee_converging", "combining_tee"], "description": "Two flows combining (use flow_split)"},
        {"name": "tee_dividing", "aliases": ["tee_diverging", "dividing_tee"], "description": "Flow dividing into two (use flow_split)"},
    ],
    "specialty": [
        {"name": "strainer", "aliases": ["y_strainer", "basket_strainer"], "description": "Y-strainer or basket strainer"},
        {"name": "screen", "aliases": ["pipe_screen"], "description": "Pipe screen (K~0.5)"},
        {"name": "reducer", "aliases": ["pipe_reducer", "concentric_reducer"], "description": "Concentric pipe reducer"},
        {"name": "eccentric_reducer", "aliases": [], "description": "Eccentric reducer"},
        {"name": "union", "aliases": ["coupling", "pipe_coupling"], "description": "Union or coupling (K~0.1)"},
        {"name": "threaded_coupling", "aliases": ["screwed_coupling"], "description": "Threaded coupling (K~0.2)"},
        {"name": "venturi", "aliases": ["venturi_meter"], "description": "Venturi meter (low loss)"},
        {"name": "flow_nozzle", "aliases": ["nozzle"], "description": "Flow nozzle"},
    ],
}

# Gas flow calculation methods with descriptions
GAS_METHODS = [
    {
        "name": "Weymouth",
        "description": "Classic high-pressure gas transmission formula. Best for large diameter pipes (>8 inch), high pressures.",
        "accuracy": "Moderate - suitable for screening calculations",
    },
    {
        "name": "Panhandle_A",
        "description": "Modified Panhandle equation. Good for medium to large diameter, clean gas.",
        "accuracy": "Better than Weymouth for moderate Reynolds numbers",
    },
    {
        "name": "Panhandle_B",
        "description": "Improved Panhandle formula for fully turbulent flow.",
        "accuracy": "Good for very high Reynolds numbers",
    },
    {
        "name": "IGT",
        "description": "Institute of Gas Technology equation. Includes viscosity dependence.",
        "accuracy": "More accurate - recommended for general use",
    },
    {
        "name": "Oliphant",
        "description": "For low pressure distribution systems (<100 kPa gauge).",
        "accuracy": "Good for low pressure applications",
    },
    {
        "name": "Spitzglass_low",
        "description": "Spitzglass formula for low pressure (<10 kPa gauge).",
        "accuracy": "Good for very low pressure distribution",
    },
    {
        "name": "Spitzglass_high",
        "description": "Spitzglass formula for medium pressure (up to 100 kPa gauge).",
        "accuracy": "Good for medium pressure distribution",
    },
    {
        "name": "isothermal_darcy",
        "description": "Rigorous isothermal compressible flow with friction. Uses Darcy-Weisbach.",
        "accuracy": "Most accurate - uses real fluid properties",
    },
]


def get_available_materials() -> list:
    """Get list of available pipe materials for roughness lookup."""
    try:
        from fluids.friction import _roughness
        materials = []
        for name, roughness in sorted(_roughness.items()):
            materials.append({
                "name": name,
                "roughness_m": roughness,
                "roughness_mm": roughness * 1000,
            })
        return materials
    except ImportError:
        return [{"error": "fluids library not available for material lookup"}]


def get_available_fluids() -> list:
    """Get list of available fluids from CoolProp."""
    try:
        import CoolProp.CoolProp as CP
        fluids = []
        for fluid in CP.FluidsList():
            try:
                # Get critical properties for reference
                Tc = CP.PropsSI('TCRIT', '', 0, '', 0, fluid) - 273.15  # to Celsius
                Pc = CP.PropsSI('PCRIT', '', 0, '', 0, fluid) / 1e5  # to bar
                fluids.append({
                    "name": fluid,
                    "critical_temp_c": round(Tc, 1) if Tc > -1e6 else None,
                    "critical_pressure_bar": round(Pc, 2) if Pc > 0 else None,
                })
            except Exception:
                fluids.append({"name": fluid})
        return fluids
    except ImportError:
        return [{"error": "CoolProp not available for fluid lookup"}]


def help_resources(
    resource_type: Literal["fittings", "fluids", "materials", "methods", "all"] = "all",
    category: Optional[str] = None,
) -> str:
    """
    List available resources for fluid dynamics calculations.

    Provides comprehensive lists of:
    - Fitting types with aliases and descriptions
    - Available fluids from CoolProp
    - Pipe materials with roughness values
    - Gas flow calculation methods

    Args:
        resource_type: Type of resources to list
            - "fittings": All available fitting types for pressure drop calculations
            - "fluids": All fluids available in CoolProp
            - "materials": Pipe materials with roughness values
            - "methods": Gas flow calculation methods
            - "all": All resources (summary)
        category: Optional category filter for fittings (e.g., "valves", "elbows_and_bends")

    Returns:
        JSON string with requested resource information

    Examples:
        List all fittings:
        >>> help_resources(resource_type="fittings")

        List only valve fittings:
        >>> help_resources(resource_type="fittings", category="valves")

        List gas flow methods:
        >>> help_resources(resource_type="methods")

        List all pipe materials:
        >>> help_resources(resource_type="materials")
    """
    result = {}

    if resource_type == "fittings" or resource_type == "all":
        if category:
            if category in FITTING_TYPES:
                result["fittings"] = {category: FITTING_TYPES[category]}
            else:
                result["fittings"] = {
                    "error": f"Unknown category: {category}",
                    "available_categories": list(FITTING_TYPES.keys()),
                }
        else:
            result["fittings"] = FITTING_TYPES

    if resource_type == "fluids" or resource_type == "all":
        fluids = get_available_fluids()
        if resource_type == "all":
            # Summary only for "all" mode
            result["fluids"] = {
                "count": len(fluids),
                "note": "Use resource_type='fluids' for full list",
                "examples": [f["name"] for f in fluids[:10]] if len(fluids) > 10 else fluids,
            }
        else:
            result["fluids"] = fluids

    if resource_type == "materials" or resource_type == "all":
        materials = get_available_materials()
        if resource_type == "all":
            result["materials"] = {
                "count": len(materials),
                "note": "Use resource_type='materials' for full list with roughness values",
                "common_materials": [
                    {"name": "Carbon steel, bare", "roughness_mm": 0.046},
                    {"name": "Stainless steel, polished", "roughness_mm": 0.015},
                    {"name": "PVC plastic", "roughness_mm": 0.0015},
                    {"name": "Cast iron", "roughness_mm": 0.26},
                ],
            }
        else:
            result["materials"] = materials

    if resource_type == "methods" or resource_type == "all":
        result["gas_methods"] = GAS_METHODS

    if resource_type == "all":
        result["notes"] = {
            "fittings_usage": "Use fitting 'type' field in fittings list: [{'type': '90_elbow', 'quantity': 2}]",
            "custom_k_value": "For unlisted fittings, provide 'K_value' directly: [{'type': 'custom', 'K_value': 0.5}]",
            "flow_split": "For tee_combining/tee_dividing, use flow_split parameter (0-1)",
            "glycol_fluids": "Glycol fluids require INCOMP:: format: 'INCOMP::MEG[0.3]' for 30% MEG",
        }

    return json.dumps(result, indent=2)
