"""
Helper functions for Fluids MCP server.

This module provides shared utility functions used by multiple tools in the server.
"""

import logging
import math
from typing import Any

import fluids
import fluids.friction
import fluids.fittings

from fluids_mcp.errors import FittingKError, PipeMaterialError

logger = logging.getLogger("fluids-mcp.helpers")


MATERIAL_ALIASES = {
    "carbon steel": "Steel",
    "carbon-steel": "Steel",
    "carbon_steel": "Steel",
    "carbonsteel": "Steel",
    "carbon pipe": "Steel",
    "CS": "Steel",
    "A53": "Steel",
    "ASTM A53": "Steel",
    "SA-53": "Steel",
    "A106": "Steel",
    "ASTM A106": "Steel",
    "SA-106": "Steel",
    "A106B": "Steel",
    "ASTM A106B": "Steel",
    "SA-106B": "Steel",
    "A106 Grade B": "Steel",
    "galvanized": "Galvanized iron",
    "galvanized steel": "Galvanized iron",
    "galvanized iron": "Galvanized iron",
    "GI": "Galvanized iron",
    "cast iron": "Cast iron",
    "CI": "Cast iron",
    "asphalted cast iron": "Asphalted cast iron",
    "riveted steel": "Riveted steel",
    "rough riveted steel": "Rough riveted steel",
}


VALID_FITTING_TYPES = frozenset({
    "30_degree_elbow", "30_elbow", "45_degree_elbow", "45_elbow", "45_miter_bend",
    "60_degree_elbow", "60_elbow", "90_degree_elbow", "90_elbow", "90_lr_elbow",
    "90_miter_bend", "angle_globe_valve", "angled_entrance", "ball_valve",
    "ball_valve_full_open", "ball_valve_half", "ball_valve_half_open", "ball_valve_open",
    "ball_valve_quarter", "ball_valve_quarter_open", "basket_strainer", "beveled_contraction",
    "beveled_entrance", "butterfly_valve", "butterfly_valve_flow_coeff", "butterfly_valve_kv",
    "butterfly_valve_open", "check_valve", "check_valve_lift", "check_valve_stop",
    "check_valve_swing", "check_valve_tilting", "coil", "coil_bend", "combining_tee",
    "concentric_reducer", "contraction", "contraction_beveled", "contraction_gradual",
    "contraction_round", "contraction_sudden", "coupling", "diaphragm_valve",
    "diaphragm_valve_open", "dividing_tee", "eccentric_reducer", "elbow_30", "elbow_45",
    "elbow_60", "elbow_90", "elbow_90_lr", "entrance_angled", "entrance_beveled",
    "entrance_distance", "entrance_rounded", "entrance_sharp", "entrance_square",
    "entrance_with_distance", "exit", "exit_normal", "exit_rounded", "expansion",
    "expansion_gradual", "expansion_pipe", "expansion_sudden", "flow_nozzle", "gate_valve",
    "gate_valve_25_open", "gate_valve_50_open", "gate_valve_full_open", "gate_valve_half",
    "gate_valve_half_open", "gate_valve_open", "gate_valve_quarter",
    "gate_valve_quarter_open", "globe_valve", "globe_valve_angle", "globe_valve_open",
    "gradual_contraction", "gradual_expansion", "helical_coil", "lift_check_valve", "miter_45",
    "miter_90", "miter_bend_45", "miter_bend_90", "needle_valve", "needle_valve_open",
    "nozzle", "orifice", "orifice_plate", "orifice_sharp", "pipe_coupling", "pipe_exit",
    "pipe_expansion", "pipe_reducer", "pipe_screen", "plug_valve", "plug_valve_half",
    "plug_valve_half_open", "plug_valve_open", "reducer", "rounded_contraction",
    "rounded_entrance", "rounded_exit", "screen", "screwed_coupling", "sharp_entrance",
    "sharp_orifice", "stop_check_valve", "strainer", "sudden_contraction", "sudden_expansion",
    "swing_check_valve", "tee_branch", "tee_branch_flow", "tee_combining", "tee_converging",
    "tee_diverging", "tee_dividing", "tee_run", "tee_run_through", "tee_side_outlet",
    "tee_straight_through", "threaded_coupling", "tilting_disk_check_valve", "union",
    "venturi", "venturi_meter", "y_strainer",
})


def validate_fitting_type(fitting_type: str) -> str:
    """Return a normalized supported fitting name without evaluating its correlation."""

    if not isinstance(fitting_type, str) or not fitting_type.strip():
        normalized = None
    else:
        normalized = fitting_type.lower().strip()
    if normalized not in VALID_FITTING_TYPES:
        raise FittingKError(
            "UNKNOWN_FITTING_TYPE",
            f"Unknown fitting type: {fitting_type!r}. Valid names include: "
            f"{', '.join(sorted(VALID_FITTING_TYPES)[:20])}...; see help_resources for the full list. "
            "Alternatively, provide a custom K_value in the fitting dict.",
            details={
                "fitting_type": fitting_type,
                "valid_fitting_types": sorted(VALID_FITTING_TYPES),
            },
        )
    return normalized


def resolve_fitting_k_basis(
    fitting: dict[str, Any],
    *,
    diameter: float,
    reynolds_number: float,
    flow_rate: float,
    darcy_friction_factor: float,
    zero_flow: bool = False,
) -> tuple[float, dict[str, Any]]:
    """Resolve one fitting from explicit K, equivalent length, or a named correlation.

    Equivalent-length inputs follow the Darcy-Weisbach/Crane TP-410 identity
    ``K = f_D (L_e/D)`` where ``f_D`` is the Darcy friction factor. Accepted
    fields are ``L_over_D`` (also ``L/D`` or ``l_over_d``) and
    ``equivalent_length_m`` (also ``Le_m`` or ``le_m``). Exactly one loss basis
    may be supplied; named fittings are evaluated by :func:`get_fitting_K` only
    when no explicit K or equivalent-length basis is present.

    No default K or equivalent length is fabricated. Invalid or ambiguous
    inputs raise :class:`FittingKError` with a public-safe typed error.
    """

    if not isinstance(fitting, dict):
        raise FittingKError(
            "INVALID_FITTING_INPUT",
            "Each fitting must be an object containing a named type, K_value, L_over_D, or equivalent_length_m.",
            details={"fitting_type": type(fitting).__name__},
        )

    ld_keys = [key for key in ("L_over_D", "L/D", "l_over_d") if fitting.get(key) is not None]
    le_keys = [
        key
        for key in ("equivalent_length_m", "Le_m", "le_m")
        if fitting.get(key) is not None
    ]
    explicit_k = fitting.get("K_value")
    bases = [
        name
        for name, present in (
            ("K_value", explicit_k is not None),
            ("L_over_D", bool(ld_keys)),
            ("equivalent_length_m", bool(le_keys)),
        )
        if present
    ]
    if len(ld_keys) > 1 or len(le_keys) > 1 or len(bases) > 1:
        raise FittingKError(
            "AMBIGUOUS_FITTING_BASIS",
            "Provide exactly one fitting loss basis: K_value, L_over_D, equivalent_length_m, or a named type.",
            details={"provided_bases": bases, "L_over_D_fields": ld_keys, "equivalent_length_fields": le_keys},
        )

    if fitting.get("type"):
        fitting_label = str(fitting["type"])
    elif explicit_k is not None:
        fitting_label = "custom_k"
    elif ld_keys or le_keys:
        fitting_label = "equivalent_length"
    else:
        fitting_label = "unnamed"
    if not math.isfinite(float(diameter)) or diameter <= 0:
        raise FittingKError(
            "INVALID_FITTING_INPUT",
            "Pipe diameter must be positive and finite when resolving fitting losses.",
            details={"diameter_m": diameter, "fitting_type": fitting_label},
        )
    if not math.isfinite(float(darcy_friction_factor)) or darcy_friction_factor < 0:
        raise FittingKError(
            "INVALID_FITTING_INPUT",
            "Darcy friction factor must be non-negative and finite when resolving fitting losses.",
            details={"darcy_friction_factor": darcy_friction_factor, "fitting_type": fitting_label},
        )

    detail: dict[str, Any] = {"type": fitting_label}
    if explicit_k is not None:
        k_value = float(explicit_k)
        source = "Provided K_value"
    elif ld_keys:
        l_over_d = float(fitting[ld_keys[0]])
        if not math.isfinite(l_over_d) or l_over_d < 0:
            raise FittingKError(
                "INVALID_FITTING_INPUT",
                "L_over_D must be non-negative and finite.",
                details={"field": ld_keys[0], "value": fitting[ld_keys[0]], "fitting_type": fitting_label},
            )
        k_value = 0.0 if zero_flow else darcy_friction_factor * l_over_d
        source = "Equivalent length L/D"
        detail.update({"L_over_D": l_over_d, "equivalent_length_m": l_over_d * diameter})
    elif le_keys:
        equivalent_length_m = float(fitting[le_keys[0]])
        if not math.isfinite(equivalent_length_m) or equivalent_length_m < 0:
            raise FittingKError(
                "INVALID_FITTING_INPUT",
                "equivalent_length_m must be non-negative and finite.",
                details={"field": le_keys[0], "value": fitting[le_keys[0]], "fitting_type": fitting_label},
            )
        l_over_d = equivalent_length_m / diameter
        k_value = 0.0 if zero_flow else darcy_friction_factor * l_over_d
        source = "Equivalent length Le"
        detail.update({"L_over_D": l_over_d, "equivalent_length_m": equivalent_length_m})
    else:
        normalized = validate_fitting_type(fitting.get("type"))
        k_value = (
            0.0
            if zero_flow
            else get_fitting_K(normalized, diameter, reynolds_number, flow_rate)
        )
        source = "Named Crane correlation"
        detail["type"] = normalized

    if not math.isfinite(k_value) or k_value < 0:
        raise FittingKError(
            "INVALID_FITTING_INPUT",
            "Resolved fitting K must be non-negative and finite.",
            details={"K_value": k_value, "fitting_type": fitting_label, "source": source},
        )
    detail.update({"K_source": source, "K_individual": k_value})
    return k_value, detail


def normalize_material_key(material: str) -> str:
    """Case and separator-insensitive material lookup key."""

    return "".join(ch for ch in material.lower() if ch.isalnum())


def resolve_pipe_material(material: str = None):
    """Return the canonical fluids.friction._roughness key for a material alias."""

    if not material:
        return None
    from fluids.friction import _roughness

    normalized = normalize_material_key(material)
    exact = {normalize_material_key(key): key for key in _roughness.keys()}
    if normalized in exact:
        return exact[normalized]

    aliases = {normalize_material_key(alias): canonical for alias, canonical in MATERIAL_ALIASES.items()}
    canonical = aliases.get(normalized)
    if canonical in _roughness:
        return canonical
    return None


def pipe_material_aliases_by_canonical() -> dict:
    """Return public material aliases grouped by canonical roughness key."""

    grouped = {}
    for alias, canonical in MATERIAL_ALIASES.items():
        grouped.setdefault(canonical, []).append(alias)
    return {key: sorted(values) for key, values in sorted(grouped.items())}


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
        FittingKError: If fitting_type is not recognized or the library calculation fails

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
    fitting_type = validate_fitting_type(fitting_type)

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
        elif fitting_type in ["butterfly_valve_kv", "butterfly_valve_flow_coeff"]:
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

        # === DEFENSIVE DEFAULT CASE ===
        else:
            raise FittingKError(
                "FITTING_K_CALCULATION_FAILED",
                f"No K-value implementation is registered for validated fitting '{fitting_type}'. "
                "Provide an audited custom K_value in the fitting dict.",
                details={"fitting_type": fitting_type},
            )

    except FittingKError:
        raise
    except Exception as e:
        logger.exception("Error calculating K-value for fitting '%s'", fitting_type)
        raise FittingKError(
            "FITTING_K_CALCULATION_FAILED",
            f"Failed to calculate K-value for fitting '{fitting_type}' using the fluids library: {e}. "
            "Verify the fitting inputs or provide an audited custom K_value in the fitting dict.",
            details={
                "fitting_type": fitting_type,
                "diameter_m": diameter,
                "reynolds_number": Re,
                "flow_rate_m3_s": flow_rate,
            },
        ) from e

def get_pipe_roughness(
    material: str = None,
    pipe_roughness: float = None,
    warnings: list[str] = None,
) -> tuple:
    """Get pipe roughness based on material or default.

    Args:
        material: Optional pipe material name
        pipe_roughness: Optional explicit roughness value
        warnings: Optional result warnings list. Unknown-material fallback is only
            permitted when this list is supplied so the fallback is auditable.

    Returns:
        Tuple of (roughness value, source description)
    """
    if pipe_roughness is not None:
        return pipe_roughness, "Provided"

    if material is not None:
        try:
            from fluids.friction import _roughness
            actual_key = resolve_pipe_material(material)
            if actual_key is not None:
                return _roughness[actual_key], f"Material Lookup '{actual_key}'"
            from fluids_mcp.utils.constants import DEFAULT_ROUGHNESS
            message = (
                f"Unknown pipe material '{material}'; using fallback roughness "
                f"{DEFAULT_ROUGHNESS:.6g} m. Provide pipe_roughness explicitly for auditable service."
            )
            if warnings is None:
                raise PipeMaterialError(
                    "UNKNOWN_PIPE_MATERIAL",
                    f"Unknown pipe material '{material}'. Provide a recognized material name or "
                    "an explicit pipe_roughness value.",
                    details={"material": material},
                )
            logger.warning(message)
            warnings.append(message)
            return DEFAULT_ROUGHNESS, f"Fallback {DEFAULT_ROUGHNESS:.6g} m (Unknown Material '{material}')"
        except PipeMaterialError:
            raise
        except Exception as e:
            from fluids_mcp.utils.constants import DEFAULT_ROUGHNESS
            message = (
                f"Pipe roughness lookup failed for material '{material}': {e}; using fallback "
                f"roughness {DEFAULT_ROUGHNESS:.6g} m. Provide pipe_roughness explicitly for auditable service."
            )
            if warnings is None:
                raise PipeMaterialError(
                    "PIPE_ROUGHNESS_LOOKUP_FAILED",
                    f"Pipe roughness lookup failed for material '{material}': {e}. "
                    "Provide an explicit pipe_roughness value.",
                    details={"material": material},
                ) from e
            logger.warning(message)
            warnings.append(message)
            return DEFAULT_ROUGHNESS, f"Fallback {DEFAULT_ROUGHNESS:.6g} m (Lookup Failed for '{material}')"

    # Default case
    from fluids_mcp.utils.constants import DEFAULT_ROUGHNESS
    return DEFAULT_ROUGHNESS, "Default"
