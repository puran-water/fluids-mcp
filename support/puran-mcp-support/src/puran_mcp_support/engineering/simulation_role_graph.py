from __future__ import annotations
from typing import Iterable
DESIGNED_GRAPH_FAMILY_ASM2D: str = 'asm2d'
DESIGNED_GRAPH_FAMILY_ADM1: str = 'adm1'
ASM2D_MLE_RECYCLE_REQUIRED_ROLE_VALUES: tuple[str, ...] = ('feed', 'mixer', 'pre_anoxic', 'aerobic', 'internal_recycle_splitter', 'sludge_splitter')
ASM2D_MBR_RETENTATE_RECYCLE_REQUIRED_ROLE_VALUES: tuple[str, ...] = ('feed', 'mixer', 'pre_anoxic', 'aerobic', 'mbr_separator', 'sludge_splitter', 'ras_pressure')
ASM2D_RECYCLE_SEPARATOR_ROLE_VALUES: tuple[str, ...] = ('mbr_separator', 'secondary_clarifier')
ADM1_TWO_STAGE_REQUIRED_ROLE_VALUES: tuple[str, ...] = ('feed', 'anaerobic_acidogenic', 'anaerobic_methanogenic')
ADM1_SINGLE_STAGE_REQUIRED_ROLE_VALUES: tuple[str, ...] = ('feed', 'anaerobic_methanogenic')
ADM1_OPTIONAL_SURFACE_ROLE_VALUES: tuple[str, ...] = ('anaerobic_gas_phase', 'anaerobic_recuperation', 'anaerobic_thickener')
ADM1_RECYCLE_REQUIRED_ROLE_VALUE: str = 'anaerobic_recuperation'

def _normalize(roles: Iterable[str | None]) -> set[str]:
    """Collapse an iterable of role tags (or a role->x mapping's keys) to a
    lowercase set of present canonical role strings, dropping ``None``/empty."""
    out: set[str] = set()
    for role in roles:
        if not role:
            continue
        out.add(str(role).strip().lower())
    return out

def role_set_forms_asm2d_recycle_family(roles: Iterable[str | None]) -> bool:
    """True when the present role tags constitute a complete ASM2D biological
    recycle family, independent of unit/topology names.

    This is the structural, name-free completeness rule ``validate_role_graph``
    enforces and the Stage-7b dispatch routes on. ``roles`` may be any iterable of
    role strings (or a ``role -> binding`` mapping, whose keys are the present
    roles). The property-package confirmation (that the units are actually
    ASM2D-backed) is enforced downstream in the solve pipeline, exactly as before.
    """
    tagged = _normalize(roles)
    if set(ASM2D_MBR_RETENTATE_RECYCLE_REQUIRED_ROLE_VALUES) <= tagged:
        return True
    return set(ASM2D_MLE_RECYCLE_REQUIRED_ROLE_VALUES) <= tagged and bool(set(ASM2D_RECYCLE_SEPARATOR_ROLE_VALUES) & tagged)

def role_set_forms_adm1_digestion_family(roles: Iterable[str | None]) -> bool:
    """True when the present role tags form an ADM1/mADM1 anaerobic digestion train
    (two-stage acidogenic->methanogenic, or a single combined methanogenic
    reactor), independent of unit/topology names. Roles-only classification."""
    tagged = _normalize(roles)
    if set(ADM1_TWO_STAGE_REQUIRED_ROLE_VALUES) <= tagged:
        return True
    return set(ADM1_SINGLE_STAGE_REQUIRED_ROLE_VALUES) <= tagged

def role_set_forms_adm1_recycle_family(roles: Iterable[str | None]) -> bool:
    """True when the present role tags form an ADM1 digestion train carrying a
    closed liquid recycle — the digestion role set PLUS the recuperation/return
    splitter. A bare recycle mixer alone is deliberately NOT sufficient."""
    tagged = _normalize(roles)
    if not role_set_forms_adm1_digestion_family(tagged):
        return False
    return ADM1_RECYCLE_REQUIRED_ROLE_VALUE in tagged

def designed_graph_family(roles: Iterable[str | None]) -> str | None:
    """Return which COMPLETE designed-graph biological family the role set forms.

    The single family-detection rule the watertap designed-graph compiler and
    ``validate_role_graph`` key on.

    Returns:
      * ``"asm2d"`` (:data:`DESIGNED_GRAPH_FAMILY_ASM2D`) when the role tags form a
        complete ASM2D biological recycle family (MLE/MBR recycle pattern).
      * ``"adm1"`` (:data:`DESIGNED_GRAPH_FAMILY_ADM1`) when they form a complete
        ADM1/mADM1 anaerobic digestion family — a single methanogenic reactor
        (UASB), a phase-separated acidogenic->methanogenic two-stage cascade, or
        either of those carrying a closed liquid recycle (EGSB / AnMBR / contact).
      * ``None`` when neither — an incomplete or non-biological role set — OR when
        BOTH families are simultaneously complete (an AMBIGUOUS mixed-family role
        set). In every ``None`` case the designed-graph route falls back to the
        existing per-template / named-preset path rather than routing to a single
        designed build.

    The ASM2D and ADM1 families are individually disjoint per train (an aerobic
    recycle train carries none of the anaerobic roles and vice-versa). But a single
    PLANT can carry BOTH — an aerobic mainline AND an anaerobic digester — and the
    role set passed here is frequently the UNION across every sizing row. A union
    that satisfies BOTH completeness rules must NOT be collapsed into one designed
    family: building it as a single ASM2D (or single ADM1) graph would weld two
    independent trains into one flowsheet. So a mixed/ambiguous union returns
    ``None`` (ambiguous) and the caller falls back to the per-template dispatch,
    which solves the aerobic and anaerobic trains as SEPARATE ``solve_flowsheet``
    calls. The single-family designed route is therefore taken only when EXACTLY one
    family is complete.
    """
    tagged = _normalize(roles)
    is_asm2d = role_set_forms_asm2d_recycle_family(tagged)
    is_adm1 = role_set_forms_adm1_digestion_family(tagged)
    if is_asm2d and is_adm1:
        return None
    if is_asm2d:
        return DESIGNED_GRAPH_FAMILY_ASM2D
    if is_adm1:
        return DESIGNED_GRAPH_FAMILY_ADM1
    return None
