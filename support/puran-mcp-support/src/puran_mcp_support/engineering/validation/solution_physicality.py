from __future__ import annotations
from dataclasses import dataclass
from typing import Any
NEGATIVE_STATE_TOLERANCE_MG_L = 1e-06
AUTOTROPH_VIABILITY_FLOOR_MG_L = 1.0
AEROBIC_DO_FLOOR_MG_L = 0.5
AMMONIA_BREAKTHROUGH_FRACTION = 0.8
MASS_BALANCE_ERROR_TOLERANCE_PCT = 5.0

@dataclass(frozen=True)
class PhysicalityViolation:
    """One reason a converged solve is not a usable physical solution."""
    code: str
    message: str
    detail: dict[str, Any]

@dataclass(frozen=True)
class PhysicalityVerdict:
    """Verdict of the post-solve physicality gate."""
    is_physical: bool
    violations: tuple[PhysicalityViolation, ...] = ()
    checks_run: tuple[str, ...] = ()

    @property
    def codes(self) -> tuple[str, ...]:
        return tuple((violation.code for violation in self.violations))

    @property
    def reason(self) -> str:
        """A single-line summary suitable for a status message."""
        if self.is_physical:
            return 'solution passed the post-solve physicality gate'
        return '; '.join((violation.message for violation in self.violations))

    def remediation(self) -> list[str]:
        """Concrete next actions, one per violation, for the caller."""
        steps: list[str] = []
        for violation in self.violations:
            steps.extend(_REMEDIATION_BY_CODE.get(violation.code, []))
        steps.append('Do not consume this solve as a design basis: it converged but is not a physical steady state.')
        return list(dict.fromkeys(steps))

    def as_dict(self) -> dict[str, Any]:
        return {'is_physical': self.is_physical, 'checks_run': list(self.checks_run), 'violations': [{'code': violation.code, 'message': violation.message, 'detail': violation.detail} for violation in self.violations]}
_REMEDIATION_BY_CODE: dict[str, list[str]] = {'negative_state': ['A state variable went negative: the iterate is not a physical solution. Re-solve from a cold initialization and add a non-negativity bound on the offending component before trusting any result from this basin.'], 'nitrifier_washout': ['Autotrophs washed out. Check that SRT is above the nitrifier washout SRT at the design temperature, that WAS flow is not over-wasting, and that aeration (KLa / DO setpoint) can hold DO above 2 mg/L in the aerobic zone.', 'If washout persists with a correct SRT and DO, the flowsheet cannot meet an ammonia guarantee as configured — escalate the design, do not re-solve.'], 'do_collapse': ['Dissolved oxygen collapsed in the aerobic zone. Raise the aeration capacity (KLa) or switch the aerobic unit to a DO-setpoint control so the solver sizes aeration to the oxygen demand instead of starving it.'], 'ammonia_breakthrough': ['Effluent ammonia is at the influent reduced-N level: the plant is not nitrifying. Verify SRT, aerobic volume, DO, and alkalinity sufficiency before re-solving.'], 'mass_balance_non_closure': ['Mass does not close across the solved flowsheet. Inspect the recycle tear streams and per-unit residuals; a converged-but-non-closing solve means the flowsheet equations do not conserve mass as built.']}

def _negative_components(components_mg_l: dict[str, Any], *, tolerance_mg_l: float) -> dict[str, float]:
    negatives: dict[str, float] = {}
    for name, raw in (components_mg_l or {}).items():
        try:
            value = float(raw)
        except (TypeError, ValueError):
            continue
        if value < -abs(tolerance_mg_l):
            negatives[name] = value
    return negatives

def _first_float(source: dict[str, Any] | None, *names: str) -> float | None:
    if not isinstance(source, dict):
        return None
    for name in names:
        if name in source and source[name] is not None:
            try:
                return float(source[name])
            except (TypeError, ValueError):
                continue
    return None

def check_asm2d_solution_physicality(*, effluent: dict[str, Any] | None=None, design_parameters_settled: dict[str, Any] | None=None, influent_components_mg_l: dict[str, Any] | None=None, mass_balance_error_pct: dict[str, Any] | None=None, nitrification_expected: bool=False, aerobic_do_mg_l: float | None=None, negative_state_tolerance_mg_l: float=NEGATIVE_STATE_TOLERANCE_MG_L, autotroph_floor_mg_l: float=AUTOTROPH_VIABILITY_FLOOR_MG_L, do_floor_mg_l: float=AEROBIC_DO_FLOOR_MG_L, ammonia_breakthrough_fraction: float=AMMONIA_BREAKTHROUGH_FRACTION, mass_balance_tolerance_pct: float=MASS_BALANCE_ERROR_TOLERANCE_PCT) -> PhysicalityVerdict:
    """Judge whether a converged ASM2D-family solve is a physical solution.

    ``nitrification_expected`` gates the biology checks: only an aerobic
    (ASM1/ASM2D) flowsheet is designed to hold a nitrifier population, so the
    washout / DO / ammonia checks are skipped for anaerobic and non-biological
    families where a zero autotroph population is the correct answer.

    ``aerobic_do_mg_l`` must be the AERATED ZONE's DO, which only the caller can
    identify. The DO check is skipped when it is absent rather than falling back
    to the effluent, whose DO is legitimately ~0 on a post-anoxic flowsheet.

    Returns a verdict; it never raises on a malformed payload — a check with no
    data to read is skipped rather than guessed at.
    """
    violations: list[PhysicalityViolation] = []
    checks_run: list[str] = []
    effluent = effluent or {}
    settled = design_parameters_settled or {}
    effluent_components = effluent.get('components_mg_L') or {}
    checks_run.append('negative_state')
    negatives = _negative_components(effluent_components, tolerance_mg_l=negative_state_tolerance_mg_l)
    settled_biomass = {name: value for name, value in settled.items() if name.startswith(('X_', 'MLSS', 'MLVSS')) and isinstance(value, (int, float))}
    negatives.update(_negative_components(settled_biomass, tolerance_mg_l=negative_state_tolerance_mg_l))
    if negatives:
        worst = min(negatives.items(), key=lambda item: item[1])
        violations.append(PhysicalityViolation(code='negative_state', message=f'negative state variable: {worst[0]} = {worst[1]:.4g} mg/L (tolerance {-abs(negative_state_tolerance_mg_l):.1g} mg/L)', detail={'negative_components_mg_L': negatives, 'tolerance_mg_L': negative_state_tolerance_mg_l}))
    if nitrification_expected:
        x_aut = _first_float(settled, 'X_AUT_mg_L', 'X_AUT')
        if x_aut is None:
            x_aut = _first_float(effluent_components, 'X_AUT')
        if x_aut is not None:
            checks_run.append('nitrifier_washout')
            if x_aut < autotroph_floor_mg_l:
                violations.append(PhysicalityViolation(code='nitrifier_washout', message=f'autotroph washout: X_AUT = {x_aut:.4g} mg/L is below the {autotroph_floor_mg_l:g} mg/L viability floor, so the flowsheet carries no nitrifier population despite an aerobic (nitrifying) design intent', detail={'X_AUT_mg_L': x_aut, 'viability_floor_mg_L': autotroph_floor_mg_l}))
        if aerobic_do_mg_l is not None:
            checks_run.append('do_collapse')
            if aerobic_do_mg_l < do_floor_mg_l:
                violations.append(PhysicalityViolation(code='do_collapse', message=f'dissolved-oxygen collapse: aerobic-zone DO = {aerobic_do_mg_l:.4g} mg/L is below the {do_floor_mg_l:g} mg/L floor (AOB oxygen half-saturation), so nitrification cannot proceed in this basin', detail={'aerobic_DO_mg_L': aerobic_do_mg_l, 'do_floor_mg_L': do_floor_mg_l}))
        effluent_nh4 = _first_float(effluent, 'NH4_N_mg_L', 'NH4_N')
        if effluent_nh4 is None:
            effluent_nh4 = _first_float(effluent_components, 'S_NH4')
        influent_reduced_n = _influent_reduced_n_mg_l(influent_components_mg_l)
        if effluent_nh4 is not None and influent_reduced_n:
            checks_run.append('ammonia_breakthrough')
            retained = effluent_nh4 / influent_reduced_n
            if retained >= ammonia_breakthrough_fraction:
                violations.append(PhysicalityViolation(code='ammonia_breakthrough', message=f'no nitrification: effluent NH4-N {effluent_nh4:.4g} mg/L is {retained:.0%} of the influent reduced N ({influent_reduced_n:.4g} mg/L), at or above the {ammonia_breakthrough_fraction:.0%} breakthrough threshold for an aerobic design', detail={'effluent_NH4_N_mg_L': effluent_nh4, 'influent_reduced_N_mg_L': influent_reduced_n, 'retained_fraction': retained, 'threshold': ammonia_breakthrough_fraction}))
    if mass_balance_error_pct:
        failed: dict[str, float] = {}
        for name, raw in mass_balance_error_pct.items():
            try:
                value = float(raw)
            except (TypeError, ValueError):
                continue
            if abs(value) > abs(mass_balance_tolerance_pct):
                failed[name] = value
        if failed or mass_balance_error_pct:
            checks_run.append('mass_balance_non_closure')
        if failed:
            violations.append(PhysicalityViolation(code='mass_balance_non_closure', message='mass balance does not close across the solved flowsheet: ' + ', '.join((f'{name} error {value:+.3g}%' for name, value in sorted(failed.items()))) + f' (tolerance +/-{mass_balance_tolerance_pct:g}%)', detail={'errors_pct': failed, 'tolerance_pct': mass_balance_tolerance_pct}))
    return PhysicalityVerdict(is_physical=not violations, violations=tuple(violations), checks_run=tuple(checks_run))

def _influent_reduced_n_mg_l(components_mg_l: dict[str, Any] | None) -> float | None:
    """Return influent reduced N (mg/L) available to nitrify.

    Reduced N = soluble ammonia + organic-bound N in biomass and particulates,
    matching `asm2d_tkn_kg_m3` minus its S_N2 closure term (S_N2 is ~0 at the
    plant boundary, and dissolved N2 is not nitrifiable anyway).
    """
    if not isinstance(components_mg_l, dict) or not components_mg_l:
        return None
    try:
        from puran_mcp_support.engineering.heuristics.asm2d_influent_fractionation import ASM2D_I_NBM, ASM2D_I_NXI, ASM2D_I_NXS
    except ImportError:
        return None

    def value(name: str) -> float:
        try:
            return float(components_mg_l.get(name, 0.0) or 0.0)
        except (TypeError, ValueError):
            return 0.0
    reduced_n = value('S_NH4')
    reduced_n += ASM2D_I_NBM * (value('X_H') + value('X_PAO') + value('X_AUT'))
    reduced_n += ASM2D_I_NXI * value('X_I')
    reduced_n += ASM2D_I_NXS * value('X_S')
    return reduced_n or None
