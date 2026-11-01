from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Callable, Mapping, Optional
P_MW = 30.973762
FE_OH3_MW = 106.867
FECL3_MW = 162.204
CACO3_EQ_MW = 50.04
STOICH_FEOH3_EQ_PER_P_1TO1 = 3.45
FECL3_PER_FEOH3_EQ = FECL3_MW / FE_OH3_MW
ALK_CACO3_PER_MOL_FE = 3.0 * CACO3_EQ_MW
PH_CRASH_THRESHOLD = 6.0

@dataclass
class ChemicalPDoseSeed:
    """Structured metal-dose convergence prior for an ASM2D chemical-P unit.

    ``metal_dose_kg_feoh3_eq_m3`` is the engine-basis seed (the value to set on
    the freed ``metal_dose`` Var). The ``fecl3_*`` / ``alkalinity_*`` fields are
    the real-reagent reporting numbers; they are populated from PHREEQC when a
    result is available and from the 1:1 stoichiometric floor otherwise.
    """
    metal_dose_kg_feoh3_eq_m3: float
    delta_p_removed_mg_l: float
    post_bio_soluble_p_mg_l: float
    target_p_mg_l: float
    dose_lower: float
    dose_upper: float
    clipped_to_bound: bool
    reagent: str = 'FeCl3'
    fe_to_p_molar: float = 1.0
    throughput_m3_d: Optional[float] = None
    feoh3_eq_kg_per_d: Optional[float] = None
    fecl3_kg_per_d_stoich_1to1: Optional[float] = None
    fecl3_kg_per_d: Optional[float] = None
    achievable_residual_p_mg_l: Optional[float] = None
    final_ph: Optional[float] = None
    ph_crash: bool = False
    alkalinity_consumed_mg_caco3_l: Optional[float] = None
    alkalinity_codose_kg_caco3eq_per_d: Optional[float] = None
    precipitated_phases: Optional[Mapping[str, float]] = None
    provenance: str = 'stoichiometric_floor'
    notes: list[str] = field(default_factory=list)
    screen_qualification: Optional[Mapping[str, Any]] = None

    def as_dict(self) -> dict[str, Any]:
        return {'metal_dose_kg_feoh3_eq_m3': self.metal_dose_kg_feoh3_eq_m3, 'delta_p_removed_mg_l': self.delta_p_removed_mg_l, 'post_bio_soluble_p_mg_l': self.post_bio_soluble_p_mg_l, 'target_p_mg_l': self.target_p_mg_l, 'dose_bounds': [self.dose_lower, self.dose_upper], 'clipped_to_bound': self.clipped_to_bound, 'reagent': self.reagent, 'fe_to_p_molar': self.fe_to_p_molar, 'throughput_m3_d': self.throughput_m3_d, 'screen_qualification': self.screen_qualification, 'feoh3_eq_kg_per_d': self.feoh3_eq_kg_per_d, 'fecl3_kg_per_d_stoich_1to1': self.fecl3_kg_per_d_stoich_1to1, 'fecl3_kg_per_d': self.fecl3_kg_per_d, 'achievable_residual_p_mg_l': self.achievable_residual_p_mg_l, 'final_ph': self.final_ph, 'ph_crash': self.ph_crash, 'alkalinity_consumed_mg_caco3_l': self.alkalinity_consumed_mg_caco3_l, 'alkalinity_codose_kg_caco3eq_per_d': self.alkalinity_codose_kg_caco3eq_per_d, 'precipitated_phases': dict(self.precipitated_phases) if self.precipitated_phases else None, 'provenance': self.provenance, 'notes': list(self.notes)}

def engine_metal_dose_seed(soluble_p_mg_l: float, target_p_mg_l: float, *, dose_lower: float=0.0, dose_upper: float=0.5) -> float:
    """Engine-basis ``metal_dose`` seed (kg Fe(OH)3-equivalent / m^3).

    The dose that, applied to a stream carrying ``soluble_p_mg_l`` of S_PO4,
    drops it to ``target_p_mg_l`` in one pass through the 1:1 ASM2D chemical-P
    surrogate. Clipped to ``[dose_lower, dose_upper]``. This is the value to set
    on the freed ``metal_dose`` Var so the stiff recycle starts at the
    precipitation knee instead of drifting toward it.
    """
    delta_p_kg_m3 = max(float(soluble_p_mg_l) - float(target_p_mg_l), 0.0) / 1000.0
    seed = STOICH_FEOH3_EQ_PER_P_1TO1 * delta_p_kg_m3
    return min(max(seed, float(dose_lower)), float(dose_upper))

def parse_phreeqc_dose_result(result: Mapping[str, Any]) -> dict[str, Any]:
    """Normalize a ``calculate_phosphorus_removal_dose`` result into seed inputs.

    Tolerant of partial payloads; returns ``{}`` when the result is missing or
    not a success. Keys mirror the water-chemistry MCP tool's response.
    """
    if not isinstance(result, Mapping) or not result:
        return {}
    status = str(result.get('status', '')).lower()
    if status and status not in {'success', 'ok', 'converged', 'infeasible', 'converged_target_met'}:
        return {}
    out: dict[str, Any] = {}
    for key in ('metal_to_p_ratio', 'optimal_dose_mmol', 'optimal_dose_mg_l', 'achieved_p_mg_l', 'final_ph'):
        value = result.get(key)
        if value is not None:
            try:
                out[key] = float(value)
            except (TypeError, ValueError):
                pass
    phases = result.get('precipitated_phases')
    if isinstance(phases, Mapping):
        out['precipitated_phases'] = {str(k): v for k, v in phases.items()}
    if not out:
        return {}
    out['status'] = status or 'unknown'
    qualification = result.get('qualification')
    if isinstance(qualification, Mapping):
        out['qualification'] = dict(qualification)
    return out

def chemical_p_dose_seed(post_bio_soluble_p_mg_l: float, target_tp_mg_l: float, *, water_matrix: Optional[Mapping[str, Any]]=None, reagent: str='FeCl3', flow_m3_d: Optional[float]=None, throughput_m3_d: Optional[float]=None, phreeqc_result: Optional[Mapping[str, Any]]=None, phreeqc_dose_provider: Optional[Callable[..., Mapping[str, Any]]]=None, dose_lower: float=0.0, dose_upper: float=0.5) -> ChemicalPDoseSeed:
    """Full metal-dose convergence prior: engine seed + real reagent feed.

    Parameters
    ----------
    post_bio_soluble_p_mg_l:
        Soluble ortho-P entering the chemical-P unit (the recycle-coupled
        post-biological value when known; the open-loop post-bio value
        otherwise). Sets ΔP for the engine seed.
    target_tp_mg_l:
        Target effluent total P. Used directly as the soluble-P target for the
        engine seed (permeate X_MeP is rejected, so effluent TP ≈ permeate
        soluble P).
    water_matrix:
        ``{pH, alkalinity (mg/L CaCO3 or "alkalinity_mg_caco3_l"), Ca, Mg,
        temperature_c}`` — passed to ``phreeqc_dose_provider`` and used for the
        alkalinity co-dose estimate. Optional.
    reagent:
        Metal-salt reagent. FeCl3 is the only basis with a delivered-mass
        conversion here; other reagents still get the engine seed (Fe(OH)3-eq
        basis) and a note.
    flow_m3_d / throughput_m3_d:
        Flow the dose acts on for the per-day reagent numbers. ``throughput``
        (mixed-liquor flow through the chemical-P unit) is preferred; falls back
        to ``flow_m3_d``.
    phreeqc_result:
        A pre-computed ``calculate_phosphorus_removal_dose`` result dict. Wins
        over ``phreeqc_dose_provider``.
    phreeqc_dose_provider:
        ``provider(post_bio_soluble_p_mg_l, target_tp_mg_l, water_matrix,
        reagent) -> result dict``. Called only when ``phreeqc_result`` is absent.
        Any exception degrades to the stoichiometric floor.
    dose_lower / dose_upper:
        Engine dose bounds; the seed is clipped to this window.
    """
    post_bio = float(post_bio_soluble_p_mg_l)
    target = float(target_tp_mg_l)
    delta_p_mg_l = max(post_bio - target, 0.0)
    raw_seed = STOICH_FEOH3_EQ_PER_P_1TO1 * (delta_p_mg_l / 1000.0)
    seed = min(max(raw_seed, float(dose_lower)), float(dose_upper))
    clipped = seed != raw_seed
    flow = throughput_m3_d if throughput_m3_d is not None else flow_m3_d
    notes: list[str] = []
    phreeqc: dict[str, Any] = {}
    if phreeqc_result is not None:
        phreeqc = parse_phreeqc_dose_result(phreeqc_result)
    elif phreeqc_dose_provider is not None:
        try:
            raw = phreeqc_dose_provider(post_bio_soluble_p_mg_l=post_bio, target_tp_mg_l=target, water_matrix=dict(water_matrix) if water_matrix else None, reagent=reagent)
            phreeqc = parse_phreeqc_dose_result(raw)
        except Exception as exc:
            notes.append(f'PHREEQC provider failed ({exc!r}); using stoichiometric floor')
    out = ChemicalPDoseSeed(metal_dose_kg_feoh3_eq_m3=seed, delta_p_removed_mg_l=round(delta_p_mg_l, 4), post_bio_soluble_p_mg_l=post_bio, target_p_mg_l=target, dose_lower=float(dose_lower), dose_upper=float(dose_upper), clipped_to_bound=clipped, reagent=reagent, throughput_m3_d=flow, notes=notes)
    if delta_p_mg_l <= 0.0:
        out.fe_to_p_molar = 0.0
        out.provenance = 'no_dose_required'
        out.notes.append('post-bio soluble P already at/below target; zero dose seed')
        _fill_reagent_numbers(out, flow)
        return out
    if phreeqc:
        out.provenance = 'water_chemistry_phreeqc'
        qualification = phreeqc.get('qualification')
        if qualification is not None:
            out.screen_qualification = qualification
            authorizes = qualification.get('authorizes') or {}
            if not any(authorizes.values()):
                out.notes.append(f"the dose is an equilibrium screen's ({qualification.get('standing', 'equilibrium_screen')}'): it seeds the engine and the reagent figures are screening figures; it authorizes no reactor size, guaranteed dose or PO clearance (review M17, PUBLIC-RULE-8226e460e1d8)'")
        out.fe_to_p_molar = float(phreeqc.get('metal_to_p_ratio', 1.0) or 1.0)
        out.final_ph = phreeqc.get('final_ph')
        out.achievable_residual_p_mg_l = phreeqc.get('achieved_p_mg_l')
        out.precipitated_phases = phreeqc.get('precipitated_phases')
        if out.final_ph is not None and out.final_ph < PH_CRASH_THRESHOLD:
            out.ph_crash = True
            out.notes.append(f'ferric acidity crashes pH to {out.final_ph:.2f}; caustic/lime co-dose mandatory')
        achieved = out.achievable_residual_p_mg_l
        if achieved is not None and achieved > target + 1e-06:
            out.notes.append(f'target {target:.2f} mg/L is below PHREEQC-achievable residual {achieved:.2f} mg/L at Fe:P={out.fe_to_p_molar:.2f}; engine 1:1 surrogate still pins the target (idealized complete precipitation)')
        dose_mg_l = phreeqc.get('optimal_dose_mg_l')
        if dose_mg_l is not None and flow is not None and (reagent.upper() == 'FECL3'):
            out.fecl3_kg_per_d = round(float(dose_mg_l) * flow / 1000.0, 2)
    else:
        out.provenance = 'stoichiometric_floor'
        out.fe_to_p_molar = 1.0
        given = phreeqc_result if phreeqc_result is not None else None
        reason = f'the water-chemistry result carried no usable converged dose (status {str(given.get('status', 'unknown'))})' if isinstance(given, Mapping) and given else 'water-chemistry unavailable'
        out.notes.append(f'{reason}; stoichiometric-floor seed at molar Fe:P=1 (sane prior, verify the real Fe:P excess + alkalinity co-dose with calculate_phosphorus_removal_dose)')
    _fill_reagent_numbers(out, flow)
    _fill_alkalinity_codose(out, water_matrix, flow)
    return out

def _fill_reagent_numbers(out: ChemicalPDoseSeed, flow: Optional[float]) -> None:
    """Populate per-day reagent masses from the engine dose, flow and Fe:P."""
    if flow is None:
        return
    feoh3_eq_kg_d = out.metal_dose_kg_feoh3_eq_m3 * flow
    out.feoh3_eq_kg_per_d = round(feoh3_eq_kg_d, 2)
    if out.reagent.upper() == 'FECL3':
        stoich_1to1 = feoh3_eq_kg_d * FECL3_PER_FEOH3_EQ
        out.fecl3_kg_per_d_stoich_1to1 = round(stoich_1to1, 2)
        if out.fecl3_kg_per_d is None:
            out.fecl3_kg_per_d = round(stoich_1to1 * max(out.fe_to_p_molar, 1e-09), 2)
    else:
        out.notes.append(f'reagent {out.reagent!r}: only Fe(OH)3-eq dose computed; no delivered-mass conversion for non-FeCl3 reagents')

def _fill_alkalinity_codose(out: ChemicalPDoseSeed, water_matrix: Optional[Mapping[str, Any]], flow: Optional[float]) -> None:
    """Estimate alkalinity consumed by the (real, excess) ferric dose."""
    delta_p_kg_m3 = out.delta_p_removed_mg_l / 1000.0
    mol_p_per_m3 = delta_p_kg_m3 / (P_MW / 1000.0)
    mol_fe_per_m3 = mol_p_per_m3 * max(out.fe_to_p_molar, 0.0)
    consumed_mg_caco3_l = mol_fe_per_m3 * ALK_CACO3_PER_MOL_FE
    out.alkalinity_consumed_mg_caco3_l = round(consumed_mg_caco3_l, 1)
    if water_matrix:
        alk = water_matrix.get('alkalinity_mg_caco3_l') or water_matrix.get('alkalinity') or water_matrix.get('alkalinity_mg_l_as_caco3')
        if alk is not None:
            try:
                alk_val = float(alk)
                deficit = consumed_mg_caco3_l - alk_val
                if deficit > 0:
                    out.notes.append(f'ferric consumes ~{consumed_mg_caco3_l:.0f} mg/L CaCO3 vs ~{alk_val:.0f} mg/L available; ~{deficit:.0f} mg/L CaCO3-eq alkalinity co-dose needed to hold pH')
            except (TypeError, ValueError):
                pass
    if flow is not None:
        out.alkalinity_codose_kg_caco3eq_per_d = round(consumed_mg_caco3_l * flow / 1000.0, 2)
