from __future__ import annotations
from dataclasses import dataclass, field
import math
from typing import Any, Dict, List, Literal, Optional
from puran_mcp_support.engineering.models.stream_state import MASM2D_IDS, MADM1_IDS, MCAS_SOLUTES
Vocabulary = Literal['mASM2d', 'mADM1', 'MCAS', 'explicit', 'unknown']
_MAP: Dict[str, Dict[str, str]] = {'Cl-': {'mASM2d': 'S_Cl', 'mADM1': 'S_Cl', 'MCAS': 'Cl_-'}, 'SO4-2': {'mADM1': 'S_SO4', 'MCAS': 'SO4_2-'}, 'Ca2+': {'mASM2d': 'S_Ca', 'mADM1': 'S_Ca', 'MCAS': 'Ca_2+'}, 'Mg2+': {'mASM2d': 'S_Mg', 'mADM1': 'S_Mg', 'MCAS': 'Mg_2+'}, 'Na+': {'mASM2d': 'S_Na', 'mADM1': 'S_Na', 'MCAS': 'Na_+'}, 'K+': {'mASM2d': 'S_K', 'mADM1': 'S_K', 'MCAS': 'K_+'}, 'Fe2+': {'mADM1': 'S_Fe2', 'MCAS': 'Fe_2+'}, 'Fe3+': {'mADM1': 'S_Fe3', 'MCAS': 'Fe_3+'}, 'total_sulfide': {'mADM1': 'S_IS'}, 'H2S': {'mADM1': 'S_H2S'}, 'H2S(g)': {'mADM1': 'S_h2s_gas'}, 'NH4+': {'mASM2d': 'S_NH4', 'mADM1': 'S_IN', 'MCAS': 'NH4_+'}, 'NO3-': {'mASM2d': 'S_NO3', 'MCAS': 'NO3_-'}, 'HCO3-': {'mASM2d': 'S_IC', 'mADM1': 'S_IC', 'MCAS': 'HCO3_-'}, 'CH3COO-': {'mASM2d': 'S_A', 'mADM1': 'S_ac'}}
_DOSE_ONLY = {'free_chlorine', 'ClO2', 'O3', 'peracetic_acid'}
_SUPPORTED_DOSE_CHEMICALS = frozenset({'fecl3', 'ferric', 'ferric_chloride', 'naocl', 'hypochlorite', 'sodium_hypochlorite', 'h2so4', 'sulfuric', 'sulphuric', 'hcl', 'hydrochloric'})
_BASIS = {'HCO3-': 'as C (S_IC)', 'total_sulfide': 'as S (S_IS, TOTAL inorganic sulfur)', 'H2S(g)': 'as S, GAS PHASE', 'NH4+': 'as N'}
_NEEDS_DOWNSTREAM_SPECIATION = {'total_sulfide': 'mADM1 S_IS is TOTAL inorganic sulfur as S, not the HS- fraction. Splitting it needs pH (pKa1 ~ 7.0, mid-band) - pass it to corrosion_screen_sulfide_exposure as total_sulfide_mg_L rather than treating it as an ionic concentration.', 'H2S(g)': "Gas-phase H2S is carried here in mg/L, the stream vocabulary's unit. corrosion_screen_sour_service needs ppmv AND total pressure, because partial pressure - not concentration - decides ISO 15156 scope. Convert before use; do not pass this number as ppmv."}

@dataclass
class CorrosionWaterRecord:
    """Species plus conditions, with every value carrying its own status."""
    species: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    temperature_C: Optional[float] = None
    pH: Optional[float] = None
    source_vocabulary: Vocabulary = 'unknown'
    charge_balance: str = 'not_enforced'
    charge_balance_error_pct: Optional[float] = None
    notes: List[str] = field(default_factory=list)

    def value(self, name: str) -> Optional[float]:
        """Concentration if usable, else None. Never 0.0-for-absent."""
        s = self.species.get(name)
        return s.get('value') if s and s.get('value') is not None else None
    _NO_VALUE_STATUSES = frozenset({'not_reported', 'below_detection', 'not_reported_dose_only'})

    def gaps(self, required: List[str]) -> List[str]:
        """Required species that are absent or unknown, and so block a verdict."""
        out = []
        for n in required:
            s = self.species.get(n)
            if s is None or s.get('status') in self._NO_VALUE_STATUSES:
                out.append(n)
        return out

def detect_vocabulary(component_ids: List[str]) -> Vocabulary:
    """Identify the vocabulary by set membership, not by `meta.source_model`.

    `meta` is optional and frequently unset, and a stream that has been through a
    converter may carry a stale label. The component list cannot lie about itself.
    """
    ids = set(component_ids)
    for name, ref in (('mADM1', MADM1_IDS), ('mASM2d', MASM2D_IDS), ('MCAS', MCAS_SOLUTES)):
        ref_set = set(ref)
        if ids and len(ids & ref_set) / len(ref_set) > 0.5:
            return name
    return 'unknown'

def _lookup(state, component_id: str) -> Optional[float]:
    """Membership-checked lookup. The whole point of this module.

    `StreamState.concentration()` raises KeyError for an absent component
    (it used to return 0.0 — a measured zero manufactured from a gap). This
    adapter wants an absent species as `None` feeding a `not_reported` record,
    not an exception, so it keeps its own membership-checked lookup.
    """
    ids = getattr(state, 'component_ids', None) or []
    if component_id not in ids:
        return None
    values = getattr(state, 'values_mg_L', None) or []
    idx = list(ids).index(component_id)
    return values[idx] if idx < len(values) else None

def stream_state_to_corrosion_water(state, *, dose_overlay: Optional[List[Dict[str, Any]]]=None, estimate_missing: bool=False, water_type: Optional[str]=None) -> CorrosionWaterRecord:
    """Build a corrosion input record from a StreamState.

    Args:
        state: a StreamState (or anything with component_ids / values_mg_L).
        dose_overlay: chemical additions, each
            {"chemical": "fecl3"|"naocl"|..., "dose_mg_L": float, ...}.
            Converted to added species by stoichiometry.
        estimate_missing: opt-in only. When False (default) an untracked species
            stays `not_reported` rather than being correlated from TDS.
        water_type: required when estimate_missing is True.

    Returns:
        CorrosionWaterRecord. Absent species are present as entries with
        status `not_reported` and NO value — never as 0.0.
    """
    vocab = detect_vocabulary(getattr(state, 'component_ids', []) or [])
    rec = CorrosionWaterRecord(source_vocabulary=vocab)
    temp_k = getattr(state, 'temperature_k', None)
    if temp_k is not None:
        rec.temperature_C = temp_k - 273.15
    rec.pH = getattr(state, 'ph', None)
    for canonical, per_vocab in _MAP.items():
        component_id = per_vocab.get(vocab)
        if component_id is None:
            rec.species[canonical] = {'value': None, 'units': 'mg/L', 'status': 'not_reported', 'source': f'{canonical} has no representation in the {vocab} vocabulary'}
            continue
        raw = _lookup(state, component_id)
        if raw is None:
            rec.species[canonical] = {'value': None, 'units': 'mg/L', 'status': 'not_reported', 'source': f"{component_id} not present in this stream's component array"}
            continue
        entry: Dict[str, Any] = {'value': float(raw), 'units': 'mg/L', 'status': 'converted' if canonical in _BASIS else 'reported', 'source': f'{vocab} {component_id}'}
        if canonical in _BASIS:
            entry['basis'] = _BASIS[canonical]
        if canonical in _NEEDS_DOWNSTREAM_SPECIATION:
            entry['needs_downstream_conversion'] = _NEEDS_DOWNSTREAM_SPECIATION[canonical]
        rec.species[canonical] = entry
    for name in _DOSE_ONLY:
        rec.species.setdefault(name, {'value': None, 'units': 'mg/L', 'status': 'not_reported_dose_only', 'source': 'dosing output, not a stream state; supply via dose_overlay'})
    if dose_overlay:
        _apply_dose_overlay(rec, dose_overlay)
    if estimate_missing:
        if not water_type:
            raise ValueError('WATER_TYPE_REQUIRED: estimate_missing=True needs an explicit water_type so the correlation used can be named')
        raise NotImplementedError('ESTIMATION_NOT_IMPLEMENTED: estimate_missing=True has no correlations behind it. It previously accepted the flag, logged a note promising estimated values, and left every species `not_reported` - a caller would believe the gaps were filled. Supply the missing species explicitly, or read `gaps()` and decide what to measure. Filling them with invented correlations is not a service this converter provides.')
    rec.notes.append('Chloride was NOT used to close any charge balance. masm2d_to_mcas._enforce_charge_balance adds Cl- to close a cation excess and hardcodes SO4 = 0.0; chaining through it would manufacture the governing corrosion variable from a residual.')
    return rec

def _apply_dose_overlay(rec: CorrosionWaterRecord, overlay: List[Dict[str, Any]]) -> None:
    """Convert chemical additions to added species by stoichiometry."""
    M_FE, M_CL = (55.845, 35.453)
    for add in overlay:
        if not isinstance(add, dict):
            raise ValueError('Dose overlay entries must be mappings')
        raw_chemical = add.get('chemical')
        chem = raw_chemical.strip().lower() if isinstance(raw_chemical, str) else ''
        if chem not in _SUPPORTED_DOSE_CHEMICALS:
            raise ValueError('Unknown dose-overlay chemical')
        raw_dose = add.get('dose_mg_L')
        if raw_dose is None or isinstance(raw_dose, bool):
            raise ValueError('dose_mg_L is required and must be numeric')
        try:
            dose = float(raw_dose)
        except (TypeError, ValueError) as exc:
            raise ValueError('dose_mg_L must be numeric') from exc
        if not math.isfinite(dose) or dose < 0.0:
            raise ValueError('dose_mg_L must be finite and non-negative')
        if chem in {'fecl3', 'ferric', 'ferric_chloride'}:
            fe = dose * (M_FE / 162.2)
            cl = fe * (3 * M_CL / M_FE)
            _add(rec, 'Cl-', cl, f'FeCl3 dosing at {dose:g} mg/L (FeCl3 -> Fe3+ + 3Cl-)')
            _add(rec, 'Fe3+', fe, f'FeCl3 dosing at {dose:g} mg/L')
            rec.notes.append('Ferric dosing adds BOTH chloride and an oxidant. Fe3+ is the ASTM G48 reagent — do not treat this as a chloride-only addition.')
        elif chem in {'naocl', 'hypochlorite', 'sodium_hypochlorite'}:
            raw_residual = add.get('residual_mg_L_as_Cl2')
            residual = None
            if raw_residual is not None:
                if isinstance(raw_residual, bool):
                    raise ValueError('INVALID_RESIDUAL_MAGNITUDE: residual_mg_L_as_Cl2 must be numeric, finite, and non-negative')
                try:
                    residual = float(raw_residual)
                except (TypeError, ValueError) as exc:
                    raise ValueError('INVALID_RESIDUAL_MAGNITUDE: residual_mg_L_as_Cl2 must be numeric, finite, and non-negative') from exc
                if not math.isfinite(residual) or residual < 0.0:
                    raise ValueError('INVALID_RESIDUAL_MAGNITUDE: residual_mg_L_as_Cl2 must be finite and non-negative')
            rec.species['free_chlorine'] = {'value': residual, 'units': 'mg/L', 'status': 'reported' if residual is not None else 'not_reported', 'basis': 'as Cl2', 'source': f'NaOCl dosing at {dose:g} mg/L as Cl2' + (f', residual {residual:g} mg/L' if residual is not None else '; RESIDUAL NOT SUPPLIED — dose is not residual, demand consumes most of it, and the residual is what acts on the metal')}
        elif chem in {'h2so4', 'sulfuric', 'sulphuric'}:
            _add(rec, 'SO4-2', dose * (96.06 / 98.08), f'H2SO4 dosing at {dose:g} mg/L')
        elif chem in {'hcl', 'hydrochloric'}:
            _add(rec, 'Cl-', dose * (M_CL / 36.46), f'HCl dosing at {dose:g} mg/L')

def _add(rec: CorrosionWaterRecord, canonical: str, delta: float, why: str) -> None:
    """Add a dosed contribution, preserving whether the base was known."""
    if not math.isfinite(delta):
        raise ValueError(f'NON_FINITE_ARITHMETIC: dose conversion produced a non-finite {canonical} contribution; no species value was emitted')
    cur = rec.species.get(canonical)
    if cur and cur.get('value') is not None:
        current = cur['value']
        combined = current + delta
        if not math.isfinite(combined):
            raise ValueError(f'NON_FINITE_ARITHMETIC: finite background and dose contribution overflowed while adding {canonical}; no species value was emitted')
        cur['value'] = combined
        cur['source'] = f'{cur['source']} + {why}'
    else:
        rec.species[canonical] = {'value': delta, 'units': 'mg/L', 'status': 'estimated', 'source': f'{why} (LOWER BOUND: background {canonical} not reported)'}
