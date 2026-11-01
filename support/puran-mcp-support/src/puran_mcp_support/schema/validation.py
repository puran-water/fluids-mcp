"""Independent-validation and parameter-contact contracts.

Cross-engine agreement establishes software parity, not physical validity.
This module keeps parameter sources, calibration data, experimental checks,
and temporally separated plant holdouts distinct and machine-readable.
"""
from __future__ import annotations
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Iterable, Literal, Mapping, Sequence
ValidationTier = Literal['analytical_limit', 'independent_equilibrium', 'published_experiment', 'plant_calibration', 'plant_holdout']
ValidationEvidenceRole = Literal['validation', 'parameter_contact', 'model_form_challenge']
CampaignRole = Literal['calibration', 'holdout']
FESP_CALIBRATION_OBSERVABLES = frozenset({'delivered_Fe', 'orthophosphate', 'total_P', 'dissolved_Fe_II', 'total_dissolved_Fe', 'sulfate', 'dissolved_sulfide', 'gas_H2S_and_flow', 'pH_and_TIC', 'solid_Fe_P_S', 'SRT_and_wasting'})

def _utc_timestamp(value: str) -> datetime:
    try:
        parsed = datetime.fromisoformat(str(value).replace('Z', '+00:00'))
    except ValueError as exc:
        raise ValueError(f'invalid ISO-8601 campaign timestamp {value!r}') from exc
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        raise ValueError('campaign timestamps must include a UTC offset')
    return parsed.astimezone(timezone.utc)

@dataclass(frozen=True)
class PlantObservationCampaign:
    """One immutable, time-bounded plant calibration or holdout campaign."""
    dataset_id: str
    site_id: str
    role: CampaignRole
    start_time_utc: str
    end_time_utc: str
    sample_count: int
    synchronized_observables: tuple[str, ...]
    source_hash: str

    def __post_init__(self) -> None:
        if not all((self.dataset_id, self.site_id, self.source_hash)):
            raise ValueError('campaign dataset, site, and source hash are required')
        if self.sample_count < 2:
            raise ValueError('campaign must contain at least two synchronized samples')
        if not self.synchronized_observables:
            raise ValueError('campaign synchronized observables are required')
        if len(set(self.synchronized_observables)) != len(self.synchronized_observables):
            raise ValueError('campaign synchronized observables contain duplicates')
        if _utc_timestamp(self.end_time_utc) <= _utc_timestamp(self.start_time_utc):
            raise ValueError('campaign end time must be later than start time')

def validate_temporal_holdout_campaigns(campaigns: Sequence[PlantObservationCampaign]) -> None:
    """Require later, nonoverlapping holdout campaigns for every site."""
    rows = tuple(campaigns)
    if not rows:
        raise ValueError('plant calibration and holdout campaigns are required')
    ids = [row.dataset_id for row in rows]
    if len(set(ids)) != len(ids):
        raise ValueError('plant observation campaign dataset IDs must be unique')
    for site_id in sorted({row.site_id for row in rows}):
        site_rows = [row for row in rows if row.site_id == site_id]
        calibration = [row for row in site_rows if row.role == 'calibration']
        holdout = [row for row in site_rows if row.role == 'holdout']
        if not calibration or not holdout:
            raise ValueError(f'site {site_id!r} requires both calibration and holdout campaigns')
        latest_calibration_end = max((_utc_timestamp(row.end_time_utc) for row in calibration))
        earliest_holdout_start = min((_utc_timestamp(row.start_time_utc) for row in holdout))
        if earliest_holdout_start <= latest_calibration_end:
            raise ValueError(f'site {site_id!r} holdout must start after all calibration data end')

@dataclass(frozen=True)
class ValidationDataset:
    dataset_id: str
    tier: ValidationTier
    source: str
    matrix_scope: str
    observables: tuple[str, ...]
    intended_use: str
    empirical: bool
    parameter_source_for: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not all((self.dataset_id, self.source, self.matrix_scope, self.intended_use)):
            raise ValueError('validation dataset identity, source, matrix, and use are required')
        if not self.observables:
            raise ValueError('validation dataset must identify at least one observable')

@dataclass(frozen=True)
class ValidationCase:
    """One scalar truth comparison with an explicit acceptance tolerance."""
    case_id: str
    dataset_id: str
    observable_id: str
    observed_value: float
    unit: str
    absolute_tolerance: float = 0.0
    relative_tolerance: float = 0.0
    parameter_ids_fitted_from_same_dataset: tuple[str, ...] = ()
    evidence_role: ValidationEvidenceRole = 'validation'

    def __post_init__(self) -> None:
        if not all((self.case_id, self.dataset_id, self.observable_id, self.unit)):
            raise ValueError('validation case identity, observable, and unit are required')
        values = (self.observed_value, self.absolute_tolerance, self.relative_tolerance)
        if any((not math.isfinite(float(value)) for value in values)):
            raise ValueError('validation values and tolerances must be finite')
        if self.absolute_tolerance < 0.0 or self.relative_tolerance < 0.0:
            raise ValueError('validation tolerances must be nonnegative')
        if len(set(self.parameter_ids_fitted_from_same_dataset)) != len(self.parameter_ids_fitted_from_same_dataset):
            raise ValueError('same-dataset fitted parameter IDs contain duplicates')
        if self.evidence_role == 'validation' and self.parameter_ids_fitted_from_same_dataset:
            raise ValueError('a same-dataset fitted case must be labeled parameter_contact, not validation')

@dataclass(frozen=True)
class ParameterContactAssessment:
    delivered_fe_mmol_l: float
    observed_p_removed_mg_l: float
    observed_p_per_fe_mol_ratio: float
    assumed_capacity_mol_p_per_mol_fe: float
    capacity_implied_p_removed_mg_l: float
    relative_capacity_gap: float
    parameter_contact: bool
    mechanism_identified: bool = False

    def public_dict(self) -> dict[str, float | bool]:
        return {'delivered_fe_mmol_l': self.delivered_fe_mmol_l, 'observed_p_removed_mg_l': self.observed_p_removed_mg_l, 'observed_p_per_fe_mol_ratio': self.observed_p_per_fe_mol_ratio, 'assumed_capacity_mol_p_per_mol_fe': self.assumed_capacity_mol_p_per_mol_fe, 'capacity_implied_p_removed_mg_l': self.capacity_implied_p_removed_mg_l, 'relative_capacity_gap': self.relative_capacity_gap, 'parameter_contact': self.parameter_contact, 'mechanism_identified': self.mechanism_identified}

@dataclass(frozen=True)
class IdentifiabilityReport:
    method: str
    estimated_parameter_ids: tuple[str, ...]
    practically_identifiable_parameter_ids: tuple[str, ...]
    sensitivity_rank: int
    confounded_parameter_groups: tuple[tuple[str, ...], ...] = ()

    @property
    def passed(self) -> bool:
        estimated = set(self.estimated_parameter_ids)
        identifiable = set(self.practically_identifiable_parameter_ids)
        return bool(self.method and estimated and (identifiable == estimated) and (self.sensitivity_rank >= len(estimated)) and (not self.confounded_parameter_groups))

@dataclass(frozen=True)
class CalibrationValidationEvidence:
    calibration_dataset_ids: tuple[str, ...]
    holdout_dataset_ids: tuple[str, ...]
    synchronized_observables: tuple[str, ...]
    campaigns: tuple[PlantObservationCampaign, ...]
    identifiability: IdentifiabilityReport
    posterior_predictive_validation_passed: bool
    target_operating_domain_covered: bool

    @property
    def release_ready(self) -> bool:
        try:
            validate_calibration_holdout_separation(calibration_dataset_ids=self.calibration_dataset_ids, holdout_dataset_ids=self.holdout_dataset_ids)
            validate_temporal_holdout_campaigns(self.campaigns)
        except ValueError:
            return False
        calibration_ids = {row.dataset_id for row in self.campaigns if row.role == 'calibration'}
        holdout_ids = {row.dataset_id for row in self.campaigns if row.role == 'holdout'}
        campaigns_complete = all((FESP_CALIBRATION_OBSERVABLES <= set(row.synchronized_observables) for row in self.campaigns))
        return FESP_CALIBRATION_OBSERVABLES <= set(self.synchronized_observables) and campaigns_complete and (calibration_ids == set(self.calibration_dataset_ids)) and (holdout_ids == set(self.holdout_dataset_ids)) and self.identifiability.passed and self.posterior_predictive_validation_passed and self.target_operating_domain_covered

def _vector_norm(values: Sequence[float]) -> float:
    return math.sqrt(sum((value * value for value in values)))

def assess_local_identifiability(*, parameter_ids: Sequence[str], sensitivity_rows: Iterable[Sequence[float]], rank_tolerance: float=1e-08, confounding_cosine_threshold: float=0.999) -> IdentifiabilityReport:
    """Screen a scaled local-sensitivity matrix for rank and confounding.

    Columns are normalized before a deterministic modified Gram--Schmidt rank
    test, so parameter units do not decide the answer. This software-side
    screen does not replace profile likelihoods, posterior diagnostics, or
    temporally separated plant validation.
    """
    ids = tuple((str(parameter_id) for parameter_id in parameter_ids))
    if not ids or any((not parameter_id for parameter_id in ids)):
        raise ValueError('parameter_ids must contain nonempty identifiers')
    if len(set(ids)) != len(ids):
        raise ValueError('parameter_ids contain duplicates')
    if not math.isfinite(rank_tolerance) or rank_tolerance <= 0:
        raise ValueError('rank_tolerance must be finite and positive')
    if not math.isfinite(confounding_cosine_threshold) or not 0 < confounding_cosine_threshold <= 1:
        raise ValueError('confounding cosine threshold must be in (0, 1]')
    rows = tuple((tuple((float(value) for value in row)) for row in sensitivity_rows))
    if not rows:
        raise ValueError('sensitivity matrix must contain at least one row')
    if any((len(row) != len(ids) for row in rows)):
        raise ValueError('sensitivity-matrix column count must match parameter_ids')
    if any((not math.isfinite(value) for row in rows for value in row)):
        raise ValueError('sensitivity matrix must contain only finite values')
    raw_columns = tuple((tuple((row[column_index] for row in rows)) for column_index in range(len(ids))))
    norms = tuple((_vector_norm(column) for column in raw_columns))
    normalized = tuple((tuple((value / norm for value in column)) if norm > rank_tolerance else tuple((0.0 for _ in column)) for column, norm in zip(raw_columns, norms)))
    orthonormal_basis: list[tuple[float, ...]] = []
    independent_ids: list[str] = []
    for parameter_id, column in zip(ids, normalized):
        residual = list(column)
        for basis_vector in orthonormal_basis:
            projection = sum((left * right for left, right in zip(residual, basis_vector)))
            residual = [value - projection * basis_value for value, basis_value in zip(residual, basis_vector)]
        residual_norm = _vector_norm(residual)
        if residual_norm > rank_tolerance:
            orthonormal_basis.append(tuple((value / residual_norm for value in residual)))
            independent_ids.append(parameter_id)
    confounded_pairs: list[tuple[str, ...]] = []
    for left_index, left_id in enumerate(ids):
        if norms[left_index] <= rank_tolerance:
            confounded_pairs.append((left_id, 'unobservable_zero_sensitivity'))
            continue
        for right_index in range(left_index + 1, len(ids)):
            if norms[right_index] <= rank_tolerance:
                continue
            cosine = abs(sum((left * right for left, right in zip(normalized[left_index], normalized[right_index]))))
            if cosine >= confounding_cosine_threshold:
                confounded_pairs.append((left_id, ids[right_index]))
    return IdentifiabilityReport(method='normalized_local_sensitivity_rank_screen', estimated_parameter_ids=ids, practically_identifiable_parameter_ids=tuple(independent_ids), sensitivity_rank=len(orthonormal_basis), confounded_parameter_groups=tuple(confounded_pairs))

def assess_hfo_capacity_parameter_contact(*, delivered_fe_mg_l: float, influent_p_mg_l: float, residual_p_mg_l: float, assumed_capacity_mol_p_per_mol_fe: float, relative_contact_tolerance: float=0.1) -> ParameterContactAssessment:
    """Identify when a net-removal match is effectively a capacity identity."""
    values = (delivered_fe_mg_l, influent_p_mg_l, residual_p_mg_l, assumed_capacity_mol_p_per_mol_fe, relative_contact_tolerance)
    if any((not math.isfinite(float(item)) for item in values)):
        raise ValueError('parameter-contact inputs must be finite')
    if delivered_fe_mg_l <= 0 or assumed_capacity_mol_p_per_mol_fe <= 0:
        raise ValueError('delivered Fe and assumed capacity must be positive')
    if influent_p_mg_l < residual_p_mg_l or residual_p_mg_l < 0:
        raise ValueError('phosphorus concentrations must imply nonnegative removal')
    if relative_contact_tolerance < 0:
        raise ValueError('relative contact tolerance must be nonnegative')
    fe_mmol_l = delivered_fe_mg_l / 55.845
    removed = influent_p_mg_l - residual_p_mg_l
    removed_p_mmol_l = removed / 30.973761998
    observed_ratio = removed_p_mmol_l / fe_mmol_l
    implied_removal = fe_mmol_l * assumed_capacity_mol_p_per_mol_fe * 30.973761998
    gap = abs(observed_ratio - assumed_capacity_mol_p_per_mol_fe) / assumed_capacity_mol_p_per_mol_fe
    return ParameterContactAssessment(delivered_fe_mmol_l=fe_mmol_l, observed_p_removed_mg_l=removed, observed_p_per_fe_mol_ratio=observed_ratio, assumed_capacity_mol_p_per_mol_fe=assumed_capacity_mol_p_per_mol_fe, capacity_implied_p_removed_mg_l=implied_removal, relative_capacity_gap=gap, parameter_contact=gap <= relative_contact_tolerance)

def validate_calibration_holdout_separation(*, calibration_dataset_ids: tuple[str, ...], holdout_dataset_ids: tuple[str, ...]) -> None:
    overlap = sorted(set(calibration_dataset_ids) & set(holdout_dataset_ids))
    if overlap:
        raise ValueError(f'calibration and holdout datasets overlap: {overlap}')
    if not calibration_dataset_ids or not holdout_dataset_ids:
        raise ValueError('both calibration and temporally separated holdout datasets are required')
HFO_VALIDATION_CORPUS: tuple[ValidationDataset, ...] = (ValidationDataset(dataset_id='usgs_phreeqc_hfo_surface_definition', tier='independent_equilibrium', source='https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqc3-html/phreeqc3-52.htm', matrix_scope='canonical Dzombak-Morel HFO surface definition', observables=('surface_site_inventory', 'surface_area', 'surface_complexation'), intended_use='thermodynamic-definition parity; not plant-mechanism validation', empirical=False, parameter_source_for=('hfo_dzombak_morel_preformed_v1',)), ValidationDataset(dataset_id='mao_yue_2016_in_situ_hfo', tier='published_experiment', source='https://doi.org/10.1038/srep35292; OpenAlex https://openalex.org/W2531332928', matrix_scope='simple inorganic phosphate matrix; preformed versus in-situ HFO', observables=('phosphate_time_series', 'formation_pH', 'effective_site_capacity'), intended_use='fresh-HFO prior and model-form challenge', empirical=True, parameter_source_for=('hfo_mao_yue_in_situ_effective_v1',)), ValidationDataset(dataset_id='ferrihydrite_crystallite_size_2013', tier='published_experiment', source='https://doi.org/10.1021/es401301z; OpenAlex https://openalex.org/W2163458761', matrix_scope='ferrihydrite crystallite-size phosphate adsorption experiments', observables=('phosphate_capacity', 'surface_area', 'multistage_adsorption_rate'), intended_use='independent morphology and rate validation', empirical=True), ValidationDataset(dataset_id='iron_sulfide_reoxidation_p_capture_2010', tier='published_experiment', source='https://doi.org/10.1016/j.watres.2010.03.023; OpenAlex https://openalex.org/W2094295724', matrix_scope='sewer sulfide-control iron entering activated sludge', observables=('phosphate_removed_per_iron', 'iron_history'), intended_use='challenge automatic aged-HFO classification after FeS reoxidation', empirical=True), ValidationDataset(dataset_id='aged_ferric_sludge_reactivity_2020', tier='published_experiment', source='https://doi.org/10.1016/j.watres.2020.116179; OpenAlex https://openalex.org/W3042085360', matrix_scope='aged ferric drinking-water sludge exposed to sulfide and phosphate', observables=('sulfide_uptake', 'phosphate_uptake', 'solid_age'), intended_use='phase-history and reoxidation validation', empirical=True), ValidationDataset(dataset_id='smp_competition_fresh_hfo_2016', tier='published_experiment', source='https://doi.org/10.2166/wst.2016.431; OpenAlex https://openalex.org/W2519116963', matrix_scope='MBR soluble microbial product competing with phosphate on fresh HFO', observables=('phosphate_adsorption', 'SMP', 'pH'), intended_use='ligand/surface-competition validation', empirical=True))

def validation_corpus_by_id() -> Mapping[str, ValidationDataset]:
    return {dataset.dataset_id: dataset for dataset in HFO_VALIDATION_CORPUS}
__all__ = ['CampaignRole', 'FESP_CALIBRATION_OBSERVABLES', 'HFO_VALIDATION_CORPUS', 'CalibrationValidationEvidence', 'IdentifiabilityReport', 'ParameterContactAssessment', 'PlantObservationCampaign', 'ValidationCase', 'ValidationDataset', 'ValidationEvidenceRole', 'ValidationTier', 'assess_local_identifiability', 'assess_hfo_capacity_parameter_contact', 'validate_calibration_holdout_separation', 'validate_temporal_holdout_campaigns', 'validation_corpus_by_id']
