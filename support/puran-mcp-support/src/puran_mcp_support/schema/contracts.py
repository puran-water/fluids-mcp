"""Safety, observability, provenance, activity-domain, and port accounting contracts."""
from __future__ import annotations
import math
from dataclasses import asdict, dataclass, field
from typing import Literal, Mapping, Optional
PRODUCT_DOSE_UNIT_FACTORS_TO_KG_PER_DAY: Mapping[str, Optional[float]] = {'kg_product_d': 1.0, 'kg_product_h': 24.0, 'l_product_d': None, 'l_product_h': None, 'us_gal_product_d': None, 'us_gal_product_h': None}

def normalized_product_dose_unit(unit: str) -> str:
    """Return the canonical unit ID for a supported as-delivered dose."""
    aliases = {'kg/d': 'kg_product_d', 'kg/day': 'kg_product_d', 'kg_product_day': 'kg_product_d', 'kg/h': 'kg_product_h', 'kg/hr': 'kg_product_h', 'l/d': 'l_product_d', 'l/day': 'l_product_d', 'l/h': 'l_product_h', 'l/hr': 'l_product_h', 'gal/d': 'us_gal_product_d', 'gpd': 'us_gal_product_d', 'gal_product_d': 'us_gal_product_d', 'gal/h': 'us_gal_product_h', 'gal/hr': 'us_gal_product_h', 'gal_product_h': 'us_gal_product_h'}
    normalized = str(unit).strip().lower()
    normalized = aliases.get(normalized, normalized)
    if normalized not in PRODUCT_DOSE_UNIT_FACTORS_TO_KG_PER_DAY:
        raise ValueError(f'unsupported as-delivered product dose unit {unit!r}')
    return normalized

def product_dose_is_volume_based(unit: str) -> bool:
    """Return whether a supported product-dose unit requires product density."""
    normalized = normalized_product_dose_unit(unit)
    return normalized.startswith(('l_', 'us_gal_'))

@dataclass(frozen=True)
class ActivityModelDomain:
    activity_model_id: str
    maximum_ionic_strength_molal: float
    minimum_temperature_c: float
    maximum_temperature_c: float

    def check(self, *, ionic_strength_molal: float, temperature_c: float, species_coverage_fraction: float) -> dict[str, object]:
        in_domain = math.isfinite(ionic_strength_molal) and math.isfinite(temperature_c) and math.isfinite(species_coverage_fraction) and (0 <= ionic_strength_molal <= self.maximum_ionic_strength_molal) and (self.minimum_temperature_c <= temperature_c <= self.maximum_temperature_c) and (0 <= species_coverage_fraction <= 1.0) and (species_coverage_fraction >= 1.0)
        return {'activity_model_id': self.activity_model_id, 'calculated_ionic_strength': ionic_strength_molal, 'activity_model_validity_limit': self.maximum_ionic_strength_molal, 'temperature_c': temperature_c, 'temperature_validity': [self.minimum_temperature_c, self.maximum_temperature_c], 'species_coverage_fraction': species_coverage_fraction, 'activity_model_in_domain': in_domain}

@dataclass(frozen=True)
class ObservabilityRequirements:
    required_measurements: tuple[str, ...]

    def missing(self, observations: Mapping[str, object]) -> tuple[str, ...]:
        return tuple((name for name in self.required_measurements if name not in observations or observations[name] is None))
FERRIC_RECOMMENDATION_OBSERVABILITY = ObservabilityRequirements(('flow', 'ferric_product_assay', 'ferric_delivery_calibration', 'pH', 'temperature', 'TIC_or_alkalinity', 'influent_TP', 'target_TP', 'orthophosphate', 'sulfate', 'dissolved_sulfide', 'gas_flow', 'gas_H2S', 'dissolved_Fe_II', 'total_dissolved_Fe', 'SRT', 'sludge_wasting', 'organic_particulate_P_characterization', 'SMP_or_DOC', 'permeate_total_Fe', 'solid_Fe_P_S_characterization', 'pipe_or_deposit_Fe_inventory'))

@dataclass(frozen=True)
class UncertainQuantity:
    mean: float
    standard_uncertainty: float

    def __post_init__(self) -> None:
        if not math.isfinite(float(self.mean)):
            raise ValueError('quantity mean must be finite')
        if not math.isfinite(float(self.standard_uncertainty)) or self.standard_uncertainty < 0:
            raise ValueError('standard_uncertainty must be finite and nonnegative')

@dataclass(frozen=True)
class OxideSurfaceState:
    """Runtime oxide-history and surface-capacity identification contract."""
    history_class: Literal['fresh_in_situ', 'preformed_aged', 'fes_reoxidized', 'plant_specific_effective', 'unresolved']
    surface_dataset_id: str
    parameter_interpretation: str
    formation_ph: Optional[float]
    formation_ph_basis: str
    active_metal_fraction: Optional[UncertainQuantity] = None
    mechanism_identified: bool = False
    mechanism_evidence_id: Optional[str] = None

    def __post_init__(self) -> None:
        if not self.surface_dataset_id or not self.parameter_interpretation:
            raise ValueError('surface dataset and parameter interpretation are required')
        if self.formation_ph is not None and (not math.isfinite(self.formation_ph)):
            raise ValueError('formation pH must be finite')
        if self.history_class == 'fresh_in_situ':
            if self.formation_ph is None or self.formation_ph_basis == 'unresolved':
                raise ValueError('fresh in-situ oxide requires a resolved formation pH')
        if self.active_metal_fraction is not None and (not 0 <= self.active_metal_fraction.mean <= 1):
            raise ValueError('active metal fraction must be in [0, 1]')
        if self.mechanism_identified and (not self.mechanism_evidence_id):
            raise ValueError('identified oxide mechanism requires an evidence identifier')

@dataclass(frozen=True)
class LigandCompetitionDefinition:
    """Declared dissolved-organic/surface-competition treatment."""
    ligand_model_id: str
    observation_basis: str
    calibrated: bool
    source: str
    smp_or_doc: Optional[UncertainQuantity] = None

    def __post_init__(self) -> None:
        if not all((self.ligand_model_id, self.observation_basis, self.source)):
            raise ValueError('ligand model, observation basis, and source are required')
        if self.calibrated and self.smp_or_doc is None:
            raise ValueError('a calibrated ligand model requires an SMP/DOC observation')

@dataclass(frozen=True)
class CommercialProductDefinition:
    product_id: str
    lot_id: Optional[str]
    assay_mass_fraction: UncertainQuantity
    density_kg_l: Optional[UncertainQuantity]
    free_acidity_eq_per_kg: Optional[UncertainQuantity]
    free_acidity_status: str
    assay_basis: Literal['active_reagent_mass_fraction'] = 'active_reagent_mass_fraction'
    free_acidity_counterion_stoichiometry: Mapping[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.product_id.strip():
            raise ValueError('commercial product ID is required')
        if not 0 < self.assay_mass_fraction.mean <= 1:
            raise ValueError('assay mass fraction must be in (0, 1]')
        if self.density_kg_l is not None and self.density_kg_l.mean <= 0:
            raise ValueError('product density must be positive')
        if self.free_acidity_status not in {'measured', 'not_measured', 'not_applicable'}:
            raise ValueError('free_acidity_status is invalid')
        if self.free_acidity_status == 'measured' and self.free_acidity_eq_per_kg is None:
            raise ValueError('measured free acidity requires a value and uncertainty')
        if self.free_acidity_eq_per_kg is not None and self.free_acidity_eq_per_kg.mean < 0:
            raise ValueError('free acidity must be nonnegative')
        if self.free_acidity_status == 'measured' and (not self.free_acidity_counterion_stoichiometry):
            raise ValueError('measured free acidity requires explicit counterion stoichiometry')
        if self.free_acidity_status != 'measured' and self.free_acidity_eq_per_kg is not None:
            raise ValueError('unmeasured or inapplicable free acidity must be null')
        if self.free_acidity_status != 'measured' and self.free_acidity_counterion_stoichiometry:
            raise ValueError('unmeasured or inapplicable free acidity cannot imply a counterion')
        if any((not component_id or not math.isfinite(float(coefficient)) or coefficient < 0 for component_id, coefficient in self.free_acidity_counterion_stoichiometry.items())):
            raise ValueError('free-acidity counterion coefficients must be finite and nonnegative')
        if self.assay_basis != 'active_reagent_mass_fraction':
            raise ValueError("the first shared dosing slice requires assay_basis='active_reagent_mass_fraction'")

@dataclass(frozen=True)
class ProductDose:
    value: float
    unit: str
    basis: str = 'as_delivered_product'

    def __post_init__(self) -> None:
        if not math.isfinite(float(self.value)) or self.value < 0:
            raise ValueError('dose must be finite and nonnegative')
        if not self.unit.strip():
            raise ValueError('dose unit is required')
        normalized_product_dose_unit(self.unit)
        if self.basis != 'as_delivered_product':
            raise ValueError("product dose basis must be 'as_delivered_product'")

@dataclass(frozen=True)
class InjectionZoneDefinition:
    mixing_model: str
    zone_volume_m3: float
    hydraulic_model: str
    dose_split: tuple[float, ...]
    injection_ports: tuple[str, ...]
    local_ph_model: str
    dispersion_or_bypass: Optional[float] = None

    def __post_init__(self) -> None:
        if not all((self.mixing_model, self.hydraulic_model, self.local_ph_model)):
            raise ValueError('mixing, hydraulic, and local-pH models are required')
        if self.zone_volume_m3 <= 0 or not math.isfinite(self.zone_volume_m3):
            raise ValueError('injection-zone volume must be positive')
        if not self.injection_ports or len(self.injection_ports) != len(self.dose_split):
            raise ValueError('each dose split requires one injection port')
        if any((value < 0 or not math.isfinite(value) for value in self.dose_split)):
            raise ValueError('dose split fractions must be finite and nonnegative')
        if abs(sum(self.dose_split) - 1.0) > 1e-09:
            raise ValueError('dose split fractions must sum to one')
        if self.dispersion_or_bypass is not None and (not (math.isfinite(self.dispersion_or_bypass) and 0 <= self.dispersion_or_bypass <= 1)):
            raise ValueError('dispersion_or_bypass must be a finite fraction in [0, 1]')
SerialZoneKind = Literal['injection_micromixing', 'anaerobic_contact', 'oxidation_contact', 'aerobic_contact', 'mbr_separation']

@dataclass(frozen=True)
class ReactiveZoneDefinition:
    """One ordered hydraulic/chemical zone in the Fe-S-P history model."""
    zone_id: str
    zone_kind: SerialZoneKind
    hydraulic_model: Literal['cstr', 'tanks_in_series']
    hydraulic_residence_time_d: float
    solids_retention_time_d: float
    volume_m3: float
    tanks_in_series: int = 1
    dose_fraction: float = 0.0
    oxygen_source_mol_m3_d: float = 0.0
    gas_boundary_mode: Literal['closed', 'finite_headspace', 'fixed_partial_pressure'] = 'closed'
    local_ph_model: Literal['fixed', 'relative_charge_balance', 'external_property_package'] = 'relative_charge_balance'
    particulate_rejection_fraction: float = 1.0

    def __post_init__(self) -> None:
        if not self.zone_id.strip():
            raise ValueError('serial chemistry zone_id is required')
        positive = (self.hydraulic_residence_time_d, self.solids_retention_time_d, self.volume_m3)
        if any((not math.isfinite(float(value)) or value <= 0.0 for value in positive)):
            raise ValueError('zone HRT, SRT, and volume must be finite and positive')
        if isinstance(self.tanks_in_series, bool) or int(self.tanks_in_series) != self.tanks_in_series:
            raise ValueError('tanks_in_series must be an integer')
        if self.tanks_in_series < 1:
            raise ValueError('tanks_in_series must be at least one')
        if self.hydraulic_model == 'cstr' and self.tanks_in_series != 1:
            raise ValueError('a CSTR zone must use tanks_in_series=1')
        bounded = (self.dose_fraction, self.particulate_rejection_fraction)
        if any((not math.isfinite(float(value)) or not 0.0 <= value <= 1.0 for value in bounded)):
            raise ValueError('dose and rejection fractions must be finite in [0, 1]')
        if not math.isfinite(float(self.oxygen_source_mol_m3_d)) or self.oxygen_source_mol_m3_d < 0.0:
            raise ValueError('zone oxygen source must be finite and nonnegative')
        if self.zone_kind not in {'oxidation_contact', 'aerobic_contact'} and self.oxygen_source_mol_m3_d:
            raise ValueError('oxygen source is valid only in oxidation/aerobic zones')
        if self.zone_kind != 'injection_micromixing' and self.dose_fraction:
            raise ValueError('dose_fraction is valid only in injection zones')
        if self.zone_kind != 'mbr_separation' and self.particulate_rejection_fraction != 1.0:
            raise ValueError('particulate rejection is valid only in the MBR zone')

@dataclass(frozen=True)
class SerialChemistryTrainDefinition:
    """Ordered injection/anaerobic/oxidation/aerobic/MBR history contract."""
    train_id: str
    zones: tuple[ReactiveZoneDefinition, ...]

    def __post_init__(self) -> None:
        if not self.train_id.strip() or not self.zones:
            raise ValueError('serial chemistry train ID and zones are required')
        zone_ids = [zone.zone_id for zone in self.zones]
        duplicates = sorted((zone_id for zone_id in set(zone_ids) if zone_ids.count(zone_id) > 1))
        if duplicates:
            raise ValueError(f'duplicate serial chemistry zone IDs: {duplicates}')
        ordering = {'injection_micromixing': 0, 'anaerobic_contact': 1, 'oxidation_contact': 2, 'aerobic_contact': 3, 'mbr_separation': 4}
        ranks = [ordering[zone.zone_kind] for zone in self.zones]
        if ranks != sorted(ranks):
            raise ValueError('serial chemistry zones are not in physical process order')
        present = {zone.zone_kind for zone in self.zones}
        required = set(ordering)
        if not required <= present:
            raise ValueError('serial chemistry train requires injection, anaerobic, oxidation, aerobic, and MBR zones')
        if self.zones[-1].zone_kind != 'mbr_separation':
            raise ValueError('the final serial chemistry zone must be MBR separation')
        dose_total = sum((zone.dose_fraction for zone in self.zones))
        if abs(dose_total - 1.0) > 1e-09:
            raise ValueError('serial injection-zone dose fractions must sum to one')
        hydraulic_flows = tuple((zone.volume_m3 / zone.hydraulic_residence_time_d for zone in self.zones))
        flow_scale = max(hydraulic_flows)
        if any((abs(flow - hydraulic_flows[0]) > 1e-06 * flow_scale for flow in hydraulic_flows[1:])):
            raise ValueError('serial chemistry zone volume/HRT pairs imply inconsistent hydraulic flows')

    @property
    def hydraulic_flow_m3_d(self) -> float:
        return self.zones[0].volume_m3 / self.zones[0].hydraulic_residence_time_d

def coerce_serial_chemistry_train(value: SerialChemistryTrainDefinition | Mapping[str, object]) -> SerialChemistryTrainDefinition:
    """Validate a JSON-shaped serial train at an engine/API boundary."""
    if isinstance(value, SerialChemistryTrainDefinition):
        return value
    if not isinstance(value, Mapping):
        raise TypeError('serial chemistry train must be a mapping')
    unknown = sorted(set(value) - {'train_id', 'zones'})
    if unknown:
        raise ValueError(f'serial chemistry train has unknown fields {unknown}')
    raw_zones = value.get('zones')
    if not isinstance(raw_zones, (list, tuple)):
        raise TypeError('serial chemistry train zones must be a list or tuple')
    zones: list[ReactiveZoneDefinition] = []
    for index, raw_zone in enumerate(raw_zones):
        if isinstance(raw_zone, ReactiveZoneDefinition):
            zones.append(raw_zone)
            continue
        if not isinstance(raw_zone, Mapping):
            raise TypeError(f'serial chemistry zone {index} must be a mapping')
        try:
            zones.append(ReactiveZoneDefinition(**dict(raw_zone)))
        except TypeError as exc:
            raise ValueError(f'invalid serial chemistry zone {index}: {exc}') from exc
    return SerialChemistryTrainDefinition(train_id=str(value.get('train_id') or ''), zones=tuple(zones))

@dataclass(frozen=True)
class ComplianceTarget:
    metric: str
    port: str
    basis: str
    value: float
    statistic: str
    window: str
    horizon: str

    def __post_init__(self) -> None:
        if self.metric not in {'orthophosphate', 'dissolved_p', 'total_p'}:
            raise ValueError('unsupported phosphorus target metric')
        if not all((self.port, self.basis, self.statistic, self.window, self.horizon)):
            raise ValueError('target port, basis, statistic, window, and horizon are required')
        if not math.isfinite(float(self.value)) or self.value < 0:
            raise ValueError('target value must be finite and nonnegative')

@dataclass(frozen=True)
class AssumptionLedger:
    target_variable: str
    thermodynamic_dataset_id: str
    database_hash: str
    phase_definition_set_id: str
    activity_model_id: str
    redox_method: str
    sulfide_source: str
    iron_sulfide_phase_case: str
    chemistry_treatment: str
    solids_inventory_id: Optional[str]
    mbr_retention_profile_id: Optional[str]
    calibration_domain_id: Optional[str]
    ligand_model_id: str
    model_commit: str
    input_hash: str
    surface_dataset_id: Optional[str] = None
    surface_parameter_interpretation: Optional[str] = None
    oxide_history_class: str = 'unresolved'
    oxide_history_id: str = 'unresolved'
    oxide_formation_ph: Optional[float] = None
    oxide_formation_ph_basis: str = 'unresolved'
    active_metal_fraction_basis: str = 'unresolved'
    effective_active_fe_fraction: Optional[UncertainQuantity] = None
    parameter_contact: Optional[bool] = None
    parameter_contact_evidence_id: Optional[str] = None
    mechanism_identified: bool = False
    mechanism_evidence_id: Optional[str] = None
    smp_or_doc_basis: str = 'unresolved'
    extrapolation_status: Literal['not_assessed', 'within_declared_domain', 'outside_declared_domain', 'not_applicable'] = 'not_assessed'
    extrapolation_reasons: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        required = (self.target_variable, self.thermodynamic_dataset_id, self.database_hash, self.phase_definition_set_id, self.activity_model_id, self.redox_method, self.sulfide_source, self.iron_sulfide_phase_case, self.chemistry_treatment, self.ligand_model_id, self.model_commit, self.input_hash, self.oxide_history_class, self.oxide_history_id, self.oxide_formation_ph_basis, self.active_metal_fraction_basis, self.smp_or_doc_basis)
        if any((not str(value).strip() for value in required)):
            raise ValueError('assumption-ledger identifiers and bases must be nonempty')
        if self.oxide_formation_ph is not None and (not math.isfinite(self.oxide_formation_ph)):
            raise ValueError('oxide formation pH must be finite')
        if self.effective_active_fe_fraction is not None and (not 0 <= self.effective_active_fe_fraction.mean <= 1):
            raise ValueError('effective active-Fe fraction must be in [0, 1]')
        if self.parameter_contact is True and (not self.parameter_contact_evidence_id):
            raise ValueError('claimed parameter contact requires an evidence identifier')
        if self.mechanism_identified and (not self.mechanism_evidence_id):
            raise ValueError('identified mechanism requires an evidence identifier')
        valid_extrapolation = {'not_assessed', 'within_declared_domain', 'outside_declared_domain', 'not_applicable'}
        if self.extrapolation_status not in valid_extrapolation:
            raise ValueError('invalid extrapolation status')
        if self.extrapolation_status == 'outside_declared_domain' and (not self.extrapolation_reasons):
            raise ValueError('out-of-domain ledger requires an extrapolation reason')

    def public_dict(self) -> dict[str, object]:
        """Return a JSON-ready ledger without dropping unresolved fields."""
        return asdict(self)

def mechanism_interpretation_from_ledger(ledger: AssumptionLedger | Mapping[str, object]) -> dict[str, object]:
    """Project one canonical oxide/mechanism record from an assumption ledger.

    Execution engines historically exposed this subset at different nesting
    levels and used both ``formation_ph_basis`` and
    ``oxide_formation_ph_basis``.  Keep the complete ledger as the source of
    truth and emit one stable, engine-neutral view for agent-facing results.
    """
    values = ledger.public_dict() if isinstance(ledger, AssumptionLedger) else dict(ledger)
    required = ('oxide_history_class', 'oxide_history_id', 'oxide_formation_ph_basis', 'active_metal_fraction_basis', 'mechanism_identified', 'extrapolation_status')
    missing = [key for key in required if key not in values]
    if missing:
        raise ValueError('assumption ledger is missing mechanism fields: ' + ', '.join(missing))
    return {'mechanism_identified': bool(values['mechanism_identified']), 'mechanism_evidence_id': values.get('mechanism_evidence_id'), 'surface_dataset_id': values.get('surface_dataset_id'), 'surface_parameter_interpretation': values.get('surface_parameter_interpretation'), 'oxide_history_class': str(values['oxide_history_class']), 'oxide_history_id': str(values['oxide_history_id']), 'formation_ph': values.get('oxide_formation_ph'), 'formation_ph_basis': str(values['oxide_formation_ph_basis']), 'active_metal_fraction_basis': str(values['active_metal_fraction_basis']), 'effective_active_fe_fraction': values.get('effective_active_fe_fraction'), 'parameter_contact': values.get('parameter_contact'), 'parameter_contact_evidence_id': values.get('parameter_contact_evidence_id'), 'extrapolation_status': str(values['extrapolation_status']), 'extrapolation_reasons': list(values.get('extrapolation_reasons') or ())}

@dataclass(frozen=True)
class SteadyStateBranchRecord:
    branch_id: str
    initial_inventory_id: str
    continuation_path: tuple[float, ...]
    dynamic_stability: Literal['stable', 'unstable', 'not_assessed']
    steady_residual_norm: float
    active_phase_set: tuple[str, ...]

    def __post_init__(self) -> None:
        if not self.branch_id or not self.initial_inventory_id:
            raise ValueError('branch and initial-inventory identifiers are required')
        if not self.continuation_path or any((not math.isfinite(value) for value in self.continuation_path)):
            raise ValueError('continuation path must contain finite values')
        if not math.isfinite(self.steady_residual_norm) or self.steady_residual_norm < 0:
            raise ValueError('steady residual norm must be finite and nonnegative')
        if len(set(self.active_phase_set)) != len(self.active_phase_set):
            raise ValueError('active phase set contains duplicate identifiers')

@dataclass(frozen=True)
class RecommendationGate:
    solver_converged: bool
    target_met: bool
    mass_balance_passed: bool
    bracket_resolved: bool
    activity_model_in_domain: bool
    observations_complete: bool
    calibrated: bool
    holdout_validation_passed: bool
    within_calibration_domain: bool
    ph_charge_solve_passed: bool
    model_disagreement_within_limit: bool
    commercial_product_characterized: bool
    target_contract_complete: bool
    injection_zone_defined: bool
    audit_translation_passed: bool
    hard_audit_findings_absent: bool
    uncertainty_quantified: bool
    surface_model_identified_or_bounded: bool = False
    ligand_competition_characterized: bool = False
    active_metal_fraction_identified: bool = False
    reasons: tuple[str, ...] = field(default_factory=tuple)

    @property
    def recommendation_eligible(self) -> bool:
        return all((self.solver_converged, self.target_met, self.mass_balance_passed, self.bracket_resolved, self.activity_model_in_domain, self.observations_complete, self.calibrated, self.holdout_validation_passed, self.within_calibration_domain, self.ph_charge_solve_passed, self.model_disagreement_within_limit, self.commercial_product_characterized, self.target_contract_complete, self.injection_zone_defined, self.audit_translation_passed, self.hard_audit_findings_absent, self.uncertainty_quantified, self.surface_model_identified_or_bounded, self.ligand_competition_characterized, self.active_metal_fraction_identified))

    def public_dict(self) -> dict[str, object]:
        return {'recommendation_eligible': self.recommendation_eligible, 'solver_converged': self.solver_converged, 'target_met': self.target_met, 'mass_balance_passed': self.mass_balance_passed, 'bracket_resolved': self.bracket_resolved, 'activity_model_in_domain': self.activity_model_in_domain, 'observations_complete': self.observations_complete, 'calibrated': self.calibrated, 'holdout_validation_passed': self.holdout_validation_passed, 'within_calibration_domain': self.within_calibration_domain, 'ph_charge_solve_passed': self.ph_charge_solve_passed, 'model_disagreement_within_limit': self.model_disagreement_within_limit, 'commercial_product_characterized': self.commercial_product_characterized, 'target_contract_complete': self.target_contract_complete, 'injection_zone_defined': self.injection_zone_defined, 'audit_translation_passed': self.audit_translation_passed, 'hard_audit_findings_absent': self.hard_audit_findings_absent, 'uncertainty_quantified': self.uncertainty_quantified, 'surface_model_identified_or_bounded': self.surface_model_identified_or_bounded, 'ligand_competition_characterized': self.ligand_competition_characterized, 'active_metal_fraction_identified': self.active_metal_fraction_identified, 'reasons': list(self.reasons)}

def commercial_product_is_characterized(product: Optional[CommercialProductDefinition], dose: Optional[ProductDose]) -> bool:
    """Return whether product-to-active-reagent conversion is fully observable."""
    if product is None or dose is None:
        return False
    if product.free_acidity_status == 'not_measured':
        return False
    try:
        volume_based = product_dose_is_volume_based(dose.unit)
    except ValueError:
        return False
    return not volume_based or product.density_kg_l is not None

def build_recommendation_gate(*, observations: Mapping[str, object], solver_converged: bool, target_met: bool, mass_balance_passed: bool, bracket_resolved: bool, activity_model_in_domain: bool, calibrated: bool, holdout_validation_passed: bool, within_calibration_domain: bool, ph_charge_solve_passed: bool, model_disagreement_within_limit: bool, product: Optional[CommercialProductDefinition], product_dose: Optional[ProductDose], target: Optional[ComplianceTarget], injection_zone: Optional[InjectionZoneDefinition], audit_translation_passed: bool, hard_audit_findings_absent: bool, uncertainty_quantified: bool, surface_model_identified_or_bounded: bool=False, ligand_competition_characterized: bool=False, active_metal_fraction_identified: bool=False) -> RecommendationGate:
    """Build a fail-closed gate and enumerate every failed condition."""
    missing = FERRIC_RECOMMENDATION_OBSERVABILITY.missing(observations)
    product_characterized = commercial_product_is_characterized(product, product_dose)
    conditions = {'solver_not_converged': solver_converged, 'target_not_met': target_met, 'mass_balance_failed': mass_balance_passed, 'dose_branch_unresolved': bracket_resolved, 'activity_model_out_of_domain': activity_model_in_domain, 'required_observations_missing': not missing, 'model_not_calibrated': calibrated, 'holdout_validation_failed': holdout_validation_passed, 'outside_calibration_domain': within_calibration_domain, 'ph_or_charge_solve_failed': ph_charge_solve_passed, 'model_disagreement_exceeds_limit': model_disagreement_within_limit, 'commercial_product_not_characterized': product_characterized, 'compliance_target_incomplete': target is not None, 'injection_zone_undefined': injection_zone is not None, 'audit_translation_failed': audit_translation_passed, 'hard_thermodynamic_audit_finding': hard_audit_findings_absent, 'uncertainty_not_quantified': uncertainty_quantified, 'surface_model_unresolved': surface_model_identified_or_bounded, 'ligand_competition_uncharacterized': ligand_competition_characterized, 'active_metal_fraction_unidentified': active_metal_fraction_identified}
    reasons = [code for code, passed in conditions.items() if not passed]
    if missing:
        reasons.append('missing_measurements:' + ','.join(missing))
    return RecommendationGate(solver_converged=solver_converged, target_met=target_met, mass_balance_passed=mass_balance_passed, bracket_resolved=bracket_resolved, activity_model_in_domain=activity_model_in_domain, observations_complete=not missing, calibrated=calibrated, holdout_validation_passed=holdout_validation_passed, within_calibration_domain=within_calibration_domain, ph_charge_solve_passed=ph_charge_solve_passed, model_disagreement_within_limit=model_disagreement_within_limit, commercial_product_characterized=product_characterized, target_contract_complete=target is not None, injection_zone_defined=injection_zone is not None, audit_translation_passed=audit_translation_passed, hard_audit_findings_absent=hard_audit_findings_absent, uncertainty_quantified=uncertainty_quantified, surface_model_identified_or_bounded=surface_model_identified_or_bounded, ligand_competition_characterized=ligand_competition_characterized, active_metal_fraction_identified=active_metal_fraction_identified, reasons=tuple(reasons))

@dataclass(frozen=True)
class SimulationProvenance:
    model_commit: str
    database_hash: str
    input_hash: str
    expression_hash: str
    site_calibrated: bool = False
    data_class: str = 'simulated'
    not_empirical: bool = True
__all__ = ['ActivityModelDomain', 'AssumptionLedger', 'CommercialProductDefinition', 'ComplianceTarget', 'FERRIC_RECOMMENDATION_OBSERVABILITY', 'InjectionZoneDefinition', 'ReactiveZoneDefinition', 'SerialChemistryTrainDefinition', 'SerialZoneKind', 'LigandCompetitionDefinition', 'ObservabilityRequirements', 'OxideSurfaceState', 'ProductDose', 'PRODUCT_DOSE_UNIT_FACTORS_TO_KG_PER_DAY', 'RecommendationGate', 'SimulationProvenance', 'SteadyStateBranchRecord', 'UncertainQuantity', 'build_recommendation_gate', 'coerce_serial_chemistry_train', 'commercial_product_is_characterized', 'mechanism_interpretation_from_ledger', 'normalized_product_dose_unit', 'product_dose_is_volume_based']
