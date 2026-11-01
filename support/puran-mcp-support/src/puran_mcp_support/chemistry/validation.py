"""Executable, leakage-aware physical-validation assessment."""
from __future__ import annotations
import math
from typing import Iterable, Mapping, Sequence
from puran_mcp_support.schema import FESP_CALIBRATION_OBSERVABLES, CalibrationValidationEvidence, ValidationCase, ValidationDataset, ValidationTier, validate_calibration_holdout_separation, validate_temporal_holdout_campaigns

def assess_calibration_validation_evidence(evidence: CalibrationValidationEvidence) -> dict[str, object]:
    """Explain the calibration/holdout release gate without a bare boolean."""
    blockers: list[str] = []
    try:
        validate_calibration_holdout_separation(calibration_dataset_ids=evidence.calibration_dataset_ids, holdout_dataset_ids=evidence.holdout_dataset_ids)
    except ValueError:
        blockers.append('calibration_and_holdout_dataset_separation_invalid')
    try:
        validate_temporal_holdout_campaigns(evidence.campaigns)
    except ValueError:
        blockers.append('temporal_holdout_campaign_separation_invalid')
    missing_summary = sorted(FESP_CALIBRATION_OBSERVABLES - set(evidence.synchronized_observables))
    if missing_summary:
        blockers.append('synchronized_observables_incomplete')
    campaign_observation_gaps = {row.dataset_id: sorted(FESP_CALIBRATION_OBSERVABLES - set(row.synchronized_observables)) for row in evidence.campaigns if FESP_CALIBRATION_OBSERVABLES - set(row.synchronized_observables)}
    if campaign_observation_gaps:
        blockers.append('campaign_synchronized_observables_incomplete')
    calibration_campaign_ids = {row.dataset_id for row in evidence.campaigns if row.role == 'calibration'}
    holdout_campaign_ids = {row.dataset_id for row in evidence.campaigns if row.role == 'holdout'}
    if calibration_campaign_ids != set(evidence.calibration_dataset_ids):
        blockers.append('calibration_campaign_ids_do_not_match_evidence')
    if holdout_campaign_ids != set(evidence.holdout_dataset_ids):
        blockers.append('holdout_campaign_ids_do_not_match_evidence')
    if not evidence.identifiability.passed:
        blockers.append('estimated_parameters_not_practically_identifiable')
    if not evidence.posterior_predictive_validation_passed:
        blockers.append('posterior_predictive_holdout_validation_failed')
    if not evidence.target_operating_domain_covered:
        blockers.append('target_operating_domain_not_covered')
    release_ready = not blockers and evidence.release_ready
    return {'schema_version': 'calibration_validation_readiness/v1', 'calibration_release_ready': release_ready, 'recommendation_eligible': False, 'recommendation_blockers': sorted(set(blockers) | {'calibration_readiness_alone_does_not_satisfy_recommendation_gate'}), 'required_synchronized_observables': sorted(FESP_CALIBRATION_OBSERVABLES), 'missing_synchronized_observables': missing_summary, 'campaign_observation_gaps': campaign_observation_gaps, 'calibration_dataset_ids': list(evidence.calibration_dataset_ids), 'holdout_dataset_ids': list(evidence.holdout_dataset_ids), 'campaigns': [{'dataset_id': row.dataset_id, 'site_id': row.site_id, 'role': row.role, 'start_time_utc': row.start_time_utc, 'end_time_utc': row.end_time_utc, 'sample_count': row.sample_count, 'source_hash': row.source_hash} for row in evidence.campaigns], 'identifiability_passed': evidence.identifiability.passed, 'posterior_predictive_validation_passed': evidence.posterior_predictive_validation_passed, 'target_operating_domain_covered': evidence.target_operating_domain_covered}

def assess_validation_suite(*, datasets: Iterable[ValidationDataset], cases: Iterable[ValidationCase], predictions: Mapping[str, float], required_tiers: Sequence[ValidationTier]=('analytical_limit', 'independent_equilibrium', 'published_experiment', 'plant_holdout')) -> dict[str, object]:
    """Compare predictions to truth without laundering fit data as validation.

    A case whose dataset also supplied any fitted parameter remains visible as
    a parameter-contact check, but it cannot satisfy independent tier coverage.
    A nonidentical-protocol model-form challenge is likewise visible but
    excluded. Missing and nonfinite predictions are explicit failed cases.
    Release readiness requires every requested tier to contain at least one
    independent passing case and every independent case in the submitted suite
    to pass.
    """
    dataset_rows = tuple(datasets)
    dataset_by_id = {dataset.dataset_id: dataset for dataset in dataset_rows}
    if len(dataset_by_id) != len(dataset_rows):
        raise ValueError('validation dataset IDs must be unique')
    case_rows = tuple(cases)
    if not case_rows:
        raise ValueError('at least one validation case is required')
    case_ids = [case.case_id for case in case_rows]
    if len(set(case_ids)) != len(case_ids):
        raise ValueError('validation case IDs must be unique')
    unknown = sorted({case.dataset_id for case in case_rows if case.dataset_id not in dataset_by_id})
    if unknown:
        raise ValueError('validation cases reference unknown datasets: ' + ', '.join(unknown))
    normalized_predictions: dict[str, float | None] = {}
    for case in case_rows:
        raw = predictions.get(case.case_id)
        normalized_predictions[case.case_id] = None if raw is None else float(raw)
    results: list[dict[str, object]] = []
    for case in case_rows:
        dataset = dataset_by_id[case.dataset_id]
        prediction = normalized_predictions[case.case_id]
        finite_prediction = prediction is not None and math.isfinite(prediction)
        tolerance = case.absolute_tolerance + case.relative_tolerance * abs(case.observed_value)
        absolute_error = abs(prediction - case.observed_value) if finite_prediction else None
        fitted_from_same_dataset = bool(case.parameter_ids_fitted_from_same_dataset)
        independent_evidence = case.evidence_role == 'validation' and (not fitted_from_same_dataset)
        passed = bool(finite_prediction and absolute_error <= tolerance)
        results.append({'case_id': case.case_id, 'dataset_id': case.dataset_id, 'tier': dataset.tier, 'observable_id': case.observable_id, 'unit': case.unit, 'observed_value': case.observed_value, 'predicted_value': prediction, 'absolute_error': absolute_error, 'acceptance_tolerance': tolerance, 'prediction_finite': finite_prediction, 'passed': passed, 'independent_evidence': independent_evidence, 'evidence_role': case.evidence_role, 'parameter_ids_fitted_from_same_dataset': list(case.parameter_ids_fitted_from_same_dataset)})
    requested = tuple(dict.fromkeys((str(tier) for tier in required_tiers)))
    tier_coverage: dict[str, dict[str, object]] = {}
    for tier in requested:
        independent = [row for row in results if row['tier'] == tier and row['independent_evidence']]
        tier_coverage[tier] = {'independent_case_count': len(independent), 'passing_independent_case_count': sum((bool(row['passed']) for row in independent)), 'passed': bool(independent) and all((bool(row['passed']) for row in independent))}
    blockers: list[str] = []
    for tier, coverage in tier_coverage.items():
        if not coverage['independent_case_count']:
            blockers.append(f'missing_independent_validation_tier:{tier}')
        elif not coverage['passed']:
            blockers.append(f'validation_tier_failed:{tier}')
    if any((not row['prediction_finite'] for row in results)):
        blockers.append('missing_or_nonfinite_validation_predictions')
    if any((row['parameter_ids_fitted_from_same_dataset'] for row in results)):
        blockers.append('parameter_source_cases_excluded_from_independent_evidence')
    if any((row['evidence_role'] == 'model_form_challenge' for row in results)):
        blockers.append('model_form_challenge_cases_excluded_from_independent_evidence')
    if any((row['independent_evidence'] and (not row['passed']) for row in results)):
        blockers.append('independent_validation_case_failed')
    return {'schema_version': 'physical_validation_suite/v1', 'case_count': len(results), 'cases': results, 'required_tiers': list(requested), 'tier_coverage': tier_coverage, 'physical_validation_passed': not any((blocker.startswith(('missing_independent_validation_tier:', 'validation_tier_failed:')) or blocker in {'missing_or_nonfinite_validation_predictions', 'independent_validation_case_failed'} for blocker in blockers)), 'recommendation_eligible': False, 'recommendation_blockers': sorted(set(blockers + ['validation_pass_does_not_replace_calibration_and_recommendation_gate']))}
__all__ = ['assess_calibration_validation_evidence', 'assess_validation_suite']
