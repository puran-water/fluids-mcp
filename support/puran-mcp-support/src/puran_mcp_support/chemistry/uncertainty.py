"""Branch-aware uncertainty summaries for fixed-dose forward ensembles."""
from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Callable, Iterable, Mapping, Sequence
from puran_mcp_support.schema import RecommendationGate

@dataclass(frozen=True)
class DoseSample:
    sample_id: str
    branch_id: str
    dose: float
    solver_converged: bool
    mass_balance_passed: bool
    target_met: bool
    branch_resolved: bool
    within_calibration_domain: bool

    def __post_init__(self) -> None:
        if not self.sample_id:
            raise ValueError('dose sample sample_id is required')
        if not self.branch_id:
            raise ValueError('dose sample branch_id is required')
        if not math.isfinite(float(self.dose)) or self.dose < 0:
            raise ValueError('dose sample must be finite and nonnegative')

@dataclass(frozen=True)
class UncertaintyRealization:
    """One reproducible uncertain-input realization for a dose scan."""
    realization_id: str
    parameters: Mapping[str, float]

    def __post_init__(self) -> None:
        if not self.realization_id:
            raise ValueError('realization_id is required')
        normalized = {str(key): float(value) for key, value in self.parameters.items()}
        if not normalized:
            raise ValueError('uncertainty realization parameters are required')
        invalid = sorted((key for key, value in normalized.items() if not math.isfinite(value)))
        if invalid:
            raise ValueError('uncertainty realization contains nonfinite parameters: ' + ', '.join(invalid))
        object.__setattr__(self, 'parameters', normalized)

def _quantile(values: list[float], probability: float) -> float:
    if not values:
        raise ValueError('cannot calculate a quantile of an empty sample')
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction

def build_branch_dose_envelope(samples: Iterable[DoseSample], *, random_seed: int, input_hash: str, model_hash: str, recommendation_gate: RecommendationGate | None=None) -> dict[str, object]:
    """Summarize feasible dose samples without pooling distinct branches."""
    rows = tuple(samples)
    if not rows:
        raise ValueError('at least one dose sample is required')
    if not input_hash or not model_hash:
        raise ValueError('input_hash and model_hash are required for reproducibility')
    branch_ids = sorted({row.branch_id for row in rows})
    branches: list[dict[str, object]] = []
    for branch_id in branch_ids:
        branch_rows = [row for row in rows if row.branch_id == branch_id]
        converged = [row for row in branch_rows if row.solver_converged]
        balanced = [row for row in converged if row.mass_balance_passed]
        resolved = [row for row in balanced if row.branch_resolved]
        resolved_target = [row.dose for row in resolved if row.target_met]
        in_domain = [row for row in resolved if row.within_calibration_domain]
        feasible = [row.dose for row in in_domain if row.target_met]
        branches.append({'branch_id': branch_id, 'sample_count': len(branch_rows), 'converged_count': len(converged), 'mass_balance_passed_count': len(balanced), 'branch_resolved_count': len(resolved), 'within_calibration_domain_count': len(in_domain), 'resolved_target_crossing_count': len(resolved_target), 'feasible_count': len(feasible), 'solver_success_probability': len(converged) / len(branch_rows), 'mass_balance_success_probability': len(balanced) / len(branch_rows), 'branch_resolution_probability': len(resolved) / len(branch_rows), 'calibration_domain_probability': len(in_domain) / len(branch_rows), 'target_probability_given_sample': len(feasible) / len(branch_rows), 'screening_target_probability_given_sample': len(resolved_target) / len(branch_rows), 'screening_dose_d10': _quantile(resolved_target, 0.1) if resolved_target else None, 'screening_dose_d50': _quantile(resolved_target, 0.5) if resolved_target else None, 'screening_dose_d90': _quantile(resolved_target, 0.9) if resolved_target else None, 'dose_d10': _quantile(feasible, 0.1) if feasible else None, 'dose_d50': _quantile(feasible, 0.5) if feasible else None, 'dose_d90': _quantile(feasible, 0.9) if feasible else None})
    feasible_branches = [branch for branch in branches if branch['feasible_count']]
    screening_branches = [branch for branch in branches if branch['resolved_target_crossing_count']]
    single_resolved_branch = len(feasible_branches) == 1
    gate_eligible = bool(recommendation_gate is not None and recommendation_gate.recommendation_eligible)
    recommendation_eligible = gate_eligible and single_resolved_branch
    blockers: list[str] = []
    if recommendation_gate is None:
        blockers.append('recommendation_gate_not_supplied')
    elif not recommendation_gate.recommendation_eligible:
        blockers.extend(recommendation_gate.reasons)
    if not feasible_branches:
        blockers.append('no_feasible_dose_samples')
        if screening_branches:
            blockers.append('screening_target_crossings_outside_validated_calibration_domain')
    elif not single_resolved_branch:
        blockers.append('multiple_feasible_branches_not_resolved')
    return {'schema_version': 'branch_dose_envelope/v1', 'random_seed': int(random_seed), 'input_hash': input_hash, 'model_hash': model_hash, 'branch_count': len(branches), 'branches': branches, 'selected_branch_id': feasible_branches[0]['branch_id'] if single_resolved_branch else None, 'selected_screening_branch_id': screening_branches[0]['branch_id'] if len(screening_branches) == 1 else None, 'recommendation_eligible': recommendation_eligible, 'recommendation_blockers': sorted(set(blockers))}

def run_branch_aware_dose_grid_ensemble(realizations: Iterable[UncertaintyRealization], *, dose_grid: Sequence[float], evaluator: Callable[[float, Mapping[str, float]], Mapping[str, object]], random_seed: int, input_hash: str, model_hash: str, recommendation_gate: RecommendationGate | None=None) -> dict[str, object]:
    """Execute fixed-dose scans and summarize only resolved local crossings.

    The evaluator must return ``branch_id``, ``solver_converged``,
    ``mass_balance_passed``, ``target_met`` and
    ``within_calibration_domain``. Persistent errors split a branch segment;
    they never imply a dose direction. A target crossing is usable only when
    adjacent successful points on the same branch bracket false-to-true target
    attainment and the segment is monotonic in that Boolean target response.
    The reported dose is the conservative upper grid bound, not an interpolated
    or optimized setpoint.
    """
    rows = tuple(realizations)
    if not rows:
        raise ValueError('at least one uncertainty realization is required')
    ids = [row.realization_id for row in rows]
    if len(set(ids)) != len(ids):
        raise ValueError('uncertainty realization IDs must be unique')
    doses = tuple(sorted({float(dose) for dose in dose_grid}))
    if not doses or any((not math.isfinite(dose) or dose < 0.0 for dose in doses)):
        raise ValueError('dose_grid must contain finite nonnegative doses')
    if doses[0] != 0.0:
        raise ValueError('dose_grid must start at zero to resolve a zero-dose target')
    evaluation_rows: list[dict[str, object]] = []
    samples: list[DoseSample] = []
    ensemble_blockers: set[str] = set()
    for realization in rows:
        scan: list[dict[str, object]] = []
        for dose in doses:
            try:
                output = dict(evaluator(dose, realization.parameters))
                required = ('branch_id', 'solver_converged', 'mass_balance_passed', 'target_met', 'within_calibration_domain')
                missing = [field for field in required if field not in output]
                if missing:
                    raise KeyError('forward evaluator omitted fields: ' + ', '.join(missing))
                branch_id = str(output['branch_id']).strip()
                if not branch_id:
                    raise ValueError('forward evaluator returned an empty branch_id')
                point = {'realization_id': realization.realization_id, 'dose': dose, 'branch_id': branch_id, 'solver_converged': bool(output['solver_converged']), 'mass_balance_passed': bool(output['mass_balance_passed']), 'target_met': bool(output['target_met']), 'within_calibration_domain': bool(output['within_calibration_domain']), 'evaluation_status': 'completed', 'error_type': None, 'forward_diagnostics': output.get('forward_diagnostics')}
            except Exception as exc:
                point = {'realization_id': realization.realization_id, 'dose': dose, 'branch_id': None, 'solver_converged': False, 'mass_balance_passed': False, 'target_met': False, 'within_calibration_domain': False, 'evaluation_status': 'solver_or_contract_error', 'error_type': type(exc).__name__, 'forward_diagnostics': None}
                ensemble_blockers.add('dose_scan_contains_solver_or_contract_errors')
            scan.append(point)
            evaluation_rows.append(point)
        segments: list[list[dict[str, object]]] = []
        current: list[dict[str, object]] = []
        current_branch: str | None = None
        for point in scan:
            usable = bool(point['evaluation_status'] == 'completed' and point['solver_converged'] and point['mass_balance_passed'] and point['branch_id'])
            if not usable:
                if current:
                    segments.append(current)
                    current = []
                    current_branch = None
                continue
            branch_id = str(point['branch_id'])
            if current and branch_id != current_branch:
                segments.append(current)
                current = []
            current.append(point)
            current_branch = branch_id
        if current:
            segments.append(current)
        if not segments:
            ensemble_blockers.add('realization_has_no_successful_dose_segment')
            samples.append(DoseSample(sample_id=realization.realization_id, branch_id='unresolved_solver_failure', dose=doses[-1], solver_converged=False, mass_balance_passed=False, target_met=False, branch_resolved=False, within_calibration_domain=False))
            continue
        branch_occurrence: dict[str, int] = {}
        for segment in segments:
            base_branch = str(segment[0]['branch_id'])
            branch_occurrence[base_branch] = branch_occurrence.get(base_branch, 0) + 1
            segment_branch = base_branch if branch_occurrence[base_branch] == 1 else f'{base_branch}::segment_{branch_occurrence[base_branch]}'
            flags = [bool(point['target_met']) for point in segment]
            nonmonotonic = any((left and (not right) for left, right in zip(flags, flags[1:])))
            if nonmonotonic:
                ensemble_blockers.add('nonmonotonic_dose_response_segment')
            bracket: tuple[float, float] | None = None
            if flags[0] and float(segment[0]['dose']) == 0.0:
                bracket = (0.0, 0.0)
            else:
                for left, right in zip(segment, segment[1:]):
                    if not bool(left['target_met']) and bool(right['target_met']):
                        bracket = (float(left['dose']), float(right['dose']))
                        break
            target_met = any(flags)
            branch_resolved = bool(bracket is not None and (not nonmonotonic))
            if target_met and (not branch_resolved):
                ensemble_blockers.add('target_met_without_resolved_local_bracket')
            dose = bracket[1] if bracket is not None else float(segment[-1]['dose'])
            samples.append(DoseSample(sample_id=realization.realization_id, branch_id=segment_branch, dose=dose, solver_converged=True, mass_balance_passed=True, target_met=target_met, branch_resolved=branch_resolved, within_calibration_domain=all((bool(point['within_calibration_domain']) for point in segment))))
    envelope = build_branch_dose_envelope(samples, random_seed=random_seed, input_hash=input_hash, model_hash=model_hash, recommendation_gate=recommendation_gate)
    combined_blockers = set(envelope['recommendation_blockers'])
    combined_blockers.update(ensemble_blockers)
    if ensemble_blockers:
        envelope['recommendation_eligible'] = False
    envelope.update({'schema_version': 'branch_dose_grid_ensemble/v1', 'dose_semantics': 'upper_bound_of_resolved_local_grid_crossing', 'dose_grid': list(doses), 'realization_count': len(rows), 'evaluation_rows': evaluation_rows, 'recommendation_blockers': sorted(combined_blockers)})
    return envelope
__all__ = ['DoseSample', 'UncertaintyRealization', 'build_branch_dose_envelope', 'run_branch_aware_dose_grid_ensemble']
