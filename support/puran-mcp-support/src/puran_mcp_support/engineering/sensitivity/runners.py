from __future__ import annotations
import hashlib
import warnings
from typing import Any
import numpy as np
from SALib.analyze.morris import analyze as _morris_analyze
from SALib.analyze.sobol import analyze as _sobol_analyze
from SALib.sample.latin import sample as _latin_sample
from SALib.sample.morris import sample as _morris_sample
from SALib.sample.sobol import sample as _sobol_sample
from .problem import GSAProblem
from .protocols import EvaluateFn, TransformFn
from .result import GSAResult

def _design_hash(design: np.ndarray) -> str:
    """Stable hash of a design matrix -- proves seed -> design determinism."""
    h = hashlib.sha256()
    h.update(str(design.shape).encode())
    h.update(np.ascontiguousarray(design, dtype=np.float64).tobytes())
    return h.hexdigest()

def _evaluate_design(problem: GSAProblem, design: np.ndarray, evaluate: EvaluateFn, transform: TransformFn | None) -> dict[str, Any]:
    """Evaluate every design row; return masked ``Y`` plus failure bookkeeping.

    ``Y`` has shape ``(n_runs, n_kpis)`` with ``NaN`` wherever a KPI could not be
    obtained (evaluate raised / returned None / yielded a missing or non-finite
    value). Nothing is imputed.
    """
    names = problem.factor_names
    log_mask = problem.log_mask
    outputs = problem.outputs
    n_runs = design.shape[0]
    Y = np.full((n_runs, len(outputs)), np.nan, dtype=np.float64)
    failed_rows: set[int] = set()
    reasons: dict[int, str] = {}
    per_kpi_failed = {kpi: 0 for kpi in outputs}
    for i in range(n_runs):
        raw: dict[str, float] = {}
        for j, name in enumerate(names):
            v = float(design[i, j])
            raw[name] = 10.0 ** v if log_mask[j] else v
        try:
            inputs = transform(raw) if transform is not None else raw
        except Exception as exc:
            failed_rows.add(i)
            reasons[i] = f'transform_error: {type(exc).__name__}: {exc}'
            for kpi in outputs:
                per_kpi_failed[kpi] += 1
            continue
        try:
            out = evaluate(inputs)
        except Exception as exc:
            failed_rows.add(i)
            reasons[i] = f'evaluate_error: {type(exc).__name__}: {exc}'
            for kpi in outputs:
                per_kpi_failed[kpi] += 1
            continue
        if out is None:
            failed_rows.add(i)
            reasons[i] = 'evaluate_returned_none'
            for kpi in outputs:
                per_kpi_failed[kpi] += 1
            continue
        for k, kpi in enumerate(outputs):
            val = out.get(kpi) if hasattr(out, 'get') else None
            if val is None:
                per_kpi_failed[kpi] += 1
                failed_rows.add(i)
                reasons.setdefault(i, f'missing_kpi: {kpi}')
                continue
            fval = float(val)
            if not np.isfinite(fval):
                per_kpi_failed[kpi] += 1
                failed_rows.add(i)
                reasons.setdefault(i, f'nonfinite_kpi: {kpi}')
                continue
            Y[i, k] = fval
    return {'Y': Y, 'failed_samples': sorted(failed_rows), 'failure_reasons': reasons, 'per_kpi_failed': per_kpi_failed, 'n_evaluated': n_runs}

def _base_result(method: str, problem: GSAProblem, N: int, seed: int, design: np.ndarray, evaldata: dict[str, Any], extra_diag: dict[str, Any]) -> GSAResult:
    failed = evaldata['failed_samples']
    diagnostics = {'design_sha256': _design_hash(design), 'per_kpi_failed': evaldata['per_kpi_failed'], 'failure_reasons': evaldata['failure_reasons'], **extra_diag}
    return GSAResult(method=method, n=N, seed=seed, factor_names=problem.factor_names, kpi_names=list(problem.outputs), n_runs=int(design.shape[0]), n_evaluated=int(evaldata['n_evaluated']), n_failed=len(failed), failed_samples=failed, metadata=problem.metadata.model_dump(), diagnostics=diagnostics)

def run_morris(problem: GSAProblem, evaluate: EvaluateFn, *, N: int, seed: int, num_levels: int=4, transform: TransformFn | None=None, num_resamples: int=100, conf_level: float=0.95) -> GSAResult:
    """Morris elementary-effects screening (mu*, mu, sigma per factor per KPI).

    ``N`` trajectories of ``num_vars + 1`` runs each. Cheap relative to Sobol --
    use it to rank/triage factors before committing to a variance decomposition.

    Failed-solve masking: a Morris trajectory is unusable if *any* of its runs
    failed for the KPI (elementary effects are finite differences along the
    trajectory). Such trajectories are dropped per-KPI before analysis; the count
    is reported in ``diagnostics``. Individual failed samples are never imputed.
    """
    salib = problem.to_salib_problem()
    design = _morris_sample(salib, N, num_levels=num_levels, seed=seed)
    evaldata = _evaluate_design(problem, design, evaluate, transform)
    Y = evaldata['Y']
    block = salib['num_vars'] + 1
    n_traj = design.shape[0] // block
    result = _base_result('morris', problem, N, seed, design, evaldata, extra_diag={'num_levels': num_levels, 'n_trajectories': n_traj})
    traj_dropped: dict[str, int] = {}
    for k, kpi in enumerate(problem.outputs):
        ycol = Y[:, k]
        valid_rows = []
        kept = 0
        for t in range(n_traj):
            sl = slice(t * block, (t + 1) * block)
            if np.all(np.isfinite(ycol[sl])):
                valid_rows.extend(range(t * block, (t + 1) * block))
                kept += 1
        traj_dropped[kpi] = n_traj - kept
        if kept == 0:
            result.kpis_with_failures.append(kpi)
            result.indices[kpi] = {m: {f: float('nan') for f in problem.factor_names} for m in ('mu_star', 'mu', 'sigma', 'mu_star_conf')}
            continue
        if kept < n_traj:
            result.kpis_with_failures.append(kpi)
        rows = np.array(valid_rows, dtype=int)
        with warnings.catch_warnings(), np.errstate(all='ignore'):
            warnings.simplefilter('ignore')
            si = _morris_analyze(salib, design[rows], ycol[rows], num_levels=num_levels, num_resamples=num_resamples, conf_level=conf_level, seed=seed)
        names = [str(x) for x in si['names']]
        result.indices[kpi] = {'mu_star': dict(zip(names, [float(x) for x in si['mu_star']])), 'mu': dict(zip(names, [float(x) for x in si['mu']])), 'sigma': dict(zip(names, [float(x) for x in si['sigma']])), 'mu_star_conf': dict(zip(names, [float(x) for x in si['mu_star_conf']]))}
    result.diagnostics['trajectories_dropped'] = traj_dropped
    return result

def run_sobol(problem: GSAProblem, evaluate: EvaluateFn, *, N: int, seed: int, calc_second_order: bool=False, transform: TransformFn | None=None, num_resamples: int=100, conf_level: float=0.95) -> GSAResult:
    """Sobol first/total-order indices (+ optional second order) with CIs.

    Saltelli design: ``N * (2D + 2)`` runs with second order, ``N * (D + 2)``
    without. ``N`` should be a power of two for the low-discrepancy sequence.

    Failed-solve masking: the Saltelli estimator is matrix-coupled, so individual
    failed runs cannot be dropped without biasing the decomposition. We therefore
    do NOT impute; if a KPI's outputs contain any failure, its indices are
    returned as ``NaN`` and the KPI is flagged in ``kpis_with_failures`` with
    ``n_failed`` recorded. Screen with Morris and ensure the design converges so
    Sobol runs on a clean output set.
    """
    salib = problem.to_salib_problem()
    design = _sobol_sample(salib, N, calc_second_order=calc_second_order, seed=seed)
    evaldata = _evaluate_design(problem, design, evaluate, transform)
    Y = evaldata['Y']
    result = _base_result('sobol', problem, N, seed, design, evaldata, extra_diag={'calc_second_order': calc_second_order})
    names = problem.factor_names
    for k, kpi in enumerate(problem.outputs):
        ycol = np.ascontiguousarray(Y[:, k], dtype=np.float64)
        has_failure = bool(np.any(~np.isfinite(ycol)))
        if has_failure:
            result.kpis_with_failures.append(kpi)
        si = None
        with warnings.catch_warnings(), np.errstate(all='ignore'):
            warnings.simplefilter('ignore')
            try:
                si = _sobol_analyze(salib, ycol, calc_second_order=calc_second_order, num_resamples=num_resamples, conf_level=conf_level, seed=seed)
            except Exception:
                si = None
        if si is None:
            result.indices[kpi] = {m: {f: float('nan') for f in names} for m in ('S1', 'S1_conf', 'ST', 'ST_conf')}
            if kpi not in result.kpis_with_failures:
                result.kpis_with_failures.append(kpi)
            continue
        indices = {'S1': dict(zip(names, [float(x) for x in si['S1']])), 'S1_conf': dict(zip(names, [float(x) for x in si['S1_conf']])), 'ST': dict(zip(names, [float(x) for x in si['ST']])), 'ST_conf': dict(zip(names, [float(x) for x in si['ST_conf']]))}
        if calc_second_order and 'S2' in si:
            s2 = np.asarray(si['S2'], dtype=float)
            s2c = np.asarray(si['S2_conf'], dtype=float)
            pair2, pair2c = ({}, {})
            for a in range(len(names)):
                for b in range(a + 1, len(names)):
                    key = f'{names[a]}::{names[b]}'
                    pair2[key] = float(s2[a, b])
                    pair2c[key] = float(s2c[a, b])
            indices['S2'] = pair2
            indices['S2_conf'] = pair2c
        result.indices[kpi] = indices
    return result

def run_lhs_mc(problem: GSAProblem, evaluate: EvaluateFn, *, N: int, seed: int, transform: TransformFn | None=None, quantiles: tuple[float, ...]=(0.05, 0.25, 0.5, 0.75, 0.95)) -> GSAResult:
    """Latin-hypercube Monte-Carlo uncertainty propagation.

    Returns, per KPI, the output distribution stats (mean/std/min/max + the
    requested quantiles) in ``uncertainty``, and Spearman rank correlations of
    each factor against the KPI in ``indices[kpi]['spearman']`` (a robust,
    monotone sensitivity measure).

    Failed-solve masking: failed rows are dropped from the statistics and
    correlations per KPI (NaN-aware); ``n_failed`` is reported, nothing imputed.
    """
    from scipy.stats import spearmanr
    salib = problem.to_salib_problem()
    design = _latin_sample(salib, N, seed=seed)
    evaldata = _evaluate_design(problem, design, evaluate, transform)
    Y = evaldata['Y']
    result = _base_result('lhs_mc', problem, N, seed, design, evaldata, extra_diag={'quantiles': list(quantiles)})
    names = problem.factor_names
    for k, kpi in enumerate(problem.outputs):
        ycol = Y[:, k]
        finite = np.isfinite(ycol)
        yv = ycol[finite]
        if yv.size < finite.size:
            result.kpis_with_failures.append(kpi)
        if yv.size == 0:
            result.uncertainty[kpi] = {}
            result.indices[kpi] = {'spearman': {f: float('nan') for f in names}}
            continue
        unc = {'mean': float(np.mean(yv)), 'std': float(np.std(yv, ddof=1)) if yv.size > 1 else 0.0, 'min': float(np.min(yv)), 'max': float(np.max(yv)), 'n': int(yv.size)}
        for q in quantiles:
            unc[f'p{int(round(q * 100)):02d}'] = float(np.quantile(yv, q))
        result.uncertainty[kpi] = unc
        spearman: dict[str, float] = {}
        for j, name in enumerate(names):
            xv = design[finite, j]
            if yv.size < 3 or np.allclose(xv, xv[0]) or np.allclose(yv, yv[0]):
                spearman[name] = float('nan')
                continue
            with warnings.catch_warnings(), np.errstate(all='ignore'):
                warnings.simplefilter('ignore')
                res = spearmanr(xv, yv)
            spearman[name] = float(getattr(res, 'statistic', getattr(res, 'correlation', np.nan)))
        result.indices[kpi] = {'spearman': spearman}
    return result

def screen_then_rank(problem: GSAProblem, evaluate: EvaluateFn, *, kpi: str, n_morris: int, n_sobol: int, seed: int, top_k: int | None=None, mu_star_threshold: float | None=None, transform: TransformFn | None=None, num_levels: int=4, calc_second_order: bool=False) -> dict[str, Any]:
    """Run Morris screening, then Sobol on the factors it flags as influential.

    Factor selection (against ``kpi``'s ``mu_star``): the top ``top_k`` factors,
    and/or those above ``mu_star_threshold``. Non-selected factors are pinned at
    their nominal value (via a wrapped transform) for the Sobol pass, so the
    variance decomposition runs over a reduced, cheaper design.

    Returns ``{"morris", "sobol", "selected_factors", "fixed_factors"}``.
    """
    if kpi not in problem.outputs:
        raise KeyError(f'kpi {kpi!r} not in problem.outputs {problem.outputs}')
    morris = run_morris(problem, evaluate, N=n_morris, seed=seed, num_levels=num_levels, transform=transform)
    ranked = morris.ranking(kpi, 'mu_star')
    selected = [name for name, _ in ranked]
    if mu_star_threshold is not None:
        selected = [name for name, v in ranked if v is not None and (not (isinstance(v, float) and np.isnan(v))) and (float(v) >= mu_star_threshold)]
    if top_k is not None:
        selected = selected[:top_k]
    if not selected:
        selected = list(problem.factor_names)
    selected_set = set(selected)
    fixed = {f.name: f.nominal_value for f in problem.factors if f.name not in selected_set}
    reduced = GSAProblem(factors=[f for f in problem.factors if f.name in selected_set], outputs=list(problem.outputs), metadata=problem.metadata)

    def screen_transform(raw: dict[str, float]) -> dict[str, float]:
        merged = {**fixed, **raw}
        return transform(merged) if transform is not None else merged
    sobol = run_sobol(reduced, evaluate, N=n_sobol, seed=seed, calc_second_order=calc_second_order, transform=screen_transform)
    return {'morris': morris, 'sobol': sobol, 'selected_factors': selected, 'fixed_factors': fixed}
