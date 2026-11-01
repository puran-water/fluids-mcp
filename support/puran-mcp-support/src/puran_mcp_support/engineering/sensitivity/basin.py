from __future__ import annotations
import math
from collections.abc import Mapping, Sequence
from typing import Any
DEFAULT_LOCAL_BASIN_THRESHOLD = 0.34

def _is_log(factor: Any) -> bool:
    is_log = getattr(factor, 'is_log', None)
    if is_log is not None:
        return bool(is_log)
    dist = getattr(factor, 'distribution', None)
    return str(getattr(dist, 'value', dist)).lower() == 'loguniform'

def _nominal(factor: Any) -> float:
    nominal = getattr(factor, 'nominal_value', None)
    if nominal is not None:
        return float(nominal)
    nominal = getattr(factor, 'nominal', None)
    if nominal is not None:
        return float(nominal)
    low = float(getattr(factor, 'low'))
    high = float(getattr(factor, 'high'))
    if _is_log(factor):
        return float(math.sqrt(low * high))
    return float(0.5 * (low + high))

def factor_relative_displacement(factor: Any, value: float) -> float:
    """Normalised displacement of ``value`` from ``factor``'s nominal in [0, ~1].

    Linear factors: ``|value - nominal| / (high - low)``. Loguniform factors:
    the same ratio computed in log10 space. Returns 0.0 for a degenerate span.
    """
    low = float(getattr(factor, 'low'))
    high = float(getattr(factor, 'high'))
    nominal = _nominal(factor)
    val = float(value)
    if _is_log(factor) and low > 0 and (high > 0) and (val > 0) and (nominal > 0):
        span = math.log10(high) - math.log10(low)
        if span <= 0:
            return 0.0
        return abs(math.log10(val) - math.log10(nominal)) / span
    span = high - low
    if span <= 0:
        return 0.0
    return abs(val - nominal) / span

def relative_displacement(inputs: Mapping[str, float], factors: Sequence[Any]) -> float:
    """Max normalised displacement-from-nominal over the sample's factors.

    Factors with no corresponding key in ``inputs`` are skipped (an annotation-
    only factor cannot move the basin). Returns 0.0 if nothing is comparable.
    """
    worst = 0.0
    for factor in factors:
        name = getattr(factor, 'name', None)
        if name is None or name not in inputs:
            continue
        value = inputs[name]
        if value is None:
            continue
        try:
            rel = factor_relative_displacement(factor, float(value))
        except (TypeError, ValueError):
            continue
        if rel > worst:
            worst = rel
    return worst

def exceeds_local_basin(inputs: Mapping[str, float], factors: Sequence[Any], threshold: float=DEFAULT_LOCAL_BASIN_THRESHOLD) -> bool:
    """True when the sample leaves the warm-start basin (force a cold solve)."""
    return relative_displacement(inputs, factors) > float(threshold)
