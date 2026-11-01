from __future__ import annotations
import math
from dataclasses import asdict, dataclass, field
from typing import Any

def _json_safe(obj: Any) -> Any:
    """Recursively convert NaN/inf floats to None for JSON serialisation."""
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj

@dataclass
class GSAResult:
    """One global-sensitivity / uncertainty study, uniform across methods."""
    method: str
    n: int
    seed: int
    factor_names: list[str]
    kpi_names: list[str]
    n_runs: int
    n_evaluated: int
    n_failed: int
    failed_samples: list[int]
    indices: dict[str, dict[str, dict[str, float]]] = field(default_factory=dict)
    uncertainty: dict[str, dict[str, float]] = field(default_factory=dict)
    kpis_with_failures: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
    diagnostics: dict[str, Any] = field(default_factory=dict)

    @property
    def failure_fraction(self) -> float:
        """Fraction of design rows with >=1 failed/NaN KPI (0.0 if none ran).

        Surfaced so a sweep reports its non-convergence rate explicitly rather
        than hiding it behind imputed values — a high failure fraction means the
        warm/cold solve selection or the factor ranges need revisiting.
        """
        if not self.n_runs:
            return 0.0
        return self.n_failed / self.n_runs

    def ranking(self, kpi: str, metric: str) -> list[tuple[str, float]]:
        """Factors ranked by ``indices[kpi][metric]`` descending (NaN last)."""
        if kpi not in self.indices or metric not in self.indices[kpi]:
            raise KeyError(f'no metric {metric!r} for KPI {kpi!r}')
        items = list(self.indices[kpi][metric].items())

        def _key(kv: tuple[str, float]) -> tuple[int, float]:
            v = kv[1]
            if v is None or (isinstance(v, float) and math.isnan(v)):
                return (1, 0.0)
            return (0, -float(v))
        return sorted(items, key=_key)

    def top_factors(self, kpi: str, metric: str, k: int) -> list[str]:
        """Names of the ``k`` highest-``metric`` factors for ``kpi``."""
        return [name for name, _ in self.ranking(kpi, metric)[:k]]

    def to_dict(self) -> dict[str, Any]:
        """JSON-safe dict (NaN -> None) for persistence / transport.

        ``failure_fraction`` is a derived property (not a dataclass field), so it
        is injected explicitly so the serialised payload carries the sweep's
        non-convergence rate.
        """
        payload = _json_safe(asdict(self))
        payload['failure_fraction'] = self.failure_fraction
        return payload
