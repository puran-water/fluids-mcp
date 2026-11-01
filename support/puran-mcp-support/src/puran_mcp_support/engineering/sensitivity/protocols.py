from __future__ import annotations
from collections.abc import Mapping
from typing import Callable, Optional, Protocol, runtime_checkable

@runtime_checkable
class EvaluateFn(Protocol):
    """Engine-agnostic KPI evaluation callable (one GSA sample in, KPIs out)."""

    def __call__(self, inputs: Mapping[str, float]) -> Optional[Mapping[str, float]]:
        ...
TransformFn = Callable[[dict[str, float]], dict[str, float]]
