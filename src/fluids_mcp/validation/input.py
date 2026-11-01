"""Public input contract validation helpers."""

from __future__ import annotations

import math
from typing import Any, Iterable, Mapping, Optional

from fluids_mcp.errors import CalculationError


def require_known_fields(payload: Mapping[str, Any], allowed_fields: Iterable[str]) -> None:
    unknown = sorted(set(payload) - set(allowed_fields))
    if unknown:
        raise CalculationError(
            "UNKNOWN_INPUT_FIELD",
            f"Unsupported input field(s): {', '.join(unknown)}",
            details={"unknown_fields": unknown},
        )


def require_one_of(values: Mapping[str, Any], *, group_name: str = "input") -> str:
    present = [name for name, value in values.items() if value is not None]
    if not present:
        raise CalculationError(
            "MISSING_REQUIRED_INPUT",
            f"Provide at least one {group_name}: {', '.join(values)}",
            details={"fields": list(values)},
        )
    return present[0]


def require_exactly_one_of(values: Mapping[str, Any], *, group_name: str = "input") -> str:
    present = [name for name, value in values.items() if value is not None]
    if len(present) != 1:
        raise CalculationError(
            "CONFLICTING_INPUTS",
            f"Provide exactly one {group_name}: {', '.join(values)}",
            details={"fields": list(values), "present_fields": present},
        )
    return present[0]


def require_at_least(values: Mapping[str, Any], count: int, *, group_name: str = "input") -> None:
    present = [name for name, value in values.items() if value is not None]
    if len(present) < count:
        raise CalculationError(
            "MISSING_REQUIRED_INPUT",
            f"Provide at least {count} {group_name} values",
            details={"fields": list(values), "present_fields": present},
        )


def require_finite(name: str, value: Optional[float]) -> None:
    if value is None:
        return
    if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
        raise CalculationError(
            "INVALID_NUMERIC_INPUT",
            f"{name} must be a finite number",
            details={"field": name, "value": repr(value)},
        )


def require_positive(name: str, value: Optional[float], *, allow_zero: bool = False) -> None:
    if value is None:
        return
    require_finite(name, value)
    limit_ok = value >= 0 if allow_zero else value > 0
    if not limit_ok:
        qualifier = "non-negative" if allow_zero else "positive"
        raise CalculationError(
            "INVALID_NUMERIC_INPUT",
            f"{name} must be {qualifier}",
            details={"field": name, "value": value},
        )


def require_literal(name: str, value: Any, allowed: Iterable[str]) -> None:
    allowed_values = list(allowed)
    if value not in allowed_values:
        raise CalculationError(
            "UNSUPPORTED_VALUE",
            f"{name} must be one of: {', '.join(allowed_values)}",
            details={"field": name, "value": value, "allowed": allowed_values},
        )
