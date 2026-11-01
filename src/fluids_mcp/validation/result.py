"""Post-calculation result validation and strict serialization."""

from __future__ import annotations

import json
import math
import re
from typing import Any, Mapping

from fluids_mcp.errors import CalculationError

_COMPLEX_TEXT = re.compile(
    r"^\(?\s*[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?\s*[-+]\s*"
    r"\d+(?:\.\d+)?(?:[eE][-+]?\d+)?j\s*\)?$"
)


def validate_result_payload(payload: Any, *, path: str = "result") -> Any:
    """Validate that a result payload contains only finite, JSON-safe values."""

    if isinstance(payload, bool) or payload is None:
        return payload
    if isinstance(payload, complex):
        raise CalculationError(
            "INVALID_RESULT_COMPLEX",
            f"{path} must not be complex",
            details={"path": path, "value": repr(payload)},
        )
    if isinstance(payload, (int, float)):
        if not math.isfinite(float(payload)):
            raise CalculationError(
                "INVALID_RESULT_NUMERIC",
                f"{path} must be finite",
                details={"path": path, "value": repr(payload)},
            )
        return payload
    if isinstance(payload, str):
        if _COMPLEX_TEXT.match(payload.strip()):
            raise CalculationError(
                "INVALID_RESULT_COMPLEX",
                f"{path} must not be stringified complex data",
                details={"path": path, "value": payload},
            )
        return payload
    if isinstance(payload, list):
        for idx, value in enumerate(payload):
            validate_result_payload(value, path=f"{path}[{idx}]")
        return payload
    if isinstance(payload, tuple):
        for idx, value in enumerate(payload):
            validate_result_payload(value, path=f"{path}[{idx}]")
        return payload
    if isinstance(payload, Mapping):
        for key, value in payload.items():
            validate_result_payload(key, path=f"{path}.<key>")
            validate_result_payload(value, path=f"{path}[{key!r}]")
        return payload
    return payload


def assert_finite_result(payload: Any, *, path: str = "result") -> Any:
    """PR 2 public alias for finite, real, JSON-safe result validation."""

    return validate_result_payload(payload, path=path)


def validate_forward_pressure_result(
    payload: Mapping[str, Any],
    *,
    inlet_key: str = "inlet_pressure_pa",
    outlet_key: str = "outlet_pressure_pa",
    pressure_drop_key: str = "pressure_drop_pa",
) -> None:
    inlet = payload.get(inlet_key)
    outlet = payload.get(outlet_key)
    pressure_drop = payload.get(pressure_drop_key)
    if inlet is not None and outlet is not None and outlet >= inlet:
        raise CalculationError(
            "INVALID_RESULT_PRESSURE_DIRECTION",
            "Forward-flow result has outlet pressure greater than or equal to inlet pressure",
            details={inlet_key: inlet, outlet_key: outlet},
        )
    if pressure_drop is not None and pressure_drop < 0:
        raise CalculationError(
            "INVALID_RESULT_SIGN",
            "Forward-flow pressure drop result must be non-negative",
            details={pressure_drop_key: pressure_drop},
        )


def strict_json_dumps(payload: Any, **kwargs: Any) -> str:
    """Serialize validated payload with strict JSON numeric behavior."""

    validate_result_payload(payload)
    return json.dumps(payload, allow_nan=False, **kwargs)
