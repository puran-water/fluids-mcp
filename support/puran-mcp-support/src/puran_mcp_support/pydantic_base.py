"""Shared strict Pydantic base classes for engineering MCP servers."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from decimal import Decimal
from typing import Any

from pydantic import BaseModel, ConfigDict, model_validator

__all__ = ["StrictFiniteBaseModel", "validate_finite_recursive"]


def validate_finite_recursive(value: Any) -> Any:
    """Return *value* if all contained numeric values are finite."""

    _validate_finite_recursive(value, "$")
    return value


def _validate_finite_recursive(value: Any, path: str) -> None:
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"non-finite numeric value at {path}: {value!r}")
        return

    if isinstance(value, int):
        return

    if isinstance(value, Decimal):
        if not value.is_finite():
            raise ValueError(f"non-finite decimal value at {path}: {value!r}")
        return

    if isinstance(value, BaseModel):
        _validate_finite_recursive(value.model_dump(), path)
        return

    if isinstance(value, Mapping):
        for key, item in value.items():
            _validate_finite_recursive(key, f"{path}.<key>")
            _validate_finite_recursive(item, f"{path}[{key!r}]")
        return

    if isinstance(value, (str, bytes, bytearray)):
        return

    if isinstance(value, Sequence):
        for index, item in enumerate(value):
            _validate_finite_recursive(item, f"{path}[{index}]")


class StrictFiniteBaseModel(BaseModel):
    """Base model that forbids extras, assignment drift, and non-finite numbers."""

    model_config = ConfigDict(
        allow_inf_nan=False,
        extra="forbid",
        validate_assignment=True,
    )

    @model_validator(mode="after")
    def _validate_finite_values(self) -> "StrictFiniteBaseModel":
        validate_finite_recursive(self)
        return self
