"""Typed calculation errors and public-message redaction."""

from __future__ import annotations

from typing import Any, Dict, Mapping, Optional
from uuid import uuid4

from puran_mcp_support.safe_paths import redact_paths


class CalculationError(Exception):
    """A typed, user-safe calculation error."""

    def __init__(
        self,
        code: str,
        message: str,
        *,
        details: Optional[Mapping[str, Any]] = None,
        severity: str = "error",
        correlation_id: Optional[str] = None,
    ) -> None:
        self.code = code
        self.message = redact_public_text(message)
        self.details = redact_public_payload(dict(details or {}))
        self.severity = severity
        self.correlation_id = correlation_id or uuid4().hex
        super().__init__(self.message)

    def to_public_dict(self) -> Dict[str, Any]:
        payload: Dict[str, Any] = {
            "code": self.code,
            "message": self.message,
            "details": self.details,
            "correlation_id": self.correlation_id,
            "severity": self.severity,
        }
        return payload


class FittingKError(CalculationError):
    """A typed error raised when a fitting K-value cannot be resolved."""


class PipeMaterialError(CalculationError):
    """A typed error raised when pipe roughness cannot be resolved safely."""


def redact_public_text(value: Any) -> Any:
    """Redact filesystem paths from strings returned to callers."""

    return redact_paths(value) if isinstance(value, str) else value


def redact_public_payload(payload: Any) -> Any:
    """Recursively redact filesystem paths from JSON-like payloads."""

    return redact_paths(payload)
