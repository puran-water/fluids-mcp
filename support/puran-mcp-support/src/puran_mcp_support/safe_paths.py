"""Path-redaction utility for public error messages."""

from __future__ import annotations

import re
from collections.abc import Mapping
from typing import Any

__all__ = ["redact_paths", "redact_mapping"]


_PATH_TOKEN = r"[^\s'\"<>]+"
_PATH_TOKEN_NO_COLON = r"[^\s'\"<>:]+"

# Patterns intentionally run most-specific to least-specific.  Venv and Python
# internals are scrubbed first so package/module tails remain useful even when

_VENV_SITE_PACKAGES = re.compile(
    rf"(?:{_PATH_TOKEN}[/\\])?\.venv[/\\]lib[/\\]python\d+(?:\.\d+)?[/\\]site-packages[/\\]"
)
_VENV_BIN_EXEC = re.compile(
    rf"(?:{_PATH_TOKEN}[/\\])?\.venv[/\\]bin[/\\]([^/\s'\"<>:]+)"
)
_VENV_PATH = re.compile(
    rf"(?:{_PATH_TOKEN}[/\\])?\.venv[/\\](?:{_PATH_TOKEN_NO_COLON}[/\\])*([^/\s'\"<>:]+)"
)
_SITE_PACKAGES = re.compile(rf"(?:{_PATH_TOKEN_NO_COLON}[/\\])*site-packages[/\\]")
_USR_LIB_PYTHON = re.compile(
    rf"/usr/lib/python\d+(?:\.\d+)?/(?:{_PATH_TOKEN_NO_COLON}/)*([^/\s'\"<>:]+)"
)


_HOME_PATH = re.compile(r"/(?:home|Users)/[^/\s'\"<>]+/(?:[^/\s'\"<>]+/)*([^/\s'\"<>:]+)")
_HOME_PATH_BARE = re.compile(r"/(?:home|Users)/[^/\s'\"<>]+/?")

# /tmp/<workspace-tag>/... only when the tag looks like one of ours.  Avoid
# blanket /tmp redaction because users may legitimately pass /tmp paths.
_TMP_WORKSPACE = re.compile(
    r"/tmp/(?:qsdsan|watertap|mcp|claude|viphreeqc)[^\s'\"]*/(?=[^\s'\"]+)"
)


def redact_paths(payload: Any) -> Any:
    """Redact developer paths from strings and nested JSON-like structures.

    Exception objects passed in (e.g. inside an error envelope's ``details``)
    are converted to a redacted string that follows the cause/context chain.
    Returning a plain str avoids leaving an unserialisable ``BaseException``
    in the public payload, which would crash downstream JSON encoders and
    typically carry the same raw paths that prompted redaction in the first
    place.
    """

    if isinstance(payload, str):
        return _redact_text(payload)
    if isinstance(payload, BaseException):
        return _redact_exception(payload)
    if isinstance(payload, Mapping):
        return {
            _redact_mapping_key(key): redact_paths(value)
            for key, value in payload.items()
        }
    if isinstance(payload, list):
        return [redact_paths(value) for value in payload]
    if isinstance(payload, tuple):
        return tuple(redact_paths(value) for value in payload)
    return payload


def _redact_mapping_key(key: Any) -> str:
    """Return a redacted JSON-safe string for mapping keys."""

    redacted = redact_paths(key)
    return _redact_text(str(redacted))


def _redact_exception(exc: BaseException) -> str:
    parts: list[str] = []
    seen: set[int] = set()
    current: BaseException | None = exc
    while current is not None and id(current) not in seen:
        seen.add(id(current))
        parts.append(f"{type(current).__name__}: {_redact_text(str(current))}")
        current = current.__cause__ or current.__context__
    return " | caused by ".join(parts)


def redact_mapping(payload: Any) -> Any:
    """Recursively redact paths in any JSON-like structure."""

    return redact_paths(payload)


def _redact_text(text: str) -> str:
    out = _VENV_SITE_PACKAGES.sub("<venv>/", text)
    out = _VENV_BIN_EXEC.sub(r"<venv>/bin/\1", out)
    out = _VENV_PATH.sub(r"<venv>/\1", out)
    out = _SITE_PACKAGES.sub("<packages>/", out)
    out = _USR_LIB_PYTHON.sub(r"<stdlib>/\1", out)
    out = _HOME_PATH.sub(r"<path>/\1", out)
    out = _HOME_PATH_BARE.sub("<path>", out)
    out = _TMP_WORKSPACE.sub("<tmp>/", out)
    return out
