"""Typed contract and validation layer for fluids-mcp.

PR 1 intentionally does not wire existing omnitools through this package.
The modules here provide schemas, unit conversion, validation, and sanitized
error primitives for later pilot conversions.
"""

from fluids_mcp.engineering_utils_bridge import ensure_engineering_utils_importable

ensure_engineering_utils_importable()

__all__ = ["__version__"]

__version__ = "0.1.0"
