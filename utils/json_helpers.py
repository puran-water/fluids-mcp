"""
JSON serialization helpers for fluids-mcp.

This module provides utilities for safe JSON serialization, particularly
handling special float values (inf, nan) that are not valid in JSON per RFC 7159.
"""

import math
import json
from typing import Any, Dict, List, Union


def sanitize_for_json(obj: Any) -> Any:
    """
    Recursively sanitize an object for JSON serialization.

    Replaces inf and nan float values with None, which serializes to null.
    This prevents json.dumps from producing invalid JSON.

    Also handles numpy scalar types if numpy is available.

    Args:
        obj: Any Python object to sanitize

    Returns:
        Sanitized object safe for json.dumps

    Examples:
        >>> sanitize_for_json({'value': float('inf')})
        {'value': None}
        >>> sanitize_for_json([1.0, float('nan'), 3.0])
        [1.0, None, 3.0]
    """
    if obj is None:
        return None

    if isinstance(obj, bool):
        return obj

    if isinstance(obj, (int, str)):
        return obj

    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return None
        return obj

    # Handle numpy types if available
    try:
        import numpy as np
        if isinstance(obj, (np.integer, np.floating)):
            val = float(obj)
            if math.isnan(val) or math.isinf(val):
                return None
            return val
        if isinstance(obj, np.ndarray):
            return sanitize_for_json(obj.tolist())
    except ImportError:
        pass

    if isinstance(obj, dict):
        return {k: sanitize_for_json(v) for k, v in obj.items()}

    if isinstance(obj, (list, tuple)):
        return [sanitize_for_json(v) for v in obj]

    # For other types, try to convert to a basic type
    try:
        # Handle objects with __dict__
        if hasattr(obj, '__dict__'):
            return sanitize_for_json(obj.__dict__)
    except Exception:
        pass

    # Fallback: convert to string
    return str(obj)


def safe_json_dumps(obj: Any, **kwargs) -> str:
    """
    Safely serialize an object to JSON string.

    Applies sanitization before serialization to handle inf/nan values.

    Args:
        obj: Object to serialize
        **kwargs: Additional arguments passed to json.dumps

    Returns:
        JSON string

    Examples:
        >>> safe_json_dumps({'value': float('inf')})
        '{"value": null}'
    """
    sanitized = sanitize_for_json(obj)
    return json.dumps(sanitized, **kwargs)


def is_valid_number(value: float) -> bool:
    """
    Check if a float value is valid (not inf, nan, or CoolProp _HUGE).

    CoolProp returns _HUGE (approximately 1e308) for invalid states.
    This function checks for all invalid float conditions.

    Args:
        value: Float value to check

    Returns:
        True if value is a valid, finite number
    """
    if not isinstance(value, (int, float)):
        return False
    if math.isnan(value) or math.isinf(value):
        return False
    # CoolProp _HUGE is approximately 1e308
    if abs(value) > 1e300:
        return False
    return True


def validate_and_sanitize_properties(
    properties: Dict[str, float],
    property_names: List[str] = None
) -> Dict[str, Any]:
    """
    Validate and sanitize a dictionary of fluid properties.

    Checks each property for valid numeric values and replaces invalid
    ones with None while logging warnings.

    Args:
        properties: Dictionary of property name -> value
        property_names: Optional list of required property names

    Returns:
        Sanitized properties dictionary with validity info
    """
    result = {
        'properties': {},
        'valid': True,
        'invalid_properties': []
    }

    for name, value in properties.items():
        if value is None:
            result['properties'][name] = None
            continue

        if is_valid_number(value):
            result['properties'][name] = value
        else:
            result['properties'][name] = None
            result['invalid_properties'].append(name)
            result['valid'] = False

    if property_names:
        missing = [n for n in property_names if n not in result['properties'] or result['properties'][n] is None]
        if missing:
            result['valid'] = False
            result['missing_properties'] = missing

    return result
