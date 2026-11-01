"""Unified property lookup for fluids and pipes."""

import json
from typing import Optional, Literal

from pydantic import ValidationError

from fluids_mcp.errors import CalculationError
from fluids_mcp.schemas.omnitools import CalculationMessage, PropertiesInput, PropertiesOutput
from fluids_mcp.validation.result import assert_finite_result
from fluids_mcp.tools.fluid_properties import get_fluid_properties, list_available_fluids
from fluids_mcp.tools.pipe_properties import get_pipe_properties


VALID_LOOKUP_TYPES = {"fluid", "pipe", "list_fluids"}
DEFAULT_LOOKUP_TYPE = "fluid"
MIN_ABSOLUTE_PROPERTY_PRESSURE_BAR = 1.0e-6


def _message_from_error(error: CalculationError) -> CalculationMessage:
    return CalculationMessage(
        code=error.code,
        message=error.message,
        severity=error.severity,  # type: ignore[arg-type]
        details=error.details,
        correlation_id=error.correlation_id,
    )


def _error_output(
    lookup_type: str,
    error: CalculationError,
    *,
    inputs_resolved: Optional[list[str]] = None,
) -> dict:
    normalized = lookup_type if lookup_type in VALID_LOOKUP_TYPES else DEFAULT_LOOKUP_TYPE
    return PropertiesOutput(
        lookup_type=normalized,  # type: ignore[arg-type]
        errors=[_message_from_error(error)],
        inputs_resolved=inputs_resolved or [],
    ).model_dump(exclude_none=True)


def _decode_tool_json(payload: str, lookup_type: str) -> dict:
    try:
        decoded = json.loads(payload)
    except json.JSONDecodeError as exc:
        raise CalculationError(
            "INVALID_TOOL_RESPONSE",
            f"Property backend returned invalid JSON: {exc}",
            details={"lookup_type": lookup_type},
        ) from exc
    if not isinstance(decoded, dict):
        raise CalculationError(
            "INVALID_TOOL_RESPONSE",
            "Property backend returned a non-object response",
            details={"lookup_type": lookup_type, "response_type": type(decoded).__name__},
        )
    return decoded


def _raise_if_backend_error(payload: dict, lookup_type: str) -> None:
    if payload.get("error"):
        raise CalculationError(
            "PROPERTY_LOOKUP_FAILED",
            str(payload["error"]),
            details={"lookup_type": lookup_type},
        )
    if payload.get("errors"):
        raise CalculationError(
            "PROPERTY_LOOKUP_FAILED",
            "; ".join(str(item) for item in payload["errors"]),
            details={"lookup_type": lookup_type},
        )
    if payload.get("successful") is False:
        raise CalculationError(
            "PROPERTY_LOOKUP_FAILED",
            str(payload.get("message") or "Property lookup failed"),
            details={"lookup_type": lookup_type},
        )


def properties(
    lookup_type: Literal["fluid", "pipe", "list_fluids"] = "fluid",

    # Fluid properties parameters
    fluid_name: Optional[str] = None,
    temperature_c: Optional[float] = None,
    pressure_bar: Optional[float] = 1.0,

    # Pipe properties parameters
    nominal_size: Optional[float] = None,
    schedule: str = "40",
    inner_diameter: Optional[float] = None,
    outer_diameter: Optional[float] = None,
    material: Optional[str] = None,
) -> dict:
    """Unified property lookup for fluids and pipes.

    Provides access to thermodynamic properties and dimensions:
    - lookup_type='fluid': Get fluid properties at specified conditions
    - lookup_type='pipe': Get pipe dimensions and roughness
    - lookup_type='list_fluids': List all available fluids

    Args:
        lookup_type: Type of property lookup

        Fluid parameters:
            fluid_name: Name of the fluid (e.g., "Water", "Air", "Methane")
            temperature_c: Temperature in Celsius
            pressure_bar: Absolute pressure in bar (default: 1.0). Values at or
                below 1e-6 bar are rejected as a zero/near-vacuum state before
                calling the thermophysical-property backends.

        Pipe parameters:
            nominal_size: Nominal pipe size in inches
            schedule: Pipe schedule designation
            inner_diameter: Inner diameter in m (to find closest)
            outer_diameter: Outer diameter in m (to find closest)
            material: Pipe material name

    Returns:
        Structured property lookup result matching PropertiesOutput

    Examples:
        Fluid properties:
        >>> properties(lookup_type="fluid", fluid_name="Water",
        ...           temperature_c=25, pressure_bar=1)

        Pipe properties:
        >>> properties(lookup_type="pipe", nominal_size=2, schedule="40")

        List fluids:
        >>> properties(lookup_type="list_fluids")
    """
    try:
        request = PropertiesInput(
            lookup_type=lookup_type,
            fluid_name=fluid_name,
            temperature_c=temperature_c,
            pressure_bar=pressure_bar,
            nominal_size=nominal_size,
            schedule=schedule,
            inner_diameter=inner_diameter,
            outer_diameter=outer_diameter,
            material=material,
        )

        inputs_resolved = [f"lookup_type={request.lookup_type}"]

        if request.lookup_type == "fluid":
            if not request.fluid_name or request.temperature_c is None:
                raise CalculationError(
                    "MISSING_REQUIRED_INPUT",
                    "fluid_name and temperature_c are required for fluid lookup",
                    details={"lookup_type": request.lookup_type},
                )
            resolved_pressure_bar = (
                request.pressure_bar if request.pressure_bar is not None else 1.0
            )
            if resolved_pressure_bar <= MIN_ABSOLUTE_PROPERTY_PRESSURE_BAR:
                raise CalculationError(
                    "INVALID_ABSOLUTE_PRESSURE",
                    "pressure_bar is an absolute pressure and must be greater than "
                    f"{MIN_ABSOLUTE_PROPERTY_PRESSURE_BAR:g} bar; zero and near-zero "
                    "states are outside this public property lookup's valid domain.",
                    details={
                        "pressure_bar": resolved_pressure_bar,
                        "minimum_exclusive_bar": MIN_ABSOLUTE_PROPERTY_PRESSURE_BAR,
                    },
                )
            raw = get_fluid_properties(
                fluid_name=request.fluid_name,
                temperature_c=request.temperature_c,
                pressure_bar=resolved_pressure_bar,
            )
            payload = _decode_tool_json(raw, request.lookup_type)
            _raise_if_backend_error(payload, request.lookup_type)
            assert_finite_result(payload)
            return PropertiesOutput(
                lookup_type=request.lookup_type,
                properties=payload,
                inputs_resolved=inputs_resolved,
            ).model_dump(exclude_none=True)

        if request.lookup_type == "pipe":
            raw = get_pipe_properties(
                nominal_size=request.nominal_size,
                schedule=request.schedule,
                inner_diameter=request.inner_diameter,
                outer_diameter=request.outer_diameter,
                material=request.material,
            )
            payload = _decode_tool_json(raw, request.lookup_type)
            _raise_if_backend_error(payload, request.lookup_type)
            assert_finite_result(payload)
            return PropertiesOutput(
                lookup_type=request.lookup_type,
                properties=payload,
                inputs_resolved=inputs_resolved,
            ).model_dump(exclude_none=True)

        raw = list_available_fluids()
        payload = _decode_tool_json(raw, request.lookup_type)
        _raise_if_backend_error(payload, request.lookup_type)
        assert_finite_result(payload)
        return PropertiesOutput(
            lookup_type=request.lookup_type,
            resources=payload,
            inputs_resolved=inputs_resolved,
        ).model_dump(exclude_none=True)

    except CalculationError as exc:
        return _error_output(str(lookup_type), exc)
    except ValidationError as exc:
        return _error_output(
            str(lookup_type),
            CalculationError(
                "INVALID_PROPERTIES_INPUT",
                str(exc),
                details={"lookup_type": lookup_type},
            ),
        )
    except Exception as exc:
        return _error_output(
            str(lookup_type),
            CalculationError(
                "PROPERTY_LOOKUP_FAILED",
                str(exc),
                details={"lookup_type": lookup_type},
            ),
        )
