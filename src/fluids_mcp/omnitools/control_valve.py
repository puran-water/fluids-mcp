"""Unified control valve sizing for liquid and gas service."""

from __future__ import annotations

import inspect
import json
from typing import Literal, Optional

from pydantic import ValidationError

from fluids_mcp.errors import CalculationError
from fluids_mcp.schemas.omnitools import CalculationMessage, ControlValveInput, ControlValveOutput
from fluids_mcp.validation.result import assert_finite_result
from fluids_mcp.tools.gas_control_valve import calculate_gas_control_valve
from fluids_mcp.tools.liquid_control_valve import calculate_liquid_control_valve


def _message_from_error(error: CalculationError) -> CalculationMessage:
    return CalculationMessage(
        code=error.code,
        message=error.message,
        severity=error.severity,  # type: ignore[arg-type]
        details=error.details,
        correlation_id=error.correlation_id,
    )


def _error_output(phase: str, error: CalculationError) -> dict:
    normalized = phase if phase in {"liquid", "gas"} else None
    output = ControlValveOutput(
        phase=normalized,  # type: ignore[arg-type]
        errors=[_message_from_error(error)],
    )
    payload = output.model_dump(exclude_none=True)
    assert_finite_result(payload)
    return payload


def _warning_messages(raw_warnings) -> list[CalculationMessage]:
    if not raw_warnings:
        return []
    warnings = raw_warnings if isinstance(raw_warnings, list) else [raw_warnings]
    return [
        CalculationMessage(code="CONTROL_VALVE_WARNING", message=str(message), severity="warning")
        for message in warnings
    ]


def _decode_tool_json(payload: str, phase: str) -> dict:
    try:
        decoded = json.loads(payload)
    except json.JSONDecodeError as exc:
        raise CalculationError(
            "INVALID_TOOL_RESPONSE",
            f"Control-valve backend returned invalid JSON: {exc}",
            details={"phase": phase},
        ) from exc
    if not isinstance(decoded, dict):
        raise CalculationError(
            "INVALID_TOOL_RESPONSE",
            "Control-valve backend returned a non-object response",
            details={"phase": phase, "response_type": type(decoded).__name__},
        )
    return decoded


def _payload_error(payload: dict, phase: str) -> Optional[CalculationError]:
    if payload.get("error_code") or payload.get("error"):
        details = {"phase": phase}
        errors = payload.get("errors")
        if isinstance(errors, list) and errors and isinstance(errors[0], dict):
            raw_details = errors[0].get("details", {})
            if isinstance(raw_details, dict):
                details.update(raw_details)
        return CalculationError(
            str(payload.get("error_code") or "CONTROL_VALVE_CALCULATION_FAILED"),
            str(payload.get("error") or "Control-valve calculation failed"),
            details=details,
        )
    errors = payload.get("errors")
    if not errors:
        return None
    if isinstance(errors, list) and errors:
        first = errors[0]
        if isinstance(first, dict):
            return CalculationError(
                str(first.get("code") or "CONTROL_VALVE_CALCULATION_FAILED"),
                str(first.get("message") or first),
                details=first.get("details", {}) if isinstance(first.get("details", {}), dict) else {},
            )
        return CalculationError(
            "CONTROL_VALVE_CALCULATION_FAILED",
            "; ".join(str(item) for item in errors),
            details={"phase": phase},
        )
    return CalculationError("CONTROL_VALVE_CALCULATION_FAILED", str(errors), details={"phase": phase})


def _backend_context(payload: dict) -> dict:
    return {
        key: value
        for key, value in payload.items()
        if key not in {"error", "error_code", "errors"} and value is not None
    }


def _map_output(payload: dict, request: ControlValveInput) -> dict:
    error = _payload_error(payload, request.phase)
    if error is not None:
        output = ControlValveOutput(
            phase=request.phase,
            errors=[_message_from_error(error)],
            warnings=_warning_messages(payload.get("warnings")),
            inputs_resolved=payload.get("inputs_resolved") or payload.get("log") or [],
            details=_backend_context(payload),
        )
    else:
        output = ControlValveOutput(
            phase=request.phase,
            flow_rate_kg_s=payload.get("flow_rate_kg_s"),
            inlet_pressure_pa=payload.get("inlet_pressure_pa"),
            outlet_pressure_pa=payload.get("outlet_pressure_pa"),
            pressure_drop_pa=payload.get("pressure_drop_pa"),
            recommended_nps_inches=payload.get("recommended_nps_inches"),
            required_cv=payload.get("required_cv"),
            rated_cv=payload.get("rated_cv"),
            cv_utilization=payload.get("cv_utilization"),
            cavitation=payload.get("cavitation"),
            flashing=payload.get("flashing"),
            choked=payload.get("choked"),
            expansion_factor=payload.get("expansion_factor"),
            pressure_ratio=payload.get("pressure_ratio"),
            critical_pressure_ratio=payload.get("critical_pressure_ratio"),
            mach_outlet=payload.get("mach_outlet"),
            warnings=_warning_messages(payload.get("warnings")),
            inputs_resolved=payload.get("inputs_resolved") or payload.get("log") or [],
            details=payload,
        )
    result = output.model_dump(exclude_none=True)
    assert_finite_result(result)
    return result


def control_valve(
    phase: Literal["liquid", "gas"] = "liquid",

    # Flow parameters
    flow_rate: Optional[float] = None,
    flow_rate_gpm: Optional[float] = None,
    flow_rate_kg_s: Optional[float] = None,
    flow_rate_norm_m3_hr: Optional[float] = None,
    flow_rate_std_m3_hr: Optional[float] = None,

    # Pressure parameters
    inlet_pressure: Optional[float] = None,
    outlet_pressure: Optional[float] = None,
    inlet_pressure_psi: Optional[float] = None,
    outlet_pressure_psi: Optional[float] = None,

    # Fluid properties
    fluid_density: Optional[float] = None,
    fluid_density_lbft3: Optional[float] = None,
    fluid_viscosity: Optional[float] = None,
    fluid_viscosity_cp: Optional[float] = None,
    fluid_name: Optional[str] = None,
    temperature_c: Optional[float] = None,
    pressure_bar: Optional[float] = None,

    # Liquid-specific valve parameters
    fluid_saturation_pressure_pa: Optional[float] = None,
    fluid_saturation_pressure_psi: Optional[float] = None,
    fluid_critical_pressure_pa: Optional[float] = None,
    fluid_critical_pressure_psi: Optional[float] = None,
    valve_fl: Optional[float] = None,
    valve_fd: Optional[float] = None,

    # Gas-specific parameters
    inlet_temperature_c: Optional[float] = None,
    gas_mw: Optional[float] = None,
    gas_gamma: Optional[float] = None,
    gas_z_factor: Optional[float] = None,
    gas_viscosity: Optional[float] = None,
    valve_xt: Optional[float] = None,
    gas_composition_mol: Optional[dict] = None,
    allow_property_defaults: bool = False,

    # Common valve parameters
    valve_type: Literal["globe", "butterfly", "ball", "gate", "other"] = "globe",
) -> dict:
    """Unified control valve sizing calculations."""
    try:
        request = ControlValveInput(**locals())
        fn = calculate_liquid_control_valve if request.phase == "liquid" else calculate_gas_control_valve
        params = request.model_dump(exclude_none=True)
        params.pop("phase", None)
        sig = inspect.signature(fn)
        allowed = set(sig.parameters.keys())
        forwarded = {key: value for key, value in params.items() if key in allowed}
        payload = _decode_tool_json(fn(**forwarded), request.phase)
        return _map_output(payload, request)
    except CalculationError as exc:
        return _error_output(str(phase), exc)
    except ValidationError as exc:
        return _error_output(
            str(phase),
            CalculationError(
                "INVALID_CONTROL_VALVE_INPUT",
                str(exc),
                details={"phase": phase},
            ),
        )
    except Exception as exc:
        return _error_output(
            str(phase),
            CalculationError(
                "CONTROL_VALVE_CALCULATION_FAILED",
                str(exc),
                details={"phase": phase},
            ),
        )
