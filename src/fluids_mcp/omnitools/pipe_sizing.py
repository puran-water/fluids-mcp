import inspect
import json
from typing import Any, Literal, Optional

from pydantic import ValidationError

from fluids_mcp.errors import CalculationError
from fluids_mcp.schemas.omnitools import CalculationMessage, PipeSizingInput, PipeSizingOutput
from fluids_mcp.validation.result import assert_finite_result, strict_json_dumps, validate_result_payload
from fluids_mcp.tools.pipe_sizing import select_liquid_pipe_size_for_dp, select_gas_pipe_size_for_dp


VALID_PHASES = {"liquid", "gas"}
DEFAULT_PHASE = "liquid"


def _message_from_error(error: CalculationError) -> CalculationMessage:
    return CalculationMessage(
        code=error.code,
        message=error.message,
        severity=error.severity,  # type: ignore[arg-type]
        details=error.details,
        correlation_id=error.correlation_id,
    )


def _error_output(phase: str, error: CalculationError) -> dict:
    normalized = phase if phase in VALID_PHASES else DEFAULT_PHASE
    output = PipeSizingOutput(
        phase=normalized,  # type: ignore[arg-type]
        errors=[_message_from_error(error)],
    )
    payload = output.model_dump(exclude_none=True)
    assert_finite_result(payload)
    strict_json_dumps(payload)
    return payload


def _warning_messages(raw_warnings: Any) -> list[CalculationMessage]:
    if not raw_warnings:
        return []
    warnings = raw_warnings if isinstance(raw_warnings, list) else [raw_warnings]
    return [
        CalculationMessage(code="PIPE_SIZING_WARNING", message=str(message), severity="warning")
        for message in warnings
    ]


def _decode_tool_json(payload: str, phase: str) -> dict:
    try:
        decoded = json.loads(payload)
    except json.JSONDecodeError as exc:
        raise CalculationError(
            "INVALID_TOOL_RESPONSE",
            f"Pipe-sizing backend returned invalid JSON: {exc}",
            details={"phase": phase},
        ) from exc
    if not isinstance(decoded, dict):
        raise CalculationError(
            "INVALID_TOOL_RESPONSE",
            "Pipe-sizing backend returned a non-object response",
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
            str(payload.get("error_code") or "PIPE_SIZING_CALCULATION_FAILED"),
            str(payload.get("error") or "Pipe-sizing calculation failed"),
            details=details,
        )
    errors = payload.get("errors")
    if not errors:
        return None
    if isinstance(errors, list):
        if errors and isinstance(errors[0], dict):
            first = errors[0]
            return CalculationError(
                str(first.get("code") or "PIPE_SIZING_CALCULATION_FAILED"),
                str(first.get("message") or first),
                details=first.get("details", {}) if isinstance(first.get("details", {}), dict) else {},
            )
        return CalculationError(
            "PIPE_SIZING_CALCULATION_FAILED",
            "; ".join(str(item) for item in errors),
            details={"phase": phase},
        )
    return CalculationError("PIPE_SIZING_CALCULATION_FAILED", str(errors), details={"phase": phase})


def _backend_context(payload: dict) -> dict:
    return {
        key: value
        for key, value in payload.items()
        if key not in {"error", "error_code", "errors"} and value is not None
    }


def _map_output(payload: dict, request: PipeSizingInput) -> dict:
    validate_result_payload(payload)
    error = _payload_error(payload, request.phase)
    if error is not None:
        output = PipeSizingOutput(
            phase=request.phase,
            errors=[_message_from_error(error)],
            warnings=_warning_messages(payload.get("warnings")),
            inputs_resolved=payload.get("inputs_resolved") or payload.get("log") or [],
            details=_backend_context(payload),
        )
    else:
        output = PipeSizingOutput(
            phase=request.phase,
            recommended_nps_in=payload.get("recommended_nps_in"),
            recommended_details=payload.get("recommended_details") or {},
            feasible_candidates=payload.get("feasible_candidates") or [],
            all_candidates=payload.get("all_candidates") or [],
            warnings=_warning_messages(payload.get("warnings")),
            inputs_resolved=payload.get("inputs_resolved") or payload.get("log") or [],
            details=payload,
        )
    result = output.model_dump(exclude_none=True)
    validate_result_payload(result)
    assert_finite_result(result)
    strict_json_dumps(result)
    return result


def pipe_sizing(
    phase: Literal["liquid", "gas"] = "liquid",

    # Flow parameters
    flow_rate: Optional[float] = None,
    flow_rate_gpm: Optional[float] = None,
    flow_rate_kg_s: Optional[float] = None,
    flow_rate_norm_m3_hr: Optional[float] = None,
    flow_rate_std_m3_hr: Optional[float] = None,

    # Pressure constraints
    allowable_dp_pa: Optional[float] = None,
    allowable_dp_psi: Optional[float] = None,
    inlet_pressure: Optional[float] = None,
    inlet_pressure_psi: Optional[float] = None,

    # Pipe specifications
    pipe_length: Optional[float] = None,
    pipe_length_ft: Optional[float] = None,
    schedule: str = "40",
    nps_candidates: Optional[list] = None,
    material: Optional[str] = None,
    pipe_roughness: Optional[float] = None,
    fittings: Optional[list] = None,

    # Fluid properties
    fluid_density: Optional[float] = None,
    fluid_density_lbft3: Optional[float] = None,
    fluid_viscosity: Optional[float] = None,
    fluid_viscosity_cp: Optional[float] = None,
    fluid_name: Optional[str] = None,
    temperature_c: Optional[float] = None,
    pressure_bar: Optional[float] = None,

    # Gas-specific parameters
    gas_mw: Optional[float] = None,
    gas_gamma: Optional[float] = None,
    gas_z_factor: Optional[float] = None,
    gas_viscosity: Optional[float] = None,
    method: str = "isothermal_darcy",

    # Constraints
    velocity_min_m_s: Optional[float] = None,
    velocity_max_m_s: Optional[float] = None,

    # Output options
    return_all_candidates: bool = False,
) -> dict:
    """Unified pipe sizing optimization.

    Finds optimal pipe size based on pressure drop constraints:
    - phase='liquid': Uses Darcy-Weisbach for incompressible flow
    - phase='gas': Uses compressible flow with equivalent length iteration

    Returns:
        Structured pipe sizing result matching PipeSizingOutput
    """
    try:
        request = PipeSizingInput(**locals())
        fn = select_liquid_pipe_size_for_dp if request.phase == "liquid" else select_gas_pipe_size_for_dp
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
                "INVALID_PIPE_SIZING_INPUT",
                str(exc),
                details={"phase": phase},
            ),
        )
    except Exception as exc:
        return _error_output(
            str(phase),
            CalculationError(
                "PIPE_SIZING_CALCULATION_FAILED",
                str(exc),
                details={"phase": phase},
            ),
        )
