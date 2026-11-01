"""Unified parameter sweep for various calculations."""

from __future__ import annotations

import inspect
import json
import math
from typing import Annotated, Any, Dict, Literal, Optional

from pydantic import Field, ValidationError

from fluids_mcp.errors import CalculationError
from fluids_mcp.schemas.omnitools import (
    MAX_SWEEP_POINTS,
    CalculationMessage,
    ParameterSweepInput,
    ParameterSweepOutput,
    SweepRow,
)
from fluids_mcp.validation.result import (
    assert_finite_result,
    strict_json_dumps,
    validate_result_payload,
)
from fluids_mcp.tools.blower_compressor import calculate_blower_compressor_requirements
from fluids_mcp.tools.gas_pipe_pressure_drop import calculate_gas_pipe_pressure_drop
from fluids_mcp.tools.pipe_pressure_drop import calculate_pipe_pressure_drop


def _message_from_error(error: CalculationError) -> CalculationMessage:
    return CalculationMessage(
        code=error.code,
        message=error.message,
        severity=error.severity,  # type: ignore[arg-type]
        details=error.details,
        correlation_id=error.correlation_id,
    )


def _safe_range_value(value: Any) -> Any:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and not math.isfinite(float(value)):
        return "<non-finite>"
    return value


def _error_output(
    *,
    variable: str,
    start: Any,
    stop: Any,
    n: Any,
    error: CalculationError,
) -> dict:
    output = ParameterSweepOutput(
        sweep_variable=variable,
        sweep_range={
            "start": _safe_range_value(start),
            "stop": _safe_range_value(stop),
            "n": _safe_range_value(n),
        },
        rows=[],
        summary={"total_points": 0, "successful_points": 0, "failed_points": 0},
        errors=[_message_from_error(error)],
    )
    payload = output.model_dump(exclude_none=True)
    assert_finite_result(payload)
    strict_json_dumps(payload)
    return payload


def _sweep_values(start: float, stop: float, n: int) -> list[float]:
    if n == 1:
        return [float(start)]
    step = (float(stop) - float(start)) / float(n - 1)
    return [float(start) + idx * step for idx in range(n)]


def _decode_tool_json(payload: str, calculation: str) -> dict:
    try:
        decoded = json.loads(payload)
    except json.JSONDecodeError as exc:
        raise CalculationError(
            "INVALID_TOOL_RESPONSE",
            f"Parameter-sweep backend returned invalid JSON: {exc}",
            details={"calculation": calculation},
        ) from exc
    if not isinstance(decoded, dict):
        raise CalculationError(
            "INVALID_TOOL_RESPONSE",
            "Parameter-sweep backend returned a non-object response",
            details={"calculation": calculation, "response_type": type(decoded).__name__},
        )
    return decoded


def _raise_if_backend_error(payload: dict, calculation: str) -> None:
    if payload.get("error_code") or payload.get("error"):
        details = {"calculation": calculation}
        errors = payload.get("errors")
        if isinstance(errors, list) and errors and isinstance(errors[0], dict):
            raw_details = errors[0].get("details", {})
            if isinstance(raw_details, dict):
                details.update(raw_details)
        raise CalculationError(
            str(payload.get("error_code") or "CALCULATION_FAILED"),
            str(payload.get("error") or payload.get("message") or "Calculation failed"),
            details=details,
        )
    if payload.get("errors"):
        errors = payload["errors"]
        if isinstance(errors, list) and errors and isinstance(errors[0], dict):
            first = errors[0]
            raise CalculationError(
                str(first.get("code") or "CALCULATION_FAILED"),
                str(first.get("message") or first),
                details=first.get("details", {}) if isinstance(first.get("details", {}), dict) else {},
            )
        message = "; ".join(str(item) for item in errors) if isinstance(errors, list) else str(errors)
        raise CalculationError(
            "CALCULATION_FAILED",
            message,
            details={"calculation": calculation},
        )


def _calculator_for(calculation: str):
    if calculation == "pipe_liquid":
        return calculate_pipe_pressure_drop
    if calculation == "pipe_gas":
        return calculate_gas_pipe_pressure_drop
    if calculation == "blower":
        return calculate_blower_compressor_requirements
    raise CalculationError(
        "INVALID_CALCULATION",
        f"Invalid calculation: {calculation}",
        details={"calculation": calculation},
    )


def _base_params(request: ParameterSweepInput) -> dict[str, Any]:
    params = request.model_dump(exclude_none=True)
    for key in ("calculation", "variable", "start", "stop", "n"):
        params.pop(key, None)
    return params


def _validate_sweep_variable(request: ParameterSweepInput, calculator: Any) -> None:
    allowed = set(inspect.signature(calculator).parameters)
    if request.calculation == "blower" and request.variable == "pressure_ratio":
        return
    if request.variable not in allowed:
        raise CalculationError(
            "INVALID_SWEEP_VARIABLE",
            f"{request.variable} is not a supported sweep variable for {request.calculation}",
            details={"variable": request.variable, "calculation": request.calculation},
        )


def _prevalidate_row(request: ParameterSweepInput, value: float) -> None:
    if (
        request.calculation == "pipe_liquid"
        and request.variable in {"flow_rate", "flow_rate_gpm"}
        and value == 0
    ):
        raise CalculationError(
            "ZERO_FLOW",
            "Zero-flow pipe_liquid sweep rows are reported as failed rows instead of "
            "calling the pipe-liquid pressure-drop backend, which divides by zero.",
            details={"variable": request.variable, "input_value": value},
        )


def _apply_sweep_value(request: ParameterSweepInput, params: dict[str, Any], value: float) -> dict[str, Any]:
    row_params = params.copy()
    if request.calculation == "blower" and request.variable == "pressure_ratio":
        if row_params.get("inlet_pressure") is not None:
            row_params["outlet_pressure"] = row_params["inlet_pressure"] * value
        elif row_params.get("inlet_pressure_psi") is not None:
            row_params["outlet_pressure_psi"] = row_params["inlet_pressure_psi"] * value
        else:
            row_params["inlet_pressure"] = 101325.0
            row_params["outlet_pressure"] = 101325.0 * value
    else:
        row_params[request.variable] = value
    return row_params


def _filter_params(calculator: Any, params: dict[str, Any]) -> dict[str, Any]:
    allowed = set(inspect.signature(calculator).parameters)
    return {key: value for key, value in params.items() if key in allowed}


def _successful_row(request: ParameterSweepInput, value: float, result: dict[str, Any]) -> SweepRow:
    validate_result_payload(result, path=f"rows[{value}].result")
    assert_finite_result(result, path=f"rows[{value}].result")
    return SweepRow(ok=True, input_value=value, result=result)


def _failed_row(value: float, error: CalculationError) -> SweepRow:
    payload = {
        "ok": False,
        "input_value": value,
        "error_code": error.code,
        "message": error.message,
    }
    assert_finite_result(payload)
    return SweepRow(**payload)


def _run_row(
    *,
    request: ParameterSweepInput,
    calculator: Any,
    base_params: dict[str, Any],
    value: float,
) -> SweepRow:
    try:
        _prevalidate_row(request, value)
        row_params = _apply_sweep_value(request, base_params, value)
        result_json = calculator(**_filter_params(calculator, row_params))
        result = _decode_tool_json(result_json, request.calculation)
        _raise_if_backend_error(result, request.calculation)
        return _successful_row(request, value, result)
    except CalculationError as exc:
        return _failed_row(value, exc)
    except Exception as exc:
        return _failed_row(
            value,
            CalculationError(
                "CALCULATION_FAILED",
                str(exc),
                details={"calculation": request.calculation, "variable": request.variable},
            ),
        )


def parameter_sweep(
    calculation: Literal["pipe_liquid", "pipe_gas", "blower"],
    variable: str,
    start: float,
    stop: float,
    n: Annotated[int, Field(ge=1, le=MAX_SWEEP_POINTS, description="Number of sweep points, from 1 to 100.")],

    # Pipe flow parameters (liquid and gas)
    flow_rate: Optional[float] = None,
    flow_rate_gpm: Optional[float] = None,
    flow_rate_kg_s: Optional[float] = None,
    flow_rate_norm_m3_hr: Optional[float] = None,
    flow_rate_std_m3_hr: Optional[float] = None,
    pipe_diameter: Optional[float] = None,
    pipe_diameter_in: Optional[float] = None,
    pipe_length: Optional[float] = None,
    pipe_length_ft: Optional[float] = None,
    nominal_size_in: Optional[float] = None,
    schedule: str = "40",
    material: Optional[str] = None,
    pipe_roughness: Optional[float] = None,

    # Fluid properties
    fluid_density: Optional[float] = None,
    fluid_density_lbft3: Optional[float] = None,
    fluid_viscosity: Optional[float] = None,
    fluid_viscosity_cp: Optional[float] = None,
    fluid_name: Optional[str] = None,
    temperature_c: Optional[float] = None,
    pressure_bar: Optional[float] = None,

    # Gas-specific parameters
    inlet_pressure: Optional[float] = None,
    outlet_pressure: Optional[float] = None,
    gas_mw: Optional[float] = None,
    gas_gamma: Optional[float] = None,
    gas_z_factor: Optional[float] = None,
    gas_viscosity: Optional[float] = None,
    gas_composition_mol: Optional[Dict[str, float]] = None,
    method: str = "Weymouth",

    # Blower-specific parameters
    inlet_pressure_psi: Optional[float] = None,
    outlet_pressure_psi: Optional[float] = None,
    inlet_temperature_c: Optional[float] = None,
    efficiency: Optional[float] = None,
    efficiency_type: str = "isentropic",

    # Compressor calculation flags
    allow_property_defaults: bool = False,
    use_eos_calculations: Optional[bool] = None,
    solve_for: Optional[str] = None,
) -> dict:
    """Unified parameter sweep calculations.

    Performs parameter sweeps for different calculation types:
    - calculation='pipe_liquid': Liquid pipe pressure drop sweep
    - calculation='pipe_gas': Gas pipe pressure drop sweep
    - calculation='blower': Blower/compressor performance sweep

    Each request accepts 1 to 100 points. Larger studies require a separate
    batch workflow rather than automatic repeated calls around this limit.

    Returns:
        Structured sweep result matching ParameterSweepOutput.
    """
    try:
        request = ParameterSweepInput(**locals())
        calculator = _calculator_for(request.calculation)
        _validate_sweep_variable(request, calculator)

        base_params = _base_params(request)
        rows = [
            _run_row(
                request=request,
                calculator=calculator,
                base_params=base_params,
                value=value,
            )
            for value in _sweep_values(request.start, request.stop, request.n)
        ]
        success_count = sum(1 for row in rows if row.ok)
        output = ParameterSweepOutput(
            sweep_variable=request.variable,
            sweep_range={"start": request.start, "stop": request.stop, "n": request.n},
            rows=rows,
            summary={
                "total_points": len(rows),
                "successful_points": success_count,
                "failed_points": len(rows) - success_count,
            },
        )
        payload = output.model_dump(exclude_none=True)
        validate_result_payload(payload)
        assert_finite_result(payload)
        strict_json_dumps(payload)
        return payload
    except CalculationError as exc:
        return _error_output(variable=variable, start=start, stop=stop, n=n, error=exc)
    except ValidationError as exc:
        return _error_output(
            variable=variable,
            start=start,
            stop=stop,
            n=n,
            error=CalculationError(
                "INVALID_PARAMETER_SWEEP_INPUT",
                "Invalid parameter_sweep input; n must be an integer from 1 to 100 and all values must be valid and finite.",
                details={"calculation": calculation, "variable": variable},
            ),
        )
    except Exception as exc:
        return _error_output(
            variable=variable,
            start=start,
            stop=stop,
            n=n,
            error=CalculationError(
                "PARAMETER_SWEEP_FAILED",
                str(exc),
                details={"calculation": calculation, "variable": variable},
            ),
        )
