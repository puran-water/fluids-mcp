"""
Unified pipe flow calculations for both liquid and gas phases.
"""

import inspect
import json
from typing import Any, Dict, Literal, Optional

from pydantic import ValidationError

from fluids_mcp.errors import CalculationError
from fluids_mcp.schemas.omnitools import CalculationMessage, PipeFlowInput, PipeFlowOutput
from fluids_mcp.validation.input import require_finite, require_one_of, require_positive
from fluids_mcp.validation.physical import (
    check_choked_flow,
    validate_dimension,
    validate_flow_rate,
    validate_forward_pressure,
    validate_pressure_drop,
)
from fluids_mcp.validation.result import (
    assert_finite_result,
    validate_forward_pressure_result,
)
from fluids_mcp.tools.pipe_pressure_drop import calculate_pipe_pressure_drop
from fluids_mcp.tools.gas_pipe_pressure_drop import calculate_gas_pipe_pressure_drop


VALID_PHASES = {"liquid", "gas"}
DEFAULT_PHASE = "liquid"
LOW_MACH_LIMIT = 0.3

GAS_METHOD_APPLICABILITY = {
    "Weymouth": (
        "empirical long-distance, high-pressure gas transmission; generally large-diameter pipe (>8 in)"
    ),
    "Panhandle_A": "long gas transmission pipelines at moderate turbulent Reynolds numbers",
    "Panhandle_B": "long gas transmission pipelines in fully turbulent, high-Reynolds-number flow",
    "IGT": "gas pipelines where the viscosity-dependent empirical correlation has been validated",
    "Oliphant": "low-pressure gas distribution below about 100 kPa gauge",
    "Spitzglass_low": "very-low-pressure gas distribution below about 10 kPa gauge",
    "Spitzglass_high": "low-pressure gas distribution up to about 100 kPa gauge",
    "isothermal_darcy": "single-phase, isothermal, unchoked low-Mach flow using Darcy friction",
}


def _message_from_error(error: CalculationError) -> CalculationMessage:
    return CalculationMessage(
        code=error.code,
        message=error.message,
        severity=error.severity,  # type: ignore[arg-type]
        details=error.details,
        correlation_id=error.correlation_id,
    )


def _error_output(phase: str, error: CalculationError) -> dict:
    normalized_phase = phase if phase in VALID_PHASES else DEFAULT_PHASE
    return PipeFlowOutput(
        phase=normalized_phase,  # type: ignore[arg-type]
        errors=[_message_from_error(error)],
    ).model_dump(exclude_none=True)


def _warning_messages(raw_warnings: Any) -> list[CalculationMessage]:
    if not raw_warnings:
        return []
    warnings = raw_warnings if isinstance(raw_warnings, list) else [raw_warnings]
    return [
        CalculationMessage(code="BACKEND_WARNING", message=str(message), severity="warning")
        for message in warnings
    ]


def _gas_method_warning(method: str) -> CalculationMessage:
    applicability = GAS_METHOD_APPLICABILITY[method]
    recommendation = (
        " For compact low-pressure process or aeration headers, use isothermal_darcy unless "
        "the selected empirical correlation is demonstrably applicable."
        if method != "isothermal_darcy"
        else " Verify the isothermal, single-phase, and low-Mach assumptions for the installation."
    )
    return CalculationMessage(
        code="GAS_METHOD_APPLICABILITY",
        message=f"Gas pipe_flow used {method}. Applicability: {applicability}.{recommendation}",
        severity="warning",
        details={"method": method, "applicability": applicability},
    )


def _decode_tool_json(payload: str, phase: str) -> dict:
    try:
        decoded = json.loads(payload)
    except json.JSONDecodeError as exc:
        raise CalculationError(
            "INVALID_TOOL_RESPONSE",
            f"Pipe-flow backend returned invalid JSON: {exc}",
            details={"phase": phase},
        ) from exc
    if not isinstance(decoded, dict):
        raise CalculationError(
            "INVALID_TOOL_RESPONSE",
            "Pipe-flow backend returned a non-object response",
            details={"phase": phase, "response_type": type(decoded).__name__},
        )
    return decoded


def _raise_if_backend_error(payload: dict, phase: str) -> None:
    if payload.get("error_code") or payload.get("error"):
        details = {"phase": phase}
        errors = payload.get("errors")
        if isinstance(errors, list) and errors and isinstance(errors[0], dict):
            raw_details = errors[0].get("details", {})
            if isinstance(raw_details, dict):
                details.update(raw_details)
        raise CalculationError(
            str(payload.get("error_code") or "PIPE_FLOW_CALCULATION_FAILED"),
            str(payload.get("error") or "Pipe-flow calculation failed"),
            details=details,
        )
    if payload.get("errors"):
        errors = payload["errors"]
        if isinstance(errors, list) and errors and isinstance(errors[0], dict):
            first = errors[0]
            raise CalculationError(
                str(first.get("code") or "PIPE_FLOW_CALCULATION_FAILED"),
                str(first.get("message") or first),
                details=first.get("details", {}) if isinstance(first.get("details", {}), dict) else {},
            )
        message = "; ".join(str(item) for item in errors) if isinstance(errors, list) else str(errors)
        raise CalculationError(
            "PIPE_FLOW_CALCULATION_FAILED",
            message,
            details={"phase": phase},
        )


def _is_zero_liquid_flow(request: PipeFlowInput) -> bool:
    if request.phase != "liquid":
        return False
    return request.flow_rate == 0 or request.flow_rate_gpm == 0


def _is_zero_gas_flow(request: PipeFlowInput) -> bool:
    if request.phase != "gas":
        return False
    return (
        request.flow_rate_kg_s == 0
        or request.flow_rate_norm_m3_hr == 0
        or request.flow_rate_std_m3_hr == 0
    )


def _reject_zero_flow(request: PipeFlowInput) -> None:
    raise CalculationError(
        "ZERO_FLOW",
        f"Zero {request.phase} flow is not supported by pipe_flow; provide a positive flow rate.",
        details={"phase": request.phase},
    )


def _gas_gamma_from_payload(payload: dict, request: PipeFlowInput) -> Optional[float]:
    if request.gas_gamma is not None:
        return request.gas_gamma
    details = payload.get("gas_property_details") or {}
    gamma = details.get("gamma")
    return gamma if isinstance(gamma, (int, float)) else None


def _gas_expansion_factor(pressure_ratio: float, gamma: float) -> float:
    # Isentropic density ratio; useful as a bounded indicator for expansion severity.
    return pressure_ratio ** (1.0 / gamma)


def _gas_choked_fields(payload: dict, request: PipeFlowInput) -> dict:
    p1 = payload.get("inlet_pressure_pa", request.inlet_pressure)
    p2 = payload.get("outlet_pressure_pa", request.outlet_pressure)
    gamma = _gas_gamma_from_payload(payload, request)
    if p1 is None or p2 is None or gamma is None:
        return {}
    check = check_choked_flow(float(p1), float(p2), float(gamma))
    return {
        "choked": check.is_choked,
        "pressure_ratio": check.pressure_ratio,
        "critical_pressure_ratio": check.critical_pressure_ratio,
        "expansion_factor": _gas_expansion_factor(check.pressure_ratio, float(gamma)),
    }


def _raise_if_gas_flow_out_of_scope(payload: dict, request: PipeFlowInput) -> None:
    if request.phase != "gas":
        return

    scope_fields = _gas_choked_fields(payload, request)
    mach = payload.get("average_mach_number")
    details = {
        "phase": "gas",
        "method": request.method,
        "max_mach": LOW_MACH_LIMIT,
        **scope_fields,
    }
    reasons = []

    if mach is not None:
        mach_value = float(mach)
        details["average_mach_number"] = mach_value
        if mach_value > LOW_MACH_LIMIT:
            reasons.append(f"average Mach number {mach_value:.3f} exceeds {LOW_MACH_LIMIT:.1f}")

    if scope_fields.get("choked") is True:
        reasons.append("pressure ratio indicates choked flow")

    if reasons:
        raise CalculationError(
            "OUT_OF_SCOPE",
            "Gas pipe_flow is limited to unchoked low-Mach flow (Mach <= 0.3); "
            f"{'; '.join(reasons)}.",
            details=details,
        )


def _first_present(payload: dict, *keys: str) -> Any:
    for key in keys:
        if payload.get(key) is not None:
            return payload[key]
    return None


def _map_backend_output(payload: dict, request: PipeFlowInput) -> dict:
    warnings = _warning_messages(payload.get("warnings"))
    if request.phase == "gas":
        warnings.insert(0, _gas_method_warning(request.method))
    common = {
        "phase": request.phase,
        "solved_variable": payload.get("solved_variable"),
        "flow_rate_m3s": payload.get("flow_rate_m3s"),
        "flow_rate_gpm": payload.get("flow_rate_gpm"),
        "flow_rate_kg_s": payload.get("flow_rate_kg_s"),
        "flow_rate_norm_m3_hr": payload.get("flow_rate_norm_m3_hr"),
        "inlet_pressure_pa": payload.get("inlet_pressure_pa"),
        "outlet_pressure_pa": payload.get("outlet_pressure_pa"),
        "pressure_drop_pa": _first_present(payload, "pressure_drop_pa", "pressure_drop_total_pa"),
        "pressure_drop_total_pa": _first_present(payload, "pressure_drop_total_pa", "pressure_drop_pa"),
        "pipe_diameter_m": payload.get("pipe_diameter_m"),
        "pipe_length_m": payload.get("pipe_length_m"),
        "reynolds_number": payload.get("reynolds_number"),
        "friction_factor": payload.get("friction_factor"),
        "flow_velocity_m_s": payload.get("flow_velocity_m_s"),
        "average_velocity_m_s": payload.get("average_velocity_m_s"),
        "average_mach_number": payload.get("average_mach_number"),
        "head_loss_m": payload.get("head_loss_m"),
        "head_loss_ft": payload.get("head_loss_ft"),
        "warnings": warnings,
        "inputs_resolved": payload.get("inputs_resolved", []),
        "details": payload,
    }
    if request.phase == "gas":
        common.update(_gas_choked_fields(payload, request))
    output = PipeFlowOutput(**common)
    result = output.model_dump(exclude_none=True)
    assert_finite_result(result)
    if request.phase == "gas":
        validate_forward_pressure_result(result)
    return result


def _validate_common_inputs(request: PipeFlowInput) -> None:
    validate_dimension("pipe_diameter", request.pipe_diameter)
    validate_dimension("pipe_diameter_in", request.pipe_diameter_in)
    validate_dimension("nominal_size_in", request.nominal_size_in)
    validate_dimension("pipe_length", request.pipe_length)
    validate_dimension("pipe_length_ft", request.pipe_length_ft)
    validate_pressure_drop(request.pressure_drop, allow_zero=True)
    validate_pressure_drop(request.pressure_drop_psi, allow_zero=True)
    require_finite("temperature_c", request.temperature_c)


def _validate_liquid_inputs(request: PipeFlowInput) -> None:
    _validate_common_inputs(request)
    require_one_of(
        {
            "flow_rate": request.flow_rate,
            "flow_rate_gpm": request.flow_rate_gpm,
            "pressure_drop": request.pressure_drop,
            "pressure_drop_psi": request.pressure_drop_psi,
        },
        group_name="liquid flow or pressure-drop input",
    )
    validate_flow_rate("flow_rate", request.flow_rate, allow_zero=True)
    validate_flow_rate("flow_rate_gpm", request.flow_rate_gpm, allow_zero=True)
    require_positive("fluid_density", request.fluid_density)
    require_positive("fluid_density_lbft3", request.fluid_density_lbft3)
    require_positive("fluid_viscosity", request.fluid_viscosity)
    require_positive("fluid_viscosity_cp", request.fluid_viscosity_cp)


def _validate_gas_inputs(request: PipeFlowInput) -> None:
    _validate_common_inputs(request)
    if request.temperature_c is None:
        raise CalculationError(
            "MISSING_REQUIRED_INPUT",
            "temperature_c is required for gas flow calculations",
        )
    require_one_of(
        {
            "flow_rate_kg_s": request.flow_rate_kg_s,
            "flow_rate_norm_m3_hr": request.flow_rate_norm_m3_hr,
            "flow_rate_std_m3_hr": request.flow_rate_std_m3_hr,
        },
        group_name="gas flow input",
    )
    validate_flow_rate("flow_rate_kg_s", request.flow_rate_kg_s, allow_zero=True)
    validate_flow_rate("flow_rate_norm_m3_hr", request.flow_rate_norm_m3_hr, allow_zero=True)
    validate_flow_rate("flow_rate_std_m3_hr", request.flow_rate_std_m3_hr, allow_zero=True)
    if request.inlet_pressure is not None and request.outlet_pressure is not None:
        validate_forward_pressure(
            request.inlet_pressure,
            request.outlet_pressure,
            context="gas pipe_flow",
        )


def pipe_flow(
    phase: Literal["liquid", "gas"] = "liquid",

    # Common parameters
    pipe_length: Optional[float] = None,
    pipe_diameter: Optional[float] = None,
    pipe_diameter_in: Optional[float] = None,
    pipe_length_ft: Optional[float] = None,
    nominal_size_in: Optional[float] = None,
    schedule: str = "40",
    material: Optional[str] = None,
    pipe_roughness: Optional[float] = None,
    fittings: Optional[list] = None,
    temperature_c: Optional[float] = None,
    fluid_name: Optional[str] = None,

    # Liquid-specific parameters
    flow_rate: Optional[float] = None,
    flow_rate_gpm: Optional[float] = None,
    fluid_density: Optional[float] = None,
    fluid_density_lbft3: Optional[float] = None,
    fluid_viscosity: Optional[float] = None,
    fluid_viscosity_cp: Optional[float] = None,
    pressure_drop: Optional[float] = None,
    pressure_drop_psi: Optional[float] = None,
    pressure_bar: Optional[float] = None,

    # Gas-specific parameters
    inlet_pressure: Optional[float] = None,
    outlet_pressure: Optional[float] = None,
    flow_rate_kg_s: Optional[float] = None,
    flow_rate_norm_m3_hr: Optional[float] = None,
    flow_rate_std_m3_hr: Optional[float] = None,
    gas_mw: Optional[float] = None,
    gas_gamma: Optional[float] = None,
    gas_z_factor: Optional[float] = None,
    gas_viscosity: Optional[float] = None,
    gas_composition_mol: Optional[Dict[str, float]] = None,
    method: Literal[
        "Weymouth", "Panhandle_A", "Panhandle_B", "IGT", "Oliphant",
        "Spitzglass_low", "Spitzglass_high", "isothermal_darcy"
    ] = "Weymouth",
    allow_property_defaults: bool = False,
) -> dict:
    """
    Unified pipe flow calculations for liquid and gas phases.

    This omnitool consolidates pipe pressure drop calculations:
    - phase='liquid': Uses incompressible flow equations (Darcy-Weisbach)
    - phase='gas': Uses compressible flow equations (various methods)

    Args:
        phase: "liquid" or "gas" - determines calculation method

        Common parameters (both phases):
            pipe_length: Pipe length in meters
            pipe_length_ft: Pipe length in feet
            pipe_diameter: Inner diameter in meters
            pipe_diameter_in: Inner diameter in inches
            nominal_size_in: Nominal pipe size in inches
            schedule: Pipe schedule (default "40")
            material: Pipe material for roughness lookup
            pipe_roughness: Absolute roughness in meters
            fittings: List of fittings with type and quantity
            temperature_c: Temperature in Celsius
            fluid_name: Fluid name for property lookup

        Liquid-specific parameters:
            flow_rate: Flow rate in m³/s
            flow_rate_gpm: Flow rate in US GPM
            fluid_density: Density in kg/m³
            fluid_density_lbft3: Density in lb/ft³
            fluid_viscosity: Dynamic viscosity in Pa·s
            fluid_viscosity_cp: Viscosity in centipoise
            pressure_drop: Pressure drop in Pa (for solving)
            pressure_drop_psi: Pressure drop in psi
            pressure_bar: Pressure in bar (for properties)

        Gas-specific parameters:
            inlet_pressure: Inlet pressure in Pa (absolute)
            outlet_pressure: Outlet pressure in Pa (absolute)
            flow_rate_kg_s: Mass flow rate in kg/s
            flow_rate_norm_m3_hr: Normal volumetric flow in m³/hr
            flow_rate_std_m3_hr: Standard volumetric flow in m³/hr
            gas_mw: Molecular weight in kg/kmol
            gas_gamma: Specific heat ratio (Cp/Cv)
            gas_z_factor: Compressibility factor
            gas_viscosity: Dynamic viscosity in Pa·s
            gas_composition_mol: Composition as mol fractions
            method: Calculation method for gas flow

    Returns:
        Structured pipe-flow result matching PipeFlowOutput

    Examples:
        Liquid flow:
        >>> pipe_flow(phase="liquid", flow_rate_gpm=100, pipe_length_ft=100,
        ...          nominal_size_in=2, fluid_name="Water", temperature_c=25)

        Gas flow:
        >>> pipe_flow(phase="gas", flow_rate_norm_m3_hr=150, inlet_pressure=101325,
        ...          pipe_length=25, nominal_size_in=3, fluid_name="Methane")
    """
    try:
        request = PipeFlowInput(**locals())
        if request.phase == "liquid":
            _validate_liquid_inputs(request)
            if _is_zero_liquid_flow(request):
                _reject_zero_flow(request)
            fn = calculate_pipe_pressure_drop
        else:
            _validate_gas_inputs(request)
            if _is_zero_gas_flow(request):
                _reject_zero_flow(request)
            fn = calculate_gas_pipe_pressure_drop

        params = request.model_dump(exclude_none=True)
        params.pop("phase", None)
        sig = inspect.signature(fn)
        allowed = set(sig.parameters.keys())
        forwarded = {key: value for key, value in params.items() if key in allowed}
        backend_payload = _decode_tool_json(fn(**forwarded), request.phase)
        _raise_if_backend_error(backend_payload, request.phase)
        _raise_if_gas_flow_out_of_scope(backend_payload, request)
        return _map_backend_output(backend_payload, request)
    except CalculationError as exc:
        return _error_output(str(phase), exc)
    except ValidationError as exc:
        return _error_output(
            str(phase),
            CalculationError(
                "INVALID_PIPE_FLOW_INPUT",
                str(exc),
                details={"phase": phase},
            ),
        )
    except Exception as exc:
        return _error_output(
            str(phase),
            CalculationError(
                "PIPE_FLOW_CALCULATION_FAILED",
                str(exc),
                details={"phase": phase},
            ),
        )
