import inspect
import json
from typing import Any, Literal, Optional

from pydantic import ValidationError

from fluids_mcp.domain.units import BAR_TO_PA, PSI_TO_PA
from fluids_mcp.errors import CalculationError
from fluids_mcp.schemas.omnitools import (
    CalculationMessage,
    OrificeRestrictorInput,
    OrificeRestrictorOutput,
)
from fluids_mcp.validation.result import assert_finite_result, strict_json_dumps, validate_result_payload
from fluids_mcp.tools.orifice_restrictor import (
    calculate_orifice_sizing,
    calculate_orifice_analysis_dp,
    calculate_orifice_analysis_flow,
    generate_plate_kit,
)


def _message_from_error(error: CalculationError) -> CalculationMessage:
    return CalculationMessage(
        code=error.code,
        message=error.message,
        severity=error.severity,  # type: ignore[arg-type]
        details=error.details,
        correlation_id=error.correlation_id,
    )


def _error_output(mode: str, phase: str, error: CalculationError) -> dict:
    normalized_mode = mode if mode in {"sizing", "analysis_dp", "analysis_flow", "plate_kit"} else None
    normalized_phase = phase if phase in {"liquid", "gas"} else None
    output = OrificeRestrictorOutput(
        mode=normalized_mode,  # type: ignore[arg-type]
        phase=normalized_phase,  # type: ignore[arg-type]
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
        CalculationMessage(code="ORIFICE_RESTRICTOR_WARNING", message=str(message), severity="warning")
        for message in warnings
    ]


def _decode_tool_json(payload: str, mode: str, phase: str) -> dict:
    try:
        decoded = json.loads(payload)
    except json.JSONDecodeError as exc:
        raise CalculationError(
            "INVALID_TOOL_RESPONSE",
            f"Orifice-restrictor backend returned invalid JSON: {exc}",
            details={"mode": mode, "phase": phase},
        ) from exc
    if not isinstance(decoded, dict):
        raise CalculationError(
            "INVALID_TOOL_RESPONSE",
            "Orifice-restrictor backend returned a non-object response",
            details={"mode": mode, "phase": phase, "response_type": type(decoded).__name__},
        )
    return decoded


def _payload_error(payload: dict, request: OrificeRestrictorInput) -> Optional[CalculationError]:
    if payload.get("error_code") or payload.get("error"):
        return CalculationError(
            str(payload.get("error_code") or "ORIFICE_RESTRICTOR_CALCULATION_FAILED"),
            str(payload.get("error") or "Orifice-restrictor calculation failed"),
            details={"mode": request.mode, "phase": request.phase},
        )
    errors = payload.get("errors")
    if not errors:
        return None
    if isinstance(errors, list):
        if errors and isinstance(errors[0], dict):
            first = errors[0]
            return CalculationError(
                str(first.get("code") or "ORIFICE_RESTRICTOR_CALCULATION_FAILED"),
                str(first.get("message") or first),
                details=first.get("details", {}) if isinstance(first.get("details", {}), dict) else {},
            )
        return CalculationError(
            "ORIFICE_RESTRICTOR_CALCULATION_FAILED",
            "; ".join(str(item) for item in errors),
            details={"mode": request.mode, "phase": request.phase},
        )
    return CalculationError(
        "ORIFICE_RESTRICTOR_CALCULATION_FAILED",
        str(errors),
        details={"mode": request.mode, "phase": request.phase},
    )


def _pressure_drop_pa(payload: dict) -> Optional[float]:
    if payload.get("pressure_drop_pa") is not None:
        return payload["pressure_drop_pa"]
    if payload.get("pressure_drop_bar") is not None:
        return payload["pressure_drop_bar"] * BAR_TO_PA
    if payload.get("pressure_drop_psi") is not None:
        return payload["pressure_drop_psi"] * PSI_TO_PA
    return None


def _map_output(payload: dict, request: OrificeRestrictorInput) -> dict:
    validate_result_payload(payload)
    error = _payload_error(payload, request)
    if error is not None:
        output = OrificeRestrictorOutput(
            mode=request.mode,
            phase=request.phase,
            errors=[_message_from_error(error)],
            warnings=_warning_messages(payload.get("warnings")),
            inputs_resolved=payload.get("inputs_resolved") or payload.get("log") or [],
            details=payload,
        )
    else:
        output = OrificeRestrictorOutput(
            mode=request.mode,
            phase=request.phase,
            orifice_diameter_m=payload.get("orifice_diameter_m"),
            orifice_diameter_mm=payload.get("orifice_diameter_mm") or payload.get("design_orifice_mm"),
            orifice_diameter_in=payload.get("orifice_diameter_in") or payload.get("design_orifice_in"),
            beta_ratio=payload.get("beta_ratio"),
            pressure_drop_pa=_pressure_drop_pa(payload),
            flow_rate_m3_hr=payload.get("flow_rate_m3_hr") or payload.get("design_flow_m3_hr"),
            mass_flow_kg_s=payload.get("mass_flow_kg_s"),
            # Permanent loss is what a restriction plate delivers to the system; the
            # pressure_drop_* fields above are the tap-to-tap differential.
            permanent_pressure_loss_pa=payload.get("permanent_pressure_loss_pa"),
            permanent_pressure_loss_psi=payload.get("permanent_pressure_loss_psi"),
            permanent_pressure_loss_bar=payload.get("permanent_pressure_loss_bar"),
            pressure_recovery_fraction=payload.get("pressure_recovery_fraction"),
            discharge_coefficient=payload.get("discharge_coefficient"),
            cavitation=payload.get("cavitation"),
            choked=payload.get("choked"),
            plates=payload.get("plates") or payload.get("plate_kit") or [],
            warnings=_warning_messages(payload.get("warnings")),
            inputs_resolved=payload.get("inputs_resolved") or payload.get("log") or [],
            details=payload,
        )
    result = output.model_dump(exclude_none=True)
    validate_result_payload(result)
    assert_finite_result(result)
    strict_json_dumps(result)
    return result


def orifice_restrictor(
    # Mode selection
    mode: Literal["sizing", "analysis_dp", "analysis_flow", "plate_kit"] = "sizing",
    phase: Literal["liquid", "gas"] = "liquid",

    # Flow parameters
    flow_rate: Optional[float] = None,
    flow_rate_gpm: Optional[float] = None,
    flow_rate_m3_hr: Optional[float] = None,

    # Pressure parameters
    pressure_drop: Optional[float] = None,
    pressure_drop_psi: Optional[float] = None,
    pressure_drop_bar: Optional[float] = None,
    inlet_pressure: Optional[float] = None,
    inlet_pressure_psi: Optional[float] = None,
    inlet_pressure_bar: Optional[float] = None,

    # Orifice geometry
    orifice_diameter: Optional[float] = None,
    orifice_diameter_in: Optional[float] = None,
    orifice_diameter_mm: Optional[float] = None,

    # Pipe geometry
    pipe_diameter: Optional[float] = None,
    pipe_diameter_in: Optional[float] = None,
    nominal_size_in: Optional[float] = None,
    schedule: str = "40",

    # Fluid properties - direct input
    fluid_density: Optional[float] = None,
    fluid_density_lbft3: Optional[float] = None,
    fluid_viscosity: Optional[float] = None,
    fluid_viscosity_cp: Optional[float] = None,

    # Fluid properties - lookup
    fluid_name: Optional[str] = None,
    temperature_c: Optional[float] = None,
    pressure_bar: Optional[float] = None,

    # Gas-specific properties
    gas_gamma: Optional[float] = None,
    gas_mw: Optional[float] = None,
    gas_z_factor: Optional[float] = None,

    # Orifice parameters
    cd: Optional[float] = None,
    meter_type: Optional[str] = None,
    taps: str = 'D',

    # Restriction-duty parameters
    dp_basis: Literal["tap", "permanent"] = "tap",
    duty: Literal["metering", "restriction"] = "metering",
    vapor_pressure_pa: Optional[float] = None,
    bubble_point_pa: Optional[float] = None,

    # Plate kit parameters
    diameter_variability_percent: float = 20.0,
    kit_sizes: int = 5,

    # Output formatting
    size_units: Literal["m", "inch", "mm"] = "mm",
) -> dict:
    """
    Unified orifice restrictor sizing tool.

    Provides fixed orifice (restriction plate) sizing for liquid and gas service.
    Uses ISO 5167 correlations with Reader-Harris-Gallagher discharge coefficients.

    Modes:
    - sizing: Q + ΔP → d (size orifice for target conditions)
    - analysis_dp: Q + d → ΔP (predict pressure drop for existing orifice)
    - analysis_flow: d + ΔP → Q (predict flow rate for existing orifice)
    - plate_kit: Generate spare plates showing ΔP sensitivity to diameter

    Args:
        mode: Calculation mode (sizing, analysis_dp, analysis_flow, plate_kit)
        phase: Fluid phase (liquid, gas)

        # Flow parameters
        flow_rate: Volumetric flow rate in m³/s
        flow_rate_gpm: Volumetric flow rate in US GPM
        flow_rate_m3_hr: Volumetric flow rate in m³/hr

        # Pressure parameters
        pressure_drop: Target or measured pressure drop in Pa
        pressure_drop_psi: Pressure drop in psi
        pressure_drop_bar: Pressure drop in bar
        inlet_pressure: Inlet pressure in Pa (absolute) - optional for liquid
        inlet_pressure_psi: Inlet pressure in psi (absolute)
        inlet_pressure_bar: Inlet pressure in bar (absolute)

        # Orifice geometry (for analysis modes)
        orifice_diameter: Orifice bore diameter in m
        orifice_diameter_in: Orifice bore diameter in inches
        orifice_diameter_mm: Orifice bore diameter in mm

        # Pipe geometry
        pipe_diameter: Pipe inner diameter in m
        pipe_diameter_in: Pipe inner diameter in inches
        nominal_size_in: Nominal pipe size in inches (for lookup)
        schedule: Pipe schedule (default "40")

        # Fluid properties
        fluid_density: Fluid density in kg/m³
        fluid_density_lbft3: Fluid density in lb/ft³
        fluid_viscosity: Fluid viscosity in Pa·s
        fluid_viscosity_cp: Fluid viscosity in centipoise
        fluid_name: Fluid name for property lookup (e.g., "Water", "Methane")
        temperature_c: Temperature in Celsius (required with fluid_name)
        pressure_bar: Pressure context for property lookup

        # Gas-specific properties
        gas_gamma: Specific heat ratio Cp/Cv (required for gas phase)
        gas_mw: Molecular weight in kg/kmol
        gas_z_factor: Compressibility factor

        # Orifice parameters
        cd: Override discharge coefficient (default: ISO 5167 correlation)
        meter_type: ISO 5167 meter type (default: ISO_5167_ORIFICE)
        taps: Tap configuration ('D', 'flange', 'corner', default: 'D')

        # Plate kit parameters
        diameter_variability_percent: ±% range for plate kit (default: 20)
        kit_sizes: Number of plates in kit (default: 5)

        # Output formatting
        size_units: Output units for diameter ('m', 'inch', 'mm', default: 'mm')

    Returns:
        Structured orifice restrictor result matching OrificeRestrictorOutput:
        - sizing: orifice diameter, beta ratio, flow parameters
        - analysis_dp: pressure drop, head loss
        - analysis_flow: flow rate (volumetric and mass)
        - plate_kit: list of plates with ΔP at design flow
    """
    try:
        request = OrificeRestrictorInput(**locals())
        if request.mode == "sizing":
            fn = calculate_orifice_sizing
        elif request.mode == "analysis_dp":
            fn = calculate_orifice_analysis_dp
        elif request.mode == "analysis_flow":
            fn = calculate_orifice_analysis_flow
        else:
            fn = generate_plate_kit

        params = request.model_dump(exclude_none=True)
        params.pop("mode", None)
        sig = inspect.signature(fn)
        allowed = set(sig.parameters.keys())
        forwarded = {key: value for key, value in params.items() if key in allowed}
        payload = _decode_tool_json(fn(**forwarded), request.mode, request.phase)
        return _map_output(payload, request)
    except CalculationError as exc:
        return _error_output(str(mode), str(phase), exc)
    except ValidationError as exc:
        return _error_output(
            str(mode),
            str(phase),
            CalculationError(
                "INVALID_ORIFICE_RESTRICTOR_INPUT",
                str(exc),
                details={"mode": mode, "phase": phase},
            ),
        )
    except Exception as exc:
        return _error_output(
            str(mode),
            str(phase),
            CalculationError(
                "ORIFICE_RESTRICTOR_CALCULATION_FAILED",
                str(exc),
                details={"mode": mode, "phase": phase},
            ),
        )
