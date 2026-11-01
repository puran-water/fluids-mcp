"""Unified machine requirements for pumps and compressors."""

from __future__ import annotations

import inspect
import json
import math
from typing import Any, Dict, Literal, Optional

from pydantic import ValidationError

from fluids_mcp.errors import CalculationError
from fluids_mcp.schemas.omnitools import (
    CalculationMessage,
    MachineRequirementsInput,
    MachineRequirementsOutput,
)
from fluids_mcp.validation.result import assert_finite_result
from fluids_mcp.tools.blower_compressor import calculate_blower_compressor_requirements
from fluids_mcp.tools.open_channel_flow_new import calculate_open_channel_flow
from fluids_mcp.tools.pump_requirements import calculate_pump_requirements
from fluids_mcp.tools.reynolds_number import calculate_reynolds_number
from fluids_mcp.omnitools.pipe_flow import pipe_flow
from fluids_mcp.utils.constants import G_GRAVITY


VALID_MACHINE_TYPES = {"pump", "compressor", "reynolds", "open_channel"}
DEFAULT_MACHINE_TYPE = "pump"


def _message_from_error(error: CalculationError) -> CalculationMessage:
    return CalculationMessage(
        code=error.code,
        message=error.message,
        severity=error.severity,  # type: ignore[arg-type]
        details=error.details,
        correlation_id=error.correlation_id,
    )


def _error_output(machine_type: str, error: CalculationError) -> dict:
    normalized = machine_type if machine_type in VALID_MACHINE_TYPES else DEFAULT_MACHINE_TYPE
    output = MachineRequirementsOutput(
        machine_type=normalized,  # type: ignore[arg-type]
        errors=[_message_from_error(error)],
    )
    payload = output.model_dump(exclude_none=True)
    assert_finite_result(payload)
    return payload


def _warning_messages(raw_warnings: Any) -> list[CalculationMessage]:
    if not raw_warnings:
        return []
    warnings = raw_warnings if isinstance(raw_warnings, list) else [raw_warnings]
    messages = []
    for warning in warnings:
        if isinstance(warning, dict):
            messages.append(
                CalculationMessage(
                    code=str(warning.get("code", "BACKEND_WARNING")),
                    message=str(warning.get("message", warning)),
                    severity="warning",
                    details=warning.get("details", {}) if isinstance(warning.get("details", {}), dict) else {},
                )
            )
        else:
            messages.append(CalculationMessage(code="BACKEND_WARNING", message=str(warning), severity="warning"))
    return messages


def _decode_tool_json(payload: str, machine_type: str) -> dict:
    try:
        decoded = json.loads(payload)
    except json.JSONDecodeError as exc:
        raise CalculationError(
            "INVALID_TOOL_RESPONSE",
            f"Machine-requirements backend returned invalid JSON: {exc}",
            details={"machine_type": machine_type},
        ) from exc
    if not isinstance(decoded, dict):
        raise CalculationError(
            "INVALID_TOOL_RESPONSE",
            "Machine-requirements backend returned a non-object response",
            details={"machine_type": machine_type, "response_type": type(decoded).__name__},
        )
    return decoded


def _payload_error(payload: dict, machine_type: str) -> Optional[CalculationError]:
    if payload.get("error_code") or payload.get("error"):
        details = {"machine_type": machine_type}
        errors = payload.get("errors")
        if isinstance(errors, list) and errors and isinstance(errors[0], dict):
            raw_details = errors[0].get("details", {})
            if isinstance(raw_details, dict):
                details.update(raw_details)
        return CalculationError(
            str(payload.get("error_code") or "MACHINE_CALCULATION_FAILED"),
            str(payload.get("error") or payload.get("message") or "Machine calculation failed"),
            details=details,
        )
    errors = payload.get("errors")
    if not errors:
        return None
    if isinstance(errors, list) and errors:
        first = errors[0]
        if isinstance(first, dict):
            return CalculationError(
                str(first.get("code") or "MACHINE_CALCULATION_FAILED"),
                str(first.get("message") or first),
                details=first.get("details", {}) if isinstance(first.get("details", {}), dict) else {},
            )
        return CalculationError(
            "MACHINE_CALCULATION_FAILED",
            "; ".join(str(item) for item in errors),
            details={"machine_type": machine_type},
        )
    return CalculationError(
        "MACHINE_CALCULATION_FAILED",
        str(errors),
        details={"machine_type": machine_type},
    )


def _backend_context(payload: dict) -> dict:
    return {
        key: value
        for key, value in payload.items()
        if key not in {"error", "error_code", "errors"} and value is not None
    }


def _map_output(payload: dict, request: MachineRequirementsInput) -> dict:
    error = _payload_error(payload, request.machine_type)
    if error is not None:
        output = MachineRequirementsOutput(
            machine_type=request.machine_type,
            errors=[_message_from_error(error)],
            warnings=_warning_messages(payload.get("warnings")),
            inputs_resolved=payload.get("inputs_resolved") or payload.get("log") or [],
            details=_backend_context(payload),
        )
    else:
        output = MachineRequirementsOutput(
            machine_type=request.machine_type,
            pump_mode=payload.get("pump_mode") if request.machine_type == "pump" else None,
            compressor_mode=(
                payload.get("compressor_mode", request.compressor_mode)
                if request.machine_type == "compressor"
                else None
            ),
            compressor_type=(
                payload.get("compressor_type", request.compressor_type)
                if request.machine_type == "compressor"
                else None
            ),
            open_channel_mode=(
                payload.get("open_channel_mode", request.open_channel_mode)
                if request.machine_type == "open_channel"
                else None
            ),
            flow_rate_m3s=payload.get("flow_rate_m3s"),
            flow_rate_gpm=payload.get("flow_rate_gpm"),
            pipe_diameter_m=payload.get("pipe_diameter_m"),
            pipe_velocity_m_s=payload.get("pipe_velocity_m_s"),
            static_head_m=payload.get("static_head_m"),
            static_head_ft=payload.get("static_head_ft"),
            target_tdh_m=payload.get("target_tdh_m"),
            target_tdh_ft=payload.get("target_tdh_ft"),
            total_dynamic_head_m=payload.get("total_dynamic_head_m"),
            total_dynamic_head_ft=payload.get("total_dynamic_head_ft"),
            npsha_m=payload.get("npsha_m"),
            npsha_ft=payload.get("npsha_ft"),
            npshr_m=payload.get("npshr_m"),
            npshr_ft=payload.get("npshr_ft"),
            npshr_status=payload.get("npshr_status"),
            npsh_margin_m=payload.get("npsh_margin_m"),
            npsh_margin_ft=payload.get("npsh_margin_ft"),
            required_npsh_margin_m=payload.get("required_npsh_margin_m"),
            npsh_margin_pass=payload.get("npsh_margin_pass"),
            estimated_power_kw=payload.get("estimated_power_kw"),
            estimated_power_hp=payload.get("estimated_power_hp"),
            required_power_kw=payload.get("required_power_kw"),
            shaft_power_kw=payload.get("shaft_power_kw"),
            required_discharge_pressure_gauge_pa=payload.get("required_discharge_pressure_gauge_pa"),
            required_discharge_pressure_absolute_pa=payload.get("required_discharge_pressure_absolute_pa"),
            pressure_component_breakdown=payload.get("pressure_component_breakdown"),
            blower_scenarios=payload.get("blower_scenarios"),
            actual_discharge_temp_c=payload.get("actual_discharge_temp_c"),
            reynolds_number=payload.get("reynolds_number"),
            normal_depth_m=payload.get("normal_depth_m"),
            critical_depth_m=payload.get("critical_depth_m"),
            specific_energy_m=payload.get("specific_energy_m"),
            critical_specific_energy_m=payload.get("critical_specific_energy_m"),
            froude_number=payload.get("froude_number"),
            flow_regime=payload.get("flow_regime"),
            flow_state=payload.get("flow_state"),
            hgl_elevation_m=payload.get("hgl_elevation_m"),
            hydraulic_grade_above_invert_m=payload.get("hydraulic_grade_above_invert_m"),
            energy_grade_above_invert_m=payload.get("energy_grade_above_invert_m"),
            energy_grade_elevation_m=payload.get("energy_grade_elevation_m"),
            available_fall_m=payload.get("available_fall_m"),
            required_fall_m=payload.get("required_fall_m"),
            residual_fall_m=payload.get("residual_fall_m"),
            gravity_capacity_pass=payload.get("gravity_capacity_pass"),
            warnings=_warning_messages(payload.get("warnings")),
            inputs_resolved=payload.get("inputs_resolved") or payload.get("log") or [],
            details=payload,
        )
    result = output.model_dump(exclude_none=True)
    assert_finite_result(result)
    return result


def _raise_from_public_result(result: dict, *, context: str) -> None:
    """Translate a nested omnitool failure without leaking its error envelope."""

    if result.get("ok", True):
        return
    errors = result.get("errors") or []
    first = errors[0] if errors and isinstance(errors[0], dict) else {}
    details = first.get("details", {}) if isinstance(first.get("details", {}), dict) else {}
    raise CalculationError(
        str(first.get("code") or "MACHINE_CALCULATION_FAILED"),
        str(first.get("message") or f"{context} calculation failed."),
        details={**details, "context": context},
    )


def _validate_positive_system_value(name: str, value: Optional[float]) -> float:
    if value is None or not isinstance(value, (int, float)) or not math.isfinite(float(value)) or value <= 0:
        raise CalculationError(
            "INVALID_AERATION_SYSTEM_INPUT",
            f"{name} must be positive and finite for aeration system-pressure mode.",
            details={name: value},
        )
    return float(value)


def _run_compressor_backend(params: dict[str, Any]) -> dict:
    sig = inspect.signature(calculate_blower_compressor_requirements)
    forwarded = {key: value for key, value in params.items() if key in sig.parameters}
    payload = _decode_tool_json(
        calculate_blower_compressor_requirements(**forwarded),
        "compressor",
    )
    error = _payload_error(payload, "compressor")
    if error is not None:
        raise error
    return payload


def _calculate_aeration_system_pressure(request: MachineRequirementsInput) -> dict:
    """Build clean/fouled aeration blower duties through the existing gas path.

    Required diffuser-inlet gauge pressure is the hydrostatic submergence term
    ``rho_water*g*h`` plus vendor clean/fouled diffuser pressure drop. The
    upstream blower pressure is solved with ``pipe_flow(phase='gas',
    method='isothermal_darcy')`` so Darcy-Weisbach/Colebrook header and fitting
    losses use gas density at actual absolute pressure. Each solved pressure is
    then passed to the compressor backend. Roots/lobe machines use its ISO 1217
    inlet-volume power basis; turbo machines retain the compression model.
    """

    if any(
        value is not None
        for value in (request.outlet_pressure, request.outlet_pressure_psi)
    ):
        raise CalculationError(
            "CONFLICTING_OUTLET_PRESSURE",
            "Do not provide outlet pressure in aeration_system_pressure mode; it is calculated from the system components.",
            details={"compressor_mode": request.compressor_mode},
        )

    missing_diffuser_fields = [
        name
        for name, value in (
            ("diffuser_dp_clean_mbar", request.diffuser_dp_clean_mbar),
            ("diffuser_dp_fouled_mbar", request.diffuser_dp_fouled_mbar),
        )
        if value is None
    ]
    if missing_diffuser_fields:
        raise CalculationError(
            "MISSING_DIFFUSER_PRESSURE_DROP",
            "Vendor clean and fouled/end-of-run diffuser pressure drops are required; neither condition is assumed to be zero.",
            details={"missing_fields": missing_diffuser_fields},
        )

    submergence_m = _validate_positive_system_value("submergence_m", request.submergence_m)
    liquid_density = _validate_positive_system_value(
        "submergence_liquid_density_kg_m3",
        request.submergence_liquid_density_kg_m3,
    )
    atmospheric_pressure = request.atmospheric_pressure
    if request.atmospheric_pressure_psi is not None:
        if atmospheric_pressure not in (None, 101325.0):
            raise CalculationError(
                "CONFLICTING_PRESSURE_INPUTS",
                "Provide atmospheric pressure in either Pa absolute or psi absolute, not both.",
                details={"field": "atmospheric_pressure"},
            )
        atmospheric_pressure = request.atmospheric_pressure_psi * 6894.76
    atmospheric_pressure = _validate_positive_system_value(
        "atmospheric_pressure",
        101325.0 if atmospheric_pressure is None else atmospheric_pressure,
    )
    if request.inlet_temperature_c is None:
        raise CalculationError(
            "MISSING_REQUIRED_INPUT",
            "inlet_temperature_c is required for aeration system-pressure mode.",
            details={"field": "inlet_temperature_c"},
        )
    if request.pipe_length is None and request.pipe_length_ft is None:
        raise CalculationError(
            "MISSING_REQUIRED_INPUT",
            "pipe_length or pipe_length_ft is required for the aeration header loss.",
            details={"field": "pipe_length"},
        )

    static_pressure_pa = liquid_density * G_GRAVITY * submergence_m
    request_params = request.model_dump(exclude_none=True)
    scenarios: dict[str, Any] = {}
    compressor_payloads: dict[str, dict[str, Any]] = {}
    combined_warnings: list[Any] = []

    for condition, diffuser_mbar in (
        ("clean", request.diffuser_dp_clean_mbar),
        ("fouled", request.diffuser_dp_fouled_mbar),
    ):
        diffuser_dp_pa = float(diffuser_mbar) * 100.0
        if not math.isfinite(diffuser_dp_pa) or diffuser_dp_pa < 0:
            raise CalculationError(
                "INVALID_AERATION_SYSTEM_INPUT",
                f"diffuser_dp_{condition}_mbar must be non-negative and finite.",
                details={f"diffuser_dp_{condition}_mbar": diffuser_mbar},
            )
        diffuser_inlet_absolute_pa = atmospheric_pressure + static_pressure_pa + diffuser_dp_pa
        header_result = pipe_flow(
            phase="gas",
            inlet_pressure=None,
            outlet_pressure=diffuser_inlet_absolute_pa,
            pipe_length=request.pipe_length,
            pipe_length_ft=request.pipe_length_ft,
            pipe_diameter=request.pipe_diameter,
            pipe_diameter_in=request.pipe_diameter_in,
            nominal_size_in=request.nominal_size_in,
            schedule=request.schedule,
            material=request.material,
            pipe_roughness=request.pipe_roughness,
            fittings=request.fittings,
            temperature_c=request.inlet_temperature_c,
            fluid_name=request.fluid_name,
            flow_rate_kg_s=request.flow_rate_kg_s,
            flow_rate_norm_m3_hr=request.flow_rate_norm_m3_hr,
            flow_rate_std_m3_hr=request.flow_rate_std_m3_hr,
            gas_mw=request.gas_mw,
            gas_gamma=request.gas_gamma,
            gas_z_factor=request.gas_z_factor,
            gas_viscosity=request.gas_viscosity,
            gas_composition_mol=request.gas_composition_mol,
            method="isothermal_darcy",
            allow_property_defaults=request.allow_property_defaults,
        )
        _raise_from_public_result(header_result, context=f"{condition} aeration header")
        required_absolute_pa = float(header_result["inlet_pressure_pa"])
        pipe_and_fitting_dp_pa = required_absolute_pa - diffuser_inlet_absolute_pa
        if pipe_and_fitting_dp_pa < -1e-6:
            raise CalculationError(
                "INVALID_HEADER_PRESSURE_RESULT",
                "Gas header calculation returned a negative pressure loss.",
                details={"condition": condition, "pressure_drop_pa": pipe_and_fitting_dp_pa},
            )
        required_gauge_pa = required_absolute_pa - atmospheric_pressure

        compressor_params = dict(request_params)
        for field in (
            "inlet_pressure",
            "inlet_pressure_psi",
            "outlet_pressure",
            "outlet_pressure_psi",
        ):
            compressor_params.pop(field, None)
        compressor_params.update(
            {
                "inlet_pressure": atmospheric_pressure,
                "outlet_pressure": required_absolute_pa,
                "atmospheric_pressure": atmospheric_pressure,
                "compressor_type": request.compressor_type,
            }
        )
        compressor_payload = _run_compressor_backend(compressor_params)
        compressor_payloads[condition] = compressor_payload
        combined_warnings.extend(header_result.get("warnings") or [])
        combined_warnings.extend(compressor_payload.get("warnings") or [])

        scenarios[condition] = {
            "condition": condition,
            "static_submergence_pressure_pa": static_pressure_pa,
            "static_submergence_pressure_mbar": static_pressure_pa / 100.0,
            "diffuser_pressure_drop_pa": diffuser_dp_pa,
            "diffuser_pressure_drop_mbar": float(diffuser_mbar),
            "pipe_and_fitting_pressure_drop_pa": pipe_and_fitting_dp_pa,
            "pipe_and_fitting_pressure_drop_mbar": pipe_and_fitting_dp_pa / 100.0,
            "required_discharge_pressure_gauge_pa": required_gauge_pa,
            "required_discharge_pressure_gauge_bar": required_gauge_pa / 100000.0,
            "required_discharge_pressure_absolute_pa": required_absolute_pa,
            "required_discharge_pressure_absolute_bar": required_absolute_pa / 100000.0,
            "shaft_power_kw": compressor_payload.get("shaft_power_kw", compressor_payload.get("required_power_kw")),
            "header_velocity_m_s": header_result.get("average_velocity_m_s"),
            "header_mach_number": header_result.get("average_mach_number"),
            "header_calculation_method": "isothermal_darcy",
            "compressor_power_method": compressor_payload.get("power_method"),
        }

    design = scenarios["fouled"]
    payload = dict(compressor_payloads["fouled"])
    payload.update(
        {
            "compressor_mode": "aeration_system_pressure",
            "compressor_type": request.compressor_type,
            "design_condition": "fouled",
            "submergence_m": submergence_m,
            "submergence_liquid_density_kg_m3": liquid_density,
            "required_discharge_pressure_gauge_pa": design["required_discharge_pressure_gauge_pa"],
            "required_discharge_pressure_absolute_pa": design["required_discharge_pressure_absolute_pa"],
            "shaft_power_kw": design["shaft_power_kw"],
            "pressure_component_breakdown": {
                condition: {
                    key: value
                    for key, value in scenario.items()
                    if key in {
                        "static_submergence_pressure_pa",
                        "diffuser_pressure_drop_pa",
                        "pipe_and_fitting_pressure_drop_pa",
                        "required_discharge_pressure_gauge_pa",
                        "required_discharge_pressure_absolute_pa",
                    }
                }
                for condition, scenario in scenarios.items()
            },
            "blower_scenarios": scenarios,
            "warnings": combined_warnings,
            "inputs_resolved": [
                "Aeration discharge pressure computed as hydrostatic submergence + vendor diffuser dP + compressible header dP.",
                "Header loss solved at actual absolute pressure with gas pipe_flow method=isothermal_darcy.",
                *(payload.get("inputs_resolved") or payload.get("log") or []),
            ],
        }
    )
    return payload


def _reject_direct_zero_pump_flow(request: MachineRequirementsInput) -> None:
    if request.machine_type != "pump" or request.pump_mode != "flow_rate_given":
        return
    zero_fields = [
        name
        for name in ("flow_rate", "flow_rate_gpm", "flow_rate_m3s")
        if getattr(request, name) == 0
    ]
    if zero_fields:
        raise CalculationError(
            "ZERO_FLOW",
            "Zero pump flow is not supported for flow_rate_given mode; use target-TDH modes for "
            "solver bracketing or provide a positive flow rate.",
            details={"fields": zero_fields, "pump_mode": request.pump_mode},
        )


def machine_requirements(
    machine_type: Literal["pump", "compressor", "reynolds", "open_channel"] = "pump",

    # Common flow parameters
    flow_rate: Optional[float] = None,
    flow_rate_gpm: Optional[float] = None,
    flow_rate_kg_s: Optional[float] = None,
    flow_rate_norm_m3_hr: Optional[float] = None,
    flow_rate_std_m3_hr: Optional[float] = None,
    flow_rate_m3s: Optional[float] = None,

    # Pump-specific parameters
    pump_mode: Literal["flow_rate_given", "target_tdh_solve_flow", "target_tdh_solve_diameter"] = "flow_rate_given",
    endpoint_basis: Literal["tank_surface", "nozzle"] = "tank_surface",
    compressor_mode: Literal["specified_outlet_pressure", "aeration_system_pressure"] = "specified_outlet_pressure",
    compressor_type: Literal["turbo", "roots_pd"] = "turbo",
    open_channel_mode: Literal["normal_flow", "gravity_leg"] = "normal_flow",
    pipe_diameter: Optional[float] = None,
    pipe_diameter_in: Optional[float] = None,
    suction_pipe_diameter: Optional[float] = None,
    suction_pipe_diameter_in: Optional[float] = None,
    discharge_pipe_diameter: Optional[float] = None,
    discharge_pipe_diameter_in: Optional[float] = None,
    pipe_length: Optional[float] = None,
    pipe_length_ft: Optional[float] = None,
    suction_pipe_length: Optional[float] = None,
    suction_pipe_length_ft: Optional[float] = None,
    discharge_pipe_length: Optional[float] = None,
    discharge_pipe_length_ft: Optional[float] = None,
    static_suction_head: Optional[float] = None,
    static_suction_head_ft: Optional[float] = None,
    static_discharge_head: Optional[float] = None,
    static_discharge_head_ft: Optional[float] = None,
    suction_fittings: Optional[list] = None,
    discharge_fittings: Optional[list] = None,
    fittings: Optional[list] = None,
    suction_nozzle_diameter: Optional[float] = None,
    discharge_nozzle_diameter: Optional[float] = None,
    suction_nozzle_diameter_in: Optional[float] = None,
    discharge_nozzle_diameter_in: Optional[float] = None,
    atmospheric_pressure: Optional[float] = None,
    atmospheric_pressure_psi: Optional[float] = None,
    suction_tank_pressure: Optional[float] = None,
    suction_tank_pressure_psi: Optional[float] = None,
    discharge_tank_pressure: Optional[float] = None,
    discharge_tank_pressure_psi: Optional[float] = None,
    suction_nozzle_pressure: Optional[float] = None,
    suction_nozzle_pressure_psi: Optional[float] = None,
    discharge_nozzle_pressure: Optional[float] = None,
    discharge_nozzle_pressure_psi: Optional[float] = None,
    target_tdh: Optional[float] = None,
    target_tdh_ft: Optional[float] = None,
    npshr: Optional[float] = None,
    npshr_ft: Optional[float] = None,
    required_npsh_margin: Optional[float] = None,
    required_npsh_margin_ft: Optional[float] = None,

    # Compressor-specific parameters
    inlet_pressure: Optional[float] = None,
    outlet_pressure: Optional[float] = None,
    inlet_pressure_psi: Optional[float] = None,
    outlet_pressure_psi: Optional[float] = None,
    inlet_temperature_c: Optional[float] = None,
    efficiency: Optional[float] = None,
    efficiency_type: str = "isentropic",
    polytropic_n: Optional[float] = None,
    submergence_m: Optional[float] = None,
    submergence_liquid_density_kg_m3: float = 998.2,
    diffuser_dp_clean_mbar: Optional[float] = None,
    diffuser_dp_fouled_mbar: Optional[float] = None,

    # Fluid properties
    fluid_density: Optional[float] = None,
    fluid_density_lbft3: Optional[float] = None,
    fluid_viscosity: Optional[float] = None,
    fluid_viscosity_cp: Optional[float] = None,
    fluid_vapor_pressure: Optional[float] = None,
    fluid_vapor_pressure_psi: Optional[float] = None,
    fluid_name: Optional[str] = None,
    temperature_c: Optional[float] = None,
    pressure_bar: Optional[float] = None,
    pipe_roughness: Optional[float] = None,
    suction_pipe_roughness: Optional[float] = None,
    discharge_pipe_roughness: Optional[float] = None,
    nominal_size_in: Optional[float] = None,
    schedule: str = "40",
    material: Optional[str] = None,
    suction_material: Optional[str] = None,
    discharge_material: Optional[str] = None,

    # Gas properties
    gas_mw: Optional[float] = None,
    gas_gamma: Optional[float] = None,
    gas_z_factor: Optional[float] = None,
    gas_viscosity: Optional[float] = None,
    gas_composition_mol: Optional[dict] = None,

    # Compressor calculation flags
    allow_property_defaults: bool = False,
    use_eos_calculations: Optional[bool] = None,
    solve_for: Optional[str] = None,

    # Reynolds number parameters
    velocity: Optional[float] = None,
    velocity_ft_s: Optional[float] = None,
    characteristic_length: Optional[float] = None,
    characteristic_length_in: Optional[float] = None,
    fluid_kinematic_viscosity: Optional[float] = None,

    # Open channel parameters
    channel_type: Optional[str] = None,
    width_m: Optional[float] = None,
    diameter_m: Optional[float] = None,
    side_slope_z: Optional[float] = None,
    slope: Optional[float] = None,
    manning_n: Optional[float] = None,
    normal_depth_m: Optional[float] = None,
    circular_depth_branch: Literal["shallow", "near_full"] = "shallow",
    invert_elevation_m: Optional[float] = None,
    upstream_head_m: Optional[float] = None,
    downstream_head_m: Optional[float] = None,
) -> dict:
    """Unified machine requirements and hydraulic calculations.

    For pumps, ``endpoint_basis='tank_surface'`` (the default) applies Bernoulli
    between reservoir free surfaces: endpoint velocities are zero, line losses
    enter TDH, and suction-flange velocity head cancels from tank-to-flange
    NPSHa. ``endpoint_basis='nozzle'`` instead uses absolute pump-flange static
    pressures and the flange velocity-head terms, consistent with ISO 9906 / HI
    pump head conventions. Side-specific pipe diameter, material, and roughness
    override the common values; omitted side values inherit the common inputs.
    Omitted atmospheric pressure is resolved by the pump backend to 101325 Pa
    absolute and disclosed in a result warning.

    ``compressor_mode='aeration_system_pressure'`` is a mode, rather than a
    new machine type, because it derives the discharge duty for the existing
    compressor/blower calculation. Likewise, ``open_channel_mode='gravity_leg'``
    evaluates the full-bore continuation of the existing gravity-flow path.
    Public omnitool names and discriminators therefore remain stable.
    """
    try:
        request = MachineRequirementsInput(**locals())
        _reject_direct_zero_pump_flow(request)
        if request.machine_type == "pump":
            fn = calculate_pump_requirements
        elif request.machine_type == "compressor":
            if request.compressor_mode == "aeration_system_pressure":
                return _map_output(_calculate_aeration_system_pressure(request), request)
            fn = calculate_blower_compressor_requirements
        elif request.machine_type == "reynolds":
            fn = calculate_reynolds_number
        else:
            fn = calculate_open_channel_flow

        params = request.model_dump(exclude_none=True)
        params.pop("machine_type", None)
        sig = inspect.signature(fn)
        allowed = set(sig.parameters.keys())
        forwarded = {key: value for key, value in params.items() if key in allowed}
        payload = _decode_tool_json(fn(**forwarded), request.machine_type)
        return _map_output(payload, request)
    except CalculationError as exc:
        return _error_output(str(machine_type), exc)
    except ValidationError as exc:
        return _error_output(
            str(machine_type),
            CalculationError(
                "INVALID_MACHINE_REQUIREMENTS_INPUT",
                str(exc),
                details={"machine_type": machine_type},
            ),
        )
    except Exception as exc:
        return _error_output(
            str(machine_type),
            CalculationError(
                "MACHINE_CALCULATION_FAILED",
                str(exc),
                details={"machine_type": machine_type},
            ),
        )
