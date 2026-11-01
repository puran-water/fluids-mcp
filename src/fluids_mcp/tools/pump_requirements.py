"""Pump TDH, NPSHA, and target-TDH sizing calculations."""

from __future__ import annotations

import json
import logging
import math
from typing import Any, Dict, List, Literal, Optional, Union

import fluids
import fluids.core
import fluids.friction
import fluids.piping
from scipy.optimize import brentq

from fluids_mcp.errors import CalculationError
from fluids_mcp.utils.constants import (
    CENTIPOISE_to_PAS,
    FT_to_M,
    GPM_to_M3S,
    G_GRAVITY,
    INCH_to_M,
    LBFT3_to_KGM3,
    PSI_to_PA,
)
from fluids_mcp.utils.helpers import get_pipe_roughness, resolve_fitting_k_basis

logger = logging.getLogger("fluids-mcp.calculate_pump_requirements")

PumpMode = Literal["flow_rate_given", "target_tdh_solve_flow", "target_tdh_solve_diameter"]
EndpointBasis = Literal["tank_surface", "nozzle"]
DEFAULT_ATMOSPHERIC_PRESSURE_PA = 101325.0
MAX_PIPE_VELOCITY_M_S = 10.0
MAX_DIAMETER_SOLVE_VELOCITY_M_S = 30.0
MIN_DIAMETER_SOLVE_VELOCITY_M_S = 0.05
TARGET_TOL_M = 1e-7
DEFAULT_REQUIRED_NPSH_MARGIN_M = 1.0


def _json_error(
    code: str,
    message: str,
    results_log: list[str],
    *,
    details: Optional[dict[str, Any]] = None,
    warnings: Optional[list[str]] = None,
) -> str:
    error = CalculationError(code, message, details=details)
    return json.dumps(
        {
            "error_code": error.code,
            "error": error.message,
            "errors": [error.to_public_dict()],
            "warnings": warnings or [],
            "log": results_log,
        },
        allow_nan=False,
    )


def _json_calculation_error(
    error: CalculationError,
    results_log: list[str],
    *,
    warnings: Optional[list[str]] = None,
) -> str:
    return json.dumps(
        {
            "error_code": error.code,
            "error": error.message,
            "errors": [error.to_public_dict()],
            "warnings": warnings or [],
            "log": results_log,
        },
        allow_nan=False,
    )


def _resolve_length(
    si_value: Optional[float],
    ft_value: Optional[float],
    name: str,
    results_log: list[str],
    errors: list[str],
) -> Optional[float]:
    if si_value is not None:
        results_log.append(f"Used provided SI {name}.")
        return si_value
    if ft_value is not None:
        results_log.append(f"Converted {name} from {ft_value} ft.")
        return ft_value * FT_to_M
    errors.append(f"Missing required input: {name} or {name}_ft.")
    return None


def _resolve_flow(
    flow_rate: Optional[float],
    flow_rate_gpm: Optional[float],
    pump_mode: PumpMode,
    results_log: list[str],
    errors: list[str],
) -> Optional[float]:
    if flow_rate is not None:
        results_log.append("Used provided SI flow_rate.")
        return flow_rate
    if flow_rate_gpm is not None:
        results_log.append(f"Converted flow_rate from {flow_rate_gpm} GPM.")
        return flow_rate_gpm * GPM_to_M3S
    if pump_mode != "target_tdh_solve_flow":
        errors.append("Missing required input: flow_rate or flow_rate_gpm.")
    return None


def _lookup_roughness(
    material: Optional[str], results_log: list[str], warnings: list[str]
) -> tuple[float, str]:
    if material is None:
        results_log.append("Used default pipe_roughness 1.5e-5 m.")
        return 1.5e-5, "Default 1.5e-5 m (no material supplied)"
    roughness, source = get_pipe_roughness(material=material, warnings=warnings)
    results_log.append(f"Looked up roughness for material '{material}': {roughness} m ({source}).")
    return roughness, source


def _resolve_pipe(
    *,
    pipe_diameter: Optional[float],
    pipe_diameter_in: Optional[float],
    nominal_size_in: Optional[float],
    schedule: str,
    material: Optional[str],
    pipe_roughness: Optional[float],
    pump_mode: PumpMode,
    diameter_required: bool,
    results_log: list[str],
    errors: list[str],
    warnings: list[str],
) -> tuple[Optional[float], Optional[float], dict[str, Any]]:
    pipe_info: dict[str, Any] = {}
    roughness = pipe_roughness
    roughness_source = "Provided" if pipe_roughness is not None else None
    diameter = None

    if pipe_diameter is not None:
        diameter = pipe_diameter
        results_log.append("Used provided SI pipe_diameter.")
    elif pipe_diameter_in is not None:
        diameter = pipe_diameter_in * INCH_to_M
        results_log.append(f"Converted pipe_diameter from {pipe_diameter_in} inches.")
    elif nominal_size_in is not None:
        try:
            nps, di, do, thickness = fluids.piping.nearest_pipe(NPS=nominal_size_in, schedule=schedule)
            diameter = di
            pipe_info = {"NPS_in": nps, "schedule": schedule, "Do_m": do, "t_m": thickness, "Di_m": di}
            results_log.append(f"Looked up pipe dimensions for NPS {nominal_size_in} Sch {schedule}.")
        except Exception as exc:
            errors.append(f"Failed to look up pipe NPS {nominal_size_in} Sch {schedule}: {exc}")
    elif pump_mode != "target_tdh_solve_diameter" and diameter_required:
        errors.append("Missing required input: pipe_diameter, pipe_diameter_in, or nominal_size_in.")

    if roughness is None:
        roughness, roughness_source = _lookup_roughness(material, results_log, warnings)
    else:
        results_log.append("Used provided SI pipe_roughness.")

    if diameter is not None:
        pipe_info.setdefault("Di_m", diameter)
    if roughness is not None:
        pipe_info["roughness_m"] = roughness
        pipe_info["roughness_source"] = roughness_source
    return diameter, roughness, pipe_info


def _resolve_side_pipe(
    *,
    side: Literal["suction", "discharge"],
    pipe_diameter: Optional[float],
    pipe_diameter_in: Optional[float],
    material: Optional[str],
    pipe_roughness: Optional[float],
    fallback_diameter: Optional[float],
    fallback_roughness: Optional[float],
    fallback_roughness_source: Optional[str],
    diameter_required: bool,
    results_log: list[str],
    errors: list[str],
    warnings: list[str],
) -> tuple[Optional[float], Optional[float], dict[str, Any]]:
    """Resolve one side of the pump piping, falling back to legacy common inputs."""

    details: dict[str, Any] = {}
    if pipe_diameter is not None:
        diameter = pipe_diameter
        diameter_source = f"Provided {side}_pipe_diameter"
        results_log.append(f"Used provided SI {side}_pipe_diameter.")
    elif pipe_diameter_in is not None:
        diameter = pipe_diameter_in * INCH_to_M
        diameter_source = f"Converted {side}_pipe_diameter_in"
        results_log.append(f"Converted {side}_pipe_diameter from {pipe_diameter_in} inches.")
    else:
        diameter = fallback_diameter
        diameter_source = "Legacy common pipe diameter"

    if diameter is None and diameter_required:
        errors.append(
            f"Missing required {side} pipe diameter: provide {side}_pipe_diameter, "
            f"{side}_pipe_diameter_in, or a common pipe diameter."
        )

    if pipe_roughness is not None:
        roughness = pipe_roughness
        roughness_source = f"Provided {side}_pipe_roughness"
        results_log.append(f"Used provided SI {side}_pipe_roughness.")
    elif material is not None:
        roughness, roughness_source = _lookup_roughness(material, results_log, warnings)
    else:
        roughness = fallback_roughness
        roughness_source = fallback_roughness_source or "Legacy common pipe roughness"

    if roughness is None:
        errors.append(
            f"Missing required {side} pipe roughness: provide {side}_pipe_roughness, "
            f"{side}_material, or a common material/roughness."
        )

    if diameter is not None:
        details["diameter_m"] = diameter
        details["diameter_source"] = diameter_source
    if roughness is not None:
        details["roughness_m"] = roughness
        details["roughness_source"] = roughness_source
    if material is not None:
        details["material"] = material
    return diameter, roughness, details


def _resolve_fluid(
    *,
    fluid_density: Optional[float],
    fluid_density_lbft3: Optional[float],
    fluid_viscosity: Optional[float],
    fluid_viscosity_cp: Optional[float],
    fluid_vapor_pressure: Optional[float],
    fluid_vapor_pressure_psi: Optional[float],
    fluid_name: Optional[str],
    temperature_c: Optional[float],
    pressure_bar: Optional[float],
    results_log: list[str],
    errors: list[str],
) -> tuple[Optional[float], Optional[float], Optional[float], dict[str, Any]]:
    fluid_info: dict[str, Any] = {}
    density = None
    viscosity = None
    vapor_pressure = None

    if fluid_density is not None and fluid_viscosity is not None:
        density = fluid_density
        viscosity = fluid_viscosity
        results_log.append("Used provided SI fluid density and viscosity.")
    elif fluid_density_lbft3 is not None and fluid_viscosity_cp is not None:
        density = fluid_density_lbft3 * LBFT3_to_KGM3
        viscosity = fluid_viscosity_cp * CENTIPOISE_to_PAS
        results_log.append(
            f"Converted density from {fluid_density_lbft3} lb/ft3 and viscosity from {fluid_viscosity_cp} cP."
        )
    elif fluid_name is not None and temperature_c is not None:
        try:
            from fluids_mcp.utils.resolve_properties import resolve_liquid_properties

            p_bar = pressure_bar if pressure_bar is not None else 1.01325
            props = resolve_liquid_properties(fluid_name, temperature_c, p_bar)
            if props is None:
                errors.append(
                    f"Failed to look up fluid properties for {fluid_name} at {temperature_c} C: all backends failed"
                )
            else:
                density = props.density
                viscosity = props.viscosity
                vapor_pressure = props.vapor_pressure
                fluid_info = {
                    "name_used": fluid_name,
                    "temperature_c": temperature_c,
                    "pressure_bar": p_bar,
                    "density_kgm3": density,
                    "viscosity_pas": viscosity,
                    "vapor_pressure_pa": vapor_pressure,
                    "vapor_pressure_source": props.source if vapor_pressure is not None else None,
                }
                results_log.append(f"Properties via {props.source} for {fluid_name} at {temperature_c} C, {p_bar} bar.")
                results_log.extend(props.log)
        except Exception as exc:
            errors.append(f"Failed to look up base fluid properties for {fluid_name} at {temperature_c} C: {exc}")
    else:
        errors.append(
            "Missing inputs for fluid density/viscosity: Need (SI pair) OR "
            "(Imperial pair) OR (fluid_name AND temperature_c)."
        )

    if vapor_pressure is None:
        if fluid_vapor_pressure is not None:
            vapor_pressure = fluid_vapor_pressure
            fluid_info["vapor_pressure_source"] = "Provided SI"
            results_log.append("Used provided SI fluid_vapor_pressure.")
        elif fluid_vapor_pressure_psi is not None:
            vapor_pressure = fluid_vapor_pressure_psi * PSI_to_PA
            fluid_info["vapor_pressure_source"] = "Provided Imperial"
            results_log.append(f"Converted fluid_vapor_pressure from {fluid_vapor_pressure_psi} psi.")
        elif density is not None:
            errors.append(
                "Missing required input for NPSHa: fluid_vapor_pressure or "
                "fluid_vapor_pressure_psi (could not be looked up)."
            )

    return density, viscosity, vapor_pressure, {k: v for k, v in fluid_info.items() if v is not None}


def _resolve_atmospheric_pressure(
    *,
    atmospheric_pressure: Optional[float],
    atmospheric_pressure_psi: Optional[float],
    results_log: list[str],
) -> tuple[float, str, bool]:
    if atmospheric_pressure_psi is not None:
        results_log.append(f"Used provided atmospheric_pressure_psi ({atmospheric_pressure_psi} psi).")
        return atmospheric_pressure_psi * PSI_to_PA, "Provided atmospheric_pressure_psi", False
    if atmospheric_pressure is not None:
        results_log.append(f"Used provided atmospheric_pressure ({atmospheric_pressure} Pa absolute).")
        return atmospheric_pressure, "Provided atmospheric_pressure", False
    results_log.append(f"Defaulted atmospheric_pressure to {DEFAULT_ATMOSPHERIC_PRESSURE_PA} Pa absolute.")
    return DEFAULT_ATMOSPHERIC_PRESSURE_PA, "Default atmospheric pressure", True


def _resolve_endpoint_pressure(
    *,
    si_value: Optional[float],
    psi_value: Optional[float],
    name: str,
    fallback: Optional[float],
    fallback_source: Optional[str],
    results_log: list[str],
) -> tuple[Optional[float], Optional[str]]:
    if si_value is not None:
        results_log.append(f"Used provided {name} ({si_value} Pa absolute).")
        return si_value, f"Provided {name}"
    if psi_value is not None:
        results_log.append(f"Converted {name}_psi from {psi_value} psi absolute.")
        return psi_value * PSI_to_PA, f"Provided {name}_psi"
    if fallback is not None:
        results_log.append(f"Defaulted {name} to atmospheric pressure ({fallback} Pa absolute).")
        return fallback, fallback_source
    return None, None


def _resolve_target_tdh(
    target_tdh: Optional[float],
    target_tdh_ft: Optional[float],
    results_log: list[str],
) -> Optional[float]:
    if target_tdh is not None:
        results_log.append("Used provided SI target_tdh.")
        return target_tdh
    if target_tdh_ft is not None:
        results_log.append(f"Converted target_tdh from {target_tdh_ft} ft.")
        return target_tdh_ft * FT_to_M
    return None


def _infer_pump_mode(
    pump_mode: PumpMode,
    *,
    flow_rate: Optional[float],
    flow_rate_gpm: Optional[float],
    pipe_diameter: Optional[float],
    pipe_diameter_in: Optional[float],
    suction_pipe_diameter: Optional[float],
    suction_pipe_diameter_in: Optional[float],
    discharge_pipe_diameter: Optional[float],
    discharge_pipe_diameter_in: Optional[float],
    nominal_size_in: Optional[float],
    target_tdh_m: Optional[float],
) -> PumpMode:
    if pump_mode != "flow_rate_given" or target_tdh_m is None:
        return pump_mode
    has_flow = flow_rate is not None or flow_rate_gpm is not None
    has_diameter = any(
        value is not None
        for value in (
            pipe_diameter,
            pipe_diameter_in,
            nominal_size_in,
            suction_pipe_diameter,
            suction_pipe_diameter_in,
            discharge_pipe_diameter,
            discharge_pipe_diameter_in,
        )
    )
    if not has_flow and has_diameter:
        return "target_tdh_solve_flow"
    if has_flow and not has_diameter:
        return "target_tdh_solve_diameter"
    return pump_mode


def _resolve_nozzle_inputs(
    *,
    suction_nozzle_diameter: Optional[float],
    suction_nozzle_diameter_in: Optional[float],
    discharge_nozzle_diameter: Optional[float],
    discharge_nozzle_diameter_in: Optional[float],
    results_log: list[str],
) -> tuple[Optional[float], Optional[float]]:
    suction = None
    discharge = None
    if suction_nozzle_diameter is not None:
        suction = suction_nozzle_diameter
        results_log.append("Used provided SI suction_nozzle_diameter.")
    elif suction_nozzle_diameter_in is not None:
        suction = suction_nozzle_diameter_in * INCH_to_M
        results_log.append(f"Converted suction_nozzle_diameter from {suction_nozzle_diameter_in} inches.")

    if discharge_nozzle_diameter is not None:
        discharge = discharge_nozzle_diameter
        results_log.append("Used provided SI discharge_nozzle_diameter.")
    elif discharge_nozzle_diameter_in is not None:
        discharge = discharge_nozzle_diameter_in * INCH_to_M
        results_log.append(f"Converted discharge_nozzle_diameter from {discharge_nozzle_diameter_in} inches.")
    return suction, discharge


def _fitting_k_total(
    fittings: list[dict[str, Any]],
    *,
    pipe_diameter: float,
    reynolds_number: float,
    flow_rate: float,
    darcy_friction_factor: float,
    zero_flow: bool,
) -> tuple[float, list[dict[str, Any]]]:
    """Return total minor-loss K using Crane or equivalent-length input.

    Named and explicit K values follow Crane TP-410 conventions. Equivalent
    lengths are converted with the Darcy-Weisbach identity
    ``K = f_D L_e/D`` at the calculated side-specific Darcy friction factor.
    """

    total = 0.0
    details = []
    for fitting in fittings:
        quantity = float(fitting.get("quantity", 1))
        if not math.isfinite(quantity) or quantity < 0:
            raise ValueError("Fitting quantity must be non-negative and finite.")
        k_value, fitting_detail = resolve_fitting_k_basis(
            fitting,
            diameter=pipe_diameter,
            reynolds_number=reynolds_number,
            flow_rate=flow_rate,
            darcy_friction_factor=darcy_friction_factor,
            zero_flow=zero_flow,
        )
        total += k_value * quantity
        fitting_detail.update(
            {
                "quantity": quantity,
                "K_individual": round(k_value, 6),
                "K_total": round(k_value * quantity, 6),
            }
        )
        details.append(fitting_detail)
    return total, details


def _pump_hydraulics(
    *,
    flow_rate: float,
    suction_pipe_diameter: float,
    discharge_pipe_diameter: float,
    suction_pipe_length: float,
    discharge_pipe_length: float,
    static_suction_head: float,
    static_discharge_head: float,
    fluid_density: float,
    fluid_viscosity: float,
    fluid_vapor_pressure: float,
    suction_pipe_roughness: float,
    discharge_pipe_roughness: float,
    endpoint_basis: EndpointBasis,
    atmospheric_pressure: float,
    atmospheric_pressure_source: str,
    suction_surface_pressure: float,
    suction_surface_pressure_source: str,
    discharge_surface_pressure: float,
    discharge_surface_pressure_source: str,
    suction_nozzle_pressure: Optional[float],
    suction_nozzle_pressure_source: Optional[str],
    discharge_nozzle_pressure: Optional[float],
    discharge_nozzle_pressure_source: Optional[str],
    suction_fittings: list[dict[str, Union[str, int, float]]],
    discharge_fittings: list[dict[str, Union[str, int, float]]],
    suction_nozzle_diameter: Optional[float],
    discharge_nozzle_diameter: Optional[float],
) -> dict[str, Any]:
    if suction_pipe_diameter <= 0 or discharge_pipe_diameter <= 0:
        raise ValueError("Suction and discharge pipe diameters must be positive.")
    if flow_rate < 0:
        raise ValueError("Flow rate must be non-negative.")
    suction_nozzle = suction_nozzle_diameter or suction_pipe_diameter
    discharge_nozzle = discharge_nozzle_diameter or discharge_pipe_diameter
    if suction_nozzle <= 0 or discharge_nozzle <= 0:
        raise ValueError("Nozzle diameters must be positive.")

    suction_pipe_area = math.pi * (suction_pipe_diameter / 2.0) ** 2
    discharge_pipe_area = math.pi * (discharge_pipe_diameter / 2.0) ** 2
    suction_nozzle_area = math.pi * (suction_nozzle / 2.0) ** 2
    discharge_nozzle_area = math.pi * (discharge_nozzle / 2.0) ** 2
    zero_flow = flow_rate == 0

    suction_pipe_velocity = flow_rate / suction_pipe_area
    discharge_pipe_velocity = flow_rate / discharge_pipe_area
    suction_nozzle_velocity = flow_rate / suction_nozzle_area
    discharge_nozzle_velocity = flow_rate / discharge_nozzle_area

    warnings: list[str] = []
    if zero_flow:
        suction_reynolds_number = 0.0
        discharge_reynolds_number = 0.0
        suction_friction_factor = 0.0
        discharge_friction_factor = 0.0
    else:
        suction_reynolds_number = fluids.core.Reynolds(
            V=suction_pipe_velocity,
            D=suction_pipe_diameter,
            rho=fluid_density,
            mu=fluid_viscosity,
        )
        discharge_reynolds_number = fluids.core.Reynolds(
            V=discharge_pipe_velocity,
            D=discharge_pipe_diameter,
            rho=fluid_density,
            mu=fluid_viscosity,
        )
        suction_friction_factor = fluids.friction.friction_factor(
            Re=suction_reynolds_number,
            eD=suction_pipe_roughness / suction_pipe_diameter,
        )
        discharge_friction_factor = fluids.friction.friction_factor(
            Re=discharge_reynolds_number,
            eD=discharge_pipe_roughness / discharge_pipe_diameter,
        )

    suction_fitting_k, suction_fitting_details = _fitting_k_total(
        suction_fittings,
        pipe_diameter=suction_pipe_diameter,
        reynolds_number=suction_reynolds_number,
        flow_rate=flow_rate,
        darcy_friction_factor=suction_friction_factor,
        zero_flow=zero_flow,
    )
    discharge_fitting_k, discharge_fitting_details = _fitting_k_total(
        discharge_fittings,
        pipe_diameter=discharge_pipe_diameter,
        reynolds_number=discharge_reynolds_number,
        flow_rate=flow_rate,
        darcy_friction_factor=discharge_friction_factor,
        zero_flow=zero_flow,
    )

    suction_k = (
        0.0
        if zero_flow
        else suction_friction_factor * suction_pipe_length / suction_pipe_diameter + suction_fitting_k
    )
    discharge_k = (
        0.0
        if zero_flow
        else discharge_friction_factor * discharge_pipe_length / discharge_pipe_diameter + discharge_fitting_k
    )
    suction_head_loss = fluids.head_from_K(K=suction_k, V=suction_pipe_velocity) if not zero_flow else 0.0
    discharge_head_loss = fluids.head_from_K(K=discharge_k, V=discharge_pipe_velocity) if not zero_flow else 0.0

    suction_velocity_head = suction_nozzle_velocity**2 / (2.0 * G_GRAVITY)
    discharge_velocity_head = discharge_nozzle_velocity**2 / (2.0 * G_GRAVITY)
    velocity_head_differential = discharge_velocity_head - suction_velocity_head
    static_head = static_discharge_head - static_suction_head
    surface_pressure_head_differential = (
        discharge_surface_pressure - suction_surface_pressure
    ) / (fluid_density * G_GRAVITY)
    nozzle_pressure_head_differential = None
    if suction_nozzle_pressure is not None and discharge_nozzle_pressure is not None:
        nozzle_pressure_head_differential = (
            discharge_nozzle_pressure - suction_nozzle_pressure
        ) / (fluid_density * G_GRAVITY)

    vapor_pressure_head = fluid_vapor_pressure / (fluid_density * G_GRAVITY)
    if endpoint_basis == "tank_surface":
        tdh = (
            static_head
            + surface_pressure_head_differential
            + suction_head_loss
            + discharge_head_loss
        )
        # Bernoulli from the suction free surface to the flange is a total-head
        # balance; the suction-flange velocity head therefore cancels from NPSHa.
        npsha = (
            suction_surface_pressure / (fluid_density * G_GRAVITY)
            + static_suction_head
            - suction_head_loss
            - vapor_pressure_head
        )
        pressure_head_differential = surface_pressure_head_differential
        velocity_head_contribution = 0.0
        npsha_basis = "suction tank surface to pump flange (total head)"
    else:
        if suction_nozzle_pressure is None or discharge_nozzle_pressure is None:
            raise ValueError("Nozzle endpoint pressures are required for endpoint_basis='nozzle'.")
        tdh = static_head + nozzle_pressure_head_differential + velocity_head_differential
        npsha = (
            suction_nozzle_pressure / (fluid_density * G_GRAVITY)
            + suction_velocity_head
            - vapor_pressure_head
        )
        pressure_head_differential = nozzle_pressure_head_differential
        velocity_head_contribution = velocity_head_differential
        npsha_basis = "pump suction flange total head"

    same_pipe_diameter = math.isclose(
        suction_pipe_diameter, discharge_pipe_diameter, rel_tol=0.0, abs_tol=1e-15
    )
    same_pipe_roughness = math.isclose(
        suction_pipe_roughness, discharge_pipe_roughness, rel_tol=0.0, abs_tol=1e-15
    )
    legacy_pipe_diameter = suction_pipe_diameter if same_pipe_diameter else None
    legacy_pipe_velocity = suction_pipe_velocity if same_pipe_diameter else None
    legacy_reynolds_number = suction_reynolds_number if same_pipe_diameter else None
    legacy_friction_factor = suction_friction_factor if same_pipe_diameter and same_pipe_roughness else None

    pump_efficiency = 0.70
    brake_power_watts = fluid_density * G_GRAVITY * flow_rate * tdh / pump_efficiency
    return {
        "warnings": warnings,
        "endpoint_basis": endpoint_basis,
        "flow_rate_m3s": flow_rate,
        "flow_rate_gpm": flow_rate / GPM_to_M3S,
        "pipe_diameter_m": legacy_pipe_diameter,
        "pipe_diameter_in": legacy_pipe_diameter / INCH_to_M if legacy_pipe_diameter is not None else None,
        "pipe_velocity_m_s": legacy_pipe_velocity,
        "suction_pipe_diameter_m": suction_pipe_diameter,
        "suction_pipe_diameter_in": suction_pipe_diameter / INCH_to_M,
        "discharge_pipe_diameter_m": discharge_pipe_diameter,
        "discharge_pipe_diameter_in": discharge_pipe_diameter / INCH_to_M,
        "suction_pipe_roughness_m": suction_pipe_roughness,
        "discharge_pipe_roughness_m": discharge_pipe_roughness,
        "suction_pipe_velocity_m_s": suction_pipe_velocity,
        "discharge_pipe_velocity_m_s": discharge_pipe_velocity,
        "suction_nozzle_velocity_m_s": suction_nozzle_velocity,
        "discharge_nozzle_velocity_m_s": discharge_nozzle_velocity,
        "suction_nozzle_diameter_m": suction_nozzle,
        "discharge_nozzle_diameter_m": discharge_nozzle,
        "reynolds_number": legacy_reynolds_number,
        "friction_factor": legacy_friction_factor,
        "suction_reynolds_number": suction_reynolds_number,
        "discharge_reynolds_number": discharge_reynolds_number,
        "suction_friction_factor": suction_friction_factor,
        "discharge_friction_factor": discharge_friction_factor,
        "suction_velocity_head_m": suction_velocity_head,
        "discharge_velocity_head_m": discharge_velocity_head,
        "velocity_head_differential_m": velocity_head_differential,
        "velocity_head_contribution_m": velocity_head_contribution,
        "pressure_head_differential_m": pressure_head_differential,
        "surface_pressure_head_differential_m": surface_pressure_head_differential,
        "nozzle_pressure_head_differential_m": nozzle_pressure_head_differential,
        "suction_head_loss_m": suction_head_loss,
        "suction_head_loss_ft": suction_head_loss / FT_to_M,
        "suction_fitting_details": suction_fitting_details,
        "discharge_head_loss_m": discharge_head_loss,
        "discharge_head_loss_ft": discharge_head_loss / FT_to_M,
        "discharge_fitting_details": discharge_fitting_details,
        "suction_velocity_head_ft": suction_velocity_head / FT_to_M,
        "discharge_velocity_head_ft": discharge_velocity_head / FT_to_M,
        "velocity_head_differential_ft": velocity_head_differential / FT_to_M,
        "static_head_m": static_head,
        "static_head_ft": static_head / FT_to_M,
        "total_dynamic_head_m": tdh,
        "total_dynamic_head_ft": tdh / FT_to_M,
        "npsha_m": npsha,
        "npsha_ft": npsha / FT_to_M,
        "estimated_power_kw": brake_power_watts / 1000.0,
        "estimated_power_hp": brake_power_watts / 745.7,
        "assumed_pump_efficiency": pump_efficiency,
        "atmospheric_pressure_pa_used": atmospheric_pressure,
        "atmospheric_pressure_source": atmospheric_pressure_source,
        "suction_surface_pressure_pa_used": suction_surface_pressure,
        "suction_surface_pressure_source": suction_surface_pressure_source,
        "discharge_surface_pressure_pa_used": discharge_surface_pressure,
        "discharge_surface_pressure_source": discharge_surface_pressure_source,
        "suction_nozzle_pressure_pa_used": suction_nozzle_pressure,
        "suction_nozzle_pressure_source": suction_nozzle_pressure_source,
        "discharge_nozzle_pressure_pa_used": discharge_nozzle_pressure,
        "discharge_nozzle_pressure_source": discharge_nozzle_pressure_source,
        "npsha_basis": npsha_basis,
    }


def _solve_flow_for_target(
    *,
    target_tdh: float,
    suction_pipe_diameter: float,
    discharge_pipe_diameter: float,
    hydraulic_kwargs: dict[str, Any],
) -> float:
    limiting_area = min(
        math.pi * (suction_pipe_diameter / 2.0) ** 2,
        math.pi * (discharge_pipe_diameter / 2.0) ** 2,
    )
    max_flow = limiting_area * MAX_PIPE_VELOCITY_M_S

    def residual(flow: float) -> float:
        return _pump_hydraulics(
            flow_rate=flow,
            suction_pipe_diameter=suction_pipe_diameter,
            discharge_pipe_diameter=discharge_pipe_diameter,
            **hydraulic_kwargs,
        )["total_dynamic_head_m"] - target_tdh

    low = residual(0.0)
    if abs(low) <= TARGET_TOL_M:
        return 0.0
    high = residual(max_flow)
    if low * high > 0:
        raise ValueError(
            f"Target TDH is outside the flow bracket 0 to {max_flow / GPM_to_M3S:.1f} GPM "
            f"for maximum pipe velocity {MAX_PIPE_VELOCITY_M_S:.1f} m/s."
        )
    return float(brentq(residual, 0.0, max_flow, xtol=1e-9, rtol=1e-8, maxiter=100))


def _solve_diameter_for_target(
    *,
    target_tdh: float,
    flow_rate: float,
    hydraulic_kwargs: dict[str, Any],
) -> float:
    if flow_rate <= 0:
        raise ValueError("Cannot solve pipe diameter for target TDH at zero flow.")
    d_min = max(0.005, math.sqrt(4.0 * flow_rate / (math.pi * MAX_DIAMETER_SOLVE_VELOCITY_M_S)))
    d_max = max(3.0, math.sqrt(4.0 * flow_rate / (math.pi * MIN_DIAMETER_SOLVE_VELOCITY_M_S)))

    def residual(diameter: float) -> float:
        return _pump_hydraulics(
            flow_rate=flow_rate,
            suction_pipe_diameter=diameter,
            discharge_pipe_diameter=diameter,
            **hydraulic_kwargs,
        )["total_dynamic_head_m"] - target_tdh

    low = residual(d_min)
    high = residual(d_max)
    if abs(low) <= TARGET_TOL_M:
        return d_min
    if abs(high) <= TARGET_TOL_M:
        return d_max
    if low * high > 0:
        raise ValueError(
            f"Target TDH is outside the diameter bracket {d_min:.4f} m to {d_max:.4f} m."
        )
    return float(brentq(residual, d_min, d_max, xtol=1e-9, rtol=1e-8, maxiter=100))


def _round_pump_result(result: dict[str, Any]) -> dict[str, Any]:
    rounded = result.copy()
    precisions = {
        "flow_rate_m3s": 6,
        "flow_rate_gpm": 2,
        "pipe_diameter_m": 6,
        "pipe_diameter_in": 4,
        "pipe_velocity_m_s": 3,
        "suction_pipe_diameter_m": 6,
        "suction_pipe_diameter_in": 4,
        "discharge_pipe_diameter_m": 6,
        "discharge_pipe_diameter_in": 4,
        "suction_pipe_roughness_m": 8,
        "discharge_pipe_roughness_m": 8,
        "suction_pipe_velocity_m_s": 3,
        "discharge_pipe_velocity_m_s": 3,
        "suction_nozzle_velocity_m_s": 3,
        "discharge_nozzle_velocity_m_s": 3,
        "suction_nozzle_diameter_m": 4,
        "discharge_nozzle_diameter_m": 4,
        "reynolds_number": 2,
        "friction_factor": 6,
        "suction_reynolds_number": 2,
        "discharge_reynolds_number": 2,
        "suction_friction_factor": 6,
        "discharge_friction_factor": 6,
        "suction_velocity_head_m": 4,
        "discharge_velocity_head_m": 4,
        "velocity_head_differential_m": 4,
        "velocity_head_contribution_m": 4,
        "pressure_head_differential_m": 4,
        "surface_pressure_head_differential_m": 4,
        "nozzle_pressure_head_differential_m": 4,
        "suction_head_loss_m": 3,
        "suction_head_loss_ft": 3,
        "discharge_head_loss_m": 3,
        "discharge_head_loss_ft": 3,
        "suction_velocity_head_ft": 4,
        "discharge_velocity_head_ft": 4,
        "velocity_head_differential_ft": 4,
        "static_head_m": 3,
        "static_head_ft": 3,
        "total_dynamic_head_m": 3,
        "total_dynamic_head_ft": 3,
        "npsha_m": 3,
        "npsha_ft": 3,
        "npshr_m": 3,
        "npshr_ft": 3,
        "npsh_margin_m": 3,
        "npsh_margin_ft": 3,
        "required_npsh_margin_m": 3,
        "required_npsh_margin_ft": 3,
        "npsh_margin_shortfall_m": 3,
        "npsh_margin_ratio": 4,
        "estimated_power_kw": 2,
        "estimated_power_hp": 2,
        "atmospheric_pressure_pa_used": 2,
        "suction_surface_pressure_pa_used": 2,
        "discharge_surface_pressure_pa_used": 2,
        "suction_nozzle_pressure_pa_used": 2,
        "discharge_nozzle_pressure_pa_used": 2,
    }
    for key, precision in precisions.items():
        if key in rounded and isinstance(rounded[key], (int, float)):
            rounded[key] = round(float(rounded[key]), precision)
    return rounded


def _negative_pump_result_details(result: dict[str, Any]) -> Optional[dict[str, Any]]:
    tdh = result.get("total_dynamic_head_m")
    power_kw = result.get("estimated_power_kw")
    if tdh is not None and tdh < -TARGET_TOL_M:
        return {
            "total_dynamic_head_m": tdh,
            "total_dynamic_head_ft": result.get("total_dynamic_head_ft"),
            "estimated_power_kw": power_kw,
            "static_head_m": result.get("static_head_m"),
            "suction_head_loss_m": result.get("suction_head_loss_m"),
            "discharge_head_loss_m": result.get("discharge_head_loss_m"),
            "velocity_head_differential_m": result.get("velocity_head_differential_m"),
        }
    if power_kw is not None and power_kw < -1e-9:
        return {
            "total_dynamic_head_m": tdh,
            "estimated_power_kw": power_kw,
            "static_head_m": result.get("static_head_m"),
        }
    return None


def _resolve_npsh_value(
    value_m: Optional[float],
    value_ft: Optional[float],
    *,
    name: str,
    results_log: list[str],
) -> Optional[float]:
    """Resolve an optional NPSH value to metres without inventing a value."""

    if value_m is not None and value_ft is not None:
        raise CalculationError(
            "CONFLICTING_NPSH_INPUTS",
            f"Provide only one unit basis for {name}: metres or feet.",
            details={"field": name},
        )
    if value_m is not None:
        resolved = float(value_m)
        results_log.append(f"Used provided {name} in metres.")
    elif value_ft is not None:
        resolved = float(value_ft) * FT_to_M
        results_log.append(f"Converted {name} from {value_ft} ft.")
    else:
        return None
    if not math.isfinite(resolved) or resolved < 0:
        raise CalculationError(
            "INVALID_NPSH_INPUT",
            f"{name} must be non-negative and finite.",
            details={name: resolved},
        )
    return resolved


def _apply_npsh_margin_assessment(
    result: dict[str, Any],
    *,
    npshr_m: Optional[float],
    required_margin_m: Optional[float],
    margin_defaulted: bool,
) -> None:
    """Attach an ANSI/HI 9.6.1-style NPSH screening assessment.

    ANSI/HI 9.6.1-2024 defines margin as ``NPSHA - NPSHR`` and recommends
    application-specific values. This tool uses a visible 1.0 m general-service
    screening default only when the caller does not supply a required margin;
    it is not a substitute for the applicable ANSI/HI market/application table.
    Manufacturer NPSHR is never inferred.
    """

    result["npshr_provided"] = npshr_m is not None
    if npshr_m is None:
        result["npshr_status"] = "not_provided"
        result["npsh_margin_assessed"] = False
        return

    margin_m = result["npsha_m"] - npshr_m
    required = DEFAULT_REQUIRED_NPSH_MARGIN_M if required_margin_m is None else required_margin_m
    result.update(
        {
            "npshr_status": "provided",
            "npsh_margin_assessed": True,
            "npshr_m": npshr_m,
            "npshr_ft": npshr_m / FT_to_M,
            "npsh_margin_m": margin_m,
            "npsh_margin_ft": margin_m / FT_to_M,
            "required_npsh_margin_m": required,
            "required_npsh_margin_ft": required / FT_to_M,
            "required_npsh_margin_defaulted": margin_defaulted,
            "npsh_margin_pass": margin_m >= required,
            "npsh_margin_shortfall_m": max(required - margin_m, 0.0),
            "npsh_margin_ratio": result["npsha_m"] / npshr_m if npshr_m > 0 else None,
            "npsh_margin_basis": "ANSI/HI 9.6.1-2024 additive margin screening",
        }
    )


def calculate_pump_requirements(
    flow_rate: Optional[float] = None,
    pipe_diameter: Optional[float] = None,
    suction_pipe_diameter: Optional[float] = None,
    discharge_pipe_diameter: Optional[float] = None,
    suction_pipe_length: Optional[float] = None,
    discharge_pipe_length: Optional[float] = None,
    static_suction_head: Optional[float] = None,
    static_discharge_head: Optional[float] = None,
    fluid_density: Optional[float] = None,
    fluid_viscosity: Optional[float] = None,
    fluid_vapor_pressure: Optional[float] = None,
    pipe_roughness: Optional[float] = None,
    suction_pipe_roughness: Optional[float] = None,
    discharge_pipe_roughness: Optional[float] = None,
    suction_nozzle_diameter: Optional[float] = None,
    discharge_nozzle_diameter: Optional[float] = None,
    flow_rate_gpm: Optional[float] = None,
    pipe_diameter_in: Optional[float] = None,
    suction_pipe_diameter_in: Optional[float] = None,
    discharge_pipe_diameter_in: Optional[float] = None,
    suction_pipe_length_ft: Optional[float] = None,
    discharge_pipe_length_ft: Optional[float] = None,
    static_suction_head_ft: Optional[float] = None,
    static_discharge_head_ft: Optional[float] = None,
    fluid_density_lbft3: Optional[float] = None,
    fluid_viscosity_cp: Optional[float] = None,
    fluid_vapor_pressure_psi: Optional[float] = None,
    suction_nozzle_diameter_in: Optional[float] = None,
    discharge_nozzle_diameter_in: Optional[float] = None,
    fluid_name: Optional[str] = None,
    temperature_c: Optional[float] = None,
    pressure_bar: Optional[float] = None,
    nominal_size_in: Optional[float] = None,
    schedule: str = "40",
    material: Optional[str] = None,
    suction_material: Optional[str] = None,
    discharge_material: Optional[str] = None,
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
    pump_mode: PumpMode = "flow_rate_given",
    endpoint_basis: EndpointBasis = "tank_surface",
    suction_fittings: List[Dict[str, Union[str, int, float]]] = None,
    discharge_fittings: List[Dict[str, Union[str, int, float]]] = None,
) -> str:
    """Calculate pump TDH/NPSHA or solve flow/diameter for a target TDH.

    ``tank_surface`` uses free-surface elevations/pressures plus suction and
    discharge losses. ``nozzle`` uses pump-flange elevations, absolute static
    pressures, and the flange velocity-head differential per ISO 9906/HI.

    When manufacturer ``npshr`` is supplied, the additive margin and pass/fail
    screen follow ANSI/HI 9.6.1-2024. Its recommendations are application
    specific; an omitted caller requirement is surfaced as a 1.0 m
    general-service screening assumption. NPSHR is never fabricated.
    """

    results_log: list[str] = []
    warnings: list[str] = []
    errors: list[str] = []

    try:
        if pump_mode not in ("flow_rate_given", "target_tdh_solve_flow", "target_tdh_solve_diameter"):
            return _json_error(
                "INVALID_PUMP_MODE",
                "pump_mode must be one of: flow_rate_given, target_tdh_solve_flow, target_tdh_solve_diameter",
                results_log,
                details={"pump_mode": pump_mode},
            )
        if endpoint_basis not in ("tank_surface", "nozzle"):
            return _json_error(
                "INVALID_ENDPOINT_BASIS",
                "endpoint_basis must be one of: tank_surface, nozzle",
                results_log,
                details={"endpoint_basis": endpoint_basis},
            )
        pressure_pairs = {
            "atmospheric_pressure": (atmospheric_pressure, atmospheric_pressure_psi),
            "suction_tank_pressure": (suction_tank_pressure, suction_tank_pressure_psi),
            "discharge_tank_pressure": (discharge_tank_pressure, discharge_tank_pressure_psi),
            "suction_nozzle_pressure": (suction_nozzle_pressure, suction_nozzle_pressure_psi),
            "discharge_nozzle_pressure": (discharge_nozzle_pressure, discharge_nozzle_pressure_psi),
        }
        conflicts = [name for name, (si_value, psi_value) in pressure_pairs.items() if si_value is not None and psi_value is not None]
        if conflicts:
            return _json_error(
                "CONFLICTING_PRESSURE_INPUTS",
                "Provide only one unit basis for each pressure input.",
                results_log,
                details={"conflicting_fields": conflicts},
            )

        local_npshr = _resolve_npsh_value(
            npshr,
            npshr_ft,
            name="npshr",
            results_log=results_log,
        )
        local_required_npsh_margin = _resolve_npsh_value(
            required_npsh_margin,
            required_npsh_margin_ft,
            name="required_npsh_margin",
            results_log=results_log,
        )
        npsh_margin_defaulted = local_npshr is not None and local_required_npsh_margin is None
        if npsh_margin_defaulted:
            results_log.append(
                "Defaulted required NPSH margin to 1.0 m for general-service screening; "
                "override with the applicable ANSI/HI 9.6.1 application margin."
            )

        local_target_tdh = _resolve_target_tdh(target_tdh, target_tdh_ft, results_log)
        pump_mode = _infer_pump_mode(
            pump_mode,
            flow_rate=flow_rate,
            flow_rate_gpm=flow_rate_gpm,
            pipe_diameter=pipe_diameter,
            pipe_diameter_in=pipe_diameter_in,
            suction_pipe_diameter=suction_pipe_diameter,
            suction_pipe_diameter_in=suction_pipe_diameter_in,
            discharge_pipe_diameter=discharge_pipe_diameter,
            discharge_pipe_diameter_in=discharge_pipe_diameter_in,
            nominal_size_in=nominal_size_in,
            target_tdh_m=local_target_tdh,
        )
        results_log.append(f"Pump mode: {pump_mode}.")
        results_log.append(f"Endpoint basis: {endpoint_basis}.")

        if pump_mode != "flow_rate_given" and local_target_tdh is None:
            return _json_error(
                "MISSING_TARGET_TDH",
                "target_tdh or target_tdh_ft is required for target-TDH pump modes.",
                results_log,
            )

        side_diameter_supplied = any(
            value is not None
            for value in (
                suction_pipe_diameter,
                suction_pipe_diameter_in,
                discharge_pipe_diameter,
                discharge_pipe_diameter_in,
            )
        )
        if pump_mode == "target_tdh_solve_diameter" and side_diameter_supplied:
            return _json_error(
                "AMBIGUOUS_DIAMETER_SOLVE",
                "target_tdh_solve_diameter solves one common diameter; do not provide side-specific diameters.",
                results_log,
                warnings=warnings,
            )

        local_flow_rate = _resolve_flow(flow_rate, flow_rate_gpm, pump_mode, results_log, errors)
        local_pipe_diameter, local_pipe_roughness, pipe_info = _resolve_pipe(
            pipe_diameter=pipe_diameter,
            pipe_diameter_in=pipe_diameter_in,
            nominal_size_in=nominal_size_in,
            schedule=schedule,
            material=material,
            pipe_roughness=pipe_roughness,
            pump_mode=pump_mode,
            diameter_required=not side_diameter_supplied,
            results_log=results_log,
            errors=errors,
            warnings=warnings,
        )
        local_suction_pipe_diameter, local_suction_pipe_roughness, suction_pipe_info = _resolve_side_pipe(
            side="suction",
            pipe_diameter=suction_pipe_diameter,
            pipe_diameter_in=suction_pipe_diameter_in,
            material=suction_material,
            pipe_roughness=suction_pipe_roughness,
            fallback_diameter=local_pipe_diameter,
            fallback_roughness=local_pipe_roughness,
            fallback_roughness_source=pipe_info.get("roughness_source"),
            diameter_required=pump_mode != "target_tdh_solve_diameter",
            results_log=results_log,
            errors=errors,
            warnings=warnings,
        )
        local_discharge_pipe_diameter, local_discharge_pipe_roughness, discharge_pipe_info = _resolve_side_pipe(
            side="discharge",
            pipe_diameter=discharge_pipe_diameter,
            pipe_diameter_in=discharge_pipe_diameter_in,
            material=discharge_material,
            pipe_roughness=discharge_pipe_roughness,
            fallback_diameter=local_pipe_diameter,
            fallback_roughness=local_pipe_roughness,
            fallback_roughness_source=pipe_info.get("roughness_source"),
            diameter_required=pump_mode != "target_tdh_solve_diameter",
            results_log=results_log,
            errors=errors,
            warnings=warnings,
        )
        local_suction_pipe_length = _resolve_length(
            suction_pipe_length, suction_pipe_length_ft, "suction_pipe_length", results_log, errors
        )
        local_discharge_pipe_length = _resolve_length(
            discharge_pipe_length, discharge_pipe_length_ft, "discharge_pipe_length", results_log, errors
        )
        local_static_suction_head = _resolve_length(
            static_suction_head, static_suction_head_ft, "static_suction_head", results_log, errors
        )
        local_static_discharge_head = _resolve_length(
            static_discharge_head, static_discharge_head_ft, "static_discharge_head", results_log, errors
        )
        local_fluid_density, local_fluid_viscosity, local_fluid_vapor_pressure, fluid_info = _resolve_fluid(
            fluid_density=fluid_density,
            fluid_density_lbft3=fluid_density_lbft3,
            fluid_viscosity=fluid_viscosity,
            fluid_viscosity_cp=fluid_viscosity_cp,
            fluid_vapor_pressure=fluid_vapor_pressure,
            fluid_vapor_pressure_psi=fluid_vapor_pressure_psi,
            fluid_name=fluid_name,
            temperature_c=temperature_c,
            pressure_bar=pressure_bar,
            results_log=results_log,
            errors=errors,
        )
        local_atmospheric_pressure, atmospheric_pressure_source, atmospheric_defaulted = (
            _resolve_atmospheric_pressure(
            atmospheric_pressure=atmospheric_pressure,
            atmospheric_pressure_psi=atmospheric_pressure_psi,
            results_log=results_log,
        )
        )
        if atmospheric_defaulted and endpoint_basis == "tank_surface":
            warnings.append(
                "Atmospheric pressure defaulted to 101325 Pa absolute (sea-level standard); "
                "provide the site atmospheric pressure for auditable NPSHa and gauge-pressure conversion."
            )
        local_suction_surface_pressure, suction_surface_pressure_source = _resolve_endpoint_pressure(
            si_value=suction_tank_pressure,
            psi_value=suction_tank_pressure_psi,
            name="suction_tank_pressure",
            fallback=local_atmospheric_pressure,
            fallback_source=atmospheric_pressure_source,
            results_log=results_log,
        )
        local_discharge_surface_pressure, discharge_surface_pressure_source = _resolve_endpoint_pressure(
            si_value=discharge_tank_pressure,
            psi_value=discharge_tank_pressure_psi,
            name="discharge_tank_pressure",
            fallback=local_atmospheric_pressure,
            fallback_source=atmospheric_pressure_source,
            results_log=results_log,
        )
        local_suction_nozzle_pressure, suction_nozzle_pressure_source = _resolve_endpoint_pressure(
            si_value=suction_nozzle_pressure,
            psi_value=suction_nozzle_pressure_psi,
            name="suction_nozzle_pressure",
            fallback=None,
            fallback_source=None,
            results_log=results_log,
        )
        local_discharge_nozzle_pressure, discharge_nozzle_pressure_source = _resolve_endpoint_pressure(
            si_value=discharge_nozzle_pressure,
            psi_value=discharge_nozzle_pressure_psi,
            name="discharge_nozzle_pressure",
            fallback=None,
            fallback_source=None,
            results_log=results_log,
        )
        if endpoint_basis == "nozzle" and (
            local_suction_nozzle_pressure is None or local_discharge_nozzle_pressure is None
        ):
            return _json_error(
                "MISSING_NOZZLE_PRESSURE",
                "endpoint_basis='nozzle' requires absolute suction_nozzle_pressure and "
                "discharge_nozzle_pressure (or their _psi forms).",
                results_log,
                warnings=warnings,
            )
        resolved_pressures = {
            "atmospheric_pressure": local_atmospheric_pressure,
            "suction_tank_pressure": local_suction_surface_pressure,
            "discharge_tank_pressure": local_discharge_surface_pressure,
        }
        if endpoint_basis == "nozzle":
            resolved_pressures.update(
                {
                    "suction_nozzle_pressure": local_suction_nozzle_pressure,
                    "discharge_nozzle_pressure": local_discharge_nozzle_pressure,
                }
            )
        invalid_pressures = {
            name: value
            for name, value in resolved_pressures.items()
            if value is None
            or not isinstance(value, (int, float))
            or not math.isfinite(float(value))
            or value <= 0
        }
        if invalid_pressures:
            return _json_error(
                "INVALID_PRESSURE_INPUT",
                "Absolute pressures must be finite and greater than zero.",
                results_log,
                details={"invalid_pressures": invalid_pressures},
                warnings=warnings,
            )
        provided_suction_nozzle_diameter, provided_discharge_nozzle_diameter = _resolve_nozzle_inputs(
            suction_nozzle_diameter=suction_nozzle_diameter,
            suction_nozzle_diameter_in=suction_nozzle_diameter_in,
            discharge_nozzle_diameter=discharge_nozzle_diameter,
            discharge_nozzle_diameter_in=discharge_nozzle_diameter_in,
            results_log=results_log,
        )

        if errors:
            return json.dumps({"errors": errors, "warnings": warnings, "log": results_log})

        required = {
            "suction_pipe_length": local_suction_pipe_length,
            "discharge_pipe_length": local_discharge_pipe_length,
            "static_suction_head": local_static_suction_head,
            "static_discharge_head": local_static_discharge_head,
            "fluid_density": local_fluid_density,
            "fluid_viscosity": local_fluid_viscosity,
            "fluid_vapor_pressure": local_fluid_vapor_pressure,
            "suction_pipe_roughness": local_suction_pipe_roughness,
            "discharge_pipe_roughness": local_discharge_pipe_roughness,
            "atmospheric_pressure": local_atmospheric_pressure,
            "suction_surface_pressure": local_suction_surface_pressure,
            "discharge_surface_pressure": local_discharge_surface_pressure,
        }
        if pump_mode != "target_tdh_solve_flow":
            required["flow_rate"] = local_flow_rate
        if pump_mode != "target_tdh_solve_diameter":
            required["suction_pipe_diameter"] = local_suction_pipe_diameter
            required["discharge_pipe_diameter"] = local_discharge_pipe_diameter
        missing = [key for key, value in required.items() if value is None]
        if missing:
            return _json_error(
                "MISSING_REQUIRED_INPUT",
                f"Missing resolved values for: {', '.join(missing)}",
                results_log,
                warnings=warnings,
            )

        if local_flow_rate is not None:
            if not isinstance(local_flow_rate, (int, float)) or not math.isfinite(float(local_flow_rate)):
                return _json_error(
                    "INVALID_NUMERIC_INPUT",
                    "flow_rate must be finite.",
                    results_log,
                    details={"flow_rate": repr(local_flow_rate)},
                    warnings=warnings,
                )
            if local_flow_rate < 0:
                return _json_error(
                    "INVALID_NUMERIC_INPUT",
                    "flow_rate must be non-negative.",
                    results_log,
                    details={"flow_rate": local_flow_rate},
                    warnings=warnings,
                )

        static_head = local_static_discharge_head - local_static_suction_head
        if endpoint_basis == "tank_surface":
            zero_flow_endpoint_head = static_head + (
                local_discharge_surface_pressure - local_suction_surface_pressure
            ) / (local_fluid_density * G_GRAVITY)
        else:
            zero_flow_endpoint_head = static_head + (
                local_discharge_nozzle_pressure - local_suction_nozzle_pressure
            ) / (local_fluid_density * G_GRAVITY)
        if local_target_tdh is not None and local_target_tdh < zero_flow_endpoint_head - TARGET_TOL_M:
            return _json_error(
                "TARGET_BELOW_STATIC_HEAD",
                "Target TDH is below static head after endpoint-pressure correction; "
                "no non-negative-flow pump solution exists.",
                results_log,
                details={
                    "target_tdh_m": local_target_tdh,
                    "target_tdh_ft": local_target_tdh / FT_to_M,
                    "static_head_m": static_head,
                    "static_head_ft": static_head / FT_to_M,
                    "zero_flow_endpoint_head_m": zero_flow_endpoint_head,
                },
                warnings=warnings,
            )

        hydraulic_kwargs = {
            "suction_pipe_length": local_suction_pipe_length,
            "discharge_pipe_length": local_discharge_pipe_length,
            "static_suction_head": local_static_suction_head,
            "static_discharge_head": local_static_discharge_head,
            "fluid_density": local_fluid_density,
            "fluid_viscosity": local_fluid_viscosity,
            "fluid_vapor_pressure": local_fluid_vapor_pressure,
            "suction_pipe_roughness": local_suction_pipe_roughness,
            "discharge_pipe_roughness": local_discharge_pipe_roughness,
            "endpoint_basis": endpoint_basis,
            "atmospheric_pressure": local_atmospheric_pressure,
            "atmospheric_pressure_source": atmospheric_pressure_source,
            "suction_surface_pressure": local_suction_surface_pressure,
            "suction_surface_pressure_source": suction_surface_pressure_source,
            "discharge_surface_pressure": local_discharge_surface_pressure,
            "discharge_surface_pressure_source": discharge_surface_pressure_source,
            "suction_nozzle_pressure": local_suction_nozzle_pressure,
            "suction_nozzle_pressure_source": suction_nozzle_pressure_source,
            "discharge_nozzle_pressure": local_discharge_nozzle_pressure,
            "discharge_nozzle_pressure_source": discharge_nozzle_pressure_source,
            "suction_fittings": suction_fittings if suction_fittings else [],
            "discharge_fittings": discharge_fittings if discharge_fittings else [],
            "suction_nozzle_diameter": provided_suction_nozzle_diameter,
            "discharge_nozzle_diameter": provided_discharge_nozzle_diameter,
        }

        if pump_mode == "target_tdh_solve_flow":
            try:
                local_flow_rate = _solve_flow_for_target(
                    target_tdh=local_target_tdh,
                    suction_pipe_diameter=local_suction_pipe_diameter,
                    discharge_pipe_diameter=local_discharge_pipe_diameter,
                    hydraulic_kwargs=hydraulic_kwargs,
                )
                results_log.append(
                    f"Solved flow_rate for target TDH: {local_flow_rate:.8f} m3/s "
                    f"({local_flow_rate / GPM_to_M3S:.2f} GPM)."
                )
            except ValueError as exc:
                return _json_error("TARGET_TDH_NOT_BRACKETED", str(exc), results_log, warnings=warnings)
        elif pump_mode == "target_tdh_solve_diameter":
            try:
                local_pipe_diameter = _solve_diameter_for_target(
                    target_tdh=local_target_tdh,
                    flow_rate=local_flow_rate,
                    hydraulic_kwargs=hydraulic_kwargs,
                )
                results_log.append(
                    f"Solved pipe_diameter for target TDH: {local_pipe_diameter:.8f} m "
                    f"({local_pipe_diameter / INCH_to_M:.3f} in)."
                )
                local_suction_pipe_diameter = local_pipe_diameter
                local_discharge_pipe_diameter = local_pipe_diameter
            except ValueError as exc:
                return _json_error("TARGET_TDH_NOT_BRACKETED", str(exc), results_log, warnings=warnings)

        raw_result = _pump_hydraulics(
            flow_rate=local_flow_rate,
            suction_pipe_diameter=local_suction_pipe_diameter,
            discharge_pipe_diameter=local_discharge_pipe_diameter,
            **hydraulic_kwargs,
        )
        warnings.extend(raw_result.pop("warnings", []))
        _apply_npsh_margin_assessment(
            raw_result,
            npshr_m=local_npshr,
            required_margin_m=local_required_npsh_margin,
            margin_defaulted=npsh_margin_defaulted,
        )
        if local_npshr is None:
            warnings.append(
                "Manufacturer NPSHr was not provided; NPSHa is reported, but NPSH margin is not held or assessed."
            )
        negative_result = _negative_pump_result_details(raw_result)
        if negative_result is not None:
            return _json_error(
                "NEGATIVE_PUMP_TDH",
                "Pump calculation produced negative TDH or power; use a turbine/gravity-flow model for this operating point.",
                results_log,
                details=negative_result,
                warnings=warnings,
            )
        results_log.append(
            f"TDH components: Static={raw_result['static_head_m']:.3f}m, "
            f"Pressure={raw_result['pressure_head_differential_m']:.3f}m, "
            f"Friction contribution="
            f"{(raw_result['discharge_head_loss_m'] + raw_result['suction_head_loss_m']) if endpoint_basis == 'tank_surface' else 0.0:.3f}m, "
            f"Velocity contribution={raw_result['velocity_head_contribution_m']:.3f}m."
        )
        results_log.append(f"NPSHA basis: {raw_result['npsha_basis']}.")
        results_log.append(
            f"NPSHA components resolved with NPSHA={raw_result['npsha_m']:.3f}m."
        )
        if raw_result["npsh_margin_assessed"]:
            results_log.append(
                f"NPSH margin={raw_result['npsh_margin_m']:.3f} m versus required "
                f"{raw_result['required_npsh_margin_m']:.3f} m: "
                f"{'PASS' if raw_result['npsh_margin_pass'] else 'FAIL'}."
            )
        results_log.append("Power calculation based on assumed efficiency: 70.0%")

        final_result = _round_pump_result(raw_result)
        final_result.update(
            {
                "inputs_resolved": results_log,
                "pump_mode": pump_mode,
                "endpoint_basis": endpoint_basis,
                "warnings": warnings if warnings else None,
                "pipe_details": {
                    "common": pipe_info,
                    "suction": suction_pipe_info,
                    "discharge": discharge_pipe_info,
                },
                "fluid_details": fluid_info if fluid_info else None,
                "target_tdh_m": round(local_target_tdh, 3) if local_target_tdh is not None else None,
                "target_tdh_ft": round(local_target_tdh / FT_to_M, 3) if local_target_tdh is not None else None,
            }
        )
        if not final_result["warnings"]:
            final_result.pop("warnings", None)
        if not final_result["pipe_details"]:
            final_result.pop("pipe_details", None)
        if not final_result["fluid_details"]:
            final_result.pop("fluid_details", None)
        if final_result["target_tdh_m"] is None:
            final_result.pop("target_tdh_m", None)
            final_result.pop("target_tdh_ft", None)
        return json.dumps(final_result, allow_nan=False)
    except CalculationError as exc:
        logger.warning("Typed pump calculation error %s: %s", exc.code, exc.message)
        return _json_calculation_error(exc, results_log, warnings=warnings)
    except Exception as exc:
        logger.error("Error in calculate_pump_requirements: %s", exc, exc_info=True)
        return _json_error(
            "PUMP_CALCULATION_FAILED",
            f"Calculation error: {exc}",
            results_log,
            warnings=warnings,
        )
