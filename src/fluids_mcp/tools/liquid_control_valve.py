"""IEC/ISA liquid control-valve sizing."""

from __future__ import annotations

import json
import logging
from typing import Literal, Optional

import fluids.control_valve

from fluids_mcp.domain.units import GPM_TO_M3S, LBFT3_TO_KGM3, PSI_TO_PA
from fluids_mcp.domain.valves import (
    initial_liquid_nps,
    liquid_velocity_m_s,
    next_standard_nps,
    rated_cv,
    valve_default,
)
from fluids_mcp.utils.constants import CENTIPOISE_to_PAS

logger = logging.getLogger("fluids-mcp.calculate_liquid_control_valve")

KV_TO_CV = 1.0 / 0.865
POOR_CONTROLLABILITY_LIMIT = 0.15
MAX_UTILIZATION = 0.80


def _json_error(code: str, message: str, log: list[str], *, details: Optional[dict] = None) -> str:
    return json.dumps(
        {
            "error_code": code,
            "error": message,
            "errors": [{"code": code, "message": message, "details": details or {}}],
            "log": log,
        }
    )


def _resolve_liquid_properties(
    *,
    fluid_density: Optional[float],
    fluid_density_lbft3: Optional[float],
    fluid_viscosity: Optional[float],
    fluid_viscosity_cp: Optional[float],
    fluid_saturation_pressure_pa: Optional[float],
    fluid_saturation_pressure_psi: Optional[float],
    fluid_critical_pressure_pa: Optional[float],
    fluid_critical_pressure_psi: Optional[float],
    fluid_name: Optional[str],
    temperature_c: Optional[float],
    pressure_bar: Optional[float],
    log: list[str],
    warnings: list[str],
) -> tuple[Optional[float], Optional[float], Optional[float], Optional[float], dict]:
    density = None
    viscosity = None
    psat = None
    pc = None
    info = {}

    if fluid_density is not None and fluid_viscosity is not None:
        density = fluid_density
        viscosity = fluid_viscosity
        log.append("Used provided SI fluid density and viscosity.")
    elif fluid_density_lbft3 is not None and fluid_viscosity_cp is not None:
        density = fluid_density_lbft3 * LBFT3_TO_KGM3
        viscosity = fluid_viscosity_cp * CENTIPOISE_to_PAS
        log.append(f"Converted density from {fluid_density_lbft3} lb/ft3 and viscosity from {fluid_viscosity_cp} cP.")
    elif fluid_name is not None and temperature_c is not None:
        try:
            from fluids_mcp.utils.resolve_properties import resolve_liquid_properties

            p_bar = pressure_bar if pressure_bar is not None else 1.01325
            props = resolve_liquid_properties(fluid_name, temperature_c, p_bar)
            if props is not None:
                density = props.density
                viscosity = props.viscosity
                psat = props.vapor_pressure
                pc = props.critical_pressure
                info = {
                    "name_used": fluid_name,
                    "temperature_c": temperature_c,
                    "pressure_bar": p_bar,
                    "density_kgm3": density,
                    "viscosity_pas": viscosity,
                    "saturation_pressure_pa": psat,
                    "critical_pressure_pa": pc,
                }
                log.append(f"Properties via {props.source} for {fluid_name} at {temperature_c} C, {p_bar} bar.")
                log.extend(props.log)
            else:
                warnings.append(f"Fluid property lookup failed for {fluid_name}; provide density/viscosity/Psat/Pc.")
        except Exception as exc:
            warnings.append(f"Fluid property lookup failed for {fluid_name}: {exc}")

    if psat is None:
        if fluid_saturation_pressure_pa is not None:
            psat = fluid_saturation_pressure_pa
            log.append("Used provided SI fluid_saturation_pressure_pa.")
        elif fluid_saturation_pressure_psi is not None:
            psat = fluid_saturation_pressure_psi * PSI_TO_PA
            log.append(f"Converted fluid_saturation_pressure from {fluid_saturation_pressure_psi} psi.")
        elif fluid_name and fluid_name.lower() in {"water", "h2o"}:
            psat = 3170.0
            warnings.append("Using fallback water saturation pressure at about 25 C.")

    if pc is None:
        if fluid_critical_pressure_pa is not None:
            pc = fluid_critical_pressure_pa
            log.append("Used provided SI fluid_critical_pressure_pa.")
        elif fluid_critical_pressure_psi is not None:
            pc = fluid_critical_pressure_psi * PSI_TO_PA
            log.append(f"Converted fluid_critical_pressure from {fluid_critical_pressure_psi} psi.")
        elif fluid_name and fluid_name.lower() in {"water", "h2o"}:
            pc = 22_064_000.0
            warnings.append("Using fallback water critical pressure.")

    return density, viscosity, psat, pc, {key: value for key, value in info.items() if value is not None}


def _select_liquid_nps(flow_m3s: float, required_cv: float, valve_type: str, warnings: list[str]) -> tuple[float, float, float]:
    selected = initial_liquid_nps(flow_m3s)
    while True:
        rated = rated_cv(valve_type, selected)  # type: ignore[arg-type]
        utilization = required_cv / rated if rated > 0 else float("inf")
        if utilization <= MAX_UTILIZATION:
            if utilization < POOR_CONTROLLABILITY_LIMIT:
                warnings.append("Poor controllability: valve is oversized for service at less than 15% Cv utilization.")
            return selected, rated, utilization
        next_nps = next_standard_nps(selected)
        if next_nps is None:
            warnings.append("Required Cv exceeds 80% utilization at largest available valve table size.")
            return selected, rated, utilization
        selected = next_nps


def calculate_liquid_control_valve(
    flow_rate: Optional[float] = None,
    inlet_pressure: Optional[float] = None,
    outlet_pressure: Optional[float] = None,
    fluid_density: Optional[float] = None,
    fluid_viscosity: Optional[float] = None,
    fluid_saturation_pressure_pa: Optional[float] = None,
    fluid_critical_pressure_pa: Optional[float] = None,
    flow_rate_gpm: Optional[float] = None,
    inlet_pressure_psi: Optional[float] = None,
    outlet_pressure_psi: Optional[float] = None,
    fluid_density_lbft3: Optional[float] = None,
    fluid_viscosity_cp: Optional[float] = None,
    fluid_saturation_pressure_psi: Optional[float] = None,
    fluid_critical_pressure_psi: Optional[float] = None,
    fluid_name: Optional[str] = None,
    temperature_c: Optional[float] = None,
    pressure_bar: Optional[float] = None,
    valve_type: Literal["globe", "butterfly", "ball", "gate", "other"] = "globe",
    valve_fl: Optional[float] = None,
    valve_fd: Optional[float] = None,
) -> str:
    """Calculate liquid control-valve Cv and select NPS by rated Cv utilization."""

    log: list[str] = []
    warnings: list[str] = []

    try:
        flow_m3s = flow_rate if flow_rate is not None else None
        if flow_m3s is None and flow_rate_gpm is not None:
            flow_m3s = flow_rate_gpm * GPM_TO_M3S
            log.append(f"Converted flow_rate from {flow_rate_gpm} GPM.")
        if flow_m3s is None:
            return _json_error("MISSING_REQUIRED_INPUT", "flow_rate or flow_rate_gpm is required.", log)
        if flow_m3s < 0:
            return _json_error("INVALID_FLOW_RATE", "Liquid control-valve flow must be non-negative.", log)

        p1 = inlet_pressure if inlet_pressure is not None else None
        if p1 is None and inlet_pressure_psi is not None:
            p1 = inlet_pressure_psi * PSI_TO_PA
            log.append(f"Converted inlet_pressure from {inlet_pressure_psi} psi.")
        p2 = outlet_pressure if outlet_pressure is not None else None
        if p2 is None and outlet_pressure_psi is not None:
            p2 = outlet_pressure_psi * PSI_TO_PA
            log.append(f"Converted outlet_pressure from {outlet_pressure_psi} psi.")
        if p1 is None or p2 is None:
            return _json_error("MISSING_REQUIRED_INPUT", "inlet and outlet pressures are required.", log)
        if p2 >= p1:
            return _json_error("INVALID_PRESSURE_DIRECTION", "Inlet pressure must be greater than outlet pressure.", log)

        density, viscosity, psat, pc, fluid_info = _resolve_liquid_properties(
            fluid_density=fluid_density,
            fluid_density_lbft3=fluid_density_lbft3,
            fluid_viscosity=fluid_viscosity,
            fluid_viscosity_cp=fluid_viscosity_cp,
            fluid_saturation_pressure_pa=fluid_saturation_pressure_pa,
            fluid_saturation_pressure_psi=fluid_saturation_pressure_psi,
            fluid_critical_pressure_pa=fluid_critical_pressure_pa,
            fluid_critical_pressure_psi=fluid_critical_pressure_psi,
            fluid_name=fluid_name,
            temperature_c=temperature_c,
            pressure_bar=pressure_bar,
            log=log,
            warnings=warnings,
        )
        if None in (density, viscosity, psat, pc):
            return _json_error(
                "MISSING_FLUID_PROPERTIES",
                "density, viscosity, saturation pressure, and critical pressure are required for liquid valve sizing.",
                log,
                details={"density": density, "viscosity": viscosity, "Psat": psat, "Pc": pc},
            )

        fl = valve_fl if valve_fl is not None else valve_default(valve_type, "FL")
        fd = valve_fd if valve_fd is not None else valve_default(valve_type, "Fd")
        log.append(f"Using FL={fl}, Fd={fd} for {valve_type} valve.")

        if flow_m3s == 0:
            required_kv = 0.0
            liquid_output = {"Kv": 0.0, "choked": False, "laminar": False}
        else:
            liquid_output = fluids.control_valve.size_control_valve_l(
                rho=density,
                Psat=psat,
                Pc=pc,
                mu=viscosity,
                P1=p1,
                P2=p2,
                Q=flow_m3s,
                FL=fl,
                Fd=fd,
                full_output=True,
            )
            required_kv = float(liquid_output["Kv"])
        required_cv = required_kv * KV_TO_CV
        selected_nps, selected_rated_cv, utilization = _select_liquid_nps(
            flow_m3s,
            required_cv,
            valve_type,
            warnings,
        )
        flashing = p2 <= psat
        cavitation = bool(liquid_output.get("choked")) and not flashing
        velocity = liquid_velocity_m_s(flow_m3s, selected_nps)

        result = {
            "inputs_resolved": log,
            "warnings": warnings if warnings else None,
            "fluid_details": fluid_info if fluid_info else None,
            "flow_rate_m3s": round(flow_m3s, 8),
            "flow_rate_gpm": round(flow_m3s / GPM_TO_M3S, 3),
            "inlet_pressure_pa": round(p1, 3),
            "outlet_pressure_pa": round(p2, 3),
            "pressure_drop_pa": round(p1 - p2, 3),
            "pressure_drop_psi": round((p1 - p2) / PSI_TO_PA, 3),
            "valve_type": valve_type,
            "recommended_nps_inches": selected_nps,
            "required_cv": round(required_cv, 4),
            "rated_cv": round(selected_rated_cv, 4),
            "cv_utilization": round(utilization, 4),
            "liquid_velocity_m_s": round(velocity, 4),
            "FL_used": fl,
            "Fd_used": fd,
            "cavitation": cavitation,
            "flashing": flashing,
            "choked": bool(liquid_output.get("choked")),
            "laminar": bool(liquid_output.get("laminar")),
            "iec_details": liquid_output,
        }
        if not result["warnings"]:
            result.pop("warnings", None)
        if not result["fluid_details"]:
            result.pop("fluid_details", None)
        return json.dumps(result, allow_nan=False)
    except Exception as exc:
        logger.error("Error in calculate_liquid_control_valve: %s", exc, exc_info=True)
        return _json_error("LIQUID_CONTROL_VALVE_FAILED", str(exc), log)
