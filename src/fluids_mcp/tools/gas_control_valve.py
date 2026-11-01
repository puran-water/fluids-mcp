"""IEC/ISA gas control-valve sizing with explicit choked-flow output."""

from __future__ import annotations

import json
import logging
import math
from typing import Dict, Literal, Optional

import fluids.control_valve

from fluids_mcp.domain.units import PSI_TO_PA
from fluids_mcp.domain.valves import (
    MAX_GAS_OUTLET_MACH,
    initial_gas_nps,
    next_standard_nps,
    nps_inside_diameter_m,
    rated_cv,
    valve_default,
)
from fluids_mcp.utils.constants import DEG_C_to_K

logger = logging.getLogger("fluids-mcp.gas_control_valve")

KV_TO_CV = 1.0 / 0.865
MAX_UTILIZATION = 0.80
POOR_CONTROLLABILITY_LIMIT = 0.15
P_NORM = 101325.0
T_NORM = 273.15
P_STD = 101325.0
T_STD = 288.15
R_UNIV = 8314.462


def _json_error(code: str, message: str, log: list[str], *, details: Optional[dict] = None) -> str:
    return json.dumps(
        {
            "error_code": code,
            "error": message,
            "errors": [{"code": code, "message": message, "details": details or {}}],
            "log": log,
        }
    )


def _resolve_gas_props(
    *,
    gas_mw: Optional[float],
    gas_gamma: Optional[float],
    gas_z_factor: Optional[float],
    gas_viscosity: Optional[float],
    fluid_name: Optional[str],
    gas_composition_mol: Optional[Dict[str, float]],
    inlet_temperature_c: float,
    inlet_pressure_pa: float,
    allow_property_defaults: bool,
    log: list[str],
    warnings: list[str],
) -> tuple[Optional[float], Optional[float], Optional[float], Optional[float], dict]:
    mw = gas_mw
    gamma = gas_gamma
    z_factor = gas_z_factor
    viscosity = gas_viscosity
    info: dict = {}
    sources: dict = {}

    if any(value is None for value in (mw, gamma, z_factor, viscosity)) and fluid_name is not None:
        try:
            from fluids_mcp.utils.resolve_properties import resolve_gas_properties

            lookup_p_bar = inlet_pressure_pa / 100000.0
            props = resolve_gas_properties(fluid_name, inlet_temperature_c, lookup_p_bar)
            if props is not None:
                if mw is None:
                    mw = props.mw
                    sources["MW"] = props.source
                if gamma is None:
                    gamma = props.gamma
                    sources["gamma"] = props.source
                if z_factor is None:
                    z_factor = props.z_factor
                    sources["Z"] = props.source
                if viscosity is None:
                    viscosity = props.viscosity
                    sources["viscosity"] = props.source
                info.update(
                    {
                        "name_used": fluid_name,
                        "lookup_temp_c": inlet_temperature_c,
                        "lookup_p_bar": lookup_p_bar,
                    }
                )
                log.extend(props.log)
            else:
                warnings.append(f"Property lookup failed for '{fluid_name}'.")
        except Exception as exc:
            warnings.append(f"Property lookup failed for '{fluid_name}': {exc}")

    if gas_composition_mol is not None and mw is None:
        mw_map = {"CH4": 16.04, "CO2": 44.01, "N2": 28.01, "H2S": 34.08, "O2": 32.0, "H2": 2.016}
        mw = sum(mw_map.get(component, 28.96) * fraction for component, fraction in gas_composition_mol.items())
        sources["MW"] = "Composition estimate"
        info["composition_mol"] = gas_composition_mol
    if gas_composition_mol is not None:
        if gamma is None:
            gamma = 1.3
            sources["gamma"] = "Composition estimate"
        if z_factor is None:
            z_factor = 0.98
            sources["Z"] = "Composition estimate"
        if viscosity is None:
            viscosity = 1.1e-5
            sources["viscosity"] = "Composition estimate"
        warnings.append("Using estimated gas gamma/Z/viscosity for mixture where lookup was unavailable.")

    missing = []
    if mw is None:
        if allow_property_defaults:
            mw = 28.96
            sources["MW"] = "Default air"
        else:
            missing.append("gas_mw")
    if gamma is None:
        if allow_property_defaults:
            gamma = 1.4
            sources["gamma"] = "Default air"
        else:
            missing.append("gas_gamma")
    if z_factor is None:
        if allow_property_defaults:
            z_factor = 1.0
            sources["Z"] = "Default ideal gas"
        else:
            missing.append("gas_z_factor")
    if viscosity is None:
        if allow_property_defaults:
            viscosity = 1.8e-5
            sources["viscosity"] = "Default air"
        else:
            missing.append("gas_viscosity")
    if missing:
        return mw, gamma, z_factor, viscosity, {"missing": missing}

    info.update(
        {
            "MW_kg_kmol": mw,
            "gamma": gamma,
            "Z_factor": z_factor,
            "viscosity_pas": viscosity,
            "source_details": sources,
        }
    )
    return mw, gamma, z_factor, viscosity, info


def _normal_flow_and_mass(
    *,
    flow_rate_kg_s: Optional[float],
    flow_rate_norm_m3_hr: Optional[float],
    flow_rate_std_m3_hr: Optional[float],
    mw: float,
) -> tuple[Optional[float], Optional[float]]:
    rho_norm = P_NORM * mw / (R_UNIV * T_NORM)
    rho_std = P_STD * mw / (R_UNIV * T_STD)
    if flow_rate_kg_s is not None:
        return flow_rate_kg_s / rho_norm, flow_rate_kg_s
    if flow_rate_norm_m3_hr is not None:
        q_norm_m3s = flow_rate_norm_m3_hr / 3600.0
        return q_norm_m3s, q_norm_m3s * rho_norm
    if flow_rate_std_m3_hr is not None:
        q_std_m3s = flow_rate_std_m3_hr / 3600.0
        mass = q_std_m3s * rho_std
        return mass / rho_norm, mass
    return None, None


def _mach_outlet(
    *,
    mass_flow_kg_s: float,
    outlet_pressure_pa: float,
    temperature_k: float,
    mw: float,
    gamma: float,
    z_factor: float,
    nps: float,
) -> float:
    rho_out = outlet_pressure_pa * mw / (z_factor * R_UNIV * temperature_k)
    area = math.pi * (nps_inside_diameter_m(nps) / 2.0) ** 2
    velocity = mass_flow_kg_s / (rho_out * area)
    speed_of_sound = math.sqrt(gamma * z_factor * (R_UNIV / mw) * temperature_k)
    return velocity / speed_of_sound


def _select_gas_nps(
    *,
    required_cv: float,
    valve_type: str,
    actual_inlet_q_m3s: float,
    mass_flow_kg_s: float,
    outlet_pressure_pa: float,
    temperature_k: float,
    mw: float,
    gamma: float,
    z_factor: float,
    warnings: list[str],
) -> tuple[float, float, float, float]:
    selected = initial_gas_nps(actual_inlet_q_m3s)
    while True:
        rated = rated_cv(valve_type, selected)  # type: ignore[arg-type]
        utilization = required_cv / rated if rated > 0 else float("inf")
        mach = _mach_outlet(
            mass_flow_kg_s=mass_flow_kg_s,
            outlet_pressure_pa=outlet_pressure_pa,
            temperature_k=temperature_k,
            mw=mw,
            gamma=gamma,
            z_factor=z_factor,
            nps=selected,
        )
        if utilization <= MAX_UTILIZATION and mach < MAX_GAS_OUTLET_MACH:
            if utilization < POOR_CONTROLLABILITY_LIMIT:
                warnings.append("Poor controllability: valve is oversized for service at less than 15% Cv utilization.")
            return selected, rated, utilization, mach
        next_nps = next_standard_nps(selected)
        if next_nps is None:
            warnings.append("Gas valve exceeds utilization or Mach limits at largest available valve table size.")
            return selected, rated, utilization, mach
        selected = next_nps


def calculate_gas_control_valve(
    flow_rate_kg_s: Optional[float] = None,
    flow_rate_norm_m3_hr: Optional[float] = None,
    flow_rate_std_m3_hr: Optional[float] = None,
    inlet_pressure: Optional[float] = None,
    outlet_pressure: Optional[float] = None,
    inlet_pressure_psi: Optional[float] = None,
    outlet_pressure_psi: Optional[float] = None,
    inlet_temperature_c: Optional[float] = None,
    gas_mw: Optional[float] = None,
    gas_gamma: Optional[float] = None,
    gas_z_factor: Optional[float] = None,
    gas_viscosity: Optional[float] = None,
    fluid_name: Optional[str] = None,
    gas_composition_mol: Optional[Dict[str, float]] = None,
    valve_type: Literal["globe", "butterfly", "ball", "gate", "other"] = "globe",
    valve_xt: Optional[float] = None,
    valve_fl: Optional[float] = None,
    valve_fd: Optional[float] = None,
    allow_property_defaults: bool = False,
) -> str:
    """Calculate gas control-valve Cv and select NPS with choked-flow metadata."""

    log: list[str] = []
    warnings: list[str] = []

    try:
        if inlet_temperature_c is None:
            return _json_error("MISSING_REQUIRED_INPUT", "inlet_temperature_c is required.", log)
        temperature_k = inlet_temperature_c + DEG_C_to_K
        if temperature_k <= 0:
            return _json_error("INVALID_TEMPERATURE", "Gas inlet temperature must be above 0 K.", log)

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

        mw, gamma, z_factor, viscosity, gas_info = _resolve_gas_props(
            gas_mw=gas_mw,
            gas_gamma=gas_gamma,
            gas_z_factor=gas_z_factor,
            gas_viscosity=gas_viscosity,
            fluid_name=fluid_name,
            gas_composition_mol=gas_composition_mol,
            inlet_temperature_c=inlet_temperature_c,
            inlet_pressure_pa=p1,
            allow_property_defaults=allow_property_defaults,
            log=log,
            warnings=warnings,
        )
        if None in (mw, gamma, z_factor, viscosity):
            return _json_error(
                "MISSING_GAS_PROPERTIES",
                "gas_mw, gas_gamma, gas_z_factor, and gas_viscosity are required for gas valve sizing.",
                log,
                details=gas_info,
            )
        if gamma <= 1.0 or z_factor <= 0.0 or viscosity <= 0.0 or mw <= 0.0:
            return _json_error("INVALID_GAS_PROPERTIES", "Gas properties must be physically positive and gamma > 1.", log)

        q_norm_m3s, mass_flow = _normal_flow_and_mass(
            flow_rate_kg_s=flow_rate_kg_s,
            flow_rate_norm_m3_hr=flow_rate_norm_m3_hr,
            flow_rate_std_m3_hr=flow_rate_std_m3_hr,
            mw=mw,
        )
        if q_norm_m3s is None or mass_flow is None:
            return _json_error("MISSING_REQUIRED_INPUT", "gas flow rate is required.", log)
        if q_norm_m3s < 0:
            return _json_error("INVALID_FLOW_RATE", "Gas control-valve flow must be non-negative.", log)

        fl = valve_fl if valve_fl is not None else valve_default(valve_type, "FL")
        fd = valve_fd if valve_fd is not None else valve_default(valve_type, "Fd")
        xt = valve_xt if valve_xt is not None else valve_default(valve_type, "xT")
        log.append(f"Using xT={xt}, FL={fl}, Fd={fd} for {valve_type} valve.")

        if q_norm_m3s == 0:
            gas_output = {"Kv": 0.0, "choked": False, "Y": 1.0, "laminar": False}
        else:
            gas_output = fluids.control_valve.size_control_valve_g(
                T=temperature_k,
                MW=mw,
                mu=viscosity,
                gamma=gamma,
                Z=z_factor,
                P1=p1,
                P2=p2,
                Q=q_norm_m3s,
                FL=fl,
                Fd=fd,
                xT=xt,
                full_output=True,
            )
        required_kv = float(gas_output["Kv"])
        required_cv = required_kv * KV_TO_CV
        rho_in = p1 * mw / (z_factor * R_UNIV * temperature_k)
        actual_inlet_q = mass_flow / rho_in if rho_in > 0 else 0.0
        selected_nps, selected_rated_cv, utilization, mach_out = _select_gas_nps(
            required_cv=required_cv,
            valve_type=valve_type,
            actual_inlet_q_m3s=actual_inlet_q,
            mass_flow_kg_s=mass_flow,
            outlet_pressure_pa=p2,
            temperature_k=temperature_k,
            mw=mw,
            gamma=gamma,
            z_factor=z_factor,
            warnings=warnings,
        )
        f_gamma = gamma / 1.4
        critical_pressure_ratio = max(0.0, 1.0 - f_gamma * xt)
        pressure_ratio = p2 / p1

        result = {
            "inputs_resolved": log,
            "warnings": warnings if warnings else None,
            "gas_property_details": gas_info,
            "flow_rate_kg_s": round(mass_flow, 6),
            "flow_rate_norm_m3_hr": round(q_norm_m3s * 3600.0, 4),
            "inlet_pressure_pa": round(p1, 3),
            "outlet_pressure_pa": round(p2, 3),
            "inlet_temperature_c": round(inlet_temperature_c, 3),
            "pressure_drop_pa": round(p1 - p2, 3),
            "valve_type": valve_type,
            "recommended_nps_inches": selected_nps,
            "required_cv": round(required_cv, 4),
            "rated_cv": round(selected_rated_cv, 4),
            "cv_utilization": round(utilization, 4),
            "choked": bool(gas_output.get("choked")),
            "expansion_factor": round(float(gas_output.get("Y", 1.0)), 6),
            "pressure_ratio": round(pressure_ratio, 6),
            "critical_pressure_ratio": round(critical_pressure_ratio, 6),
            "mach_outlet": round(mach_out, 6),
            "xT_used": xt,
            "FL_used": fl,
            "Fd_used": fd,
            "iec_details": gas_output,
        }
        if not result["warnings"]:
            result.pop("warnings", None)
        return json.dumps(result, allow_nan=False)
    except Exception as exc:
        logger.error("Error in calculate_gas_control_valve: %s", exc, exc_info=True)
        return _json_error("GAS_CONTROL_VALVE_FAILED", str(exc), log)
