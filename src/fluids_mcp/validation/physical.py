"""Physical-sanity guards for fluid calculations."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

from fluids_mcp.errors import CalculationError
from fluids_mcp.validation.input import require_positive


@dataclass(frozen=True)
class ChokedFlowCheck:
    is_choked: bool
    pressure_ratio: float
    critical_pressure_ratio: float


def validate_flow_rate(name: str, value: Optional[float], *, allow_zero: bool = True) -> None:
    require_positive(name, value, allow_zero=allow_zero)


def validate_dimension(name: str, value: Optional[float]) -> None:
    require_positive(name, value, allow_zero=False)


def validate_pressure_drop(pressure_drop: Optional[float], *, allow_zero: bool = False) -> None:
    require_positive("pressure_drop", pressure_drop, allow_zero=allow_zero)


def validate_forward_pressure(
    inlet_pressure: Optional[float],
    outlet_pressure: Optional[float],
    *,
    context: str = "forward flow",
) -> None:
    require_positive("inlet_pressure", inlet_pressure)
    require_positive("outlet_pressure", outlet_pressure)
    if inlet_pressure is None or outlet_pressure is None:
        return
    if outlet_pressure >= inlet_pressure:
        raise CalculationError(
            "INVALID_PRESSURE_DIRECTION",
            f"Outlet pressure must be lower than inlet pressure for {context}",
            details={"inlet_pressure": inlet_pressure, "outlet_pressure": outlet_pressure},
        )


def critical_pressure_ratio(gamma: float) -> float:
    require_positive("gamma", gamma)
    if gamma <= 1.0:
        raise CalculationError(
            "INVALID_GAS_PROPERTY",
            "gas_gamma must be greater than 1.0 for choked-flow checks",
            details={"gamma": gamma},
        )
    return (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))


def check_choked_flow(
    inlet_pressure: float,
    outlet_pressure: float,
    gamma: float,
) -> ChokedFlowCheck:
    validate_forward_pressure(inlet_pressure, outlet_pressure, context="compressible flow")
    ratio = outlet_pressure / inlet_pressure
    critical_ratio = critical_pressure_ratio(gamma)
    return ChokedFlowCheck(
        is_choked=ratio <= critical_ratio,
        pressure_ratio=ratio,
        critical_pressure_ratio=critical_ratio,
    )


def validate_not_choked(
    inlet_pressure: float,
    outlet_pressure: float,
    gamma: float,
    *,
    allow_choked: bool = False,
) -> ChokedFlowCheck:
    check = check_choked_flow(inlet_pressure, outlet_pressure, gamma)
    if check.is_choked and not allow_choked:
        raise CalculationError(
            "CHOKED_FLOW",
            "Choked flow detected",
            details={
                "pressure_ratio": check.pressure_ratio,
                "critical_pressure_ratio": check.critical_pressure_ratio,
            },
        )
    return check


def mach_number(velocity: float, speed_of_sound: float) -> float:
    require_positive("speed_of_sound", speed_of_sound)
    require_positive("velocity", velocity, allow_zero=True)
    return velocity / speed_of_sound


def validate_mach(
    mach: float,
    *,
    max_mach: float = 1.0,
    code: str = "EXCESSIVE_MACH_NUMBER",
) -> None:
    require_positive("mach", mach, allow_zero=True)
    require_positive("max_mach", max_mach)
    if mach >= max_mach:
        raise CalculationError(
            code,
            f"Mach number {mach:.3f} exceeds allowed limit {max_mach:.3f}",
            details={"mach": mach, "max_mach": max_mach},
        )


def validate_reynolds_number(reynolds_number: Optional[float], *, allow_zero: bool = False) -> None:
    if reynolds_number is None:
        return
    if not math.isfinite(reynolds_number):
        raise CalculationError(
            "INVALID_REYNOLDS_NUMBER",
            "Reynolds number must be finite",
            details={"reynolds_number": reynolds_number},
        )
    if reynolds_number < 0 or (reynolds_number == 0 and not allow_zero):
        raise CalculationError(
            "INVALID_REYNOLDS_NUMBER",
            "Reynolds number must be positive",
            details={"reynolds_number": reynolds_number},
        )
