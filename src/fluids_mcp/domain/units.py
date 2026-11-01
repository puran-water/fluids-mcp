"""Explicit unit conversion helpers for public fluids-mcp contracts."""

from __future__ import annotations

from typing import Literal

from fluids_mcp.errors import CalculationError

PressureUnit = Literal["pa", "kpa", "bar", "psi"]
LengthUnit = Literal["m", "mm", "inch", "ft"]
DensityUnit = Literal["kg/m3", "lb/ft3"]
FlowUnit = Literal["m3/s", "m3/hr", "gpm"]
VelocityUnit = Literal["m/s", "ft/s"]

PSI_TO_PA = 6894.757293168
BAR_TO_PA = 100000.0
KPA_TO_PA = 1000.0
INCH_TO_M = 0.0254
FT_TO_M = 0.3048
MM_TO_M = 0.001
LBFT3_TO_KGM3 = 16.01846337396
GPM_TO_M3S = 0.0000630901964
M3HR_TO_M3S = 1.0 / 3600.0
FPS_TO_MS = 0.3048


def _factor(unit: str, factors: dict[str, float], quantity: str) -> float:
    try:
        return factors[unit]
    except KeyError as exc:
        raise CalculationError(
            "UNSUPPORTED_UNIT",
            f"Unsupported {quantity} unit: {unit}",
            details={"unit": unit, "supported_units": sorted(factors)},
        ) from exc


def convert_pressure(value: float, from_unit: PressureUnit, to_unit: PressureUnit = "pa") -> float:
    factors = {"pa": 1.0, "kpa": KPA_TO_PA, "bar": BAR_TO_PA, "psi": PSI_TO_PA}
    pa = value * _factor(from_unit, factors, "pressure")
    return pa / _factor(to_unit, factors, "pressure")


def convert_length(value: float, from_unit: LengthUnit, to_unit: LengthUnit = "m") -> float:
    factors = {"m": 1.0, "mm": MM_TO_M, "inch": INCH_TO_M, "ft": FT_TO_M}
    meters = value * _factor(from_unit, factors, "length")
    return meters / _factor(to_unit, factors, "length")


def convert_density(value: float, from_unit: DensityUnit, to_unit: DensityUnit = "kg/m3") -> float:
    factors = {"kg/m3": 1.0, "lb/ft3": LBFT3_TO_KGM3}
    kg_m3 = value * _factor(from_unit, factors, "density")
    return kg_m3 / _factor(to_unit, factors, "density")


def convert_flow(value: float, from_unit: FlowUnit, to_unit: FlowUnit = "m3/s") -> float:
    factors = {"m3/s": 1.0, "m3/hr": M3HR_TO_M3S, "gpm": GPM_TO_M3S}
    m3_s = value * _factor(from_unit, factors, "flow")
    return m3_s / _factor(to_unit, factors, "flow")


def convert_velocity(value: float, from_unit: VelocityUnit, to_unit: VelocityUnit = "m/s") -> float:
    factors = {"m/s": 1.0, "ft/s": FPS_TO_MS}
    m_s = value * _factor(from_unit, factors, "velocity")
    return m_s / _factor(to_unit, factors, "velocity")


def psi_to_pa(value: float) -> float:
    return convert_pressure(value, "psi", "pa")


def bar_to_pa(value: float) -> float:
    return convert_pressure(value, "bar", "pa")


def inch_to_m(value: float) -> float:
    return convert_length(value, "inch", "m")


def mm_to_m(value: float) -> float:
    return convert_length(value, "mm", "m")


def lbft3_to_kgm3(value: float) -> float:
    return convert_density(value, "lb/ft3", "kg/m3")


def gpm_to_m3s(value: float) -> float:
    return convert_flow(value, "gpm", "m3/s")


def fps_to_ms(value: float) -> float:
    return convert_velocity(value, "ft/s", "m/s")
