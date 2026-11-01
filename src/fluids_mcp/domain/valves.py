"""Control-valve size tables and selection helpers."""

from __future__ import annotations

import math
from typing import Literal, Optional

from fluids_mcp.domain.units import FT_TO_M, GPM_TO_M3S, INCH_TO_M
from fluids_mcp.errors import CalculationError

ValveType = Literal["globe", "butterfly", "ball", "gate", "other"]

STANDARD_NPS_INCHES = [0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0]

# Typical full-open rated Cv values for standard trim. These are intentionally
# conservative enough for screening; final selection still belongs to a vendor
# catalog for the selected body/trim.
RATED_CV_TABLES: dict[str, dict[float, float]] = {
    "globe": {
        0.5: 4.0,
        0.75: 7.5,
        1.0: 12.0,
        1.5: 25.0,
        2.0: 50.0,
        2.5: 75.0,
        3.0: 110.0,
        4.0: 190.0,
        6.0: 400.0,
        8.0: 700.0,
        10.0: 1100.0,
        12.0: 1600.0,
    },
    "butterfly": {
        0.5: 8.0,
        0.75: 15.0,
        1.0: 25.0,
        1.5: 60.0,
        2.0: 115.0,
        2.5: 180.0,
        3.0: 260.0,
        4.0: 470.0,
        6.0: 1050.0,
        8.0: 1900.0,
        10.0: 3000.0,
        12.0: 4300.0,
    },
    "ball": {
        0.5: 15.0,
        0.75: 30.0,
        1.0: 55.0,
        1.5: 120.0,
        2.0: 220.0,
        2.5: 340.0,
        3.0: 500.0,
        4.0: 900.0,
        6.0: 2000.0,
        8.0: 3600.0,
        10.0: 5600.0,
        12.0: 8000.0,
    },
    "gate": {
        0.5: 12.0,
        0.75: 25.0,
        1.0: 45.0,
        1.5: 95.0,
        2.0: 180.0,
        2.5: 280.0,
        3.0: 410.0,
        4.0: 740.0,
        6.0: 1650.0,
        8.0: 3000.0,
        10.0: 4700.0,
        12.0: 6700.0,
    },
    "other": {
        0.5: 4.0,
        0.75: 7.5,
        1.0: 12.0,
        1.5: 25.0,
        2.0: 50.0,
        2.5: 75.0,
        3.0: 110.0,
        4.0: 190.0,
        6.0: 400.0,
        8.0: 700.0,
        10.0: 1100.0,
        12.0: 1600.0,
    },
}

VALVE_DEFAULTS = {
    "globe": {"FL": 0.9, "Fd": 0.46, "xT": 0.75},
    "butterfly": {"FL": 0.7, "Fd": 0.25, "xT": 0.35},
    "ball": {"FL": 0.6, "Fd": 0.98, "xT": 0.6},
    "gate": {"FL": 0.9, "Fd": 0.8, "xT": 0.8},
    "other": {"FL": 0.9, "Fd": 1.0, "xT": 0.7},
}

LIQUID_TARGET_VELOCITY_M_S = 8.0 * FT_TO_M
GAS_TARGET_VELOCITY_M_S = 100.0 * FT_TO_M
MAX_GAS_OUTLET_MACH = 0.3


def standard_nps_at_least(min_nps: float) -> float:
    for nps in STANDARD_NPS_INCHES:
        if nps >= min_nps:
            return nps
    raise CalculationError(
        "VALVE_SIZE_OUT_OF_RANGE",
        "Required valve size exceeds available rated Cv table",
        details={"minimum_nps": min_nps, "max_table_nps": STANDARD_NPS_INCHES[-1]},
    )


def next_standard_nps(nps: float) -> Optional[float]:
    for candidate in STANDARD_NPS_INCHES:
        if candidate > nps:
            return candidate
    return None


def rated_cv(valve_type: ValveType, nps: float) -> float:
    table = RATED_CV_TABLES.get(valve_type, RATED_CV_TABLES["other"])
    try:
        return table[nps]
    except KeyError as exc:
        raise CalculationError(
            "VALVE_SIZE_NOT_RATED",
            f"No rated Cv table entry for {valve_type} NPS {nps}",
            details={"valve_type": valve_type, "nps": nps},
        ) from exc


def valve_default(valve_type: ValveType, key: str) -> float:
    return VALVE_DEFAULTS.get(valve_type, VALVE_DEFAULTS["other"])[key]


def initial_liquid_nps(flow_m3s: float, target_velocity_m_s: float = LIQUID_TARGET_VELOCITY_M_S) -> float:
    if flow_m3s <= 0:
        return STANDARD_NPS_INCHES[0]
    required_diameter_m = math.sqrt(4.0 * flow_m3s / (math.pi * target_velocity_m_s))
    return standard_nps_at_least(required_diameter_m / INCH_TO_M)


def initial_gas_nps(actual_volumetric_flow_m3s: float, target_velocity_m_s: float = GAS_TARGET_VELOCITY_M_S) -> float:
    if actual_volumetric_flow_m3s <= 0:
        return STANDARD_NPS_INCHES[0]
    required_diameter_m = math.sqrt(4.0 * actual_volumetric_flow_m3s / (math.pi * target_velocity_m_s))
    return standard_nps_at_least(required_diameter_m / INCH_TO_M)


def nps_inside_diameter_m(nps: float) -> float:
    # Screening approximation: use nominal diameter. It is conservative enough
    # for valve velocity/Mach checks and avoids schedule coupling in valve tables.
    return nps * INCH_TO_M


def liquid_velocity_m_s(flow_m3s: float, nps: float) -> float:
    area = math.pi * (nps_inside_diameter_m(nps) / 2.0) ** 2
    return flow_m3s / area if area > 0 else math.inf


def gpm_to_m3s(value: float) -> float:
    return value * GPM_TO_M3S
