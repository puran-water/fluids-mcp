"""
Fixed orifice (restriction plate) sizing calculations.

This module provides functions for sizing fixed restriction orifices used in
applications like RO permeate backpressure control.

Leverages the fluids library's differential_pressure_meter_solver() for
ISO 5167 orifice calculations with Reader-Harris-Gallagher discharge coefficients.
"""

import json
from fluids_mcp.resource_limits import MAX_PLATE_KIT_SIZES
import logging
import math
from typing import Optional, Literal, List, Dict, Any

# Third-party libraries
import fluids
import fluids.piping
from fluids.flow_meter import (
    differential_pressure_meter_solver,
    C_Reader_Harris_Gallagher,
    dP_orifice,
    ISO_5167_ORIFICE,
)

# Import shared utilities
from fluids_mcp.utils.constants import (
    GPM_to_M3S, INCH_to_M, FT_to_M, PSI_to_PA,
    LBFT3_to_KGM3, CENTIPOISE_to_PAS, DEG_C_to_K,
    DEFAULT_ATMOSPHERIC_PRESSURE
)
# Configure logging
logger = logging.getLogger("fluids-mcp.orifice_restrictor")

# Unit conversion constants
BAR_TO_PA = 100000.0
M3HR_TO_M3S = 1.0 / 3600.0
MM_TO_M = 0.001


# Tullis cavitation-index thresholds for a sharp-edged orifice,
# sigma = (P2 - P_reference) / (P1 - P2). Representative values near beta 0.4; they
# move with beta, geometry and absolute scale, so a returned regime is a screening
# flag, not proof of damage. Ref: Tullis, "Hydraulics of Pipelines", Wiley 1989, ch. 5.
TULLIS_THRESHOLDS = {
    "incipient": 1.10,
    "critical": 0.96,
    "damage": 0.45,
    "choking": 0.27,
}

# Contraction coefficient for a sharp-edged plate, used only to report the
# vena-contracta pressure. Not used in any sizing solve.
VENA_CONTRACTA_CD = 0.62


def classify_cavitation(sigma: float) -> str:
    """Classify a Tullis sigma. See TULLIS_THRESHOLDS for the caveat."""
    if sigma >= TULLIS_THRESHOLDS["incipient"]:
        return "none"
    if sigma >= TULLIS_THRESHOLDS["critical"]:
        return "incipient"
    if sigma >= TULLIS_THRESHOLDS["damage"]:
        return "critical"
    if sigma >= TULLIS_THRESHOLDS["choking"]:
        return "damage"
    return "choking"


def cavitation_assessment(P1: float, P2: float, Q: float, D2: float, rho: float,
                          vapor_pressure_pa: Optional[float] = None,
                          bubble_point_pa: Optional[float] = None) -> Optional[Dict[str, Any]]:
    """
    Liquid cavitation and dissolved-gas breakout screening for a restriction orifice.

    Vapour pressure and bubble point are deliberately separate inputs. A liquid
    carrying dissolved gas (an AnMBR permeate holds dissolved CH4 and CO2) comes out
    of solution at its bubble point, which can be more than an order of magnitude
    above the water vapour pressure. Substituting one for the other silently turns a
    degassing problem into a clean cavitation result, so both indices are reported
    and neither is inferred from the other.

    Returns None when P1/P2 are unusable, so the caller can omit the block rather
    than publish a meaningless index.
    """
    dP = P1 - P2
    if dP <= 0 or P1 <= 0:
        return None

    out: Dict[str, Any] = {"basis": "Tullis sigma = (P2 - P_ref) / (P1 - P2), P2 recovered"}

    p_vc = None
    if D2 > 0 and rho > 0 and Q > 0:
        v_vc = Q / (VENA_CONTRACTA_CD * math.pi * D2 ** 2 / 4)
        p_vc = P1 - 0.5 * rho * v_vc ** 2
        out["vena_contracta_pa_abs"] = round(p_vc, 1)
        out["velocity_vena_contracta_m_s"] = round(v_vc, 2)

    if vapor_pressure_pa is not None:
        sigma = (P2 - vapor_pressure_pa) / dP
        out["sigma_cavitation"] = round(sigma, 3)
        out["cavitation_regime"] = classify_cavitation(sigma)
        out["flashing"] = P2 <= vapor_pressure_pa

    if bubble_point_pa is not None:
        sigma_gas = (P2 - bubble_point_pa) / dP
        out["sigma_gas_breakout"] = round(sigma_gas, 3)
        out["outlet_below_bubble_point"] = P2 < bubble_point_pa
        # Gas leaves solution wherever the LOCAL pressure drops below the bubble
        # point, and the local minimum is the jet - not the recovered outlet. A
        # plate can be degassing inside the vena contracta while its outlet reads
        # comfortably above, which is exactly the case that voids Cd.
        if p_vc is not None:
            out["vena_contracta_below_bubble_point"] = p_vc < bubble_point_pa
            out["gas_breakout"] = (P2 < bubble_point_pa) or (p_vc < bubble_point_pa)
        else:
            out["gas_breakout"] = P2 < bubble_point_pa

    return out or None


def permanent_pressure_loss(D: float, D2: float, dP_tap: float, m: float,
                            rho: float, mu: float, taps: str = 'D',
                            cd: Optional[float] = None) -> Dict[str, float]:
    """
    Convert an ISO 5167 tap-to-tap differential into the permanent (unrecovered) loss.

    The solver sizes against the differential measured at the specified taps, part of
    which is recovered downstream of the vena contracta. For a metering application
    that is the quantity of interest; for a restriction plate it is not - the system
    only ever sees the permanent loss. At beta 0.19 the recovery is about 4% of the
    tap differential, at beta 0.38 about 16%, so sizing a schedule of plates against
    tap differential under-burns by an amount that varies across the schedule.

    Ref: ISO 5167-2:2022; implemented by fluids.flow_meter.dP_orifice.
    """
    C = cd if cd is not None else C_Reader_Harris_Gallagher(
        D=D, Do=D2, rho=rho, mu=mu, m=m, taps=taps
    )
    # dP_orifice depends on P1 and P2 only through their difference, so any positive
    # reference works; use one comfortably above dP so P2 stays positive.
    P1_ref = 10.0 * dP_tap + DEFAULT_ATMOSPHERIC_PRESSURE
    dP_perm = dP_orifice(D=D, Do=D2, P1=P1_ref, P2=P1_ref - dP_tap, C=C)
    return {
        "discharge_coefficient": round(C, 4),
        "permanent_pressure_loss_pa": round(dP_perm, 1),
        "permanent_pressure_loss_psi": round(dP_perm / PSI_to_PA, 2),
        "permanent_pressure_loss_bar": round(dP_perm / BAR_TO_PA, 4),
        "pressure_recovery_fraction": round(1.0 - dP_perm / dP_tap, 4) if dP_tap > 0 else 0.0,
    }


def solve_tap_dp_for_permanent_loss(D: float, target_perm: float, m: float, P1: float,
                                    rho: float, mu: float, meter_type,
                                    taps: str = 'D', cd: Optional[float] = None) -> float:
    """
    Invert the recovery: find the tap differential whose permanent loss hits target.

    Bracketed on [target, 3x target] because permanent loss is always below the tap
    differential and the ratio stays above 1/3 for any beta this tool can size.
    """
    from scipy.optimize import brentq

    def residual(dP_tap: float) -> float:
        kwargs = {
            'D': D, 'm': m, 'P1': P1, 'P2': P1 - dP_tap, 'rho': rho, 'mu': mu,
            'meter_type': meter_type, 'taps': taps, 'epsilon_specified': 1.0,
        }
        if cd is not None:
            kwargs['C_specified'] = cd
        D2 = differential_pressure_meter_solver(**kwargs)
        return permanent_pressure_loss(
            D, D2, dP_tap, m, rho, mu, taps, cd
        )["permanent_pressure_loss_pa"] - target_perm

    hi = min(target_perm * 3.0, 0.95 * P1)
    if residual(hi) < 0:
        raise ValueError(
            f"cannot reach {target_perm/PSI_to_PA:.2f} psi permanent loss at "
            f"P1={P1/PSI_to_PA:.1f} psia; stage the let-down or raise P1"
        )
    return brentq(residual, target_perm * 1.0001, hi, xtol=0.5)


def check_iso5167_validity(D: float, D2: float, dP: float, Re_D: Optional[float] = None,
                           taps: str = 'D', P1: Optional[float] = None,
                           P2: Optional[float] = None, phase: str = "liquid",
                           duty: str = "metering") -> List[str]:
    """
    Check ISO 5167 validity limits and return warnings.

    Args:
        D: Pipe inner diameter in meters
        D2: Orifice diameter in meters
        dP: Pressure drop in Pa
        Re_D: Reynolds number based on pipe diameter (optional)
        taps: Tap configuration ('D', 'flange', 'corner')
        P1: Inlet pressure in Pa (optional, for gas expansibility check)
        P2: Outlet pressure in Pa (optional, for gas expansibility check)
        phase: 'liquid' or 'gas' (for gas expansibility check)
        duty: 'metering' or 'restriction'. A restriction plate is routinely bored
            below the metering minimum in a line below the metering minimum; the
            envelope warnings still fire because Reader-Harris/Gallagher genuinely
            becomes an extrapolation there, but a restriction duty appends what that
            costs so the reader is not left to guess.

    Returns:
        List of warning strings for any limit violations
    """
    warnings = []
    beta = D2 / D if D > 0 else 0

    # Beta ratio limits (ISO 5167-2)
    if not (0.1 <= beta <= 0.75):
        warnings.append(f"Beta ratio {beta:.3f} outside ISO 5167 range [0.1, 0.75]")

    # Pressure drop limit
    if dP > 250000:  # 250 kPa
        warnings.append(f"Pressure drop {dP/1000:.1f} kPa exceeds ISO 5167 limit of 250 kPa")

    # Orifice bore minimum
    if D2 < 0.0125:  # 12.5 mm
        warnings.append(f"Orifice bore {D2*1000:.1f} mm below ISO 5167 minimum of 12.5 mm")

    # Pipe diameter range
    if D < 0.050:  # 50 mm
        warnings.append(f"Pipe diameter {D*1000:.1f} mm below ISO 5167 minimum of 50 mm")
    elif D > 1.0:  # 1000 mm
        warnings.append(f"Pipe diameter {D*1000:.1f} mm exceeds ISO 5167 maximum of 1000 mm")

    # Reynolds number limits (vary by tap type and beta)
    if Re_D is not None:
        if taps in ['D', 'corner']:
            if beta <= 0.56:
                if Re_D < 5000:
                    warnings.append(f"Re_D={Re_D:.0f} below minimum 5000 for beta<0.56")
            else:
                min_Re = 16000 * beta**2
                if Re_D < min_Re:
                    warnings.append(f"Re_D={Re_D:.0f} below minimum {min_Re:.0f} for beta={beta:.2f}")
        elif taps == 'flange':
            min_Re = max(5000, 170000 * beta**2 * D)
            if Re_D < min_Re:
                warnings.append(f"Re_D={Re_D:.0f} below minimum {min_Re:.0f} for flange taps")

    # Gas expansibility validity check (ISO 5167 orifice_expansibility)
    # The ISO 5167 expansibility formula was determined for P2/P1 >= 0.80
    if phase == "gas" and P1 is not None and P2 is not None and P1 > 0:
        pressure_ratio = P2 / P1
        if pressure_ratio < 0.80:
            warnings.append(
                f"Pressure ratio P2/P1={pressure_ratio:.3f} is below 0.80. "
                f"ISO 5167 expansibility correlation may be inaccurate. "
                f"Consider using a different method for high pressure drops."
            )

    # A restriction plate legitimately lives outside the metering envelope. Keep every
    # warning above - the correlation really is extrapolated - but say what that means
    # so the reader does not either dismiss it or treat it as disqualifying.
    if duty == "restriction" and any(
        "below ISO 5167 minimum" in w or "exceeds ISO 5167" in w for w in warnings
    ):
        warnings.append(
            "Restriction duty: the limits above are the ISO 5167 *metering* envelope. "
            "Reader-Harris/Gallagher is extrapolated here (Cd ~0.61 for a sharp-edged "
            "plate at high bore Reynolds number), so treat the bore as +/-5% before "
            "machining tolerance. ASME MFC-14M covers 6-40 mm lines; for a "
            "non-metering restriction, calibrate or provide field-adjustable trim."
        )

    return warnings


def resolve_fluid_properties(
    fluid_name: Optional[str],
    temperature_c: Optional[float],
    pressure_bar: Optional[float],
    fluid_density: Optional[float],
    fluid_viscosity: Optional[float],
    fluid_density_lbft3: Optional[float] = None,
    fluid_viscosity_cp: Optional[float] = None,
    results_log: Optional[List[str]] = None,
    error_log: Optional[List[str]] = None
) -> tuple:
    """
    Resolve fluid density and viscosity from various input options.

    Returns:
        Tuple of (density_kg_m3, viscosity_pa_s, fluid_info_dict)
    """
    if results_log is None:
        results_log = []
    if error_log is None:
        error_log = []

    local_density = None
    local_viscosity = None
    fluid_info = {}

    # Priority 1: Direct SI inputs
    if fluid_density is not None:
        local_density = fluid_density
        results_log.append("Used provided SI fluid_density.")

    if fluid_viscosity is not None:
        local_viscosity = fluid_viscosity
        results_log.append("Used provided SI fluid_viscosity.")

    # Priority 2: Imperial inputs
    if local_density is None and fluid_density_lbft3 is not None:
        local_density = fluid_density_lbft3 * LBFT3_to_KGM3
        results_log.append(f"Converted density from {fluid_density_lbft3} lb/ft³.")

    if local_viscosity is None and fluid_viscosity_cp is not None:
        local_viscosity = fluid_viscosity_cp * CENTIPOISE_to_PAS
        results_log.append(f"Converted viscosity from {fluid_viscosity_cp} cP.")

    # Priority 3: Property lookup via centralized resolver (CoolProp → fluidprop)
    if (local_density is None or local_viscosity is None) and fluid_name is not None and temperature_c is not None:
        try:
            from fluids_mcp.utils.resolve_properties import resolve_liquid_properties
            p_bar = pressure_bar if pressure_bar is not None else 1.01325
            props = resolve_liquid_properties(fluid_name, temperature_c, p_bar)
            if props is not None:
                if local_density is None:
                    local_density = props.density
                if local_viscosity is None:
                    local_viscosity = props.viscosity
                fluid_info = {
                    "name_used": fluid_name,
                    "temperature_c": temperature_c,
                    "pressure_bar": p_bar,
                    "density_kg_m3": local_density,
                    "viscosity_pa_s": local_viscosity
                }
                results_log.append(f"Properties via {props.source} for {fluid_name} at {temperature_c}°C, {p_bar} bar.")
                results_log.extend(props.log)
            else:
                error_log.append(f"Fluid property lookup failed for {fluid_name}: all backends failed")
        except Exception as e:
            error_log.append(f"Fluid property lookup failed: {e}")

    return local_density, local_viscosity, fluid_info


def resolve_pipe_diameter(
    pipe_diameter: Optional[float],
    pipe_diameter_in: Optional[float],
    nominal_size_in: Optional[float],
    schedule: str = "40",
    results_log: Optional[List[str]] = None,
    error_log: Optional[List[str]] = None
) -> Optional[float]:
    """
    Resolve pipe inner diameter from various input options.

    Returns:
        Pipe inner diameter in meters, or None if resolution fails
    """
    if results_log is None:
        results_log = []
    if error_log is None:
        error_log = []

    local_D = None

    # Priority 1: Direct SI input
    if pipe_diameter is not None:
        local_D = pipe_diameter
        results_log.append(f"Used provided pipe_diameter: {local_D*1000:.2f} mm")

    # Priority 2: Imperial input
    elif pipe_diameter_in is not None:
        local_D = pipe_diameter_in * INCH_to_M
        results_log.append(f"Converted pipe_diameter from {pipe_diameter_in} inches: {local_D*1000:.2f} mm")

    # Priority 3: NPS lookup
    elif nominal_size_in is not None:
        try:
            NPS, Di, Do, t = fluids.piping.nearest_pipe(NPS=nominal_size_in, schedule=schedule)
            local_D = Di
            results_log.append(f"Looked up NPS {nominal_size_in}\" Sch {schedule}: ID = {local_D*1000:.2f} mm")
        except Exception as e:
            error_log.append(f"Pipe dimension lookup failed: {e}")

    return local_D


def resolve_orifice_diameter(
    orifice_diameter: Optional[float],
    orifice_diameter_in: Optional[float],
    orifice_diameter_mm: Optional[float],
    results_log: Optional[List[str]] = None
) -> Optional[float]:
    """
    Resolve orifice diameter from various input options.

    Returns:
        Orifice diameter in meters, or None if not provided
    """
    if results_log is None:
        results_log = []

    local_D2 = None

    if orifice_diameter is not None:
        local_D2 = orifice_diameter
        results_log.append(f"Used provided orifice_diameter: {local_D2*1000:.2f} mm")
    elif orifice_diameter_in is not None:
        local_D2 = orifice_diameter_in * INCH_to_M
        results_log.append(f"Converted orifice_diameter from {orifice_diameter_in} inches: {local_D2*1000:.2f} mm")
    elif orifice_diameter_mm is not None:
        local_D2 = orifice_diameter_mm * MM_TO_M
        results_log.append(f"Converted orifice_diameter from {orifice_diameter_mm} mm")

    return local_D2


def resolve_pressure_drop(
    pressure_drop: Optional[float],
    pressure_drop_psi: Optional[float],
    pressure_drop_bar: Optional[float],
    results_log: Optional[List[str]] = None
) -> Optional[float]:
    """
    Resolve pressure drop from various input options.

    Returns:
        Pressure drop in Pa, or None if not provided
    """
    if results_log is None:
        results_log = []

    local_dP = None

    if pressure_drop is not None:
        local_dP = pressure_drop
        results_log.append(f"Used provided pressure_drop: {local_dP/1000:.2f} kPa")
    elif pressure_drop_psi is not None:
        local_dP = pressure_drop_psi * PSI_to_PA
        results_log.append(f"Converted pressure_drop from {pressure_drop_psi} psi: {local_dP/1000:.2f} kPa")
    elif pressure_drop_bar is not None:
        local_dP = pressure_drop_bar * BAR_TO_PA
        results_log.append(f"Converted pressure_drop from {pressure_drop_bar} bar: {local_dP/1000:.2f} kPa")

    return local_dP


def resolve_flow_rate(
    flow_rate: Optional[float],
    flow_rate_gpm: Optional[float],
    flow_rate_m3_hr: Optional[float],
    results_log: Optional[List[str]] = None
) -> Optional[float]:
    """
    Resolve volumetric flow rate from various input options.

    Returns:
        Flow rate in m³/s, or None if not provided
    """
    if results_log is None:
        results_log = []

    local_Q = None

    if flow_rate is not None:
        local_Q = flow_rate
        results_log.append(f"Used provided flow_rate: {local_Q*3600:.2f} m³/h")
    elif flow_rate_gpm is not None:
        local_Q = flow_rate_gpm * GPM_to_M3S
        results_log.append(f"Converted flow_rate from {flow_rate_gpm} GPM: {local_Q*3600:.2f} m³/h")
    elif flow_rate_m3_hr is not None:
        local_Q = flow_rate_m3_hr * M3HR_TO_M3S
        results_log.append(f"Converted flow_rate from {flow_rate_m3_hr} m³/h")

    return local_Q


def format_output_diameter(D2: float, size_units: str) -> Dict[str, Any]:
    """Format orifice diameter for output in requested units."""
    result = {
        "orifice_diameter_m": round(D2, 6),
        "orifice_diameter_mm": round(D2 * 1000, 2),
        "orifice_diameter_in": round(D2 / INCH_to_M, 3),
    }

    if size_units == "mm":
        result["orifice_diameter"] = result["orifice_diameter_mm"]
        result["orifice_diameter_units"] = "mm"
    elif size_units == "inch":
        result["orifice_diameter"] = result["orifice_diameter_in"]
        result["orifice_diameter_units"] = "inches"
    else:  # meters
        result["orifice_diameter"] = result["orifice_diameter_m"]
        result["orifice_diameter_units"] = "m"

    return result


def calculate_orifice_sizing(
    # Flow inputs
    flow_rate: Optional[float] = None,
    flow_rate_gpm: Optional[float] = None,
    flow_rate_m3_hr: Optional[float] = None,

    # Pressure inputs
    pressure_drop: Optional[float] = None,
    pressure_drop_psi: Optional[float] = None,
    pressure_drop_bar: Optional[float] = None,
    inlet_pressure: Optional[float] = None,
    inlet_pressure_psi: Optional[float] = None,
    inlet_pressure_bar: Optional[float] = None,

    # Pipe geometry
    pipe_diameter: Optional[float] = None,
    pipe_diameter_in: Optional[float] = None,
    nominal_size_in: Optional[float] = None,
    schedule: str = "40",

    # Fluid properties
    fluid_name: Optional[str] = None,
    temperature_c: Optional[float] = None,
    pressure_bar: Optional[float] = None,
    fluid_density: Optional[float] = None,
    fluid_viscosity: Optional[float] = None,
    fluid_density_lbft3: Optional[float] = None,
    fluid_viscosity_cp: Optional[float] = None,

    # Gas properties (for compressible flow)
    gas_gamma: Optional[float] = None,
    gas_mw: Optional[float] = None,
    gas_z_factor: Optional[float] = None,

    # Orifice parameters
    cd: Optional[float] = None,
    meter_type: str = None,  # Will default to ISO_5167_ORIFICE
    taps: str = 'D',

    # Restriction-duty parameters
    dp_basis: Literal["tap", "permanent"] = "tap",
    duty: Literal["metering", "restriction"] = "metering",
    vapor_pressure_pa: Optional[float] = None,
    bubble_point_pa: Optional[float] = None,

    # Phase selection
    phase: Literal["liquid", "gas"] = "liquid",

    # Output options
    size_units: Literal["m", "inch", "mm"] = "mm",
) -> str:
    """
    Size an orifice restrictor for target flow and pressure drop.

    Given Q and ΔP, calculates the required orifice diameter using ISO 5167
    correlations with Reader-Harris-Gallagher discharge coefficients.

    Args:
        flow_rate: Volumetric flow rate in m³/s
        flow_rate_gpm: Volumetric flow rate in US GPM
        flow_rate_m3_hr: Volumetric flow rate in m³/hr
        pressure_drop: Target pressure drop in Pa
        pressure_drop_psi: Target pressure drop in psi
        pressure_drop_bar: Target pressure drop in bar
        inlet_pressure: Inlet pressure in Pa (absolute) - optional for liquid
        inlet_pressure_psi: Inlet pressure in psi (absolute)
        inlet_pressure_bar: Inlet pressure in bar (absolute)
        pipe_diameter: Pipe inner diameter in m
        pipe_diameter_in: Pipe inner diameter in inches
        nominal_size_in: Nominal pipe size in inches (for lookup)
        schedule: Pipe schedule (default "40")
        fluid_name: Fluid name for property lookup
        temperature_c: Temperature in Celsius
        pressure_bar: Pressure context for property lookup
        fluid_density: Fluid density in kg/m³
        fluid_viscosity: Fluid viscosity in Pa·s
        gas_gamma: Specific heat ratio Cp/Cv (for gas)
        gas_mw: Molecular weight in kg/kmol (for gas)
        gas_z_factor: Compressibility factor (for gas)
        cd: Override discharge coefficient
        meter_type: Meter type string (default ISO 5167 orifice)
        taps: Tap configuration ('D', 'flange', 'corner')
        phase: 'liquid' or 'gas'
        size_units: Output units ('m', 'inch', 'mm')

    Returns:
        JSON string with orifice diameter and calculation details
    """
    results_log = []
    error_log = []

    try:
        # Use constant instead of string for meter type
        local_meter_type = meter_type if meter_type else ISO_5167_ORIFICE

        # Resolve flow rate
        local_Q = resolve_flow_rate(flow_rate, flow_rate_gpm, flow_rate_m3_hr, results_log)
        if local_Q is None:
            error_log.append("Missing required input: flow_rate")

        # Resolve pressure drop
        local_dP = resolve_pressure_drop(pressure_drop, pressure_drop_psi, pressure_drop_bar, results_log)
        if local_dP is None:
            error_log.append("Missing required input: pressure_drop")

        # Resolve pipe diameter
        local_D = resolve_pipe_diameter(pipe_diameter, pipe_diameter_in, nominal_size_in,
                                        schedule, results_log, error_log)
        if local_D is None:
            error_log.append("Missing required input: pipe_diameter or nominal_size_in")

        # For gas phase, derive pressure_bar from inlet_pressure if not explicitly set
        # This ensures correct density calculation at elevated pressures
        local_pressure_bar = pressure_bar
        if local_pressure_bar is None:
            if inlet_pressure is not None:
                local_pressure_bar = inlet_pressure / 100000.0  # Pa to bar
                results_log.append(f"Derived pressure_bar={local_pressure_bar:.3f} from inlet_pressure for property lookup")
            elif inlet_pressure_psi is not None:
                local_pressure_bar = (inlet_pressure_psi * PSI_to_PA) / 100000.0
                results_log.append(f"Derived pressure_bar={local_pressure_bar:.3f} from inlet_pressure_psi for property lookup")
            elif inlet_pressure_bar is not None:
                local_pressure_bar = inlet_pressure_bar
                results_log.append(f"Using inlet_pressure_bar={local_pressure_bar:.3f} for property lookup")

        # Resolve fluid properties
        local_rho, local_mu, fluid_info = resolve_fluid_properties(
            fluid_name, temperature_c, local_pressure_bar,
            fluid_density, fluid_viscosity,
            fluid_density_lbft3, fluid_viscosity_cp,
            results_log, error_log
        )

        if local_rho is None:
            error_log.append("Missing required input: fluid_density or (fluid_name + temperature_c)")
        if local_mu is None:
            error_log.append("Missing required input: fluid_viscosity or (fluid_name + temperature_c)")

        # Resolve inlet pressure
        local_P1 = None
        if inlet_pressure is not None:
            local_P1 = inlet_pressure
            results_log.append(f"Used provided inlet_pressure: {local_P1/1000:.1f} kPa")
        elif inlet_pressure_psi is not None:
            local_P1 = inlet_pressure_psi * PSI_to_PA
            results_log.append(f"Converted inlet_pressure from {inlet_pressure_psi} psi")
        elif inlet_pressure_bar is not None:
            local_P1 = inlet_pressure_bar * BAR_TO_PA
            results_log.append(f"Converted inlet_pressure from {inlet_pressure_bar} bar")

        # For liquid phase, we can use a reference pressure if not provided
        p1_was_supplied = local_P1 is not None
        if local_P1 is None and phase == "liquid" and local_dP is not None:
            # Incompressible, so the bore does not depend on absolute pressure - but
            # the cavitation index does, and is suppressed below rather than computed
            # against an invented anchor.
            local_P1 = 200000.0 + local_dP  # 2 bar above outlet
            results_log.append(
                f"Using reference inlet pressure: {local_P1/1000:.1f} kPa (liquid shortcut). "
                f"Bore is unaffected; cavitation assessment suppressed."
            )

        # Check for critical errors
        if error_log:
            return json.dumps({"errors": error_log, "log": results_log})

        # Calculate outlet pressure from inlet and pressure drop
        local_P2 = local_P1 - local_dP

        if local_P2 <= 0:
            error_log.append(f"Calculated outlet pressure {local_P2:.0f} Pa is not positive. Check inlet pressure and pressure drop.")
            return json.dumps({"errors": error_log, "log": results_log})

        # Calculate mass flow rate
        local_m = local_Q * local_rho
        results_log.append(f"Calculated mass flow rate: {local_m:.4f} kg/s")

        # For liquid, set expansibility to 1.0 (incompressible)
        # For gas, need k (gamma) for expansibility calculation
        epsilon_specified = None
        local_k = None

        if phase == "liquid":
            epsilon_specified = 1.0  # Incompressible
            results_log.append("Using epsilon=1.0 for incompressible liquid")
        else:
            # Gas phase - need isentropic exponent
            local_k = gas_gamma
            if local_k is None:
                error_log.append("Missing required input for gas: gas_gamma (specific heat ratio)")
                return json.dumps({"errors": error_log, "log": results_log})
            results_log.append(f"Using gas_gamma (k) = {local_k}")

            # Check for choked flow
            pressure_ratio = local_P2 / local_P1
            critical_ratio = (2 / (local_k + 1)) ** (local_k / (local_k - 1))
            if pressure_ratio < critical_ratio:
                error_log.append(f"Choked flow condition: P2/P1 = {pressure_ratio:.3f} < critical ratio {critical_ratio:.3f}")
                return json.dumps({"errors": error_log, "log": results_log})

        # Build solver arguments
        solver_kwargs = {
            'D': local_D,
            'm': local_m,
            'P1': local_P1,
            'P2': local_P2,
            'rho': local_rho,
            'mu': local_mu,
            'meter_type': local_meter_type,
            'taps': taps,
        }

        if epsilon_specified is not None:
            solver_kwargs['epsilon_specified'] = epsilon_specified
        if local_k is not None:
            solver_kwargs['k'] = local_k
        if cd is not None:
            solver_kwargs['C_specified'] = cd
            results_log.append(f"Using override Cd = {cd}")

        # Call the solver
        try:
            if dp_basis == "permanent":
                if phase != "liquid":
                    error_log.append("dp_basis='permanent' is only defined for liquid phase")
                    return json.dumps({"errors": error_log, "log": results_log})
                # The requested dP is what the system must see. Invert the recovery to
                # find the tap differential that delivers it, then size against that.
                dP_tap = solve_tap_dp_for_permanent_loss(
                    local_D, local_dP, local_m, local_P1, local_rho, local_mu,
                    local_meter_type, taps, cd
                )
                results_log.append(
                    f"dp_basis='permanent': target {local_dP/PSI_to_PA:.2f} psi permanent "
                    f"requires {dP_tap/PSI_to_PA:.2f} psi tap differential"
                )
                solver_kwargs['P2'] = local_P1 - dP_tap
                local_P2 = local_P1 - dP_tap
            D2 = differential_pressure_meter_solver(**solver_kwargs)
            results_log.append(f"Solver converged: D2 = {D2*1000:.2f} mm")
        except Exception as solver_e:
            error_log.append(f"Orifice solver failed: {solver_e}")
            return json.dumps({"errors": error_log, "log": results_log})

        dP_tap_actual = local_P1 - local_P2

        # Calculate derived values
        beta = D2 / local_D
        A_orifice = math.pi * D2**2 / 4
        V_orifice = local_Q / A_orifice
        Re_D = local_rho * (local_Q / (math.pi * local_D**2 / 4)) * local_D / local_mu
        Re_orifice = local_rho * V_orifice * D2 / local_mu

        # Check ISO 5167 validity (including gas expansibility check)
        validity_warnings = check_iso5167_validity(
            local_D, D2, dP_tap_actual, Re_D, taps,
            P1=local_P1, P2=local_P2, phase=phase, duty=duty
        )
        if validity_warnings:
            error_log.extend(validity_warnings)

        # Permanent (unrecovered) loss - what the system actually sees. Always
        # reported so the tap-vs-permanent distinction cannot be missed.
        perm_block = {}
        if phase == "liquid":
            perm_block = permanent_pressure_loss(
                local_D, D2, dP_tap_actual, local_m, local_rho, local_mu, taps, cd
            )

        # Tullis sigma is referenced to the *recovered* downstream pressure - the
        # pressure in the pipe well past the vena contracta - not the low-pressure
        # tap. Using the tap value understates sigma and can cross a regime boundary.
        cav_block = None
        if phase == "liquid" and p1_was_supplied:
            P2_recovered = local_P1 - perm_block["permanent_pressure_loss_pa"]
            cav_block = cavitation_assessment(
                local_P1, P2_recovered, local_Q, D2, local_rho,
                vapor_pressure_pa=vapor_pressure_pa,
                bubble_point_pa=bubble_point_pa,
            )

        if cav_block and cav_block.get("cavitation_regime") in ("damage", "choking"):
            error_log.append(
                f"Cavitation regime '{cav_block['cavitation_regime']}' at sigma "
                f"{cav_block['sigma_cavitation']}: the sqrt(dP) orifice law this sizing "
                f"assumes is not valid there. Stage the let-down or raise back pressure."
            )
        if cav_block and cav_block.get("gas_breakout"):
            where = ("the vena contracta" if cav_block.get("vena_contracta_below_bubble_point")
                     and not cav_block.get("outlet_below_bubble_point") else "the outlet")
            error_log.append(
                f"Pressure at {where} is below the liquid's bubble point: dissolved gas "
                f"breaks out across the plate, so the discharge coefficient no longer applies."
            )

        # Format output
        output_dims = format_output_diameter(D2, size_units)

        result = {
            "mode": "sizing",
            "phase": phase,
            "duty": duty,
            "dp_basis": dp_basis,
            "inputs_resolved": results_log,
            "warnings": error_log if error_log else None,
            "fluid_details": fluid_info if fluid_info else None,

            # Core results
            **output_dims,
            "pipe_diameter_mm": round(local_D * 1000, 2),
            "beta_ratio": round(beta, 4),

            # Flow parameters
            "flow_rate_m3_hr": round(local_Q * 3600, 3),
            "flow_rate_gpm": round(local_Q / GPM_to_M3S, 2),
            "mass_flow_kg_s": round(local_m, 4),

            # Pressure parameters. pressure_drop_* is the tap-to-tap differential at
            # `taps`; permanent_pressure_loss_* is the unrecovered loss.
            "pressure_drop_bar": round(dP_tap_actual / BAR_TO_PA, 3),
            "pressure_drop_psi": round(dP_tap_actual / PSI_to_PA, 2),
            "inlet_pressure_bar": round(local_P1 / BAR_TO_PA, 3),
            "outlet_pressure_bar": round(local_P2 / BAR_TO_PA, 3),
            **perm_block,
            "cavitation": cav_block,

            # Orifice characteristics
            "velocity_orifice_m_s": round(V_orifice, 2),
            "reynolds_pipe": round(Re_D, 0),
            "reynolds_orifice": round(Re_orifice, 0),
            "meter_type": str(local_meter_type),
            "taps": taps,
        }

        # Remove null values
        result = {k: v for k, v in result.items() if v is not None}

        return json.dumps(result)

    except Exception as e:
        logger.error(f"Error in calculate_orifice_sizing: {e}", exc_info=True)
        return json.dumps({
            "error": f"Calculation error: {str(e)}",
            "log": results_log,
            "errors": error_log
        })


def calculate_orifice_analysis_dp(
    # Flow inputs
    flow_rate: Optional[float] = None,
    flow_rate_gpm: Optional[float] = None,
    flow_rate_m3_hr: Optional[float] = None,

    # Orifice geometry
    orifice_diameter: Optional[float] = None,
    orifice_diameter_in: Optional[float] = None,
    orifice_diameter_mm: Optional[float] = None,

    # Pipe geometry
    pipe_diameter: Optional[float] = None,
    pipe_diameter_in: Optional[float] = None,
    nominal_size_in: Optional[float] = None,
    schedule: str = "40",

    # Pressure inputs (for gas expansibility)
    inlet_pressure: Optional[float] = None,
    inlet_pressure_psi: Optional[float] = None,
    inlet_pressure_bar: Optional[float] = None,

    # Fluid properties
    fluid_name: Optional[str] = None,
    temperature_c: Optional[float] = None,
    pressure_bar: Optional[float] = None,
    fluid_density: Optional[float] = None,
    fluid_viscosity: Optional[float] = None,
    fluid_density_lbft3: Optional[float] = None,
    fluid_viscosity_cp: Optional[float] = None,

    # Gas properties
    gas_gamma: Optional[float] = None,

    # Orifice parameters
    cd: Optional[float] = None,
    meter_type: str = None,
    taps: str = 'D',

    # Restriction-duty parameters
    duty: Literal["metering", "restriction"] = "metering",
    vapor_pressure_pa: Optional[float] = None,
    bubble_point_pa: Optional[float] = None,

    # Phase selection
    phase: Literal["liquid", "gas"] = "liquid",

    # Output options
    size_units: Literal["m", "inch", "mm"] = "mm",
) -> str:
    """
    Calculate pressure drop for a given orifice at specified flow.

    Given Q and d, calculates the resulting ΔP using ISO 5167 correlations.

    Returns:
        JSON string with pressure drop and calculation details
    """
    results_log = []
    error_log = []

    try:
        local_meter_type = meter_type if meter_type else ISO_5167_ORIFICE

        # Resolve inputs
        local_Q = resolve_flow_rate(flow_rate, flow_rate_gpm, flow_rate_m3_hr, results_log)
        if local_Q is None:
            error_log.append("Missing required input: flow_rate")

        local_D2 = resolve_orifice_diameter(orifice_diameter, orifice_diameter_in,
                                            orifice_diameter_mm, results_log)
        if local_D2 is None:
            error_log.append("Missing required input: orifice_diameter")

        local_D = resolve_pipe_diameter(pipe_diameter, pipe_diameter_in, nominal_size_in,
                                        schedule, results_log, error_log)
        if local_D is None:
            error_log.append("Missing required input: pipe_diameter or nominal_size_in")

        # For gas phase, derive pressure_bar from inlet_pressure if not explicitly set
        # This ensures correct density calculation at elevated pressures
        local_pressure_bar = pressure_bar
        if phase == "gas" and local_pressure_bar is None:
            if inlet_pressure is not None:
                local_pressure_bar = inlet_pressure / 100000.0  # Pa to bar
                results_log.append(f"Derived pressure_bar={local_pressure_bar:.3f} from inlet_pressure for property lookup")
            elif inlet_pressure_psi is not None:
                local_pressure_bar = (inlet_pressure_psi * PSI_to_PA) / 100000.0
                results_log.append(f"Derived pressure_bar={local_pressure_bar:.3f} from inlet_pressure_psi for property lookup")
            elif inlet_pressure_bar is not None:
                local_pressure_bar = inlet_pressure_bar
                results_log.append(f"Using inlet_pressure_bar={local_pressure_bar:.3f} for property lookup")

        local_rho, local_mu, fluid_info = resolve_fluid_properties(
            fluid_name, temperature_c, local_pressure_bar,
            fluid_density, fluid_viscosity,
            fluid_density_lbft3, fluid_viscosity_cp,
            results_log, error_log
        )

        if local_rho is None:
            error_log.append("Missing required input: fluid_density")
        if local_mu is None:
            error_log.append("Missing required input: fluid_viscosity")

        if error_log:
            return json.dumps({"errors": error_log, "log": results_log})

        # Calculate mass flow rate
        local_m = local_Q * local_rho

        # Calculate beta
        beta = local_D2 / local_D

        # For liquid, use simplified calculation with epsilon=1.0
        # For gas, need inlet pressure and gamma for expansibility
        supplied_P1 = None
        if phase == "liquid":
            epsilon = 1.0
            # For incompressible flow the differential is independent of absolute
            # pressure, so a reference anchor gives the right dP. It does not give a
            # meaningful cavitation index, which needs the real P1 - so prefer what
            # the caller supplied and say so when falling back.
            supplied_P1 = None
            if inlet_pressure is not None:
                supplied_P1 = inlet_pressure
            elif inlet_pressure_psi is not None:
                supplied_P1 = inlet_pressure_psi * PSI_to_PA
            elif inlet_pressure_bar is not None:
                supplied_P1 = inlet_pressure_bar * BAR_TO_PA

            if supplied_P1 is not None:
                local_P1 = supplied_P1
                results_log.append(
                    f"Using supplied inlet pressure {local_P1/PSI_to_PA:.2f} psia"
                )
            else:
                local_P1 = 1000000.0
                results_log.append(
                    "No inlet pressure supplied; using a 10 bar reference anchor. "
                    "dP is independent of absolute pressure for incompressible liquid, "
                    "but the cavitation assessment is suppressed."
                )
        else:
            if gas_gamma is None:
                error_log.append("Missing required input for gas: gas_gamma")
                return json.dumps({"errors": error_log, "log": results_log})

            # Resolve inlet pressure
            local_P1 = None
            if inlet_pressure is not None:
                local_P1 = inlet_pressure
            elif inlet_pressure_psi is not None:
                local_P1 = inlet_pressure_psi * PSI_to_PA
            elif inlet_pressure_bar is not None:
                local_P1 = inlet_pressure_bar * BAR_TO_PA

            if local_P1 is None:
                error_log.append("Missing required input for gas: inlet_pressure")
                return json.dumps({"errors": error_log, "log": results_log})

        # Iteratively solve for P2 (pressure drop)
        # Using differential_pressure_meter_solver with P2=None
        solver_kwargs = {
            'D': local_D,
            'D2': local_D2,
            'm': local_m,
            'P1': local_P1,
            'rho': local_rho,
            'mu': local_mu,
            'meter_type': local_meter_type,
            'taps': taps,
        }

        if phase == "liquid":
            solver_kwargs['epsilon_specified'] = 1.0
        else:
            solver_kwargs['k'] = gas_gamma

        if cd is not None:
            solver_kwargs['C_specified'] = cd

        try:
            local_P2 = differential_pressure_meter_solver(**solver_kwargs)
            local_dP = local_P1 - local_P2
            results_log.append(f"Solver converged: ΔP = {local_dP/1000:.2f} kPa")
        except Exception as solver_e:
            error_log.append(f"Orifice solver failed: {solver_e}")
            return json.dumps({"errors": error_log, "log": results_log})

        # Calculate derived values
        A_orifice = math.pi * local_D2**2 / 4
        V_orifice = local_Q / A_orifice
        Re_D = local_rho * (local_Q / (math.pi * local_D**2 / 4)) * local_D / local_mu

        # Check validity (including gas expansibility check)
        validity_warnings = check_iso5167_validity(
            local_D, local_D2, local_dP, Re_D, taps,
            P1=local_P1, P2=local_P2, phase=phase, duty=duty
        )
        if validity_warnings:
            error_log.extend(validity_warnings)

        # Permanent (unrecovered) loss - what the system sees, as opposed to the
        # tap-to-tap differential reported as pressure_drop_*.
        perm_block = {}
        if phase == "liquid":
            perm_block = permanent_pressure_loss(
                local_D, local_D2, local_dP, local_m, local_rho, local_mu, taps, cd
            )

        # Only meaningful against a real inlet pressure; suppressed when the
        # incompressible reference anchor was used instead.
        cav_block = None
        if phase == "liquid" and supplied_P1 is not None:
            # Referenced to the recovered downstream pressure, not the tap.
            cav_block = cavitation_assessment(
                local_P1, local_P1 - perm_block["permanent_pressure_loss_pa"],
                local_Q, local_D2, local_rho,
                vapor_pressure_pa=vapor_pressure_pa,
                bubble_point_pa=bubble_point_pa,
            )
            if cav_block and cav_block.get("cavitation_regime") in ("damage", "choking"):
                error_log.append(
                    f"Cavitation regime '{cav_block['cavitation_regime']}' at sigma "
                    f"{cav_block['sigma_cavitation']}: the sqrt(dP) orifice law is not "
                    f"valid there, so this predicted dP is unreliable."
                )
            if cav_block and cav_block.get("gas_breakout"):
                error_log.append(
                    "Pressure falls below the liquid's bubble point across this plate: "
                    "dissolved gas breaks out and the discharge coefficient no longer applies."
                )

        # Format output
        output_dims = format_output_diameter(local_D2, size_units)

        result = {
            "mode": "analysis_dp",
            "phase": phase,
            "duty": duty,
            "inputs_resolved": results_log,
            "warnings": error_log if error_log else None,
            "fluid_details": fluid_info if fluid_info else None,

            # Calculated pressure drop. pressure_drop_* is the tap-to-tap differential
            # at `taps`; permanent_pressure_loss_* is the unrecovered loss.
            "pressure_drop_pa": round(local_dP, 1),
            "pressure_drop_bar": round(local_dP / BAR_TO_PA, 4),
            "pressure_drop_psi": round(local_dP / PSI_to_PA, 3),
            # Permanent head, not the tap differential: this is the head the
            # system actually loses across the plate.
            "head_loss_m": round(
                perm_block.get("permanent_pressure_loss_pa", local_dP)
                / (local_rho * 9.80665), 3),
            "head_loss_basis": "permanent" if perm_block else "tap",
            **perm_block,
            "cavitation": cav_block,

            # Input geometry
            **output_dims,
            "pipe_diameter_mm": round(local_D * 1000, 2),
            "beta_ratio": round(beta, 4),

            # Flow parameters
            "flow_rate_m3_hr": round(local_Q * 3600, 3),
            "flow_rate_gpm": round(local_Q / GPM_to_M3S, 2),
            "velocity_orifice_m_s": round(V_orifice, 2),
            "reynolds_pipe": round(Re_D, 0),
        }

        result = {k: v for k, v in result.items() if v is not None}
        return json.dumps(result)

    except Exception as e:
        logger.error(f"Error in calculate_orifice_analysis_dp: {e}", exc_info=True)
        return json.dumps({
            "error": f"Calculation error: {str(e)}",
            "log": results_log,
            "errors": error_log
        })


def calculate_orifice_analysis_flow(
    # Pressure inputs
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

    # Fluid properties
    fluid_name: Optional[str] = None,
    temperature_c: Optional[float] = None,
    pressure_bar: Optional[float] = None,
    fluid_density: Optional[float] = None,
    fluid_viscosity: Optional[float] = None,
    fluid_density_lbft3: Optional[float] = None,
    fluid_viscosity_cp: Optional[float] = None,

    # Gas properties
    gas_gamma: Optional[float] = None,

    # Orifice parameters
    cd: Optional[float] = None,
    meter_type: str = None,
    taps: str = 'D',

    # Restriction-duty parameters
    duty: Literal["metering", "restriction"] = "metering",
    vapor_pressure_pa: Optional[float] = None,
    bubble_point_pa: Optional[float] = None,

    # Phase selection
    phase: Literal["liquid", "gas"] = "liquid",

    # Output options
    size_units: Literal["m", "inch", "mm"] = "mm",
) -> str:
    """
    Calculate flow rate for a given orifice and measured pressure drop.

    Given d and ΔP, calculates the resulting Q using ISO 5167 correlations.

    Returns:
        JSON string with flow rate and calculation details
    """
    results_log = []
    error_log = []

    try:
        local_meter_type = meter_type if meter_type else ISO_5167_ORIFICE

        # Resolve inputs
        local_dP = resolve_pressure_drop(pressure_drop, pressure_drop_psi, pressure_drop_bar, results_log)
        if local_dP is None:
            error_log.append("Missing required input: pressure_drop")

        local_D2 = resolve_orifice_diameter(orifice_diameter, orifice_diameter_in,
                                            orifice_diameter_mm, results_log)
        if local_D2 is None:
            error_log.append("Missing required input: orifice_diameter")

        local_D = resolve_pipe_diameter(pipe_diameter, pipe_diameter_in, nominal_size_in,
                                        schedule, results_log, error_log)
        if local_D is None:
            error_log.append("Missing required input: pipe_diameter")

        # Resolve inlet pressure first (needed for gas property lookup)
        local_P1 = None
        if inlet_pressure is not None:
            local_P1 = inlet_pressure
        elif inlet_pressure_psi is not None:
            local_P1 = inlet_pressure_psi * PSI_to_PA
        elif inlet_pressure_bar is not None:
            local_P1 = inlet_pressure_bar * BAR_TO_PA

        if local_P1 is None and phase == "liquid" and local_dP is not None:
            local_P1 = 200000.0 + local_dP  # Reference pressure for liquid
            results_log.append("Using reference inlet pressure for liquid calculation")

        if local_P1 is None and phase == "gas":
            error_log.append("Missing required input for gas: inlet_pressure")

        # For gas phase, derive pressure_bar from inlet_pressure if not explicitly set
        # This ensures correct density calculation at elevated pressures
        local_pressure_bar = pressure_bar
        if phase == "gas" and local_pressure_bar is None and local_P1 is not None:
            local_pressure_bar = local_P1 / 100000.0  # Pa to bar
            results_log.append(f"Derived pressure_bar={local_pressure_bar:.3f} from inlet_pressure for property lookup")

        local_rho, local_mu, fluid_info = resolve_fluid_properties(
            fluid_name, temperature_c, local_pressure_bar,
            fluid_density, fluid_viscosity,
            fluid_density_lbft3, fluid_viscosity_cp,
            results_log, error_log
        )

        if local_rho is None:
            error_log.append("Missing required input: fluid_density")
        if local_mu is None:
            error_log.append("Missing required input: fluid_viscosity")

        if error_log:
            return json.dumps({"errors": error_log, "log": results_log})

        local_P2 = local_P1 - local_dP

        # Build solver arguments
        solver_kwargs = {
            'D': local_D,
            'D2': local_D2,
            'P1': local_P1,
            'P2': local_P2,
            'rho': local_rho,
            'mu': local_mu,
            'meter_type': local_meter_type,
            'taps': taps,
        }

        if phase == "liquid":
            solver_kwargs['epsilon_specified'] = 1.0
        else:
            if gas_gamma is None:
                error_log.append("Missing required input for gas: gas_gamma")
                return json.dumps({"errors": error_log, "log": results_log})
            solver_kwargs['k'] = gas_gamma

        if cd is not None:
            solver_kwargs['C_specified'] = cd

        try:
            local_m = differential_pressure_meter_solver(**solver_kwargs)
            local_Q = local_m / local_rho
            results_log.append(f"Solver converged: Q = {local_Q*3600:.3f} m³/h")
        except Exception as solver_e:
            error_log.append(f"Orifice solver failed: {solver_e}")
            return json.dumps({"errors": error_log, "log": results_log})

        # Calculate derived values
        beta = local_D2 / local_D
        A_orifice = math.pi * local_D2**2 / 4
        V_orifice = local_Q / A_orifice
        Re_D = local_rho * (local_Q / (math.pi * local_D**2 / 4)) * local_D / local_mu

        # Check validity (including gas expansibility check)
        validity_warnings = check_iso5167_validity(
            local_D, local_D2, local_dP, Re_D, taps,
            P1=local_P1, P2=local_P2, phase=phase, duty=duty
        )
        if validity_warnings:
            error_log.extend(validity_warnings)

        # The supplied dP is tap-to-tap; report what the system actually loses.
        perm_block = {}
        cav_block = None
        if phase == "liquid":
            perm_block = permanent_pressure_loss(
                local_D, local_D2, local_dP, local_m, local_rho, local_mu, taps, cd
            )
            if inlet_pressure is not None or inlet_pressure_psi is not None \
                    or inlet_pressure_bar is not None:
                cav_block = cavitation_assessment(
                    local_P1, local_P1 - perm_block["permanent_pressure_loss_pa"],
                    local_Q, local_D2, local_rho,
                    vapor_pressure_pa=vapor_pressure_pa,
                    bubble_point_pa=bubble_point_pa,
                )
                if cav_block and cav_block.get("cavitation_regime") in ("damage", "choking"):
                    error_log.append(
                        f"Cavitation regime '{cav_block['cavitation_regime']}' at sigma "
                        f"{cav_block['sigma_cavitation']}: the sqrt(dP) law used to infer "
                        f"this flow is not valid there, so the flow is unreliable."
                    )
                if cav_block and cav_block.get("gas_breakout"):
                    error_log.append(
                        "Pressure falls below the liquid's bubble point across this plate: "
                        "dissolved gas breaks out and the discharge coefficient no longer applies."
                    )

        # Format output
        output_dims = format_output_diameter(local_D2, size_units)

        result = {
            "mode": "analysis_flow",
            "phase": phase,
            "duty": duty,
            "inputs_resolved": results_log,
            "warnings": error_log if error_log else None,
            "fluid_details": fluid_info if fluid_info else None,
            **perm_block,
            "cavitation": cav_block,

            # Calculated flow rate
            "flow_rate_m3_s": round(local_Q, 6),
            "flow_rate_m3_hr": round(local_Q * 3600, 3),
            "flow_rate_gpm": round(local_Q / GPM_to_M3S, 2),
            "mass_flow_kg_s": round(local_m, 4),

            # Input geometry
            **output_dims,
            "pipe_diameter_mm": round(local_D * 1000, 2),
            "beta_ratio": round(beta, 4),

            # Pressure parameters
            "pressure_drop_bar": round(local_dP / BAR_TO_PA, 4),
            "pressure_drop_psi": round(local_dP / PSI_to_PA, 3),

            # Performance
            "velocity_orifice_m_s": round(V_orifice, 2),
            "reynolds_pipe": round(Re_D, 0),
        }

        result = {k: v for k, v in result.items() if v is not None}
        return json.dumps(result)

    except Exception as e:
        logger.error(f"Error in calculate_orifice_analysis_flow: {e}", exc_info=True)
        return json.dumps({
            "error": f"Calculation error: {str(e)}",
            "log": results_log,
            "errors": error_log
        })


def generate_plate_kit(
    # Flow inputs (design point)
    flow_rate: Optional[float] = None,
    flow_rate_gpm: Optional[float] = None,
    flow_rate_m3_hr: Optional[float] = None,

    # Pressure inputs (design point)
    pressure_drop: Optional[float] = None,
    pressure_drop_psi: Optional[float] = None,
    pressure_drop_bar: Optional[float] = None,
    inlet_pressure: Optional[float] = None,
    inlet_pressure_psi: Optional[float] = None,
    inlet_pressure_bar: Optional[float] = None,

    # Pipe geometry
    pipe_diameter: Optional[float] = None,
    pipe_diameter_in: Optional[float] = None,
    nominal_size_in: Optional[float] = None,
    schedule: str = "40",

    # Fluid properties
    fluid_name: Optional[str] = None,
    temperature_c: Optional[float] = None,
    pressure_bar: Optional[float] = None,
    fluid_density: Optional[float] = None,
    fluid_viscosity: Optional[float] = None,
    fluid_density_lbft3: Optional[float] = None,
    fluid_viscosity_cp: Optional[float] = None,

    # Gas properties
    gas_gamma: Optional[float] = None,

    # Orifice parameters
    meter_type: str = None,
    taps: str = 'D',

    # Kit parameters
    diameter_variability_percent: float = 20.0,
    kit_sizes: int = 5,

    # Phase selection
    phase: Literal["liquid", "gas"] = "liquid",

    # Output options
    size_units: Literal["m", "inch", "mm"] = "mm",
) -> str:
    """
    Generate a plate kit showing ΔP sensitivity to orifice diameter.

    Logic:
    1. Size design orifice for target Q and ΔP
    2. Generate plate diameters at ±X% of design diameter
    3. Calculate resulting ΔP at FIXED design flow for each plate

    This shows: "If you install the 45mm plate at design flow, expect ~3.1 bar ΔP"

    Returns:
        JSON string with plate kit and sensitivity analysis
    """
    if isinstance(kit_sizes, bool) or not isinstance(kit_sizes, int) or not 2 <= kit_sizes <= MAX_PLATE_KIT_SIZES:
        return json.dumps({"error_code": "RESOURCE_LIMIT", "errors": [
            f"kit_sizes must be an integer from 2 to {MAX_PLATE_KIT_SIZES}."
        ]})
    results_log = []
    error_log = []

    try:
        local_meter_type = meter_type if meter_type else ISO_5167_ORIFICE

        # Resolve inputs
        local_Q = resolve_flow_rate(flow_rate, flow_rate_gpm, flow_rate_m3_hr, results_log)
        if local_Q is None:
            error_log.append("Missing required input: flow_rate")

        local_dP = resolve_pressure_drop(pressure_drop, pressure_drop_psi, pressure_drop_bar, results_log)
        if local_dP is None:
            error_log.append("Missing required input: pressure_drop")

        local_D = resolve_pipe_diameter(pipe_diameter, pipe_diameter_in, nominal_size_in,
                                        schedule, results_log, error_log)
        if local_D is None:
            error_log.append("Missing required input: pipe_diameter")

        local_rho, local_mu, fluid_info = resolve_fluid_properties(
            fluid_name, temperature_c, pressure_bar,
            fluid_density, fluid_viscosity,
            fluid_density_lbft3, fluid_viscosity_cp,
            results_log, error_log
        )

        if local_rho is None or local_mu is None:
            error_log.append("Missing required input: fluid properties")

        # Resolve inlet pressure
        local_P1 = None
        if inlet_pressure is not None:
            local_P1 = inlet_pressure
        elif inlet_pressure_psi is not None:
            local_P1 = inlet_pressure_psi * PSI_to_PA
        elif inlet_pressure_bar is not None:
            local_P1 = inlet_pressure_bar * BAR_TO_PA

        if local_P1 is None and phase == "liquid" and local_dP is not None:
            local_P1 = 200000.0 + local_dP
            results_log.append("Using reference inlet pressure for liquid calculation")

        if error_log:
            return json.dumps({"errors": error_log, "log": results_log})

        local_P2 = local_P1 - local_dP
        local_m = local_Q * local_rho

        # Step 1: Size design orifice
        sizing_kwargs = {
            'D': local_D,
            'm': local_m,
            'P1': local_P1,
            'P2': local_P2,
            'rho': local_rho,
            'mu': local_mu,
            'meter_type': local_meter_type,
            'taps': taps,
        }

        if phase == "liquid":
            sizing_kwargs['epsilon_specified'] = 1.0
        else:
            if gas_gamma is None:
                error_log.append("Missing required input for gas: gas_gamma")
                return json.dumps({"errors": error_log, "log": results_log})
            sizing_kwargs['k'] = gas_gamma

        try:
            design_D2 = differential_pressure_meter_solver(**sizing_kwargs)
            results_log.append(f"Design orifice sized: {design_D2*1000:.2f} mm")
        except Exception as e:
            error_log.append(f"Failed to size design orifice: {e}")
            return json.dumps({"errors": error_log, "log": results_log})

        # Check ISO 5167 validity for design point (including gas expansibility check)
        design_Re_D = local_rho * (local_Q / (math.pi * local_D**2 / 4)) * local_D / local_mu
        validity_warnings = check_iso5167_validity(
            local_D, design_D2, local_dP, design_Re_D, taps,
            P1=local_P1, P2=local_P2, phase=phase
        )
        if validity_warnings:
            error_log.extend(validity_warnings)

        # Step 2: Generate plate diameters
        variability_fraction = diameter_variability_percent / 100.0
        D2_min = design_D2 * (1 - variability_fraction)
        D2_max = design_D2 * (1 + variability_fraction)

        # Generate evenly spaced diameters
        D2_values = []
        for i in range(kit_sizes):
            D2_values.append(D2_min + (D2_max - D2_min) * i / (kit_sizes - 1))

        # Step 3: Calculate ΔP for each plate at design flow
        plate_kit = []

        for i, D2 in enumerate(D2_values):
            # Calculate ΔP at design flow with this orifice
            analysis_kwargs = {
                'D': local_D,
                'D2': D2,
                'm': local_m,
                'P1': local_P1,
                'rho': local_rho,
                'mu': local_mu,
                'meter_type': local_meter_type,
                'taps': taps,
            }

            if phase == "liquid":
                analysis_kwargs['epsilon_specified'] = 1.0
            else:
                analysis_kwargs['k'] = gas_gamma

            try:
                P2_result = differential_pressure_meter_solver(**analysis_kwargs)
                dP_result = local_P1 - P2_result
            except Exception as e:
                # If solver fails for this plate, note it
                dP_result = None
                error_log.append(f"Plate {i+1} ({D2*1000:.1f} mm) calculation failed: {e}")

            beta = D2 / local_D
            is_design = abs(D2 - design_D2) < 0.0001

            plate_info = {
                "plate_number": i + 1,
                "orifice_diameter_mm": round(D2 * 1000, 1),
                "orifice_diameter_in": round(D2 / INCH_to_M, 3),
                "beta_ratio": round(beta, 3),
                "is_design_plate": is_design,
            }

            if dP_result is not None:
                plate_info["pressure_drop_bar"] = round(dP_result / BAR_TO_PA, 3)
                plate_info["pressure_drop_psi"] = round(dP_result / PSI_to_PA, 2)
                plate_info["pressure_drop_kpa"] = round(dP_result / 1000, 2)
            else:
                plate_info["pressure_drop_bar"] = "Error"

            plate_kit.append(plate_info)

        # Build result
        result = {
            "mode": "plate_kit",
            "phase": phase,
            "inputs_resolved": results_log,
            "warnings": error_log if error_log else None,
            "fluid_details": fluid_info if fluid_info else None,

            # Design point
            "design_flow_m3_hr": round(local_Q * 3600, 3),
            "design_flow_gpm": round(local_Q / GPM_to_M3S, 2),
            "design_pressure_drop_bar": round(local_dP / BAR_TO_PA, 3),
            "design_orifice_mm": round(design_D2 * 1000, 1),
            "design_orifice_in": round(design_D2 / INCH_to_M, 3),

            # Kit parameters
            "diameter_variability_percent": diameter_variability_percent,
            "kit_sizes": kit_sizes,
            "pipe_diameter_mm": round(local_D * 1000, 2),

            # Plate kit
            "plate_kit": plate_kit,

            # Interpretation help
            "interpretation": (
                f"At design flow of {local_Q*3600:.1f} m³/h, pressure drop varies with orifice size. "
                f"Design plate ({design_D2*1000:.1f} mm) gives {local_dP/BAR_TO_PA:.2f} bar. "
                f"Smaller plates increase ΔP, larger plates decrease ΔP."
            ),
        }

        result = {k: v for k, v in result.items() if v is not None}
        return json.dumps(result)

    except Exception as e:
        logger.error(f"Error in generate_plate_kit: {e}", exc_info=True)
        return json.dumps({
            "error": f"Calculation error: {str(e)}",
            "log": results_log,
            "errors": error_log
        })
