"""
Fixed orifice (restriction plate) sizing calculations.

This module provides functions for sizing fixed restriction orifices used in
applications like RO permeate backpressure control.

Leverages the fluids library's differential_pressure_meter_solver() for
ISO 5167 orifice calculations with Reader-Harris-Gallagher discharge coefficients.
"""

import json
import logging
import math
from typing import Optional, Literal, List, Dict, Any

# Third-party libraries
import fluids
import fluids.piping
from fluids.flow_meter import (
    differential_pressure_meter_solver,
    flow_meter_discharge,
    ISO_5167_ORIFICE,
    C_Reader_Harris_Gallagher,
    orifice_expansibility,
)

# Import shared utilities
from utils.constants import (
    GPM_to_M3S, INCH_to_M, FT_to_M, PSI_to_PA,
    LBFT3_to_KGM3, CENTIPOISE_to_PAS, DEG_C_to_K,
    DEFAULT_ATMOSPHERIC_PRESSURE
)
from utils.import_helpers import FLUIDPROP_AVAILABLE, FluidProperties, COOLPROP_AVAILABLE, CP
from utils.fluid_aliases import map_fluid_name

# Configure logging
logger = logging.getLogger("fluids-mcp.orifice_restrictor")

# Unit conversion constants
BAR_TO_PA = 100000.0
M3HR_TO_M3S = 1.0 / 3600.0
MM_TO_M = 0.001


def check_iso5167_validity(D: float, D2: float, dP: float, Re_D: Optional[float] = None,
                           taps: str = 'D') -> List[str]:
    """
    Check ISO 5167 validity limits and return warnings.

    Args:
        D: Pipe inner diameter in meters
        D2: Orifice diameter in meters
        dP: Pressure drop in Pa
        Re_D: Reynolds number based on pipe diameter (optional)
        taps: Tap configuration ('D', 'flange', 'corner')

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

    # Priority 3: Property lookup via FluidProp
    if (local_density is None or local_viscosity is None) and fluid_name is not None and temperature_c is not None:
        if FLUIDPROP_AVAILABLE and FluidProperties is not None:
            try:
                mapped_name = map_fluid_name(fluid_name)
                p_bar = pressure_bar if pressure_bar is not None else 1.01325

                fluid_props = FluidProperties(
                    coolprop_name=mapped_name,
                    T_in_deg_C=temperature_c,
                    P_in_bar=p_bar
                )

                if local_density is None:
                    local_density = float(fluid_props.rho[0])
                if local_viscosity is None:
                    local_viscosity = float(fluid_props.eta[0])

                fluid_info = {
                    "name_used": mapped_name,
                    "temperature_c": temperature_c,
                    "pressure_bar": p_bar,
                    "density_kg_m3": local_density,
                    "viscosity_pa_s": local_viscosity
                }
                results_log.append(f"Looked up properties for {mapped_name} at {temperature_c}°C, {p_bar} bar.")

            except Exception as e:
                error_log.append(f"Fluid property lookup failed: {e}")
        else:
            error_log.append("Fluid property lookup unavailable (FluidProp not installed).")

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

        # Resolve fluid properties
        local_rho, local_mu, fluid_info = resolve_fluid_properties(
            fluid_name, temperature_c, pressure_bar,
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
        if local_P1 is None and phase == "liquid" and local_dP is not None:
            # Use a safe reference pressure for incompressible liquid calculation
            local_P1 = 200000.0 + local_dP  # 2 bar above outlet
            results_log.append(f"Using reference inlet pressure: {local_P1/1000:.1f} kPa (liquid shortcut)")

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
            D2 = differential_pressure_meter_solver(**solver_kwargs)
            results_log.append(f"Solver converged: D2 = {D2*1000:.2f} mm")
        except Exception as solver_e:
            error_log.append(f"Orifice solver failed: {solver_e}")
            return json.dumps({"errors": error_log, "log": results_log})

        # Calculate derived values
        beta = D2 / local_D
        A_orifice = math.pi * D2**2 / 4
        V_orifice = local_Q / A_orifice
        Re_D = local_rho * (local_Q / (math.pi * local_D**2 / 4)) * local_D / local_mu
        Re_orifice = local_rho * V_orifice * D2 / local_mu

        # Check ISO 5167 validity
        validity_warnings = check_iso5167_validity(local_D, D2, local_dP, Re_D, taps)
        if validity_warnings:
            error_log.extend(validity_warnings)

        # Format output
        output_dims = format_output_diameter(D2, size_units)

        result = {
            "mode": "sizing",
            "phase": phase,
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

            # Pressure parameters
            "pressure_drop_bar": round(local_dP / BAR_TO_PA, 3),
            "pressure_drop_psi": round(local_dP / PSI_to_PA, 2),
            "inlet_pressure_bar": round(local_P1 / BAR_TO_PA, 3),
            "outlet_pressure_bar": round(local_P2 / BAR_TO_PA, 3),

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

        local_rho, local_mu, fluid_info = resolve_fluid_properties(
            fluid_name, temperature_c, pressure_bar,
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
        if phase == "liquid":
            epsilon = 1.0
            # Use a high reference pressure for iteration (10 bar)
            # This ensures the solver can find P2 even for high pressure drops
            local_P1 = 1000000.0  # 10 bar reference
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

        # Check validity
        validity_warnings = check_iso5167_validity(local_D, local_D2, local_dP, Re_D, taps)
        if validity_warnings:
            error_log.extend(validity_warnings)

        # Format output
        output_dims = format_output_diameter(local_D2, size_units)

        result = {
            "mode": "analysis_dp",
            "phase": phase,
            "inputs_resolved": results_log,
            "warnings": error_log if error_log else None,
            "fluid_details": fluid_info if fluid_info else None,

            # Calculated pressure drop
            "pressure_drop_pa": round(local_dP, 1),
            "pressure_drop_bar": round(local_dP / BAR_TO_PA, 4),
            "pressure_drop_psi": round(local_dP / PSI_to_PA, 3),
            "head_loss_m": round(local_dP / (local_rho * 9.80665), 3),

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

        local_rho, local_mu, fluid_info = resolve_fluid_properties(
            fluid_name, temperature_c, pressure_bar,
            fluid_density, fluid_viscosity,
            fluid_density_lbft3, fluid_viscosity_cp,
            results_log, error_log
        )

        if local_rho is None:
            error_log.append("Missing required input: fluid_density")
        if local_mu is None:
            error_log.append("Missing required input: fluid_viscosity")

        # Resolve inlet pressure
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

        # Check validity
        validity_warnings = check_iso5167_validity(local_D, local_D2, local_dP, Re_D, taps)
        if validity_warnings:
            error_log.extend(validity_warnings)

        # Format output
        output_dims = format_output_diameter(local_D2, size_units)

        result = {
            "mode": "analysis_flow",
            "phase": phase,
            "inputs_resolved": results_log,
            "warnings": error_log if error_log else None,
            "fluid_details": fluid_info if fluid_info else None,

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

        # Step 2: Generate plate diameters
        variability_fraction = diameter_variability_percent / 100.0
        D2_min = design_D2 * (1 - variability_fraction)
        D2_max = design_D2 * (1 + variability_fraction)

        if kit_sizes <= 1:
            kit_sizes = 5

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
