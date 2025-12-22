"""Unified orifice restrictor sizing for liquid and gas service."""

from typing import Optional, Literal
import inspect
from tools.orifice_restrictor import (
    calculate_orifice_sizing,
    calculate_orifice_analysis_dp,
    calculate_orifice_analysis_flow,
    generate_plate_kit,
)


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

    # Plate kit parameters
    diameter_variability_percent: float = 20.0,
    kit_sizes: int = 5,

    # Output formatting
    size_units: Literal["m", "inch", "mm"] = "mm",
) -> str:
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
        JSON string with calculation results:
        - sizing: orifice diameter, beta ratio, flow parameters
        - analysis_dp: pressure drop, head loss
        - analysis_flow: flow rate (volumetric and mass)
        - plate_kit: list of plates with ΔP at design flow
    """
    # Collect all parameters
    params = locals().copy()
    params.pop("mode")

    # Select function based on mode
    if mode == "sizing":
        fn = calculate_orifice_sizing
    elif mode == "analysis_dp":
        fn = calculate_orifice_analysis_dp
    elif mode == "analysis_flow":
        fn = calculate_orifice_analysis_flow
    elif mode == "plate_kit":
        fn = generate_plate_kit
    else:
        return f'{{"error": "Invalid mode: {mode}. Use sizing, analysis_dp, analysis_flow, or plate_kit"}}'

    # Get allowed parameters for the selected function
    sig = inspect.signature(fn)
    allowed = set(sig.parameters.keys())

    # Forward only allowed parameters that are not None
    forwarded = {k: v for k, v in params.items() if k in allowed and v is not None}

    return fn(**forwarded)
