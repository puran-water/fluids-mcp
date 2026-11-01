"""
Reynolds number calculation tool.

This module provides a tool to calculate Reynolds number for various fluid flows.
"""

import json
import logging
import math
import fluids
import fluids.core
from typing import Optional

from fluids_mcp.errors import CalculationError

# Import shared utilities
from fluids_mcp.utils.constants import (
    FT_to_M, INCH_to_M, LBFT3_to_KGM3, CENTIPOISE_to_PAS
)

# Configure logging
logger = logging.getLogger("fluids-mcp.reynolds_number")


def _json_error(code: str, message: str, log_messages: list, details: Optional[dict] = None) -> str:
    error = CalculationError(code, message, details=details)
    return json.dumps({
        "error_code": error.code,
        "error": error.message,
        "errors": [error.to_public_dict()],
        "log": log_messages,
    }, allow_nan=False)


def _numeric_detail(value):
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return value if math.isfinite(float(value)) else repr(value)
    return value


def _finite_input_error(name: str, value: Optional[float], log_messages: list) -> str:
    return _json_error(
        "INVALID_INPUT",
        f"{name} must be finite.",
        log_messages,
        details={name: _numeric_detail(value)},
    )


def _positive_input_error(name: str, value: Optional[float], log_messages: list) -> str:
    return _json_error(
        "INVALID_INPUT",
        f"{name} must be positive.",
        log_messages,
        details={name: value},
    )


def _validate_positive(name: str, value: Optional[float], log_messages: list) -> Optional[str]:
    if value is None:
        return None
    if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(float(value)):
        return _finite_input_error(name, value, log_messages)
    if value <= 0:
        return _positive_input_error(name, value, log_messages)
    return None


def calculate_reynolds_number(
    # --- Core SI inputs ---
    velocity: Optional[float] = None,        # Fluid velocity in m/s
    characteristic_length: Optional[float] = None,  # Characteristic length (pipe diameter) in m
    fluid_density: Optional[float] = None,   # Fluid density in kg/m³
    fluid_viscosity: Optional[float] = None, # Fluid dynamic viscosity in Pa·s
    fluid_kinematic_viscosity: Optional[float] = None,  # Fluid kinematic viscosity in m²/s
    
    # --- Alternative unit inputs ---
    velocity_ft_s: Optional[float] = None,   # Velocity in ft/s
    characteristic_length_in: Optional[float] = None, # Characteristic length in inches
    fluid_density_lbft3: Optional[float] = None, # Fluid density in lb/ft³
    fluid_viscosity_cp: Optional[float] = None, # Fluid viscosity in centipoise
    
    # --- Property Lookup Inputs ---
    fluid_name: Optional[str] = None,        # e.g., "Water", "Air"
    temperature_c: Optional[float] = None,   # Temperature in Celsius (for property lookup)
    pressure_bar: Optional[float] = None     # Pressure in bar (for property lookup)
) -> str:
    """Calculate Reynolds number for a flow with flexible input options.
    
    Args:
        velocity: Fluid velocity in m/s
        characteristic_length: Characteristic length (pipe diameter) in m
        fluid_density: Fluid density in kg/m³ (required if using dynamic viscosity)
        fluid_viscosity: Fluid dynamic viscosity in Pa·s (can use this OR kinematic viscosity)
        fluid_kinematic_viscosity: Fluid kinematic viscosity in m²/s (can use this OR density+viscosity)
        velocity_ft_s: Velocity in ft/s
        characteristic_length_in: Characteristic length in inches
        fluid_density_lbft3: Fluid density in lb/ft³
        fluid_viscosity_cp: Fluid viscosity in centipoise
        fluid_name: Name of the fluid (e.g., "Water", "Air")
        temperature_c: Temperature in Celsius (for property lookup)
        pressure_bar: Pressure in bar (for property lookup)
    
    Returns:
        JSON string with Reynolds number and flow regime
    """
    results_log = []
    error_log = []
    
    try:
        if temperature_c is not None:
            temperature_error = _validate_positive("temperature_k", temperature_c + 273.15, results_log)
            if temperature_error is not None:
                return temperature_error
        if pressure_bar is not None:
            pressure_error = _validate_positive("pressure_bar", pressure_bar, results_log)
            if pressure_error is not None:
                return pressure_error

        # 1. Resolve Velocity (SI)
        local_velocity = None
        if velocity is not None:
            local_velocity = velocity
            results_log.append("Used provided SI velocity.")
        elif velocity_ft_s is not None:
            local_velocity = velocity_ft_s * FT_to_M
            results_log.append(f"Converted velocity from {velocity_ft_s} ft/s.")
        else:
            error_log.append("Missing required input: velocity or velocity_ft_s.")
            
        # 2. Resolve Characteristic Length (SI)
        local_characteristic_length = None
        if characteristic_length is not None:
            local_characteristic_length = characteristic_length
            results_log.append("Used provided SI characteristic_length.")
        elif characteristic_length_in is not None:
            local_characteristic_length = characteristic_length_in * INCH_to_M
            results_log.append(f"Converted characteristic_length from {characteristic_length_in} inches.")
        else:
            error_log.append("Missing required input: characteristic_length or characteristic_length_in.")
        
        # 3. Resolve Fluid Properties (SI)
        local_fluid_density = None
        local_fluid_viscosity = None
        local_fluid_kinematic_viscosity = None

        # Priority 1: Direct kinematic viscosity (simplest case)
        if fluid_kinematic_viscosity is not None:
            local_fluid_kinematic_viscosity = fluid_kinematic_viscosity
            results_log.append("Used provided SI fluid_kinematic_viscosity.")
            calculation_method = "kinematic"
        # Priority 2: Direct density and viscosity
        elif fluid_density is not None and fluid_viscosity is not None:
            local_fluid_density = fluid_density
            local_fluid_viscosity = fluid_viscosity
            results_log.append("Used provided SI fluid density and viscosity.")
            calculation_method = "dynamic"
        # Priority 3: Imperial conversions
        elif fluid_density_lbft3 is not None and fluid_viscosity_cp is not None:
            local_fluid_density = fluid_density_lbft3 * LBFT3_to_KGM3
            local_fluid_viscosity = fluid_viscosity_cp * CENTIPOISE_to_PAS
            results_log.append(f"Converted density from {fluid_density_lbft3} lb/ft³ and viscosity from {fluid_viscosity_cp} cP.")
            calculation_method = "dynamic"
        # Priority 4: Fluid property lookup via centralized resolver (CoolProp → fluidprop)
        elif fluid_name is not None and temperature_c is not None:
            try:
                from fluids_mcp.utils.resolve_properties import resolve_liquid_properties
                p_bar = pressure_bar if pressure_bar is not None else 1.01325
                props = resolve_liquid_properties(fluid_name, temperature_c, p_bar)
                if props is not None:
                    local_fluid_kinematic_viscosity = props.kinematic_viscosity
                    local_fluid_density = props.density
                    local_fluid_viscosity = props.viscosity
                    results_log.append(f"Properties via {props.source} for {fluid_name} at {temperature_c}°C, {p_bar} bar.")
                    results_log.extend(props.log)
                    calculation_method = "kinematic"
                else:
                    error_log.append(f"Failed to look up fluid properties for {fluid_name} at {temperature_c}°C: all backends failed")
                    local_fluid_density = None
                    local_fluid_viscosity = None
                    local_fluid_kinematic_viscosity = None
                    calculation_method = None
            except Exception as fluid_lookup_e:
                error_log.append(f"Failed to look up fluid properties for {fluid_name} at {temperature_c}°C: {fluid_lookup_e}")
                local_fluid_density = None
                local_fluid_viscosity = None
                local_fluid_kinematic_viscosity = None
                calculation_method = None
        else:
            error_log.append("Missing required inputs for fluid properties: fluid_kinematic_viscosity OR (fluid_density AND fluid_viscosity) OR (fluid_density_lbft3 AND fluid_viscosity_cp) OR (fluid_name AND temperature_c).")
            calculation_method = None

        if local_velocity is None:
            return _json_error(
                "INVALID_INPUT",
                "velocity or velocity_ft_s is required.",
                results_log,
                details={"missing": "velocity"},
            )
        velocity_error = _validate_positive("velocity", local_velocity, results_log)
        if velocity_error is not None:
            return velocity_error
        if local_characteristic_length is None:
            return _json_error(
                "INVALID_INPUT",
                "characteristic_length or characteristic_length_in is required.",
                results_log,
                details={"missing": "characteristic_length"},
            )
        length_error = _validate_positive("characteristic_length", local_characteristic_length, results_log)
        if length_error is not None:
            return length_error
        if calculation_method == "kinematic" and local_fluid_kinematic_viscosity is not None:
            kinematic_error = _validate_positive(
                "fluid_kinematic_viscosity",
                local_fluid_kinematic_viscosity,
                results_log,
            )
            if kinematic_error is not None:
                return kinematic_error
        if calculation_method == "dynamic":
            density_error = _validate_positive("fluid_density", local_fluid_density, results_log)
            if density_error is not None:
                return density_error
            viscosity_error = _validate_positive("fluid_viscosity", local_fluid_viscosity, results_log)
            if viscosity_error is not None:
                return viscosity_error

        # After inputs are resolved, calculate Re based on available data
        if calculation_method == "kinematic" and local_fluid_kinematic_viscosity is not None:
            Re = fluids.core.Reynolds(V=local_velocity, D=local_characteristic_length, nu=local_fluid_kinematic_viscosity)
            results_log.append("Calculated Reynolds number using kinematic viscosity.")
        elif calculation_method == "dynamic" and local_fluid_density is not None and local_fluid_viscosity is not None:
            Re = fluids.core.Reynolds(V=local_velocity, D=local_characteristic_length, rho=local_fluid_density, mu=local_fluid_viscosity)
            results_log.append("Calculated Reynolds number using density and viscosity.")
        else:
            return json.dumps({"errors": ["Failed to calculate Reynolds number. Missing fluid properties."], "log": results_log})
        if not math.isfinite(float(Re)):
            return _json_error(
                "INVALID_RESULT_NUMERIC",
                "Reynolds number must be finite.",
                results_log,
                details={"reynolds_number": repr(Re)},
            )
        # Determine flow regime
        if Re < 2300:
            regime = "laminar"
        elif Re < 4000:
            regime = "transitional"
        else:
            regime = "turbulent"
        
        # Return the result
        result = {
            "inputs_resolved": results_log,
            "warnings": error_log if error_log else None,
            "reynolds_number": round(Re, 2),
            "flow_regime": regime
        }
        
        # Remove null warnings if empty
        if not result["warnings"]: del result["warnings"]
        
        return json.dumps(result, allow_nan=False)
        
    except Exception as e:
        logger.error(f"Error in calculate_reynolds_number: {e}", exc_info=True)
        return _json_error(
            "REYNOLDS_CALCULATION_FAILED",
            f"Calculation error: {str(e)}",
            results_log,
            details={"errors_occurred": error_log},
        )
