"""
Pipe sizing tools for automatic selection based on pressure drop constraints.
"""

import json
import logging
import math
from typing import Any, Optional
from fluids import friction_factor, Reynolds
from fluids.piping import nearest_pipe
import CoolProp.CoolProp as CP

from fluids_mcp.domain.units import FT_TO_M, GPM_TO_M3S, PSI_TO_PA
from fluids_mcp.errors import CalculationError, FittingKError
from fluids_mcp.resource_limits import MAX_NPS_CANDIDATES

# Import helpers from existing modules
from fluids_mcp.utils.helpers import get_fitting_K, get_pipe_roughness, validate_fitting_type

logger = logging.getLogger(__name__)

# Default NPS candidates for Schedule 40
DEFAULT_NPS_CANDIDATES = [0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4, 6, 8, 10, 12]
CANONICAL_GAS_PIPE_SIZING_METHOD = "isothermal_darcy"
GAS_PIPE_SIZING_METHOD_ALIASES = {
    "isothermal": CANONICAL_GAS_PIPE_SIZING_METHOD,
    CANONICAL_GAS_PIPE_SIZING_METHOD: CANONICAL_GAS_PIPE_SIZING_METHOD,
}
SUPPORTED_GAS_PIPE_SIZING_METHODS = set(GAS_PIPE_SIZING_METHOD_ALIASES)


def _json_error(code: str, message: str, log_messages: list, details: Optional[dict] = None) -> str:
    error = CalculationError(code, message, details=details)
    return json.dumps(
        {
            "error_code": error.code,
            "error": error.message,
            "errors": [error.to_public_dict()],
            "log": log_messages,
        },
        allow_nan=False,
    )


def _numeric_detail(value: Any) -> Any:
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return value if math.isfinite(float(value)) else repr(value)
    return value


def _positive_input_error(name: str, value: Any, log_messages: list) -> str:
    return _json_error(
        "INVALID_NUMERIC_INPUT",
        f"{name} must be positive.",
        log_messages,
        details={name: value},
    )


def _validate_finite_number(name: str, value: Any, log_messages: list) -> Optional[str]:
    if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(float(value)):
        return _json_error(
            "INVALID_NUMERIC_INPUT",
            f"{name} must be finite.",
            log_messages,
            details={name: _numeric_detail(value)},
        )
    return None


def _validate_positive_number(name: str, value: Any, log_messages: list) -> Optional[str]:
    finite_error = _validate_finite_number(name, value, log_messages)
    if finite_error is not None:
        return finite_error
    if value <= 0:
        return _positive_input_error(name, value, log_messages)
    return None


def _validate_nonnegative_number(name: str, value: Any, log_messages: list) -> Optional[str]:
    finite_error = _validate_finite_number(name, value, log_messages)
    if finite_error is not None:
        return finite_error
    if value < 0:
        return _json_error(
            "INVALID_NUMERIC_INPUT",
            f"{name} must be non-negative.",
            log_messages,
            details={name: value},
        )
    return None


def _validate_nps_candidates(nps_list: Any, log_messages: list) -> Optional[str]:
    if not isinstance(nps_list, list):
        return _json_error(
            "INVALID_NUMERIC_INPUT",
            "nps_candidates must be a list of positive finite values.",
            log_messages,
            details={"nps_candidates": nps_list},
        )
    if not 1 <= len(nps_list) <= MAX_NPS_CANDIDATES:
        return _json_error(
            "RESOURCE_LIMIT", f"nps_candidates must contain 1 to {MAX_NPS_CANDIDATES} values.",
            log_messages, details={"count": len(nps_list), "maximum": MAX_NPS_CANDIDATES},
        )
    for index, nps in enumerate(nps_list):
        error = _validate_positive_number(f"nps_candidates[{index}]", nps, log_messages)
        if error is not None:
            return error
    return None


def _validate_fittings(fittings: Any, log_messages: list) -> Optional[str]:
    if fittings is None:
        return None
    if not isinstance(fittings, list):
        return _json_error(
            "INVALID_NUMERIC_INPUT",
            "fittings must be a list.",
            log_messages,
            details={"fittings": fittings},
        )
    for index, fitting in enumerate(fittings):
        if not isinstance(fitting, dict):
            return _json_error(
                "INVALID_NUMERIC_INPUT",
                "each fitting must be an object.",
                log_messages,
                details={"fitting_index": index, "fitting": fitting},
            )
        quantity = fitting.get("quantity", 1)
        quantity_error = _validate_nonnegative_number(f"fittings[{index}].quantity", quantity, log_messages)
        if quantity_error is not None:
            return quantity_error
        if "K_value" in fitting and fitting.get("K_value") is not None:
            k_error = _validate_nonnegative_number(
                f"fittings[{index}].K_value",
                fitting.get("K_value"),
                log_messages,
            )
            if k_error is not None:
                return k_error
        else:
            validate_fitting_type(fitting.get("type"))
    return None


def _candidate_failure_error(candidates: list, log_messages: list) -> Optional[str]:
    if not candidates or all("error" in candidate for candidate in candidates):
        return _json_error(
            "PIPE_SIZING_NO_VALID_CANDIDATES",
            "No pipe-size candidates could be evaluated.",
            log_messages,
            details={"candidates": candidates},
        )
    return None


def select_liquid_pipe_size_for_dp(
    # Flow rate (provide one)
    flow_rate: Optional[float] = None,          # m³/s
    flow_rate_gpm: Optional[float] = None,      # US GPM

    # Pipe length (provide one)
    pipe_length: Optional[float] = None,        # m
    pipe_length_ft: Optional[float] = None,     # feet

    # Fluid properties (provide directly or via lookup)
    fluid_density: Optional[float] = None,      # kg/m³
    fluid_density_lbft3: Optional[float] = None, # lb/ft³
    fluid_viscosity: Optional[float] = None,    # Pa·s
    fluid_viscosity_cp: Optional[float] = None, # centipoise
    fluid_name: Optional[str] = None,           # For property lookup
    temperature_c: Optional[float] = None,      # °C
    pressure_bar: Optional[float] = None,       # bar

    # Pressure drop constraint (provide one)
    allowable_dp_pa: Optional[float] = None,    # Pa
    allowable_dp_psi: Optional[float] = None,   # psi

    # Pipe specifications
    schedule: str = "40",
    nps_candidates: Optional[list] = None,
    material: Optional[str] = None,
    pipe_roughness: Optional[float] = None,     # m

    # Fittings
    fittings: Optional[list] = None,

    # Constraints
    velocity_min_m_s: Optional[float] = None,
    velocity_max_m_s: Optional[float] = None,

    # Output options
    return_all_candidates: bool = False
) -> str:
    """
    Select optimal pipe size for liquid flow based on allowable pressure drop.

    Iterates through standard pipe sizes to find the smallest that meets
    the pressure drop constraint.

    Args:
        flow_rate: Flow rate in m³/s (SI)
        flow_rate_gpm: Flow rate in US GPM (alternative)
        pipe_length: Pipe length in m (SI)
        pipe_length_ft: Pipe length in feet (alternative)
        fluid_density: Fluid density in kg/m³ (SI)
        fluid_density_lbft3: Fluid density in lb/ft³ (alternative)
        fluid_viscosity: Dynamic viscosity in Pa·s (SI)
        fluid_viscosity_cp: Dynamic viscosity in centipoise (alternative)
        fluid_name: Name of fluid for property lookup
        temperature_c: Temperature in Celsius (for property lookup)
        pressure_bar: Pressure in bar (for property lookup)
        allowable_dp_pa: Maximum allowable pressure drop in Pa
        allowable_dp_psi: Maximum allowable pressure drop in psi
        schedule: Pipe schedule (default: "40")
        nps_candidates: List of NPS sizes to evaluate (default: standard sizes)
        material: Pipe material for roughness
        pipe_roughness: Absolute roughness in m
        fittings: List of fittings with type and quantity
        velocity_min_m_s: Minimum velocity constraint in m/s
        velocity_max_m_s: Maximum velocity constraint in m/s
        return_all_candidates: Return all evaluated sizes (default: False)

    Returns:
        JSON string with recommended size and evaluation results
    """
    try:
        log_messages = []
        warnings = []

        kwargs = {k: v for k, v in locals().items() if k != 'kwargs' and v is not None}

        # Get NPS candidates
        nps_list = kwargs.get('nps_candidates', DEFAULT_NPS_CANDIDATES)
        limit_error = _validate_nps_candidates(nps_list, log_messages)
        if limit_error is not None:
            return limit_error
        schedule = kwargs.get('schedule', '40')

        # Convert flow rate to SI
        flow_rate_m3s = kwargs.get('flow_rate')
        if flow_rate_m3s is None and 'flow_rate_gpm' in kwargs:
            flow_rate_m3s = kwargs['flow_rate_gpm'] * GPM_TO_M3S
            log_messages.append(f"Converted flow_rate from {kwargs['flow_rate_gpm']} GPM to {flow_rate_m3s:.6f} m³/s")

        if flow_rate_m3s is None:
            return json.dumps({"error": "Flow rate must be specified"})
        flow_error = _validate_positive_number("flow_rate", flow_rate_m3s, log_messages)
        if flow_error is not None:
            return flow_error

        # Convert pipe length to SI
        pipe_length_m = kwargs.get('pipe_length')
        if pipe_length_m is None and 'pipe_length_ft' in kwargs:
            pipe_length_m = kwargs['pipe_length_ft'] * FT_TO_M
            log_messages.append(f"Converted pipe_length from {kwargs['pipe_length_ft']} ft to {pipe_length_m:.2f} m")

        if pipe_length_m is None:
            return json.dumps({"error": "Pipe length must be specified"})
        length_error = _validate_positive_number("pipe_length", pipe_length_m, log_messages)
        if length_error is not None:
            return length_error

        # Get fluid properties
        fluid_density = kwargs.get('fluid_density')
        fluid_viscosity = kwargs.get('fluid_viscosity')

        # Try property lookup if not provided
        if (fluid_density is None or fluid_viscosity is None) and 'fluid_name' in kwargs:
            temp_c = kwargs.get('temperature_c', 25)
            pressure_bar = kwargs.get('pressure_bar', 1.01325)

            try:
                fluid_name = kwargs['fluid_name']
                if fluid_name.lower() == 'water':
                    fluid_name = 'Water'
                elif fluid_name.lower() == 'air':
                    fluid_name = 'Air'

                pressure_pa = pressure_bar * 100000
                temp_k = temp_c + 273.15

                fluid_density = CP.PropsSI('D', 'T', temp_k, 'P', pressure_pa, fluid_name)
                fluid_viscosity = CP.PropsSI('V', 'T', temp_k, 'P', pressure_pa, fluid_name)
                log_messages.append(f"Looked up properties for {fluid_name} at {temp_c}°C, {pressure_bar} bar")
            except Exception as e:
                log_messages.append(f"Property lookup failed: {str(e)}")

        # Convert from imperial if needed
        if fluid_density is None and 'fluid_density_lbft3' in kwargs:
            fluid_density = kwargs['fluid_density_lbft3'] * 16.0185
            log_messages.append(f"Converted density from {kwargs['fluid_density_lbft3']} lb/ft³")

        if fluid_viscosity is None and 'fluid_viscosity_cp' in kwargs:
            fluid_viscosity = kwargs['fluid_viscosity_cp'] * 0.001
            log_messages.append(f"Converted viscosity from {kwargs['fluid_viscosity_cp']} cp")

        if fluid_density is None or fluid_viscosity is None:
            return json.dumps({"error": "Fluid properties must be specified or derivable"})
        density_error = _validate_positive_number("fluid_density", fluid_density, log_messages)
        if density_error is not None:
            return density_error
        viscosity_error = _validate_positive_number("fluid_viscosity", fluid_viscosity, log_messages)
        if viscosity_error is not None:
            return viscosity_error

        # Get allowable pressure drop
        allowable_dp = kwargs.get('allowable_dp_pa')
        if allowable_dp is None and 'allowable_dp_psi' in kwargs:
            allowable_dp = kwargs['allowable_dp_psi'] * PSI_TO_PA
            log_messages.append(f"Converted allowable_dp from {kwargs['allowable_dp_psi']} psi")

        if allowable_dp is None:
            return json.dumps({"error": "Allowable pressure drop must be specified"})
        allowable_dp_error = _validate_positive_number("allowable_dp", allowable_dp, log_messages)
        if allowable_dp_error is not None:
            return allowable_dp_error

        # Get velocity constraints
        v_min = kwargs.get('velocity_min_m_s', 0)
        v_max_provided = 'velocity_max_m_s' in kwargs
        v_max = kwargs.get('velocity_max_m_s', float('inf'))
        v_min_error = _validate_nonnegative_number("velocity_min_m_s", v_min, log_messages)
        if v_min_error is not None:
            return v_min_error
        if v_max_provided:
            v_max_error = _validate_positive_number("velocity_max_m_s", v_max, log_messages)
            if v_max_error is not None:
                return v_max_error
        if v_min > v_max:
            return _json_error(
                "INVALID_NUMERIC_INPUT",
                "velocity_min_m_s must be less than or equal to velocity_max_m_s.",
                log_messages,
                details={"velocity_min_m_s": v_min, "velocity_max_m_s": v_max},
            )

        # Get pipe roughness
        pipe_roughness = kwargs.get('pipe_roughness')
        if pipe_roughness is None and 'material' in kwargs:
            pipe_roughness, source = get_pipe_roughness(material=kwargs['material'], warnings=warnings)
            log_messages.append(f"Looked up roughness for material '{kwargs['material']}': {pipe_roughness} m ({source})")
        if pipe_roughness is None:
            pipe_roughness = 1.5e-5  # Default for smooth pipe
            log_messages.append(f"Using default roughness: {pipe_roughness} m")
        roughness_error = _validate_nonnegative_number("pipe_roughness", pipe_roughness, log_messages)
        if roughness_error is not None:
            return roughness_error
        nps_error = _validate_nps_candidates(nps_list, log_messages)
        if nps_error is not None:
            return nps_error

        # Get fittings
        fittings = kwargs.get('fittings', [])
        fittings_error = _validate_fittings(fittings, log_messages)
        if fittings_error is not None:
            return fittings_error

        # Evaluate each NPS candidate
        candidates = []
        recommended_nps = None

        for nps in nps_list:
            try:
                # Get pipe dimensions
                pipe_data = nearest_pipe(NPS=nps, schedule=schedule)
                # Returns tuple: (NPS, Di, Do, t)
                Di = pipe_data[1]  # Inner diameter in meters
                Do = pipe_data[2]  # Outer diameter in meters

                # Calculate flow parameters
                area = 3.14159 * (Di/2)**2
                velocity = flow_rate_m3s / area
                Re = Reynolds(V=velocity, D=Di, rho=fluid_density, mu=fluid_viscosity)

                # Check velocity constraints
                velocity_ok = v_min <= velocity <= v_max

                # Calculate friction factor
                relative_roughness = pipe_roughness / Di
                f = friction_factor(Re=Re, eD=relative_roughness)

                # Calculate straight pipe pressure drop
                dp_straight = f * (pipe_length_m / Di) * (fluid_density * velocity**2 / 2)

                # Calculate fitting losses
                K_total = 0
                fitting_details = []
                for fitting in fittings:
                    fitting_type = fitting.get('type')
                    quantity = fitting.get('quantity', 1)
                    K_value = fitting.get('K_value')

                    if K_value is None:
                        # Use standard K-factor calculation
                        K_value = get_fitting_K(fitting_type, Di, Re, flow_rate_m3s)

                    K_total += K_value * quantity
                    fitting_details.append({
                        'type': fitting_type,
                        'quantity': quantity,
                        'K_individual': K_value,
                        'K_total': K_value * quantity
                    })

                dp_fittings = K_total * (fluid_density * velocity**2 / 2)
                dp_total = dp_straight + dp_fittings

                # Store candidate results
                candidate = {
                    'nps_in': nps,
                    'Di_m': round(Di, 5),
                    'Do_m': round(Do, 5),
                    'velocity_m_s': round(velocity, 3),
                    'reynolds_number': round(Re, 0),
                    'friction_factor': round(f, 6),
                    'dp_straight_pa': round(dp_straight, 1),
                    'dp_fittings_pa': round(dp_fittings, 1),
                    'dp_total_pa': round(dp_total, 1),
                    'dp_total_psi': round(dp_total / PSI_TO_PA, 3),
                    'meets_dp_constraint': dp_total <= allowable_dp,
                    'meets_velocity_constraints': velocity_ok,
                    'feasible': dp_total <= allowable_dp and velocity_ok
                }

                if fitting_details:
                    candidate['fitting_details'] = fitting_details

                candidates.append(candidate)

                # Check if this is the recommended size (smallest feasible)
                if candidate['feasible'] and recommended_nps is None:
                    recommended_nps = nps

            except FittingKError:
                raise
            except Exception as e:
                log_messages.append(f"Error evaluating NPS {nps}: {str(e)}")
                candidates.append({
                    'nps_in': nps,
                    'error': str(e)
                })

        # Prepare results
        candidate_failure = _candidate_failure_error(candidates, log_messages)
        if candidate_failure is not None:
            return candidate_failure

        result = {
            'inputs_resolved': log_messages,
            'allowable_dp_pa': allowable_dp,
            'allowable_dp_psi': allowable_dp / PSI_TO_PA,
            'flow_rate_m3s': flow_rate_m3s,
            'flow_rate_gpm': flow_rate_m3s / GPM_TO_M3S,
            'pipe_length_m': pipe_length_m,
            'schedule': schedule
        }
        if warnings:
            result['warnings'] = warnings

        if recommended_nps is not None:
            result['recommended_nps_in'] = recommended_nps
            # Find the recommended candidate details
            for c in candidates:
                if c.get('nps_in') == recommended_nps:
                    result['recommended_details'] = {
                        'Di_m': c['Di_m'],
                        'velocity_m_s': c['velocity_m_s'],
                        'dp_total_pa': c['dp_total_pa'],
                        'dp_total_psi': c['dp_total_psi'],
                        'reynolds_number': c['reynolds_number']
                    }
                    break
        else:
            result['recommended_nps_in'] = None
            result['note'] = "No pipe size meets the specified constraints"

        if kwargs.get('return_all_candidates', False):
            result['all_candidates'] = candidates
        else:
            # Return only feasible candidates
            result['feasible_candidates'] = [c for c in candidates if c.get('feasible', False)]

        return json.dumps(result, indent=2, allow_nan=False)

    except FittingKError as e:
        logger.warning("Fitting K-value resolution failed: %s", e.message)
        return _json_error(e.code, e.message, log_messages, details=e.details)
    except Exception as e:
        logger.error(f"Error in select_liquid_pipe_size_for_dp: {str(e)}")
        return _json_error("PIPE_SIZING_CALCULATION_FAILED", str(e), [])


def select_gas_pipe_size_for_dp(
    # Flow rate (provide one)
    flow_rate_kg_s: Optional[float] = None,
    flow_rate_norm_m3_hr: Optional[float] = None,
    flow_rate_std_m3_hr: Optional[float] = None,

    # Pressures
    inlet_pressure: Optional[float] = None,      # Pa (absolute)
    inlet_pressure_psi: Optional[float] = None,  # psi (absolute)

    # Pipe length
    pipe_length: Optional[float] = None,         # m
    pipe_length_ft: Optional[float] = None,      # feet

    # Temperature
    temperature_c: Optional[float] = None,       # °C

    # Gas properties
    gas_mw: Optional[float] = None,              # kg/kmol
    gas_gamma: Optional[float] = None,           # Cp/Cv
    gas_z_factor: Optional[float] = None,        # Compressibility
    gas_viscosity: Optional[float] = None,       # Pa·s
    fluid_name: Optional[str] = None,            # For property lookup

    # Pressure drop constraint
    allowable_dp_pa: Optional[float] = None,     # Pa
    allowable_dp_psi: Optional[float] = None,    # psi

    # Pipe specifications
    schedule: str = "40",
    nps_candidates: Optional[list] = None,
    material: Optional[str] = None,
    pipe_roughness: Optional[float] = None,

    # Fittings
    fittings: Optional[list] = None,

    # Options
    method: str = CANONICAL_GAS_PIPE_SIZING_METHOD,
    velocity_max_m_s: Optional[float] = None,
    return_all_candidates: bool = False
) -> str:
    """
    Select optimal pipe size for gas flow based on allowable pressure drop.

    Uses iterative equivalent length method for fittings with compressible flow.

    Args:
        flow_rate_kg_s: Mass flow rate in kg/s
        flow_rate_norm_m3_hr: Normal volumetric flow in m³/hr (0°C, 1 atm)
        flow_rate_std_m3_hr: Standard volumetric flow in m³/hr (15°C, 1 atm)
        inlet_pressure: Inlet pressure in Pa (absolute)
        inlet_pressure_psi: Inlet pressure in psi (absolute)
        pipe_length: Pipe length in m
        pipe_length_ft: Pipe length in feet
        temperature_c: Gas temperature in Celsius
        gas_mw: Molecular weight in kg/kmol
        gas_gamma: Specific heat ratio (Cp/Cv)
        gas_z_factor: Compressibility factor
        gas_viscosity: Dynamic viscosity in Pa·s
        fluid_name: Gas name for property lookup
        allowable_dp_pa: Maximum allowable pressure drop in Pa
        allowable_dp_psi: Maximum allowable pressure drop in psi
        schedule: Pipe schedule (default: "40")
        nps_candidates: List of NPS sizes to evaluate
        material: Pipe material for roughness
        pipe_roughness: Absolute roughness in m
        fittings: List of fittings with type and quantity
        method: Calculation method. ``isothermal_darcy`` is canonical; ``isothermal`` is an alias.
        velocity_max_m_s: Maximum velocity constraint
        return_all_candidates: Return all evaluated sizes

    Returns:
        JSON string with recommended size and evaluation results
    """
    try:
        log_messages = []
        warnings = []

        kwargs = {k: v for k, v in locals().items() if k != 'kwargs' and v is not None}

        # Get NPS candidates
        nps_list = kwargs.get('nps_candidates', DEFAULT_NPS_CANDIDATES)
        limit_error = _validate_nps_candidates(nps_list, log_messages)
        if limit_error is not None:
            return limit_error
        schedule = kwargs.get('schedule', '40')

        # Get temperature
        temp_c = kwargs.get('temperature_c')
        if temp_c is None:
            return json.dumps({"error": "Temperature must be specified"})
        temp_k = temp_c + 273.15
        temp_error = _validate_positive_number("temperature_k", temp_k, log_messages)
        if temp_error is not None:
            return temp_error

        # Convert inlet pressure to SI
        inlet_pressure = kwargs.get('inlet_pressure')
        if inlet_pressure is None and 'inlet_pressure_psi' in kwargs:
            inlet_pressure = kwargs['inlet_pressure_psi'] * PSI_TO_PA
            log_messages.append(f"Converted inlet_pressure from {kwargs['inlet_pressure_psi']} psi")

        if inlet_pressure is None:
            return json.dumps({"error": "Inlet pressure must be specified"})
        inlet_pressure_error = _validate_positive_number("inlet_pressure", inlet_pressure, log_messages)
        if inlet_pressure_error is not None:
            return inlet_pressure_error

        # Convert pipe length to SI
        pipe_length_m = kwargs.get('pipe_length')
        if pipe_length_m is None and 'pipe_length_ft' in kwargs:
            pipe_length_m = kwargs['pipe_length_ft'] * FT_TO_M
            log_messages.append(f"Converted pipe_length from {kwargs['pipe_length_ft']} ft")

        if pipe_length_m is None:
            return json.dumps({"error": "Pipe length must be specified"})
        pipe_length_error = _validate_positive_number("pipe_length", pipe_length_m, log_messages)
        if pipe_length_error is not None:
            return pipe_length_error

        # Get gas properties
        MW = kwargs.get('gas_mw')
        gamma = kwargs.get('gas_gamma')
        Z = kwargs.get('gas_z_factor', 1.0)
        viscosity = kwargs.get('gas_viscosity')
        fluid_name = kwargs.get('fluid_name')

        # Try property lookup if not provided
        if fluid_name is not None and (MW is None or gamma is None or viscosity is None):
            try:
                pressure_pa = inlet_pressure
                temp_k = temp_c + 273.15

                if fluid_name.lower() in ['air', 'nitrogen', 'oxygen', 'methane', 'carbondioxide', 'hydrogen']:
                    # Try CoolProp lookup
                    if fluid_name.lower() == 'air':
                        fluid_name = 'Air'
                    elif fluid_name.lower() == 'methane':
                        fluid_name = 'Methane'
                    elif fluid_name.lower() == 'nitrogen':
                        fluid_name = 'Nitrogen'

                    if MW is None:
                        MW = CP.PropsSI('M', 'T', temp_k, 'P', pressure_pa, fluid_name) * 1000  # Convert to kg/kmol
                    if viscosity is None:
                        viscosity = CP.PropsSI('V', 'T', temp_k, 'P', pressure_pa, fluid_name)
                    if gamma is None:
                        try:
                            Cp = CP.PropsSI('Cpmass', 'T', temp_k, 'P', pressure_pa, fluid_name)
                            try:
                                Cv = CP.PropsSI('Cvmass', 'T', temp_k, 'P', pressure_pa, fluid_name)
                            except Exception as e:
                                logger.debug("Cvmass lookup failed: %s, calculating from Cp-R", e)
                                # If Cvmass fails, calculate from Cp and gas constant
                                if MW is not None and MW > 0:
                                    R = 8314.46 / MW  # Specific gas constant (J/kg·K)
                                    Cv = Cp - R
                                    log_messages.append("Cv calculated from Cp-R relation")
                                else:
                                    Cv = None

                            if Cv is not None and Cv > 0:
                                gamma = Cp / Cv
                            else:
                                # Fluid-specific fallback
                                if fluid_name and fluid_name.lower() in ('methane', 'natural gas', 'biogas'):
                                    gamma = 1.3
                                    log_messages.append(f"Using methane/biogas gamma: 1.3")
                                else:
                                    gamma = 1.4
                                    log_messages.append(f"Using default gamma: 1.4")
                        except Exception as e:
                            logger.debug("Cpmass lookup failed: %s", e)
                            # If Cpmass also fails, try Cp0mass as last resort
                            try:
                                Cp0 = CP.PropsSI('Cp0mass', 'T', temp_k, 'P', pressure_pa, fluid_name)
                                if MW is not None and MW > 0:
                                    R = 8314.46 / MW
                                    Cv = Cp0 - R
                                    gamma = Cp0 / Cv if Cv > 0 else 1.4
                                    log_messages.append("Using ideal gas Cp0 with calculated Cv")
                                else:
                                    gamma = None
                            except Exception as e:
                                logger.debug("Cp0mass lookup also failed: %s", e)
                                gamma = None

                    if gamma is not None:
                        log_messages.append(f"Looked up properties for {fluid_name}")
            except Exception as e:
                log_messages.append(f"Property lookup failed: {str(e)}")

        # Use defaults if still missing
        if MW is None:
            MW = 28.97  # Air molecular weight
            log_messages.append("Using default MW for air: 28.97 kg/kmol")
        if gamma is None:
            # Better fluid-specific defaults
            if fluid_name and fluid_name.lower() in ('methane', 'natural gas', 'biogas'):
                gamma = 1.3  # More accurate for methane/biogas
                log_messages.append("Using default gamma: 1.3 (methane/biogas)")
            else:
                gamma = 1.4  # Typical for air and diatomic gases
                log_messages.append("Using default gamma: 1.4")
        if Z is None:
            Z = 1.0  # Ideal gas
            log_messages.append("Using ideal gas Z-factor: 1.0")
        mw_error = _validate_positive_number("gas_mw", MW, log_messages)
        if mw_error is not None:
            return mw_error
        gamma_error = _validate_finite_number("gas_gamma", gamma, log_messages)
        if gamma_error is not None:
            return gamma_error
        if gamma <= 1:
            return _json_error(
                "INVALID_NUMERIC_INPUT",
                "gas_gamma must be greater than 1.",
                log_messages,
                details={"gas_gamma": gamma},
            )
        z_error = _validate_positive_number("gas_z_factor", Z, log_messages)
        if z_error is not None:
            return z_error
        if viscosity is not None:
            viscosity_error = _validate_positive_number("gas_viscosity", viscosity, log_messages)
            if viscosity_error is not None:
                return viscosity_error

        # Get flow rate (convert to kg/s)
        flow_rate_kg_s = kwargs.get('flow_rate_kg_s')
        if flow_rate_kg_s is None:
            if 'flow_rate_norm_m3_hr' in kwargs:
                # Normal conditions: 0°C, 1 atm - use Z=1.0 for ideal gas at normal conditions
                rho_norm = 101325 * MW / (8314.46 * 273.15)  # Z=1 at normal conditions
                flow_rate_kg_s = kwargs['flow_rate_norm_m3_hr'] * rho_norm / 3600
                log_messages.append(f"Converted from normal m³/hr: {flow_rate_kg_s:.4f} kg/s (using Z=1 for normal conditions)")
            elif 'flow_rate_std_m3_hr' in kwargs:
                # Standard conditions: 15°C, 1 atm - use Z=1.0 for ideal gas at standard conditions
                rho_std = 101325 * MW / (8314.46 * 288.15)  # Z=1 at standard conditions
                flow_rate_kg_s = kwargs['flow_rate_std_m3_hr'] * rho_std / 3600
                log_messages.append(f"Converted from standard m³/hr: {flow_rate_kg_s:.4f} kg/s (using Z=1 for standard conditions)")

        if flow_rate_kg_s is None:
            return json.dumps({"error": "Flow rate must be specified"})
        gas_flow_error = _validate_positive_number("flow_rate_kg_s", flow_rate_kg_s, log_messages)
        if gas_flow_error is not None:
            return gas_flow_error

        # Get allowable pressure drop
        allowable_dp = kwargs.get('allowable_dp_pa')
        if allowable_dp is None and 'allowable_dp_psi' in kwargs:
            allowable_dp = kwargs['allowable_dp_psi'] * PSI_TO_PA
            log_messages.append(f"Converted allowable_dp from {kwargs['allowable_dp_psi']} psi")

        if allowable_dp is None:
            return json.dumps({"error": "Allowable pressure drop must be specified"})
        gas_dp_error = _validate_positive_number("allowable_dp", allowable_dp, log_messages)
        if gas_dp_error is not None:
            return gas_dp_error
        if allowable_dp >= inlet_pressure:
            return _json_error(
                "INVALID_NUMERIC_INPUT",
                "allowable_dp must be less than inlet_pressure.",
                log_messages,
                details={"allowable_dp": allowable_dp, "inlet_pressure": inlet_pressure},
            )

        # Get velocity constraint
        v_max = kwargs.get('velocity_max_m_s', 30)  # Default 30 m/s for gas
        gas_vmax_error = _validate_positive_number("velocity_max_m_s", v_max, log_messages)
        if gas_vmax_error is not None:
            return gas_vmax_error

        # Get pipe roughness
        pipe_roughness = kwargs.get('pipe_roughness')
        if pipe_roughness is None and 'material' in kwargs:
            pipe_roughness, source = get_pipe_roughness(material=kwargs['material'], warnings=warnings)
            log_messages.append(f"Looked up roughness for '{kwargs['material']}': {pipe_roughness} m ({source})")
        if pipe_roughness is None:
            pipe_roughness = 4.5e-5  # Default for steel
            log_messages.append(f"Using default roughness: {pipe_roughness} m")
        gas_roughness_error = _validate_nonnegative_number("pipe_roughness", pipe_roughness, log_messages)
        if gas_roughness_error is not None:
            return gas_roughness_error
        gas_nps_error = _validate_nps_candidates(nps_list, log_messages)
        if gas_nps_error is not None:
            return gas_nps_error

        # Get fittings
        fittings = kwargs.get('fittings', [])
        gas_fittings_error = _validate_fittings(fittings, log_messages)
        if gas_fittings_error is not None:
            return gas_fittings_error
        requested_method = kwargs.get('method', CANONICAL_GAS_PIPE_SIZING_METHOD)
        method = GAS_PIPE_SIZING_METHOD_ALIASES.get(requested_method)
        if method is None:
            return _json_error(
                "INVALID_INPUT",
                f"Unsupported gas pipe sizing method: {requested_method}. Supported methods: "
                f"{', '.join(sorted(SUPPORTED_GAS_PIPE_SIZING_METHODS))}.",
                log_messages,
                details={
                    "method": requested_method,
                    "supported_methods": sorted(SUPPORTED_GAS_PIPE_SIZING_METHODS),
                },
            )
        if requested_method != method:
            log_messages.append(
                f"Canonicalized gas pipe sizing method '{requested_method}' to '{method}'."
            )
        else:
            log_messages.append(f"Using gas pipe sizing method '{method}'.")

        # Evaluate each NPS candidate
        candidates = []
        recommended_nps = None

        for nps in nps_list:
            try:
                # Get pipe dimensions
                pipe_data = nearest_pipe(NPS=nps, schedule=schedule)
                # Returns tuple: (NPS, Di, Do, t)
                Di = pipe_data[1]  # Inner diameter in meters
                Do = pipe_data[2]  # Outer diameter in meters

                # Iterative calculation with fittings
                outlet_pressure = inlet_pressure - allowable_dp  # Initial guess
                L_eff = pipe_length_m  # Start with just pipe length
                convergence_log = []

                # Iteration for equivalent length
                max_iterations = 10
                tolerance = 0.001  # 0.1% tolerance

                for iteration in range(max_iterations):
                    # Calculate average conditions
                    P_avg = (inlet_pressure + outlet_pressure) / 2
                    rho_avg = P_avg * MW / (Z * 8314.46 * temp_k)

                    # Calculate velocity and Reynolds
                    area = 3.14159 * (Di/2)**2
                    Q_avg = flow_rate_kg_s / rho_avg
                    V_avg = Q_avg / area

                    if viscosity:
                        Re = Reynolds(V=V_avg, D=Di, rho=rho_avg, mu=viscosity)
                    else:
                        Re = 1e6  # Default high Reynolds for turbulent flow

                    # Calculate friction factor
                    relative_roughness = pipe_roughness / Di
                    f = friction_factor(Re=Re, eD=relative_roughness)

                    # Calculate fitting K-factors and equivalent length
                    K_total = 0
                    for fitting in fittings:
                        fitting_type = fitting.get('type')
                        quantity = fitting.get('quantity', 1)
                        K_value = fitting.get('K_value')

                        if K_value is None:
                            K_value = get_fitting_K(fitting_type, Di, Re, Q_avg)

                        K_total += K_value * quantity

                    L_eq = K_total * Di / f if f > 0 else 0
                    L_eff_new = pipe_length_m + L_eq

                    # Calculate pressure drop with effective length
                    if method == CANONICAL_GAS_PIPE_SIZING_METHOD:
                        # Use isothermal flow equation
                        from fluids import isothermal_gas
                        P2_new = isothermal_gas(
                            rho=rho_avg,
                            fd=f,
                            P1=inlet_pressure,
                            L=L_eff_new,
                            D=Di,
                            m=flow_rate_kg_s
                        )
                    else:
                        # Simple approximation for other methods
                        dp = f * (L_eff_new / Di) * (rho_avg * V_avg**2 / 2)
                        P2_new = inlet_pressure - dp

                    # Check convergence
                    if abs(L_eff_new - L_eff) / max(L_eff, 1) < tolerance:
                        convergence_log.append(f"Converged at iteration {iteration+1}")
                        L_eff = L_eff_new
                        outlet_pressure = P2_new
                        break

                    L_eff = L_eff_new
                    outlet_pressure = P2_new

                # Final calculations
                dp_total = inlet_pressure - outlet_pressure
                P_avg_final = (inlet_pressure + outlet_pressure) / 2
                rho_avg_final = P_avg_final * MW / (Z * 8314.46 * temp_k)
                V_final = flow_rate_kg_s / (rho_avg_final * area)

                # Check constraints
                velocity_ok = V_final <= v_max
                dp_ok = dp_total <= allowable_dp

                # Store candidate results
                candidate = {
                    'nps_in': nps,
                    'Di_m': round(Di, 5),
                    'Do_m': round(Do, 5),
                    'velocity_m_s': round(V_final, 2),
                    'reynolds_number': round(Re, 0),
                    'friction_factor': round(f, 6),
                    'L_pipe_m': pipe_length_m,
                    'L_equivalent_m': round(L_eq, 2),
                    'L_effective_m': round(L_eff, 2),
                    'dp_total_pa': round(dp_total, 1),
                    'dp_total_psi': round(dp_total / PSI_TO_PA, 3),
                    'outlet_pressure_pa': round(outlet_pressure, 1),
                    'pressure_ratio': round(inlet_pressure / outlet_pressure, 3),
                    'meets_dp_constraint': dp_ok,
                    'meets_velocity_constraint': velocity_ok,
                    'feasible': dp_ok and velocity_ok,
                    'iterations': len(convergence_log)
                }

                candidates.append(candidate)

                # Check if this is the recommended size
                if candidate['feasible'] and recommended_nps is None:
                    recommended_nps = nps

            except FittingKError:
                raise
            except Exception as e:
                log_messages.append(f"Error evaluating NPS {nps}: {str(e)}")
                candidates.append({
                    'nps_in': nps,
                    'error': str(e)
                })

        # Prepare results
        candidate_failure = _candidate_failure_error(candidates, log_messages)
        if candidate_failure is not None:
            return candidate_failure

        result = {
            'inputs_resolved': log_messages,
            'gas_properties': {
                'MW_kg_kmol': MW,
                'gamma': gamma,
                'Z_factor': Z,
                'viscosity_pas': viscosity
            },
            'allowable_dp_pa': allowable_dp,
            'allowable_dp_psi': allowable_dp / PSI_TO_PA,
            'flow_rate_kg_s': flow_rate_kg_s,
            'inlet_pressure_pa': inlet_pressure,
            'pipe_length_m': pipe_length_m,
            'temperature_c': temp_c,
            'schedule': schedule,
            'method': method
        }
        if warnings:
            result['warnings'] = warnings

        if recommended_nps is not None:
            result['recommended_nps_in'] = recommended_nps
            # Find recommended details
            for c in candidates:
                if c.get('nps_in') == recommended_nps:
                    result['recommended_details'] = {
                        'Di_m': c['Di_m'],
                        'velocity_m_s': c['velocity_m_s'],
                        'dp_total_pa': c['dp_total_pa'],
                        'dp_total_psi': c['dp_total_psi'],
                        'L_equivalent_m': c.get('L_equivalent_m', 0)
                    }
                    break
        else:
            result['recommended_nps_in'] = None
            result['note'] = "No pipe size meets the specified constraints"

        if kwargs.get('return_all_candidates', False):
            result['all_candidates'] = candidates
        else:
            result['feasible_candidates'] = [c for c in candidates if c.get('feasible', False)]

        return json.dumps(result, indent=2, allow_nan=False)

    except FittingKError as e:
        logger.warning("Fitting K-value resolution failed: %s", e.message)
        return _json_error(e.code, e.message, log_messages, details=e.details)
    except Exception as e:
        logger.error(f"Error in select_gas_pipe_size_for_dp: {str(e)}")
        return _json_error("PIPE_SIZING_CALCULATION_FAILED", str(e), [])
