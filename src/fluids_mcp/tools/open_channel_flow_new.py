import json
import logging
import math
from typing import Optional, Literal, Dict, Any, Tuple, Union, List

# Scientific computing
import scipy.optimize

# Import from pydraulics (now properly packaged)
try:
    from fluids_mcp.pydraulics.open_flow import (
        RectangularChannel,
        TrapezoidalChannel,
        CircularChannel,
        TriangularChannel
    )
    PYDRAULICS_AVAILABLE = True

except ImportError as import_err:
    # Define dummy classes if pydraulics is not available
    class MockChannel:
        def __init__(self, *args, **kwargs): pass
        def area(self, y): return None
        def hydraulic_radius(self, y): return None
        def calc_properties(self): pass
        a = p = rh = tw = dh = y = None

    RectangularChannel = TrapezoidalChannel = CircularChannel = TriangularChannel = MockChannel
    PYDRAULICS_AVAILABLE = False
    logging.error(f"Pydraulics library not found: {import_err}")

# Try to import fluids library
try:
    import fluids
    import fluids.open_flow
    import fluids.geometry
    FLUIDS_AVAILABLE = True
except ImportError:
    fluids = None
    FLUIDS_AVAILABLE = False
    logging.error("Fluids library not found. Open channel calculations will fail.")

# Import shared utilities
from fluids_mcp.utils.constants import FT_to_M, G_GRAVITY
from fluids_mcp.errors import CalculationError
from fluids_mcp.utils.helpers import get_pipe_roughness, resolve_fitting_k_basis

# Configure logging
logger = logging.getLogger("fluids-mcp.open_channel_flow")


def _open_channel_error(code: str, message: str, *, details: Optional[dict] = None) -> str:
    error = CalculationError(code, message, details=details)
    return json.dumps(
        {
            "error_code": error.code,
            "error": error.message,
            "errors": [error.to_public_dict()],
            "log": [],
        },
        allow_nan=False,
    )


def _calculate_gravity_leg(
    *,
    flow_rate_m3s: Optional[float],
    diameter_m: Optional[float],
    pipe_length: Optional[float],
    pipe_length_ft: Optional[float],
    upstream_head_m: Optional[float],
    downstream_head_m: Optional[float],
    fluid_density: Optional[float],
    fluid_viscosity: Optional[float],
    pipe_roughness: Optional[float],
    material: Optional[str],
    fittings: Optional[list],
) -> str:
    """Evaluate a full-bore tank-to-tank gravity leg.

    Available fall is the upstream minus downstream tank total head. Required
    fall is Darcy-Weisbach straight-pipe loss plus Crane TP-410 minor losses;
    tank-surface endpoint velocity heads are zero. A negative residual means
    the specified full-bore flow cannot be conveyed by the available fall and
    would require additional head or a different diameter.
    """

    required_inputs = {
        "flow_rate_m3s": flow_rate_m3s,
        "diameter_m": diameter_m,
        "upstream_head_m": upstream_head_m,
        "downstream_head_m": downstream_head_m,
        "fluid_density": fluid_density,
        "fluid_viscosity": fluid_viscosity,
    }
    missing = [name for name, value in required_inputs.items() if value is None]
    if pipe_length is None and pipe_length_ft is None:
        missing.append("pipe_length or pipe_length_ft")
    if missing:
        return _open_channel_error(
            "MISSING_GRAVITY_LEG_INPUT",
            f"Missing required gravity-leg inputs: {', '.join(missing)}.",
            details={"missing_fields": missing},
        )
    if pipe_length is not None and pipe_length_ft is not None:
        return _open_channel_error(
            "CONFLICTING_LENGTH_INPUTS",
            "Provide gravity-leg pipe length in metres or feet, not both.",
        )

    values = {
        "flow_rate_m3s": flow_rate_m3s,
        "diameter_m": diameter_m,
        "fluid_density": fluid_density,
        "fluid_viscosity": fluid_viscosity,
    }
    invalid = {
        name: value
        for name, value in values.items()
        if not isinstance(value, (int, float)) or not math.isfinite(float(value)) or value <= 0
    }
    if invalid:
        return _open_channel_error(
            "INVALID_GRAVITY_LEG_INPUT",
            "Flow, diameter, density, and viscosity must be positive and finite.",
            details={"invalid_fields": invalid},
        )
    length_m = float(pipe_length) if pipe_length is not None else float(pipe_length_ft) * FT_to_M
    if not math.isfinite(length_m) or length_m < 0:
        return _open_channel_error(
            "INVALID_GRAVITY_LEG_INPUT",
            "Gravity-leg pipe length must be non-negative and finite.",
            details={"pipe_length_m": length_m},
        )
    if not all(
        isinstance(value, (int, float)) and math.isfinite(float(value))
        for value in (upstream_head_m, downstream_head_m)
    ):
        return _open_channel_error(
            "INVALID_GRAVITY_LEG_INPUT",
            "Upstream and downstream heads must be finite.",
        )

    warnings: list[str] = []
    roughness_m, roughness_source = get_pipe_roughness(
        material=material,
        pipe_roughness=pipe_roughness,
        warnings=warnings,
    )
    area_m2 = math.pi * float(diameter_m) ** 2 / 4.0
    velocity_m_s = float(flow_rate_m3s) / area_m2
    reynolds_number = float(fluid_density) * velocity_m_s * float(diameter_m) / float(fluid_viscosity)
    darcy_friction_factor = fluids.friction.friction_factor(
        Re=reynolds_number,
        eD=roughness_m / float(diameter_m),
    )
    velocity_head_m = velocity_m_s**2 / (2.0 * G_GRAVITY)
    straight_pipe_head_loss_m = (
        darcy_friction_factor * length_m / float(diameter_m) * velocity_head_m
    )

    total_k = 0.0
    fitting_details = []
    for fitting in fittings or []:
        quantity = float(fitting.get("quantity", 1))
        if not math.isfinite(quantity) or quantity < 0:
            return _open_channel_error(
                "INVALID_FITTING_INPUT",
                "Fitting quantity must be non-negative and finite.",
                details={"quantity": fitting.get("quantity")},
            )
        try:
            k_value, detail = resolve_fitting_k_basis(
                fitting,
                diameter=float(diameter_m),
                reynolds_number=reynolds_number,
                flow_rate=float(flow_rate_m3s),
                darcy_friction_factor=darcy_friction_factor,
            )
        except CalculationError as exc:
            return json.dumps(
                {
                    "error_code": exc.code,
                    "error": exc.message,
                    "errors": [exc.to_public_dict()],
                    "log": [],
                },
                allow_nan=False,
            )
        total_k += k_value * quantity
        detail.update({"quantity": quantity, "K_total": k_value * quantity})
        fitting_details.append(detail)

    minor_head_loss_m = total_k * velocity_head_m
    required_fall_m = straight_pipe_head_loss_m + minor_head_loss_m
    available_fall_m = float(upstream_head_m) - float(downstream_head_m)
    residual_fall_m = available_fall_m - required_fall_m
    result = {
        "status": "Success",
        "open_channel_mode": "gravity_leg",
        "flow_state": "full_bore_surcharged",
        "flow_rate_m3s": float(flow_rate_m3s),
        "pipe_diameter_m": float(diameter_m),
        "pipe_length_m": length_m,
        "pipe_velocity_m_s": velocity_m_s,
        "reynolds_number": reynolds_number,
        "darcy_friction_factor": darcy_friction_factor,
        "roughness_m": roughness_m,
        "roughness_source": roughness_source,
        "pipe_material": material,
        "available_fall_m": available_fall_m,
        "required_fall_m": required_fall_m,
        "residual_fall_m": residual_fall_m,
        "gravity_capacity_pass": residual_fall_m >= 0,
        "straight_pipe_head_loss_m": straight_pipe_head_loss_m,
        "minor_head_loss_m": minor_head_loss_m,
        "total_fitting_k": total_k,
        "fitting_details": fitting_details,
        "endpoint_basis": "tank_surface",
        "governing_basis": "Darcy-Weisbach/Colebrook plus Crane TP-410 minor losses",
        "warnings": warnings,
        "log": [
            "Available fall = upstream tank head - downstream tank head.",
            "Required fall = Darcy straight-pipe loss + minor losses; tank endpoint velocity heads are zero.",
        ],
    }
    return json.dumps(result, allow_nan=False)

def calculate_open_channel_flow(
    # --- Channel Geometry ---
    channel_type: Literal["rectangular", "trapezoidal", "circular", "triangular"] = "circular",
    width_m: Optional[float] = None,       # Bottom width (rect/trap)
    diameter_m: Optional[float] = None,    # Diameter (circular)
    side_slope_z: Optional[float] = None,  # Side slope z:1 (trap/triangular)

    # --- Flow Conditions ---
    slope: Optional[float] = None,         # Channel bottom slope (dimensionless)
    manning_n: Optional[float] = None,     # Manning's n

    # --- Specify EITHER flow_rate OR normal_depth ---
    flow_rate_m3s: Optional[float] = None, # Flow rate in m³/s
    normal_depth_m: Optional[float] = None, # Normal flow depth in m

    # --- HGL / branch controls ---
    open_channel_mode: Literal["normal_flow", "gravity_leg"] = "normal_flow",
    circular_depth_branch: Literal["shallow", "near_full"] = "shallow",
    invert_elevation_m: Optional[float] = None,

    # --- Full-bore tank-to-tank gravity leg ---
    pipe_length: Optional[float] = None,
    pipe_length_ft: Optional[float] = None,
    upstream_head_m: Optional[float] = None,
    downstream_head_m: Optional[float] = None,
    fluid_density: Optional[float] = None,
    fluid_viscosity: Optional[float] = None,
    pipe_roughness: Optional[float] = None,
    material: Optional[str] = None,
    fittings: Optional[list] = None,

) -> str:
    """Calculate Manning normal flow/HGL or check a full-bore gravity leg.

    Open-channel geometry and Manning capacity follow ASCE MOP 60. Critical
    depth is solved from ``Q²T/(gA³)=1`` rather than the pydraulics WIP ``yc``
    implementation, and specific energy is ``E=y+V²/(2g)``. Circular normal
    depth exposes shallow and near-full roots around the non-monotonic maximum
    capacity near 0.94D. ``gravity_leg`` uses Darcy-Weisbach/Colebrook and
    Crane TP-410 losses between tank surfaces.

    Parameters
    ----------
    channel_type : str
        Type of channel: "rectangular", "trapezoidal", "circular", or "triangular"
    width_m : float, optional
        Bottom width for rectangular or trapezoidal channels in meters
    diameter_m : float, optional
        Diameter for circular channels in meters
    side_slope_z : float, optional
        Side slope ratio (horizontal:vertical) for trapezoidal and triangular channels
    slope : float
        Channel bottom slope (m/m)
    manning_n : float
        Manning's roughness coefficient
    flow_rate_m3s : float, optional
        Flow rate in m³/s (specify this OR normal_depth_m, not both)
    normal_depth_m : float, optional
        Normal flow depth in m (specify this OR flow_rate_m3s, not both)

    Returns
    -------
    str
        JSON string with calculation results or error messages
    """
    results_log: List[str] = []
    error_log: List[str] = []
    warnings: List[str] = []
    calculated_results: Dict[str, Any] = {}
    geometry_params: Dict[str, float] = {}
    g = G_GRAVITY  # m/s² (standard gravity)

    if open_channel_mode == "gravity_leg":
        return _calculate_gravity_leg(
            flow_rate_m3s=flow_rate_m3s,
            diameter_m=diameter_m,
            pipe_length=pipe_length,
            pipe_length_ft=pipe_length_ft,
            upstream_head_m=upstream_head_m,
            downstream_head_m=downstream_head_m,
            fluid_density=fluid_density,
            fluid_viscosity=fluid_viscosity,
            pipe_roughness=pipe_roughness,
            material=material,
            fittings=fittings,
        )
    if open_channel_mode != "normal_flow":
        return _open_channel_error(
            "INVALID_OPEN_CHANNEL_MODE",
            "open_channel_mode must be one of: normal_flow, gravity_leg.",
            details={"open_channel_mode": open_channel_mode},
        )
    if circular_depth_branch not in ("shallow", "near_full"):
        return _open_channel_error(
            "INVALID_CIRCULAR_DEPTH_BRANCH",
            "circular_depth_branch must be one of: shallow, near_full.",
            details={"circular_depth_branch": circular_depth_branch},
        )

    # --- Validate basic inputs ---
    if channel_type not in ["rectangular", "trapezoidal", "circular", "triangular"]:
        error_log.append("Invalid channel_type. Must be one of 'rectangular', 'trapezoidal', 'circular', 'triangular'.")

    if slope is None or slope <= 0:
        error_log.append("Missing or non-positive channel slope required.")
    if manning_n is None or manning_n <= 0:
        error_log.append("Missing or non-positive Manning's n required.")

    # Validate geometry
    if channel_type == "rectangular":
        if width_m is None or width_m <= 0:
            error_log.append("Missing or non-positive width_m for rectangular channel.")
        else:
            geometry_params['b'] = width_m

    elif channel_type == "trapezoidal":
        if width_m is None or width_m <= 0:
            error_log.append("Missing or non-positive bottom width_m for trapezoidal channel.")
        else:
            geometry_params['b'] = width_m
        if side_slope_z is None or side_slope_z < 0:
            error_log.append("Missing or negative side_slope_z for trapezoidal channel.")
        else:
            geometry_params['z'] = side_slope_z

    elif channel_type == "circular":
        if diameter_m is None or diameter_m <= 0:
            error_log.append("Missing or non-positive diameter_m for circular channel.")
        else:
            geometry_params['D'] = diameter_m

    elif channel_type == "triangular":
        if side_slope_z is None or side_slope_z < 0:
            error_log.append("Missing or negative side_slope_z for triangular channel.")
        else:
            geometry_params['z'] = side_slope_z

    # Check which variable to solve for
    solve_for = None
    if flow_rate_m3s is not None and normal_depth_m is None:
        solve_for = "depth"
        if flow_rate_m3s < 0:
            error_log.append("flow_rate_m3s cannot be negative.")
        calculated_results['target_flow_rate_m3s'] = flow_rate_m3s

    elif normal_depth_m is not None and flow_rate_m3s is None:
        solve_for = "flow_rate"
        if normal_depth_m <= 0:
            error_log.append("normal_depth_m must be positive.")
        if channel_type == "circular" and diameter_m and normal_depth_m > diameter_m:
            error_log.append("normal_depth_m cannot exceed diameter_m for circular channel.")
        calculated_results['normal_depth_m'] = normal_depth_m

    else:
        error_log.append("Provide exactly one of 'flow_rate_m3s' or 'normal_depth_m'.")

    if error_log:
        return json.dumps({"errors": error_log, "log": results_log})

    # --- Create channel object and calculate flow parameters ---
    def calculate_geometry(y: float, chan_type: str, geom_params: Dict[str, float]) -> Tuple[Optional[float], ...]:
        """Calculate channel geometry parameters based on channel type and depth.

        Returns (A, P, Rh, T, Dh) where:
        - A: Cross-sectional flow area (m²)
        - P: Wetted perimeter (m)
        - Rh: Hydraulic radius (m)
        - T: Top width (m)
        - Dh: Hydraulic depth (m)
        """
        # Ensure y is a valid positive number for calculation
        if not isinstance(y, (int, float)) or y <= 1e-9:
            return None, None, None, None, None

        # Use pydraulics for channel geometry calculations if available
        channel = None
        try:
            if chan_type == "rectangular":
                if 'b' not in geom_params:
                    return None, None, None, None, None
                b = geom_params['b']

                # Direct calculation as fallback
                A = b * y
                P = b + 2 * y
                Rh = A / P if P > 0 else None
                T = b
                Dh = A / T if T > 0 else None

                # Try using pydraulics for better accuracy
                if PYDRAULICS_AVAILABLE:
                    try:
                        channel = RectangularChannel(n=1.0, So=0.001, b=b, Q=None, y=y)
                        channel.calc_properties()
                        if channel.a is not None and channel.p is not None:
                            A = float(channel.a)
                            P = float(channel.p)
                            Rh = float(channel.rh) if channel.rh is not None else (A/P if P > 0 else None)
                            T = float(channel.tw) if channel.tw is not None else b
                            Dh = float(channel.dh) if channel.dh is not None else y
                    except Exception as e:
                        logger.debug("Rectangular channel calc failed: %s", e)
                        pass  # Fall back to direct calculations

            elif chan_type == "trapezoidal":
                if 'b' not in geom_params or 'z' not in geom_params:
                    return None, None, None, None, None
                b = geom_params['b']
                z = geom_params['z']

                # Direct calculation as fallback
                A = (b + z * y) * y
                P = b + 2 * y * (1 + z**2)**0.5
                Rh = A / P if P > 0 else None
                T = b + 2 * z * y
                Dh = A / T if T > 0 else None

                # Try using pydraulics for better accuracy
                if PYDRAULICS_AVAILABLE:
                    try:
                        channel = TrapezoidalChannel(n=1.0, So=0.001, z=z, b=b, Q=None, y=y)
                        channel.calc_properties()
                        if channel.a is not None and channel.p is not None:
                            A = float(channel.a)
                            P = float(channel.p)
                            Rh = float(channel.rh) if channel.rh is not None else (A/P if P > 0 else None)
                            T = float(channel.tw) if channel.tw is not None else (b + 2*z*y)
                            Dh = float(channel.dh) if channel.dh is not None else (A/T if T > 0 else None)
                    except Exception as e:
                        logger.debug("Trapezoidal channel calc failed: %s", e)
                        pass  # Fall back to direct calculations

            elif chan_type == "circular":
                if 'D' not in geom_params:
                    return None, None, None, None, None
                D = geom_params['D']

                # Limit depth to pipe diameter
                y_eff = min(y, D)

                # Direct calculation as fallback
                if y_eff >= D:  # Full pipe
                    A = math.pi * (D/2)**2
                    P = math.pi * D
                    Rh = D/4
                    T = 0  # No free surface
                    Dh = D/4
                else:  # Partially full pipe
                    theta = 2 * math.acos(1 - 2*y_eff/D)  # Central angle
                    A = (D**2/8) * (theta - math.sin(theta))
                    P = (D/2) * theta
                    Rh = A / P if P > 0 else None
                    T = D * math.sin(theta/2)
                    Dh = A / T if T > 0 else None

                # Try using pydraulics for better accuracy
                if PYDRAULICS_AVAILABLE:
                    try:
                        channel = CircularChannel(n=1.0, So=0.001, D=D, Q=None, y=y_eff)
                        channel.calc_properties()
                        if channel.a is not None and channel.p is not None:
                            A = float(channel.a)
                            P = float(channel.p)
                            Rh = float(channel.rh) if channel.rh is not None else (A/P if P > 0 else None)
                            T = float(channel.tw) if channel.tw is not None else T
                            Dh = float(channel.dh) if channel.dh is not None else (A/T if T > 0 else None)
                    except Exception as e:
                        logger.debug("Circular channel calc failed: %s", e)
                        pass  # Fall back to direct calculations

            elif chan_type == "triangular":
                if 'z' not in geom_params:
                    return None, None, None, None, None
                z = geom_params['z']

                # Direct calculation as fallback
                A = z * y**2
                P = 2 * y * (1 + z**2)**0.5
                Rh = A / P if P > 0 else None
                T = 2 * z * y
                Dh = A / T if T > 0 else None

                # Try using pydraulics for better accuracy
                if PYDRAULICS_AVAILABLE:
                    try:
                        channel = TriangularChannel(n=1.0, So=0.001, z=z, Q=None, y=y)
                        channel.calc_properties()
                        if channel.a is not None and channel.p is not None:
                            A = float(channel.a)
                            P = float(channel.p)
                            Rh = float(channel.rh) if channel.rh is not None else (A/P if P > 0 else None)
                            T = float(channel.tw) if channel.tw is not None else (2*z*y)
                            Dh = float(channel.dh) if channel.dh is not None else (A/T if T > 0 else None)
                    except Exception as e:
                        logger.debug("Triangular channel calc failed: %s", e)
                        pass  # Fall back to direct calculations
            else:
                return None, None, None, None, None

            return A, P, Rh, T, Dh

        except Exception as e:
            logger.error(f"Error in calculate_geometry for {chan_type}: {e}")
            return None, None, None, None, None

    def calculate_q_from_y(y: float, chan_type: str, geom_params: Dict[str, float], S: float, n: float) -> float:
        """Calculate flow rate for a given depth using Manning's equation."""
        # Calculate geometry parameters
        A, P, Rh, T, Dh = calculate_geometry(y, chan_type, geom_params)

        # Check for valid geometry
        if A is None or Rh is None or A < 1e-12 or Rh < 1e-12:
            return 0.0

        # Calculate velocity using Manning's equation
        if FLUIDS_AVAILABLE:
            try:
                V = fluids.open_flow.V_Manning(Rh=Rh, S=S, n=n)
                return V * A
            except Exception as e:
                logger.debug("Manning velocity calc failed: %s", e)
                pass

        # Direct calculation as fallback
        V = (1.0/n) * (Rh**(2/3)) * (S**0.5)
        return V * A

    # Objective function for depth calculation
    def manning_error_func(y: float, Q_target: float, chan_type: str, geom_params: Dict[str, float], S: float, n: float) -> float:
        """Calculate difference between target flow rate and calculated flow rate."""
        Q_calc = calculate_q_from_y(y, chan_type, geom_params, S, n)
        return Q_calc - Q_target

    # Robust root finding using scipy.optimize
    def find_normal_depth(error_func, args: Tuple, channel_type: str, geometry_params: Dict[str, float]) -> float:
        """Find normal depth with explicit circular shallow/near-full branches."""

        # Establish reasonable bounds for different channel types
        y_min = 1e-6  # Very small positive depth

        if channel_type == "circular" and 'D' in geometry_params:
            diameter = geometry_params['D']
            capacity_result = scipy.optimize.minimize_scalar(
                lambda depth: -calculate_q_from_y(
                    depth,
                    channel_type,
                    geometry_params,
                    slope,
                    manning_n,
                ),
                bounds=(diameter * 1e-6, diameter),
                method="bounded",
                options={"xatol": diameter * 1e-12},
            )
            peak_depth = float(capacity_result.x)
            peak_capacity = calculate_q_from_y(
                peak_depth, channel_type, geometry_params, slope, manning_n
            )
            full_capacity = calculate_q_from_y(
                diameter, channel_type, geometry_params, slope, manning_n
            )
            target_flow = float(args[0])
            calculated_results.update(
                {
                    "circular_depth_branch": circular_depth_branch,
                    "manning_capacity_max_m3s": peak_capacity,
                    "manning_capacity_full_m3s": full_capacity,
                    "capacity_peak_depth_m": peak_depth,
                    "capacity_peak_depth_ratio": peak_depth / diameter,
                }
            )
            capacity_tolerance = max(peak_capacity, 1.0) * 1e-10
            if target_flow > peak_capacity + capacity_tolerance:
                raise RuntimeError(
                    f"Target flow {target_flow:.6g} m3/s exceeds circular Manning maximum "
                    f"{peak_capacity:.6g} m3/s at y/D={peak_depth/diameter:.4f}."
                )

            if circular_depth_branch == "near_full":
                if target_flow < full_capacity - capacity_tolerance:
                    raise RuntimeError(
                        f"Near-full branch is unavailable below the full-pipe Manning capacity "
                        f"{full_capacity:.6g} m3/s; select circular_depth_branch='shallow'."
                    )
                lower = peak_depth
                upper = diameter
            else:
                lower = y_min
                upper = peak_depth

            f_lower = error_func(lower, *args)
            f_upper = error_func(upper, *args)
            if abs(f_lower) <= capacity_tolerance:
                return lower
            if abs(f_upper) <= capacity_tolerance:
                return upper
            if f_lower * f_upper > 0:
                raise RuntimeError(
                    f"Target flow is not bracketed on circular {circular_depth_branch} branch."
                )
            return float(
                scipy.optimize.brentq(
                    error_func,
                    lower,
                    upper,
                    args=args,
                    xtol=1e-10,
                    rtol=1e-10,
                    maxiter=100,
                )
            )
        else:
            y_max = 10.0  # Reasonable maximum depth for other channels

        # Try to establish brackets by testing error function at bounds
        try:
            f_min = error_func(y_min, *args)
            f_max = error_func(y_max, *args)

            # Check if we have proper bracketing (opposite signs)
            if f_min * f_max < 0:
                # Use brentq for robust bracketed root finding
                y_solution = scipy.optimize.brentq(
                    error_func, y_min, y_max, args=args,
                    xtol=1e-8, rtol=1e-8, maxiter=100
                )
                return y_solution
            else:
                # No bracketing found, use fsolve as fallback
                logger.info("No bracketing found, using fsolve as fallback")
                # Use reasonable initial guess
                if channel_type == "circular" and 'D' in geometry_params:
                    initial_guess = geometry_params['D'] * 0.5  # Half full
                else:
                    initial_guess = 0.5  # 0.5 m depth

                result = scipy.optimize.fsolve(
                    error_func, initial_guess, args=args,
                    xtol=1e-8, full_output=True
                )

                if result[2] == 1:  # Successful convergence
                    return float(result[0][0])
                else:
                    raise RuntimeError("fsolve failed to converge")

        except Exception as e:
            logger.error(f"Both brentq and fsolve failed: {e}")
            raise RuntimeError(f"Root finding failed: {e}")

    def find_critical_depth(
        flow_rate: float,
        chan_type: str,
        geom_params: Dict[str, float],
    ) -> float:
        """Solve Q²T/(gA³)=1, the exact arbitrary-section critical condition."""

        if flow_rate <= 0:
            return 0.0

        def critical_residual(depth: float) -> float:
            area, _perimeter, _radius, top_width, _hydraulic_depth = calculate_geometry(
                depth, chan_type, geom_params
            )
            if area is None or top_width is None or area <= 0 or top_width < 0:
                return math.inf
            return flow_rate**2 * top_width / (g * area**3) - 1.0

        lower = 1e-8
        if chan_type == "circular":
            upper = geom_params["D"] * (1.0 - 1e-10)
        else:
            upper = max(1.0, float(normal_depth_m or 0.0), float(calculated_results.get("normal_depth_m", 0.0)))
            while critical_residual(upper) > 0 and upper < 1e6:
                upper *= 2.0

        f_lower = critical_residual(lower)
        f_upper = critical_residual(upper)
        if not math.isfinite(f_lower) or f_lower <= 0 or f_upper >= 0:
            raise RuntimeError("Critical-depth equation could not be bracketed for the resolved section and flow.")
        return float(
            scipy.optimize.brentq(
                critical_residual,
                lower,
                upper,
                xtol=1e-10,
                rtol=1e-10,
                maxiter=100,
            )
        )

    # Main calculation
    try:
        y = None
        Q = None

        if solve_for == "depth":
            Q_target = flow_rate_m3s
            results_log.append(f"Solving for normal depth, target Q={Q_target:.4f} m³/s.")

            # Solve for normal depth using robust scipy methods
            args = (Q_target, channel_type, geometry_params, slope, manning_n)
            try:
                y_sol = find_normal_depth(manning_error_func, args, channel_type, geometry_params)
            except RuntimeError as solve_error:
                error_log.append(f"Solver failed: {solve_error}")
                return json.dumps({"errors": error_log, "log": results_log})

            if y_sol <= 0:
                error_log.append("Solver failed to find a positive depth.")
                return json.dumps({"errors": error_log, "log": results_log})

            y = y_sol
            Q = calculate_q_from_y(y, channel_type, geometry_params, slope, manning_n)
            calculated_results['normal_depth_m'] = y
            calculated_results['flow_rate_m3s'] = Q
            results_log.append(f"Solver converged: y={y:.6f} m => Q={Q:.6f} m³/s.")

        elif solve_for == "flow_rate":
            y = normal_depth_m
            results_log.append(f"Calculating flow rate at y={y:.4f} m.")
            Q = calculate_q_from_y(y, channel_type, geometry_params, slope, manning_n)
            calculated_results['flow_rate_m3s'] = Q
            results_log.append(f"Flow rate calculated: Q={Q:.6f} m³/s at y={y:.4f} m.")

        # Calculate additional parameters
        if y is not None and Q is not None:
            # Get geometry parameters
            A, P, Rh, T, Dh = calculate_geometry(y, channel_type, geometry_params)

            # Store geometry parameters
            if A is not None:
                calculated_results['flow_area_m2'] = A
            if P is not None:
                calculated_results['wetted_perimeter_m'] = P
            if Rh is not None:
                calculated_results['hydraulic_radius_m'] = Rh
            if T is not None:
                calculated_results['top_width_m'] = T
            if Dh is not None:
                calculated_results['hydraulic_depth_m'] = Dh

            # Calculate HGL, specific energy, critical depth, and Froude regime.
            if A is not None and A > 1e-12:
                V = Q / A
                calculated_results['velocity_m_s'] = V
                specific_energy_m = y + V**2 / (2.0 * g)
                calculated_results['specific_energy_m'] = specific_energy_m
                calculated_results['hydraulic_grade_above_invert_m'] = y
                calculated_results['energy_grade_above_invert_m'] = specific_energy_m
                if invert_elevation_m is not None:
                    calculated_results['hgl_elevation_m'] = invert_elevation_m + y
                    calculated_results['energy_grade_elevation_m'] = invert_elevation_m + specific_energy_m

                is_full_circular = (
                    channel_type == "circular"
                    and diameter_m is not None
                    and y >= diameter_m * (1.0 - 1e-9)
                )
                if channel_type == "circular" and diameter_m is not None:
                    depth_ratio = y / diameter_m
                    calculated_results['depth_ratio_y_over_d'] = depth_ratio
                    if "manning_capacity_max_m3s" not in calculated_results:
                        capacity_result = scipy.optimize.minimize_scalar(
                            lambda depth: -calculate_q_from_y(
                                depth, channel_type, geometry_params, slope, manning_n
                            ),
                            bounds=(diameter_m * 1e-6, diameter_m),
                            method="bounded",
                        )
                        peak_depth = float(capacity_result.x)
                        calculated_results.update(
                            {
                                "manning_capacity_max_m3s": calculate_q_from_y(
                                    peak_depth, channel_type, geometry_params, slope, manning_n
                                ),
                                "manning_capacity_full_m3s": calculate_q_from_y(
                                    diameter_m, channel_type, geometry_params, slope, manning_n
                                ),
                                "capacity_peak_depth_m": peak_depth,
                                "capacity_peak_depth_ratio": peak_depth / diameter_m,
                            }
                        )
                    if is_full_circular:
                        calculated_results['flow_state'] = "full_bore"
                        calculated_results['near_full_region'] = True
                    elif depth_ratio >= 0.94:
                        calculated_results['flow_state'] = "near_full_open_channel"
                        calculated_results['near_full_region'] = True
                        branch_label = (
                            circular_depth_branch if solve_for == "depth" else "provided-depth"
                        )
                        warnings.append(
                            "Circular depth is in the near-full non-monotonic Manning region; "
                            f"reported basis is {branch_label}."
                        )
                    else:
                        calculated_results['flow_state'] = "partially_full"
                        calculated_results['near_full_region'] = False
                else:
                    calculated_results['flow_state'] = "open_channel"
                    calculated_results['near_full_region'] = False

                if is_full_circular:
                    calculated_results['froude_number'] = None
                    calculated_results['flow_regime'] = "Full-bore/pressurized; Froude not applicable"
                elif T is not None and T > 1e-12:
                    try:
                        hydraulic_depth = A / T
                        calculated_results['hydraulic_depth_m'] = hydraulic_depth
                        Fr = V / math.sqrt(g * hydraulic_depth)
                        calculated_results['froude_number'] = Fr

                        # Determine flow regime
                        if Fr < 0.99:
                            flow_regime = "Subcritical"
                        elif Fr > 1.01:
                            flow_regime = "Supercritical"
                        else:
                            flow_regime = "Critical"
                        calculated_results['flow_regime'] = flow_regime
                    except Exception as fr_err:
                        logger.error(f"Error calculating Froude number: {fr_err}")
                        calculated_results['froude_number'] = None
                        calculated_results['flow_regime'] = "Unknown"
                else:
                    calculated_results['froude_number'] = None
                    calculated_results['flow_regime'] = "Unknown (Invalid hydraulic depth)"

                critical_depth = find_critical_depth(Q, channel_type, geometry_params)
                critical_area, _cp, _crh, _ct, _cdh = calculate_geometry(
                    critical_depth, channel_type, geometry_params
                )
                critical_velocity = Q / critical_area
                calculated_results['critical_depth_m'] = critical_depth
                calculated_results['critical_velocity_m_s'] = critical_velocity
                calculated_results['critical_specific_energy_m'] = (
                    critical_depth + critical_velocity**2 / (2.0 * g)
                )
            else:
                calculated_results['velocity_m_s'] = 0.0
                calculated_results['froude_number'] = None
                calculated_results['flow_regime'] = "Unknown (Invalid area)"

    except Exception as e:
        logger.error(f"Error during open channel calculation: {e}", exc_info=True)
        error_log.append(f"Calculation error: {str(e)}")
        return json.dumps({"errors": error_log, "log": results_log})

    # Format and return results
    final_result = {
        "status": "Success" if not error_log else "Failure",
        "open_channel_mode": "normal_flow",
        "log": results_log,
        "warnings": warnings,
        "errors": error_log if error_log else None,
        "inputs": {
            "channel_type": channel_type,
            "geometry_params": geometry_params,
            "slope": slope,
            "manning_n": manning_n,
            "provided_flow_rate_m3s": flow_rate_m3s,
            "provided_normal_depth_m": normal_depth_m,
            "circular_depth_branch": circular_depth_branch,
            "invert_elevation_m": invert_elevation_m,
        },
        "results": {
            k: round(v, 6) if isinstance(v, float) else v
            for k, v in calculated_results.items()
        } if calculated_results else None
    }

    # Clean up null/empty fields
    if not final_result["errors"]:
        del final_result["errors"]
    if not final_result["results"]:
        del final_result["results"]

    # Preserve the legacy nested ``results`` object while exposing HGL and
    # critical-flow fields directly through machine_requirements.
    final_result.update(
        {
            key: round(value, 6) if isinstance(value, float) else value
            for key, value in calculated_results.items()
        }
    )

    return json.dumps(final_result, allow_nan=False)
