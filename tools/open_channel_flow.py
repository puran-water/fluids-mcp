import json
import logging
import math
import sys  # <-- Import sys
import os   # <-- Import os
from typing import Optional, Literal

# --- Third-Party Libraries ---

# Fluids library (no change needed for this)
try:
    import fluids
    import fluids.open_flow
    import fluids.geometry
    FLUIDS_AVAILABLE = True
except ImportError:
    fluids = None
    FLUIDS_AVAILABLE = False
    logging.error("Fluids library not found. Open channel calculations will fail.")

# --- Add parent directory to sys.path to find pydraulics ---
try:
    # Get the directory containing the current script (__file__)
    current_script_dir = os.path.dirname(os.path.abspath(__file__))
    # Get the parent directory (one level up)
    parent_dir = os.path.dirname(current_script_dir)

    # Add the parent directory to the Python path if it's not already there
    if parent_dir not in sys.path:
        sys.path.insert(0, parent_dir) # Insert at the beginning

    # --- Now import from pydraulics (located in the parent directory) ---
    # Import from open_flow.py instead of the non-existent channel.py
    from pydraulics.open_flow import (
        RectangularChannel,
        TrapezoidalChannel,
        CircularChannel,
        TriangularChannel
    )
    PYDRAULICS_AVAILABLE = True

    # Optional: Clean up sys.path if desired, though often not necessary for scripts
    if parent_dir == sys.path[0]:
        sys.path.pop(0)

except ImportError as import_err:
    # Define dummy classes if pydraulics is optional or for testing
    class MockChannel:
        def __init__(self, *args, **kwargs): pass
        def area(self, y): return None
        def hydraulic_radius(self, y): return None
        def calc_properties(self): pass
        a = p = rh = tw = dh = y = None

    RectangularChannel = TrapezoidalChannel = CircularChannel = TriangularChannel = MockChannel
    PYDRAULICS_AVAILABLE = False
    logging.error(f"Pydraulics library (or channel class source) not found in parent directory or Python path: {import_err}")

# Configure logging
logger = logging.getLogger("fluids-mcp.open_channel_flow")

def calculate_open_channel_flow(
    # --- Channel Geometry ---
    channel_type: Literal["rectangular", "trapezoidal", "circular", "triangular"],
    width_m: Optional[float] = None,       # Bottom width (rect/trap)
    diameter_m: Optional[float] = None,    # Diameter (circular)
    side_slope_z: Optional[float] = None,  # Side slope z:1 (trap/triangular)

    # --- Flow Conditions ---
    slope: Optional[float] = None,         # Channel bottom slope (dimensionless)
    manning_n: Optional[float] = None,     # Manning's n

    # --- Specify EITHER flow_rate OR normal_depth ---
    flow_rate_m3s: Optional[float] = None, # Flow rate in m³/s
    normal_depth_m: Optional[float] = None # Normal flow depth in m

) -> str:
    """Calculates open channel flow rate or normal depth using Manning's equation."""

    results_log = []
    error_log = []
    calculated_results = {}
    geometry_params = {}
    g = 9.80665  # m/s² (standard gravity)

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

    # --- Helper to compute Q given y, using pydraulics + Manning ---
    def calculate_q_from_y(y, chan_type, geom_params, S, n):
        # Ensure y is a valid positive number for calculation
        if not isinstance(y, (int, float)) or y <= 1e-9:
            # logger.debug(f"calculate_q_from_y called with invalid y: {y}. Returning 0.")
            return 0.0

        channel = None
        A = None
        Rh = None
        try:
            # Instantiate the correct pydraulics channel object.
            # Provide the *numerical* depth 'y' and Q=None.
            # pydraulics calculates properties internally during __init__.
            if chan_type == "rectangular":
                if 'b' not in geom_params: raise ValueError("Missing 'b' for rectangular")
                channel = RectangularChannel(n=n, So=S, b=geom_params['b'], Q=None, y=y)
            elif chan_type == "trapezoidal":
                if 'b' not in geom_params: raise ValueError("Missing 'b' for trapezoidal")
                if 'z' not in geom_params: raise ValueError("Missing 'z' for trapezoidal")
                channel = TrapezoidalChannel(n=n, So=S, z=geom_params['z'], b=geom_params['b'], Q=None, y=y)
            elif chan_type == "circular":
                if 'D' not in geom_params: raise ValueError("Missing 'D' for circular")
                D_circ = geom_params['D']
                # Ensure depth doesn't exceed diameter for calculation stability
                y_eff = min(float(y), D_circ) # Convert y to float here
                if y_eff <= 1e-9: return 0.0 # Avoid issues with zero/negative depth in circular
                channel = CircularChannel(n=n, So=S, D=D_circ, Q=None, y=y_eff)
            elif chan_type == "triangular":
                if 'z' not in geom_params: raise ValueError("Missing 'z' for triangular")
                channel = TriangularChannel(n=n, So=S, z=geom_params['z'], Q=None, y=y)
            else:
                logger.error(f"Unknown channel type '{chan_type}' in calculate_q_from_y.")
                return 0.0

            # Access the calculated attributes directly from the pydraulics object
            A_sym = channel.a
            Rh_sym = channel.rh

            # Convert potential sympy Floats/expressions to standard Python floats
            # Handle cases where attributes might still be None if pydraulics failed internally
            if A_sym is not None:
                try:
                    A = float(A_sym)
                except (TypeError, ValueError):
                    logger.warning(f"Could not convert pydraulics area '{A_sym}' to float for y={y}. Type: {type(A_sym)}")
                    A = None
            if Rh_sym is not None:
                try:
                    Rh = float(Rh_sym)
                except (TypeError, ValueError):
                    logger.warning(f"Could not convert pydraulics hyd_radius '{Rh_sym}' to float for y={y}. Type: {type(Rh_sym)}")
                    Rh = None

            # Check if Area or Hydraulic Radius are invalid *before* calling V_Manning
            if A is None or Rh is None or A < 1e-12 or Rh < 1e-12:
                # logger.debug(f"Invalid A or Rh calculated for y={y}. A={A}, Rh={Rh}. Returning 0.")
                return 0.0

            # Calculate velocity using fluids.open_flow.V_Manning
            V_calc = fluids.open_flow.V_Manning(Rh=Rh, S=S, n=n)
            Q_calc = V_calc * A
            # logger.debug(f"y={y:.6f}, A={A:.6f}, Rh={Rh:.6f}, V={V_calc:.6f}, Q={Q_calc:.6f}")
            return Q_calc

        except ImportError:
            logger.error("Pydraulics library not found during calculate_q_from_y.")
            return 0.0 # Cannot calculate without pydraulics
        except Exception as e:
            # Catch any other error during instantiation or calculation
            logger.error(f"Error in calculate_q_from_y for y={y}, type={chan_type}: {type(e).__name__}: {e}", exc_info=False) # Set exc_info=True for full traceback if needed
            return 0.0 # Return 0 flow on any error

    # Objective for root-finding
    def manning_error_func(y, Q_target, *args):
        Q_calc = calculate_q_from_y(y, *args)
        return Q_calc - Q_target

    # Simple secant method
    def secant_method(func, y0, y1, f_args, tol=1e-7, max_iter=100):
        f0 = func(y0, *f_args)
        for i in range(max_iter):
            f1 = func(y1, *f_args)
            if abs(f1) < tol:
                return y1
            denom = (f1 - f0)
            if abs(denom) < 1e-15:
                # Attempt a small step
                y_next = y1 + (y1 - y0) * 0.01 if abs(y1) > 1e-6 else 0.01
            else:
                y_next = y1 - f1 * (y1 - y0) / denom

            if y_next <= 0:
                y_next = tol * (i + 1)
            if abs(y_next - y1) < tol * abs(y_next) + tol:
                return y_next
            y0, f0 = y1, f1
            y1 = y_next
        logger.warning("Secant method did not converge within %d iterations.", max_iter)
        return y1

    # Main calculation
    try:
        y = None
        Q = None

        if solve_for == "depth":
            Q_target = flow_rate_m3s
            results_log.append(f"Solving for normal depth, target Q={Q_target:.4f} m³/s.")
            y0 = 0.01  # Small positive guess

            if channel_type == "circular":
                D = geometry_params['D']
                # Rough second guess by scaling half-full flow
                theta_half = math.pi
                A_half = fluids.geometry.circular_A_from_theta(D=D, angle=theta_half)
                Rh_half = fluids.geometry.circular_Rh_from_theta(D=D, angle=theta_half)
                Q_half = fluids.open_flow.Q_Manning(A=A_half, Rh=Rh_half, S=slope, n=manning_n)
                y1_guess = (Q_target / Q_half) * (D / 2.0) if Q_half > 1e-9 else (D * 0.5)
                y1 = min(max(y1_guess, 0.02), D * 0.99)
                y0 = min(y0, D * 0.1)
            else:
                y1 = 1.0

            if abs(y0 - y1) < 1e-6:
                y1 = y0 + 0.1

            args = (channel_type, geometry_params, slope, manning_n)
            # Solve Q(y) = Q_target
            y_sol = secant_method(manning_error_func, y0, y1, (Q_target,) + args)
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

        # --- Calculate Derived Parameters (Velocity, Froude, etc.) ---
        if y is not None and Q is not None:
            channel = None
            A = P = T = Rh = Dh = None
            flow_regime = "Unknown"

            try:
                # 1. Instantiate pydraulics channel
                if channel_type == "rectangular":
                    channel = RectangularChannel(n=manning_n, So=slope, b=geometry_params['b'])
                elif channel_type == "trapezoidal":
                    channel = TrapezoidalChannel(n=manning_n, So=slope,
                                                 z=geometry_params['z'],
                                                 b=geometry_params['b'])
                elif channel_type == "circular":
                    D_circ = geometry_params['D']
                    channel = CircularChannel(n=manning_n, So=slope, D=D_circ)
                elif channel_type == "triangular":
                    channel = TriangularChannel(n=manning_n, So=slope, z=geometry_params['z'])

                if not channel:
                    raise ValueError("Channel object is None after instantiation.")

                # 2. Depth assignment (handle circular max)
                if channel_type == "circular":
                    y_eff = min(y, geometry_params['D'])
                    channel.y = y_eff
                else:
                    channel.y = y

                # 3. Force internal property calculations
                try:
                    channel.calc_properties()
                except Exception as calc_err:
                    logger.error(f"Error in calc_properties: {calc_err}")
                    # Fall back to direct calculations without relying on channel object properties
                    if channel_type == "rectangular":
                        b = geometry_params['b']
                        A = b * y
                        P = b + 2 * y
                        Rh = A / P if P > 0 else None
                        T = b
                        Dh = y
                    elif channel_type == "trapezoidal":
                        b = geometry_params['b']
                        z = geometry_params['z']
                        A = (b + z * y) * y
                        P = b + 2 * y * (1 + z**2)**0.5
                        Rh = A / P if P > 0 else None
                        T = b + 2 * z * y
                        Dh = A / T if T > 0 else None
                    elif channel_type == "triangular":
                        z = geometry_params['z']
                        A = z * y**2
                        P = (y * 2) * (1 + z**2)**0.5
                        Rh = A / P if P > 0 else None
                        T = 2 * z * y
                        Dh = 0.5 * y
                    elif channel_type == "circular":
                        D = geometry_params['D']
                        # Use simplified calculations for circular channel
                        if y >= D:  # Full pipe
                            A = math.pi * (D/2)**2
                            P = math.pi * D
                            Rh = D/4
                            T = 0
                            Dh = D/4
                        else:  # Partially full pipe
                            theta = 2 * math.acos(1 - 2*y/D)  # Central angle
                            A = (D**2/8) * (theta - math.sin(theta))
                            P = (D/2) * theta
                            Rh = A / P if P > 0 else None
                            T = D * math.sin(theta/2)
                            Dh = A / T if T > 0 else None
                    else:
                        A = P = Rh = T = Dh = None
                    return A, P, Rh, T, Dh
                
                # 4. Retrieve geometry
                try:
                    A = float(channel.a) if channel.a is not None else None
                    P = float(channel.p) if channel.p is not None else None
                    Rh = float(channel.rh) if channel.rh is not None else None
                    T = float(channel.tw) if channel.tw is not None else None
                    Dh = float(channel.dh) if channel.dh is not None else None
                    return A, P, Rh, T, Dh
                except (TypeError, ValueError) as convert_err:
                    logger.error(f"Error converting channel properties to floats: {convert_err}")
                    # If conversion fails, use our direct calculations
                    return None, None, None, None, None

                # --- Validation/Conversion for numeric properties ---
                numeric_props = {'A': A, 'P': P, 'Rh': Rh, 'T': T, 'Dh': Dh}
                for name, val in numeric_props.items():
                    if val is None:
                        numeric_props[name] = None
                        error_log.append(f"Warning: {name} returned None from channel object.")
                    elif isinstance(val, (int, float)):
                        # Already numeric
                        pass
                    elif hasattr(val, 'evalf'):
                        # Symbolic, attempt eval
                        try:
                            numeric_props[name] = float(val.evalf())
                        except Exception as eval_err:
                            error_log.append(f"Warning: Could not evalf() property '{name}': {eval_err}.")
                            numeric_props[name] = None
                    else:
                        # Unexpected type
                        error_log.append(f"Warning: {name} has unexpected type {type(val)}, treating as None.")
                        numeric_props[name] = None

                # Reassign validated values
                A = numeric_props['A']
                P = numeric_props['P']
                Rh = numeric_props['Rh']
                T = numeric_props['T']
                Dh = numeric_props['Dh']

                # --- Velocity, Froude number ---
                if A is not None and A > 1e-12:
                    V = Q / A
                    calculated_results['velocity_m_s'] = V
                    if Dh is not None and Dh > 1e-9:
                        try:
                            Fr = V / math.sqrt(g * Dh)
                            calculated_results['froude_number'] = Fr
                            calculated_results['hydraulic_depth_m'] = Dh
                            if Fr < 1.0:
                                flow_regime = "Subcritical"
                            elif abs(Fr - 1.0) < 1e-6:
                                flow_regime = "Critical"
                            else:
                                flow_regime = "Supercritical"
                            calculated_results['flow_regime'] = flow_regime
                        except Exception as fr_e:
                            error_log.append(f"Error calculating Froude number: {fr_e}")
                            calculated_results['froude_number'] = "Error"
                            calculated_results['flow_regime'] = "Error"
                    else:
                        results_log.append("Dh is None/invalid; skipping Froude calc.")
                        calculated_results['froude_number'] = None
                        calculated_results['hydraulic_depth_m'] = Dh
                        calculated_results['flow_regime'] = "N/A (Invalid Dh)"
                else:
                    results_log.append("Area is None/invalid; velocity = 0, skipping Froude calc.")
                    calculated_results['velocity_m_s'] = 0.0
                    calculated_results['froude_number'] = None
                    calculated_results['hydraulic_depth_m'] = None
                    calculated_results['flow_regime'] = "N/A (Invalid Area)"

                # Store geometry
                calculated_results['flow_area_m2'] = A
                calculated_results['wetted_perimeter_m'] = P
                calculated_results['hydraulic_radius_m'] = Rh
                calculated_results['top_width_m'] = T

            except Exception as geo_derived_e:
                error_log.append(f"Error getting derived geometry: {type(geo_derived_e).__name__}: {geo_derived_e}")
                logger.error(f"Error getting derived geometry: {geo_derived_e}", exc_info=True)
                calculated_results.update({
                    'velocity_m_s': "Error",
                    'froude_number': "Error",
                    'flow_regime': "Error",
                    'flow_area_m2': "Error",
                    'wetted_perimeter_m': "Error",
                    'hydraulic_radius_m': "Error",
                    'top_width_m': "Error",
                    'hydraulic_depth_m': "Error"
                })
        else:
            error_log.append("Calculation failed: no valid depth/flow rate to derive geometry.")

    except Exception as e:
        logger.error(f"Error during open channel calculation: {e}", exc_info=True)
        error_log.append(f"Calculation error: {str(e)}")
        return json.dumps({"errors": error_log, "log": results_log})

    # --- Format and return ---
    final_result = {
        "status": "Success" if not error_log else "Failure",
        "log": results_log,
        "errors": error_log if error_log else None,
        "inputs": {
            "channel_type": channel_type,
            "geometry_params": geometry_params,
            "slope": slope,
            "manning_n": manning_n,
            "provided_flow_rate_m3s": flow_rate_m3s,
            "provided_normal_depth_m": normal_depth_m
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

    return json.dumps(final_result)
