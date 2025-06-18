import json
import logging
import math
import os
from typing import Optional, Literal, Dict, Any, Tuple, Union, List

# Scientific computing
import scipy.optimize

# Import from pydraulics (now properly packaged)
try:
    from pydraulics.open_flow import (
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
from utils.constants import G_GRAVITY

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
    """Calculates open channel flow rate or normal depth using Manning's equation.
    
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
    calculated_results: Dict[str, Any] = {}
    geometry_params: Dict[str, float] = {}
    g = G_GRAVITY  # m/s² (standard gravity)

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
                    except:
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
                    except:
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
                    except:
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
                    except:
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
            except:
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
        """Find normal depth using robust scipy methods with proper bracketing."""
        
        # Establish reasonable bounds for different channel types
        y_min = 1e-6  # Very small positive depth
        
        if channel_type == "circular" and 'D' in geometry_params:
            y_max = geometry_params['D'] * 0.99  # Just under full pipe
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
                
            # Calculate velocity and Froude number
            if A is not None and A > 1e-12:
                V = Q / A
                calculated_results['velocity_m_s'] = V
                
                if Dh is not None and Dh > 1e-9:
                    try:
                        Fr = V / math.sqrt(g * Dh)
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
