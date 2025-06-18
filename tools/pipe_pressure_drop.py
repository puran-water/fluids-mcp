"""
Pipe pressure drop calculation tool.

This module provides a tool to calculate pressure drop in pipes with various fluids and fittings.
"""

import json
import math
import logging
import fluids
import fluids.core
import fluids.friction
from typing import List, Dict, Optional, Union, Any

# Import shared utilities
from utils.constants import (
    GPM_to_M3S, INCH_to_M, FT_to_M, PSI_to_PA, 
    LBFT3_to_KGM3, CENTIPOISE_to_PAS, DEFAULT_ROUGHNESS
)
from utils.helpers import get_fitting_K
from utils.import_helpers import FLUIDPROP_AVAILABLE, FluidProperties, FLUID_SELECTION

# Configure logging
logger = logging.getLogger("fluids-mcp.pipe_pressure_drop")

def calculate_total_pressure_drop(flow_rate, pipe_diameter, pipe_length, fluid_density, fluid_viscosity, pipe_roughness, fittings_list):
    """Helper function to calculate total pressure drop for given parameters"""
    area = math.pi * (pipe_diameter / 2) ** 2
    velocity = flow_rate / area
    
    Re = fluids.core.Reynolds(V=velocity, D=pipe_diameter, rho=fluid_density, mu=fluid_viscosity)
    fd = fluids.friction.friction_factor(Re=Re, eD=pipe_roughness / pipe_diameter)
    dP_straight = fd * (pipe_length / pipe_diameter) * (fluid_density * velocity**2 / 2.0)
    
    # Calculate fittings pressure drop
    K_total = 0
    for fitting in fittings_list:
        fitting_type = fitting.get("type", "").lower()
        quantity = fitting.get("quantity", 1)
        K_value_provided = fitting.get("K_value")
        
        if K_value_provided is not None:
            K_fitting = float(K_value_provided)
        else:
            K_fitting = get_fitting_K(fitting_type, pipe_diameter, Re, flow_rate)
        
        K_total += K_fitting * quantity
    
    dP_fittings = fluids.dP_from_K(K=K_total, rho=fluid_density, V=velocity)
    return dP_straight + dP_fittings

def solve_for_flow_rate(target_pressure_drop, pipe_diameter, pipe_length, fluid_density, fluid_viscosity, pipe_roughness, fittings_list):
    """Solve for flow rate given target pressure drop using numerical methods"""
    from scipy.optimize import brentq
    
    def pressure_drop_error(flow_rate):
        if flow_rate <= 0:
            return float('inf')
        try:
            calculated_dp = calculate_total_pressure_drop(
                flow_rate, pipe_diameter, pipe_length, fluid_density, fluid_viscosity, pipe_roughness, fittings_list
            )
            return calculated_dp - target_pressure_drop
        except:
            return float('inf')
    
    # Initial bounds for flow rate search
    flow_min = 1e-6  # Very small flow
    flow_max = 10.0   # Large flow rate
    
    # Find bounds where function changes sign
    try:
        result = brentq(pressure_drop_error, flow_min, flow_max, xtol=1e-8)
        return result
    except ValueError:
        # If brentq fails, try a different approach
        from scipy.optimize import minimize_scalar
        result = minimize_scalar(lambda q: abs(pressure_drop_error(q)), bounds=(flow_min, flow_max), method='bounded')
        if result.success:
            return result.x
        else:
            raise ValueError("Failed to converge on flow rate solution")

def solve_for_diameter(target_pressure_drop, flow_rate, pipe_length, fluid_density, fluid_viscosity, pipe_roughness, fittings_list):
    """Solve for pipe diameter given target pressure drop using numerical methods"""
    from scipy.optimize import brentq
    
    def pressure_drop_error(diameter):
        if diameter <= 0:
            return float('inf')
        try:
            calculated_dp = calculate_total_pressure_drop(
                flow_rate, diameter, pipe_length, fluid_density, fluid_viscosity, pipe_roughness, fittings_list
            )
            return calculated_dp - target_pressure_drop
        except:
            return float('inf')
    
    # Initial bounds for diameter search
    D_min = 0.001   # 1 mm
    D_max = 2.0     # 2 meters
    
    try:
        result = brentq(pressure_drop_error, D_min, D_max, xtol=1e-8)
        return result
    except ValueError:
        from scipy.optimize import minimize_scalar
        result = minimize_scalar(lambda d: abs(pressure_drop_error(d)), bounds=(D_min, D_max), method='bounded')
        if result.success:
            return result.x
        else:
            raise ValueError("Failed to converge on diameter solution")

def solve_for_length(target_pressure_drop, flow_rate, pipe_diameter, fluid_density, fluid_viscosity, pipe_roughness, fittings_list):
    """Solve for pipe length given target pressure drop using numerical methods"""
    from scipy.optimize import brentq
    
    def pressure_drop_error(length):
        if length <= 0:
            return float('inf')
        try:
            calculated_dp = calculate_total_pressure_drop(
                flow_rate, pipe_diameter, length, fluid_density, fluid_viscosity, pipe_roughness, fittings_list
            )
            return calculated_dp - target_pressure_drop
        except:
            return float('inf')
    
    # Initial bounds for length search
    L_min = 0.001   # 1 mm
    L_max = 10000.0 # 10 km
    
    try:
        result = brentq(pressure_drop_error, L_min, L_max, xtol=1e-8)
        return result
    except ValueError:
        from scipy.optimize import minimize_scalar
        result = minimize_scalar(lambda l: abs(pressure_drop_error(l)), bounds=(L_min, L_max), method='bounded')
        if result.success:
            return result.x
        else:
            raise ValueError("Failed to converge on length solution")

def calculate_pipe_pressure_drop(
    # --- Core SI Inputs (still supported) ---
    flow_rate: Optional[float] = None,           # Flow rate in m³/s
    pipe_diameter: Optional[float] = None,       # Pipe inner diameter in m
    pipe_length: Optional[float] = None,         # Pipe length in m
    fluid_density: Optional[float] = None,       # Fluid density in kg/m³
    fluid_viscosity: Optional[float] = None,     # Fluid dynamic viscosity in Pa·s
    pipe_roughness: Optional[float] = None,      # Pipe roughness in m

    # --- Alternative Unit Inputs ---
    flow_rate_gpm: Optional[float] = None,       # Flow rate in US GPM
    pipe_diameter_in: Optional[float] = None,    # Pipe inner diameter in inches
    pipe_length_ft: Optional[float] = None,      # Pipe length in feet
    fluid_density_lbft3: Optional[float] = None, # Fluid density in lb/ft³
    fluid_viscosity_cp: Optional[float] = None,  # Fluid dynamic viscosity in centipoise

    # --- Property Lookup Inputs ---
    fluid_name: Optional[str] = None,            # e.g., "Water", "Air"
    temperature_c: Optional[float] = None,       # Temperature in Celsius (for property lookup)
    pressure_bar: Optional[float] = None,        # Pressure in bar (for property lookup, default ~1 atm)
    nominal_size_in: Optional[float] = None,     # Nominal pipe size in inches
    schedule: str = "40",                        # Pipe schedule
    material: Optional[str] = None,              # e.g., "Steel", "PVC"

    # --- Pressure Drop / Design Parameters ---
    pressure_drop: Optional[float] = None,          # Total pressure drop in Pa
    pressure_drop_psi: Optional[float] = None,      # Total pressure drop in psi
    
    # --- Fittings ---
    fittings: List[Dict[str, Union[str, int, float]]] = None  # List of fittings with type, quantity, and optional K_value
) -> str:
    """Calculate the pressure drop through straight pipe and fittings with flexible input options.
    
    This function accepts multiple input formats for flexibility:
    1. SI units (direct inputs)
    2. Imperial units (with automatic conversion)
    3. Property lookups based on standard materials, fluids, and pipe sizes
    
    Requirements:
    - Must provide either flow_rate (SI) OR flow_rate_gpm
    - Must provide either pipe_diameter (SI) OR pipe_diameter_in OR nominal_size_in with schedule
    - Must provide either pipe_length (SI) OR pipe_length_ft
    - Must provide either (fluid_density AND fluid_viscosity) OR (fluid_density_lbft3 AND fluid_viscosity_cp) 
      OR (fluid_name AND temperature_c)
    - Roughness is optional with defaults available
    
    Args:
        flow_rate: Optional - Flow rate in m³/s
        pipe_diameter: Optional - Pipe inner diameter in m
        pipe_length: Optional - Pipe length in m
        fluid_density: Optional - Fluid density in kg/m³
        fluid_viscosity: Optional - Fluid viscosity in Pa·s
        pipe_roughness: Optional - Pipe roughness in m
        flow_rate_gpm: Optional - Flow rate in US GPM
        pipe_diameter_in: Optional - Pipe inner diameter in inches
        pipe_length_ft: Optional - Pipe length in feet
        fluid_density_lbft3: Optional - Fluid density in lb/ft³
        fluid_viscosity_cp: Optional - Fluid dynamic viscosity in centipoise
        fluid_name: Optional - Name of the fluid (e.g., "Water", "Air")
        temperature_c: Optional - Temperature in Celsius (for property lookup)
        pressure_bar: Optional - Pressure in bar (for property lookup)
        nominal_size_in: Optional - Nominal pipe size in inches
        schedule: Pipe schedule (default: "40")
        material: Optional - Pipe material name (e.g., "Steel", "PVC")
        fittings: Optional - List of fittings, each with 'type' and 'quantity', and optionally 'K_value'.
                 Types include: 90_elbow, 45_elbow, tee_run_through, tee_branch_flow,
                 gate_valve, globe_valve, check_valve_swing, ball_valve, butterfly_valve,
                 entrance_sharp, exit_normal, etc.
    
    Returns:
        Detailed pressure drop analysis including head loss with unit conversions
    """
    
    results_log = [] # To track how values were obtained
    error_log = []
    
    try:
        # 1. Resolve Flow Rate (SI)
        local_flow_rate = None
        if flow_rate is not None:
            local_flow_rate = flow_rate
            results_log.append("Used provided SI flow_rate.")
        elif flow_rate_gpm is not None:
            local_flow_rate = flow_rate_gpm * GPM_to_M3S
            results_log.append(f"Converted flow_rate from {flow_rate_gpm} GPM.")
        else:
            error_log.append("Missing required input: flow_rate or flow_rate_gpm.")
            
        # 2. Resolve Pipe Length (SI)
        local_pipe_length = None
        if pipe_length is not None:
            local_pipe_length = pipe_length
            results_log.append("Used provided SI pipe_length.")
        elif pipe_length_ft is not None:
            local_pipe_length = pipe_length_ft * FT_to_M
            results_log.append(f"Converted pipe_length from {pipe_length_ft} ft.")
        else:
            error_log.append("Missing required input: pipe_length or pipe_length_ft.")
            
        # 3. Resolve Pipe Diameter & Roughness (SI)
        local_pipe_diameter = None
        local_pipe_roughness = None
        pipe_info = {} # Store details from lookup
        
        if pipe_diameter is not None:
            local_pipe_diameter = pipe_diameter
            results_log.append("Used provided SI pipe_diameter.")
            # Roughness: Use provided, lookup by material, or default
            if pipe_roughness is not None:
                local_pipe_roughness = pipe_roughness
                results_log.append("Used provided SI pipe_roughness.")
            elif material is not None:
                # Internal lookup logic (similar to get_pipe_properties)
                try:
                    from fluids.friction import _roughness
                    mat_lower = {k.lower(): k for k in _roughness.keys()}
                    if material.lower() in mat_lower:
                        actual_key = mat_lower[material.lower()]
                        local_pipe_roughness = _roughness[actual_key]
                        results_log.append(f"Looked up roughness for material '{actual_key}'.")
                    else:
                        local_pipe_roughness = DEFAULT_ROUGHNESS # Default
                        results_log.append(f"Material '{material}' not found, using default roughness.")
                        error_log.append(f"Warning: Material '{material}' not found for roughness lookup.")
                except Exception as lookup_e:
                     local_pipe_roughness = DEFAULT_ROUGHNESS # Default
                     results_log.append(f"Roughness lookup failed ({lookup_e}), using default.")
                     error_log.append(f"Warning: Roughness lookup failed: {lookup_e}")
            else:
                local_pipe_roughness = DEFAULT_ROUGHNESS # Default if nothing else provided
                results_log.append("Used default pipe_roughness.")
            
        elif pipe_diameter_in is not None:
            local_pipe_diameter = pipe_diameter_in * INCH_to_M
            results_log.append(f"Converted pipe_diameter from {pipe_diameter_in} inches.")
            # Roughness logic (same as above)
            if pipe_roughness is not None: # Allow explicit SI roughness even if diameter is imperial
                 local_pipe_roughness = pipe_roughness
                 results_log.append("Used provided SI pipe_roughness.")
            elif material is not None:
                 # Internal lookup logic
                 try:
                    from fluids.friction import _roughness
                    mat_lower = {k.lower(): k for k in _roughness.keys()}
                    if material.lower() in mat_lower:
                        actual_key = mat_lower[material.lower()]
                        local_pipe_roughness = _roughness[actual_key]
                        results_log.append(f"Looked up roughness for material '{actual_key}'.")
                    else:
                        local_pipe_roughness = DEFAULT_ROUGHNESS
                        results_log.append(f"Material '{material}' not found, using default roughness.")
                        error_log.append(f"Warning: Material '{material}' not found for roughness lookup.")
                 except Exception as lookup_e:
                     local_pipe_roughness = DEFAULT_ROUGHNESS
                     results_log.append(f"Roughness lookup failed ({lookup_e}), using default.")
                     error_log.append(f"Warning: Roughness lookup failed: {lookup_e}")
            else:
                local_pipe_roughness = DEFAULT_ROUGHNESS
                results_log.append("Used default pipe_roughness.")
                
        elif nominal_size_in is not None:
            # Use fluids.piping.nearest_pipe to get dimensions
            # For pipe dimension lookup failure
            try:
                NPS, Di, Do, t = fluids.piping.nearest_pipe(NPS=nominal_size_in, schedule=schedule)
                local_pipe_diameter = Di
                pipe_info = {"NPS_in": NPS, "schedule": schedule, "Do_m": Do, "t_m": t}
                results_log.append(f"Looked up pipe dimensions for NPS {nominal_size_in} Sch {schedule}.")
                
                # Roughness logic (same as above)
                if pipe_roughness is not None:
                    local_pipe_roughness = pipe_roughness
                    results_log.append("Used provided SI pipe_roughness.")
                elif material is not None:
                    # Internal lookup logic
                    try:
                        from fluids.friction import _roughness
                        mat_lower = {k.lower(): k for k in _roughness.keys()}
                        if material.lower() in mat_lower:
                            actual_key = mat_lower[material.lower()]
                            local_pipe_roughness = _roughness[actual_key]
                            results_log.append(f"Looked up roughness for material '{actual_key}'.")
                        else:
                            local_pipe_roughness = DEFAULT_ROUGHNESS
                            results_log.append(f"Material '{material}' not found, using default roughness.")
                            error_log.append(f"Warning: Material '{material}' not found for roughness lookup.")
                    except Exception as lookup_e:
                        local_pipe_roughness = DEFAULT_ROUGHNESS
                        results_log.append(f"Roughness lookup failed ({lookup_e}), using default.")
                        error_log.append(f"Warning: Roughness lookup failed: {lookup_e}")
                else:
                    local_pipe_roughness = DEFAULT_ROUGHNESS # Default if material not specified
                    results_log.append("Used default pipe_roughness (material not specified).")
                    pipe_info["roughness_used_m"] = local_pipe_roughness
                    pipe_info["roughness_source"] = "Default"
            except Exception as pipe_lookup_e:
                error_log.append(f"Failed to look up pipe NPS {nominal_size_in} Sch {schedule}: {pipe_lookup_e}")
                local_pipe_diameter = None  # Explicitly set to None on failure
                local_pipe_roughness = None # Explicitly set to None on failure
        else:
            error_log.append("Missing required input: pipe_diameter, pipe_diameter_in, or nominal_size_in.")

        # 4. Resolve Fluid Properties (SI)
        local_fluid_density = None
        local_fluid_viscosity = None
        
        if fluid_density is not None and fluid_viscosity is not None:
            local_fluid_density = fluid_density
            local_fluid_viscosity = fluid_viscosity
            results_log.append("Used provided SI fluid density and viscosity.")
        elif fluid_density_lbft3 is not None and fluid_viscosity_cp is not None:
             local_fluid_density = fluid_density_lbft3 * LBFT3_to_KGM3
             local_fluid_viscosity = fluid_viscosity_cp * CENTIPOISE_to_PAS
             results_log.append(f"Converted density from {fluid_density_lbft3} lb/ft³ and viscosity from {fluid_viscosity_cp} cP.")
        elif fluid_name is not None and temperature_c is not None and FLUIDPROP_AVAILABLE:
            # Internal lookup using fluidprop or similar
            try:
                valid_fluids = [f[0] for f in FLUID_SELECTION]
                actual_fluid_name = fluid_name
                
                if fluid_name not in valid_fluids:
                    # Try case-insensitive match
                    fluid_lower = fluid_name.lower()
                    match = next((f for f in valid_fluids if f.lower() == fluid_lower), None)
                    if match:
                        actual_fluid_name = match
                    else:
                        raise ValueError(f"Fluid '{fluid_name}' not found.")
                
                p_bar = pressure_bar if pressure_bar is not None else 1.01325 # Default pressure if not given
                fluid_props = FluidProperties(coolprop_name=actual_fluid_name, T_in_deg_C=temperature_c, P_in_bar=p_bar)
                local_fluid_density = float(fluid_props.rho[0])
                local_fluid_viscosity = float(fluid_props.eta[0])
                results_log.append(f"Looked up properties for {actual_fluid_name} at {temperature_c}°C, {p_bar} bar.")
            except Exception as fluid_lookup_e:
                error_log.append(f"Failed to look up fluid properties for {fluid_name} at {temperature_c}°C: {fluid_lookup_e}")
                local_fluid_density = None  # Explicitly set to None on failure
                local_fluid_viscosity = None # Explicitly set to None on failure        
        else:
            error_log.append("Missing required inputs for fluid properties: (fluid_density AND fluid_viscosity) OR (fluid_density_lbft3 AND fluid_viscosity_cp) OR (fluid_name AND temperature_c).")

        # 5. Resolve Pressure Drop (SI)
        local_pressure_drop = None
        if pressure_drop is not None:
            local_pressure_drop = pressure_drop
            results_log.append("Used provided SI pressure_drop.")
        elif pressure_drop_psi is not None:
            local_pressure_drop = pressure_drop_psi * PSI_to_PA
            results_log.append(f"Converted pressure_drop from {pressure_drop_psi} psi.")
        # else: to be solved for or calculated

        # 6. Determine solve_for variable
        specified_vars = {
            'flow_rate': local_flow_rate,
            'pipe_diameter': local_pipe_diameter, 
            'pipe_length': local_pipe_length,
            'pressure_drop': local_pressure_drop
        }
        unknown_vars = [k for k, v in specified_vars.items() if v is None]
        
        if len(unknown_vars) == 0:
            # All variables specified - just calculate and verify
            solve_for = None
            results_log.append("All primary variables specified - calculating for verification.")
        elif len(unknown_vars) == 1:
            solve_for = unknown_vars[0]
            results_log.append(f"Solving for: {solve_for}")
        else:
            error_log.append(f"Too many unknown variables. Specify 3 of 4 primary variables: {list(specified_vars.keys())}. Unknown: {unknown_vars}")
            
        # --- Check for any errors before proceeding ---
        if error_log:
             # If critical inputs are missing, return error early
             missing_critical = any(e.startswith("Missing required input") or e.startswith("Too many unknown variables") for e in error_log)
             if missing_critical:
                 return json.dumps({"errors": error_log, "log": results_log})
             # Otherwise, proceed but include warnings

        # --- Core Calculation Logic (using resolved local_ SI variables) ---
        # Check critical variables (excluding the one we're solving for)
        critical_base_vars = [local_fluid_density, local_fluid_viscosity, local_pipe_roughness]
        if solve_for != 'flow_rate' and local_flow_rate is None:
            critical_base_vars.append(local_flow_rate)
        if solve_for != 'pipe_diameter' and local_pipe_diameter is None:
            critical_base_vars.append(local_pipe_diameter)
        if solve_for != 'pipe_length' and local_pipe_length is None:
            critical_base_vars.append(local_pipe_length)
        if solve_for != 'pressure_drop' and local_pressure_drop is None and solve_for is not None:
            critical_base_vars.append(local_pressure_drop)
            
        if any(v is None for v in critical_base_vars):
            return json.dumps({"errors": ["Critical input parameter resolution failed."] + error_log, "log": results_log})

        # --- Solving Logic ---
        if solve_for is not None:
            try:
                if solve_for == 'pressure_drop':
                    # Standard calculation - no solving needed
                    pass
                elif solve_for == 'flow_rate':
                    # Solve for flow rate given pressure drop
                    local_flow_rate = solve_for_flow_rate(
                        local_pressure_drop, local_pipe_diameter, local_pipe_length,
                        local_fluid_density, local_fluid_viscosity, local_pipe_roughness, fittings_list
                    )
                    results_log.append(f"Solved flow_rate = {local_flow_rate:.6f} m³/s")
                elif solve_for == 'pipe_diameter':
                    # Solve for diameter given pressure drop and flow rate
                    local_pipe_diameter = solve_for_diameter(
                        local_pressure_drop, local_flow_rate, local_pipe_length,
                        local_fluid_density, local_fluid_viscosity, local_pipe_roughness, fittings_list
                    )
                    results_log.append(f"Solved pipe_diameter = {local_pipe_diameter:.6f} m")
                elif solve_for == 'pipe_length':
                    # Solve for length given pressure drop, flow rate, and diameter
                    local_pipe_length = solve_for_length(
                        local_pressure_drop, local_flow_rate, local_pipe_diameter,
                        local_fluid_density, local_fluid_viscosity, local_pipe_roughness, fittings_list
                    )
                    results_log.append(f"Solved pipe_length = {local_pipe_length:.3f} m")
            except Exception as solve_error:
                error_log.append(f"Failed to solve for {solve_for}: {solve_error}")
                return json.dumps({"errors": error_log, "log": results_log})
             
        area = math.pi * (local_pipe_diameter/2)**2
        if area == 0:
             return json.dumps({"error": "Pipe diameter resolved to zero.", "log": results_log})
        velocity = local_flow_rate / area
        
        Re = fluids.core.Reynolds(V=velocity, D=local_pipe_diameter, rho=local_fluid_density, mu=local_fluid_viscosity)
        fd = fluids.friction.friction_factor(Re=Re, eD=local_pipe_roughness/local_pipe_diameter)
        dP_straight = fd * (local_pipe_length / local_pipe_diameter) * (local_fluid_density * velocity**2 / 2.0)
        
        # --- Fitting Calculation (using helper function) ---
        K_total = 0
        fitting_details = []
        fittings_list = fittings if fittings else []
        
        for fitting in fittings_list:
             fitting_type = fitting.get("type", "").lower()
             quantity = fitting.get("quantity", 1)
             K_value_provided = fitting.get("K_value") # Explicit K
             K_fitting = 0
             
             if K_value_provided is not None:
                 K_fitting = float(K_value_provided)
                 fitting_details.append({"type": fitting_type, "quantity": quantity, "K_source": "Provided", "K_total": round(K_fitting * quantity, 4)})
             else:
                 # Use simplified fitting mapping via helper function
                 K_fitting = get_fitting_K(fitting_type, local_pipe_diameter, Re, local_flow_rate)
                 fitting_details.append({"type": fitting_type, "quantity": quantity, "K_source": "Calculated", "K_individual": round(K_fitting, 4), "K_total": round(K_fitting * quantity, 4)})
             
             K_total += K_fitting * quantity
        
        dP_fittings = fluids.dP_from_K(K=K_total, rho=local_fluid_density, V=velocity)
        dP_total = dP_straight + dP_fittings
        head_loss = dP_total / (local_fluid_density * 9.81) if local_fluid_density > 0 else 0

        # --- Consolidate Results ---
        final_result = {
            "inputs_resolved": results_log,
            "warnings": error_log if error_log else None, # Include warnings if any
            "pipe_details": pipe_info if pipe_info else None,
            # Solving Information
            "solved_variable": solve_for if solve_for else None,
            # Core Results
            "flow_rate_m3s": round(local_flow_rate, 6),
            "flow_rate_gpm": round(local_flow_rate / GPM_to_M3S, 2),
            "pipe_diameter_m": round(local_pipe_diameter, 6),
            "pipe_diameter_in": round(local_pipe_diameter / INCH_to_M, 4),
            "pipe_length_m": round(local_pipe_length, 3),
            "pipe_length_ft": round(local_pipe_length / FT_to_M, 1),
            "reynolds_number": round(Re, 2),
            "friction_factor": round(fd, 6),
            "flow_velocity_m_s": round(velocity, 3),
            "pressure_drop_straight_pipe_pa": round(dP_straight, 2),
            "pressure_drop_fittings_pa": round(dP_fittings, 2),
            "pressure_drop_total_pa": round(dP_total, 2),
            "pressure_drop_total_psi": round(dP_total / PSI_to_PA, 3), # Example output conversion
            "head_loss_m": round(head_loss, 3),
            "head_loss_ft": round(head_loss / FT_to_M, 3), # Example output conversion
            "fitting_details": fitting_details
        }
        # Remove null warnings/pipe_details if empty
        if not final_result["warnings"]: del final_result["warnings"]
        if not final_result["pipe_details"]: del final_result["pipe_details"]
        
        return json.dumps(final_result)
    
    except Exception as e:
        logger.error(f"Error in calculate_pipe_pressure_drop: {e}", exc_info=True)
        return json.dumps({"error": f"Calculation error: {str(e)}", "log": results_log, "errors_occurred": error_log})

def pipe_pressure_drop_sweep(
    variable: str,
    start: float, 
    stop: float, 
    n: int,
    # Base calculation parameters
    flow_rate: Optional[float] = None,
    pipe_diameter: Optional[float] = None,
    pipe_length: Optional[float] = None,
    fluid_density: Optional[float] = None,
    fluid_viscosity: Optional[float] = None,
    pipe_roughness: Optional[float] = None,
    flow_rate_gpm: Optional[float] = None,
    pipe_diameter_in: Optional[float] = None,
    pipe_length_ft: Optional[float] = None,
    fluid_density_lbft3: Optional[float] = None,
    fluid_viscosity_cp: Optional[float] = None,
    pressure_drop: Optional[float] = None,
    pressure_drop_psi: Optional[float] = None,
    fluid_name: Optional[str] = None,
    temperature_c: Optional[float] = None,
    pressure_bar: Optional[float] = None,
    nominal_size_in: Optional[float] = None,
    schedule: str = "40",
    material: Optional[str] = None,
    fittings: List[Dict[str, Union[str, int, float]]] = None
) -> str:
    """Parameter sweep for liquid pipe pressure drop analysis
    
    Enables design optimization and parameter studies for pipe sizing.
    Sweeps one variable while keeping others constant to analyze system behavior.
    
    Args:
        variable: Variable to sweep ('flow_rate', 'pipe_diameter', 'pipe_length', 'pressure_drop')
        start: Start value for sweep
        stop: Stop value for sweep  
        n: Number of points in sweep
        **other_params: All other parameters from calculate_pipe_pressure_drop
    
    Returns:
        JSON string with sweep results as list of dictionaries
    """
    try:
        import numpy as np
    except ImportError:
        return json.dumps({"error": "numpy required for sweep functionality"})
    
    # Generate sweep values
    sweep_values = np.linspace(start, stop, n)
    results = []
    
    # Build base parameters dictionary
    base_kwargs = {
        'flow_rate': flow_rate,
        'pipe_diameter': pipe_diameter,
        'pipe_length': pipe_length,
        'fluid_density': fluid_density,
        'fluid_viscosity': fluid_viscosity,
        'pipe_roughness': pipe_roughness,
        'flow_rate_gpm': flow_rate_gpm,
        'pipe_diameter_in': pipe_diameter_in,
        'pipe_length_ft': pipe_length_ft,
        'fluid_density_lbft3': fluid_density_lbft3,
        'fluid_viscosity_cp': fluid_viscosity_cp,
        'pressure_drop': pressure_drop,
        'pressure_drop_psi': pressure_drop_psi,
        'fluid_name': fluid_name,
        'temperature_c': temperature_c,
        'pressure_bar': pressure_bar,
        'nominal_size_in': nominal_size_in,
        'schedule': schedule,
        'material': material,
        'fittings': fittings
    }
    
    # Remove None values
    base_kwargs = {k: v for k, v in base_kwargs.items() if v is not None}
    
    for value in sweep_values:
        try:
            # Set up parameters for this sweep point
            kwargs = base_kwargs.copy()
            kwargs[variable] = value
            
            # Calculate
            result_json = calculate_pipe_pressure_drop(**kwargs)
            result = json.loads(result_json)
            
            if 'errors' not in result or not result['errors']:
                # Extract key results
                row = {
                    variable: round(value, 6),
                    'flow_rate_m3s': result.get('flow_rate_m3s', None),
                    'pipe_diameter_m': result.get('pipe_diameter_m', None),
                    'pipe_length_m': result.get('pipe_length_m', None),
                    'pressure_drop_total_pa': result.get('pressure_drop_total_pa', None),
                    'pressure_drop_total_psi': result.get('pressure_drop_total_psi', None),
                    'flow_velocity_m_s': result.get('flow_velocity_m_s', None),
                    'reynolds_number': result.get('reynolds_number', None),
                    'head_loss_m': result.get('head_loss_m', None),
                    'solved_variable': result.get('solved_variable', None)
                }
                results.append(row)
            else:
                # Add row with errors
                row = {variable: round(value, 6), 'error': str(result.get('errors', 'Unknown error'))}
                results.append(row)
                
        except Exception as e:
            row = {variable: round(value, 6), 'error': str(e)}
            results.append(row)
    
    return json.dumps({
        "sweep_variable": variable,
        "sweep_range": {"start": start, "stop": stop, "n": n},
        "results": results,
        "summary": {
            "total_points": len(results),
            "successful_points": len([r for r in results if 'error' not in r]),
            "failed_points": len([r for r in results if 'error' in r])
        }
    })
