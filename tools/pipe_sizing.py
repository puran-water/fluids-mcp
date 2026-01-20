"""
Pipe sizing tools for automatic selection based on pressure drop constraints.
"""

import json
import logging
from typing import Dict, Any, List, Optional, Union, Tuple
import fluids
from fluids import friction_factor, Reynolds
from fluids.piping import nearest_pipe
import CoolProp.CoolProp as CP
from CoolProp.HumidAirProp import HAPropsSI

# Import helpers from existing modules
from utils.helpers import get_fitting_K, get_pipe_roughness
from tools.pipe_pressure_drop import calculate_total_pressure_drop

logger = logging.getLogger(__name__)

# Default NPS candidates for Schedule 40
DEFAULT_NPS_CANDIDATES = [0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4, 6, 8, 10, 12]

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
        
        # Collect all parameters into kwargs for backwards compatibility
        kwargs = {k: v for k, v in locals().items() if k != 'kwargs' and v is not None}
        
        # Get NPS candidates
        nps_list = kwargs.get('nps_candidates', DEFAULT_NPS_CANDIDATES)
        schedule = kwargs.get('schedule', '40')
        
        # Convert flow rate to SI
        flow_rate_m3s = kwargs.get('flow_rate')
        if flow_rate_m3s is None and 'flow_rate_gpm' in kwargs:
            flow_rate_m3s = kwargs['flow_rate_gpm'] * 6.30902e-5
            log_messages.append(f"Converted flow_rate from {kwargs['flow_rate_gpm']} GPM to {flow_rate_m3s:.6f} m³/s")
        
        if flow_rate_m3s is None:
            return json.dumps({"error": "Flow rate must be specified"})
        
        # Convert pipe length to SI
        pipe_length_m = kwargs.get('pipe_length')
        if pipe_length_m is None and 'pipe_length_ft' in kwargs:
            pipe_length_m = kwargs['pipe_length_ft'] * 0.3048
            log_messages.append(f"Converted pipe_length from {kwargs['pipe_length_ft']} ft to {pipe_length_m:.2f} m")
        
        if pipe_length_m is None:
            return json.dumps({"error": "Pipe length must be specified"})
        
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
        
        # Get allowable pressure drop
        allowable_dp = kwargs.get('allowable_dp_pa')
        if allowable_dp is None and 'allowable_dp_psi' in kwargs:
            allowable_dp = kwargs['allowable_dp_psi'] * 6894.76
            log_messages.append(f"Converted allowable_dp from {kwargs['allowable_dp_psi']} psi")
        
        if allowable_dp is None:
            return json.dumps({"error": "Allowable pressure drop must be specified"})
        
        # Get velocity constraints
        v_min = kwargs.get('velocity_min_m_s', 0)
        v_max = kwargs.get('velocity_max_m_s', float('inf'))
        
        # Get pipe roughness
        pipe_roughness = kwargs.get('pipe_roughness')
        if pipe_roughness is None and 'material' in kwargs:
            pipe_roughness, source = get_pipe_roughness(material=kwargs['material'])
            log_messages.append(f"Looked up roughness for material '{kwargs['material']}': {pipe_roughness} m ({source})")
        if pipe_roughness is None:
            pipe_roughness = 1.5e-5  # Default for smooth pipe
            log_messages.append(f"Using default roughness: {pipe_roughness} m")
        
        # Get fittings
        fittings = kwargs.get('fittings', [])
        
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
                        if K_value is None:
                            K_value = 0.5  # Default if unknown
                    
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
                    'dp_total_psi': round(dp_total / 6894.76, 3),
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
                    
            except Exception as e:
                log_messages.append(f"Error evaluating NPS {nps}: {str(e)}")
                candidates.append({
                    'nps_in': nps,
                    'error': str(e)
                })
        
        # Prepare results
        result = {
            'inputs_resolved': log_messages,
            'allowable_dp_pa': allowable_dp,
            'allowable_dp_psi': allowable_dp / 6894.76,
            'flow_rate_m3s': flow_rate_m3s,
            'flow_rate_gpm': flow_rate_m3s / 6.30902e-5,
            'pipe_length_m': pipe_length_m,
            'schedule': schedule
        }
        
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
            
        return json.dumps(result, indent=2)
        
    except Exception as e:
        logger.error(f"Error in select_liquid_pipe_size_for_dp: {str(e)}")
        return json.dumps({"error": str(e)})


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
    method: str = "isothermal",
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
        method: Calculation method (default: "isothermal")
        velocity_max_m_s: Maximum velocity constraint
        return_all_candidates: Return all evaluated sizes
        
    Returns:
        JSON string with recommended size and evaluation results
    """
    try:
        log_messages = []
        
        # Collect all parameters into kwargs for backwards compatibility
        kwargs = {k: v for k, v in locals().items() if k != 'kwargs' and v is not None}
        
        # Get NPS candidates
        nps_list = kwargs.get('nps_candidates', DEFAULT_NPS_CANDIDATES)
        schedule = kwargs.get('schedule', '40')
        
        # Get temperature
        temp_c = kwargs.get('temperature_c')
        if temp_c is None:
            return json.dumps({"error": "Temperature must be specified"})
        temp_k = temp_c + 273.15
        
        # Convert inlet pressure to SI
        inlet_pressure = kwargs.get('inlet_pressure')
        if inlet_pressure is None and 'inlet_pressure_psi' in kwargs:
            inlet_pressure = kwargs['inlet_pressure_psi'] * 6894.76
            log_messages.append(f"Converted inlet_pressure from {kwargs['inlet_pressure_psi']} psi")
        
        if inlet_pressure is None:
            return json.dumps({"error": "Inlet pressure must be specified"})
        
        # Convert pipe length to SI
        pipe_length_m = kwargs.get('pipe_length')
        if pipe_length_m is None and 'pipe_length_ft' in kwargs:
            pipe_length_m = kwargs['pipe_length_ft'] * 0.3048
            log_messages.append(f"Converted pipe_length from {kwargs['pipe_length_ft']} ft")
        
        if pipe_length_m is None:
            return json.dumps({"error": "Pipe length must be specified"})
        
        # Get gas properties
        MW = kwargs.get('gas_mw')
        gamma = kwargs.get('gas_gamma')
        Z = kwargs.get('gas_z_factor', 1.0)
        viscosity = kwargs.get('gas_viscosity')
        
        # Try property lookup if not provided
        if 'fluid_name' in kwargs and (MW is None or gamma is None or viscosity is None):
            try:
                fluid_name = kwargs['fluid_name']
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
        
        # Get allowable pressure drop
        allowable_dp = kwargs.get('allowable_dp_pa')
        if allowable_dp is None and 'allowable_dp_psi' in kwargs:
            allowable_dp = kwargs['allowable_dp_psi'] * 6894.76
            log_messages.append(f"Converted allowable_dp from {kwargs['allowable_dp_psi']} psi")
        
        if allowable_dp is None:
            return json.dumps({"error": "Allowable pressure drop must be specified"})
        
        # Get velocity constraint
        v_max = kwargs.get('velocity_max_m_s', 30)  # Default 30 m/s for gas
        
        # Get pipe roughness
        pipe_roughness = kwargs.get('pipe_roughness')
        if pipe_roughness is None and 'material' in kwargs:
            pipe_roughness, source = get_pipe_roughness(material=kwargs['material'])
            log_messages.append(f"Looked up roughness for '{kwargs['material']}': {pipe_roughness} m ({source})")
        if pipe_roughness is None:
            pipe_roughness = 4.5e-5  # Default for steel
            log_messages.append(f"Using default roughness: {pipe_roughness} m")
        
        # Get fittings
        fittings = kwargs.get('fittings', [])
        method = kwargs.get('method', 'isothermal')
        
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
                            if K_value is None:
                                K_value = 0.5  # Default if unknown
                        
                        K_total += K_value * quantity
                    
                    L_eq = K_total * Di / f if f > 0 else 0
                    L_eff_new = pipe_length_m + L_eq
                    
                    # Calculate pressure drop with effective length
                    if method == 'isothermal':
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
                    'dp_total_psi': round(dp_total / 6894.76, 3),
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
                    
            except Exception as e:
                log_messages.append(f"Error evaluating NPS {nps}: {str(e)}")
                candidates.append({
                    'nps_in': nps,
                    'error': str(e)
                })
        
        # Prepare results
        result = {
            'inputs_resolved': log_messages,
            'gas_properties': {
                'MW_kg_kmol': MW,
                'gamma': gamma,
                'Z_factor': Z,
                'viscosity_pas': viscosity
            },
            'allowable_dp_pa': allowable_dp,
            'allowable_dp_psi': allowable_dp / 6894.76,
            'flow_rate_kg_s': flow_rate_kg_s,
            'inlet_pressure_pa': inlet_pressure,
            'pipe_length_m': pipe_length_m,
            'temperature_c': temp_c,
            'schedule': schedule,
            'method': method
        }
        
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
            
        return json.dumps(result, indent=2)
        
    except Exception as e:
        logger.error(f"Error in select_gas_pipe_size_for_dp: {str(e)}")
        return json.dumps({"error": str(e)})