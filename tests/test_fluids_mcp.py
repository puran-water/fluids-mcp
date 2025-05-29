import json
import sys
import os

# Add the parent directory to the path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import all tools
from tools.gas_pipe_pressure_drop import calculate_gas_pipe_pressure_drop
from tools.gas_control_valve import calculate_gas_control_valve
from tools.blower_compressor import calculate_blower_compressor_requirements
from tools.fluid_properties import get_fluid_properties, list_available_fluids

def test_gas_tools():
    """Test all gas-related tools with typical wastewater treatment scenarios"""
    
    print("="*50)
    print("TESTING GAS TOOLS FOR WASTEWATER TREATMENT")
    print("="*50)
    
    # Test 1: Biogas pipe pressure drop
    print("\n1. Testing biogas pipe pressure drop (digester to boiler):")
    result = calculate_gas_pipe_pressure_drop(
        flow_rate_norm_m3_hr=500,  # Typical biogas flow
        fluid_name="Methane",       # Biogas approximation
        inlet_pressure=103000,      # ~1.5 psig
        pipe_length=200,            # 200m pipe run
        nominal_size_in=6,          # 6" pipe
        schedule="40",
        material="Steel",
        temperature_c=35,           # Digester temperature
        method="Weymouth"
    )
    print(f"Result: {json.loads(result).get('errors', 'SUCCESS')}")
    if 'pressure_drop_psi' in json.loads(result):
        print(f"Pressure drop: {json.loads(result)['pressure_drop_psi']:.2f} psi")
    
    # Test 2: Aeration blower requirements
    print("\n2. Testing aeration blower requirements:")
    result = calculate_blower_compressor_requirements(
        flow_rate_norm_m3_hr=5000,  # Typical aeration flow
        fluid_name="Air",
        inlet_pressure_psi=14.7,    # Atmospheric
        outlet_pressure_psi=22.0,   # ~7.3 psig for 5m water depth
        inlet_temperature_c=20,
        efficiency=0.75
    )
    result_dict = json.loads(result)
    print(f"Result: {result_dict.get('errors', 'SUCCESS')}")
    if 'required_power_kw' in result_dict:
        print(f"Required power: {result_dict['required_power_kw']:.1f} kW")
        print(f"Discharge temp: {result_dict['actual_discharge_temp_c']:.1f} °C")
    
    # Test 3: Gas control valve for digester pressure control
    print("\n3. Testing gas control valve sizing:")
    result = calculate_gas_control_valve(
        flow_rate_norm_m3_hr=200,   # Small biogas flow
        fluid_name="Methane",
        inlet_pressure_psi=2.0,     # 2 psig
        outlet_pressure_psi=0.5,    # 0.5 psig
        inlet_temperature_c=35,
        valve_type="butterfly"
    )
    result_dict = json.loads(result)
    print(f"Result: {result_dict.get('errors', 'SUCCESS')}")
    if 'required_cv_us' in result_dict:
        print(f"Required Cv: {result_dict['required_cv_us']:.1f}")
        print(f"Estimated size: {result_dict.get('estimated_valve_size', 'N/A')}")

def test_liquid_tools():
    """Test liquid tools for wastewater applications"""
    
    print("\n" + "="*50)
    print("TESTING LIQUID TOOLS FOR WASTEWATER TREATMENT")
    print("="*50)
    
    # Add liquid tool tests here...

if __name__ == "__main__":
    test_gas_tools()
    # test_liquid_tools()