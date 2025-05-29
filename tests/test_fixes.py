#!/usr/bin/env python3
"""
Test script to verify fixes for fluids MCP server issues.
Run this to test the fixes for gamma attribute access, density calculations, and valve sizing.
"""

import json
import sys
import os

# Add the parent directory to the path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the fixed tools
from tools.gas_pipe_pressure_drop import calculate_gas_pipe_pressure_drop
from tools.gas_control_valve import calculate_gas_control_valve
from tools.blower_compressor import calculate_blower_compressor_requirements

def test_gas_tools_fixed():
    """Test the fixed gas tools with scenarios that previously failed"""
    
    print("="*60)
    print("TESTING FIXED GAS TOOLS")
    print("="*60)
    
    # Test 1: Gas pipe pressure drop with Methane (previously failed with property lookup)
    print("\n1. Testing Methane pipe pressure drop (should work with property lookup now):")
    try:
        result = calculate_gas_pipe_pressure_drop(
            flow_rate_norm_m3_hr=1000,
            fluid_name="Methane",
            inlet_pressure=500000,  # 5 bar
            pipe_length=1000,
            nominal_size_in=4,
            schedule="40",
            material="Steel",
            temperature_c=20
        )
        result_dict = json.loads(result)
        if "errors" in result_dict:
            print(f"  WARNING: {result_dict['errors']}")
        else:
            print("  SUCCESS: No errors")
        if "pressure_drop_total_pa" in result_dict:
            print(f"  Pressure drop: {result_dict['pressure_drop_total_pa']:.0f} Pa")
    except Exception as e:
        print(f"  ERROR: {e}")
    
    # Test 2: Blower compressor with Air (check density calculation)
    print("\n2. Testing Air blower requirements (check flow rate conversion):")
    try:
        result = calculate_blower_compressor_requirements(
            flow_rate_norm_m3_hr=1000,
            fluid_name="Air",
            inlet_pressure_psi=14.7,
            outlet_pressure_psi=29.4,
            inlet_temperature_c=20,
            efficiency=0.75
        )
        result_dict = json.loads(result)
        if "errors" in result_dict:
            print(f"  WARNING: {result_dict['errors']}")
        else:
            print("  SUCCESS: No errors")
        
        if "flow_rate_kg_s" in result_dict:
            kg_s = result_dict["flow_rate_kg_s"]
            print(f"  Flow rate: {kg_s:.3f} kg/s (should be ~0.36, not 0.0004)")
            if 0.3 < kg_s < 0.4:
                print("  SUCCESS: Flow rate is reasonable now!")
            else:
                print("  ERROR: Flow rate still seems wrong")
        
        if "required_power_kw" in result_dict:
            print(f"  Power: {result_dict['required_power_kw']:.1f} kW")
    except Exception as e:
        print(f"  ERROR: {e}")
    
    # Test 3: Gas control valve sizing (check Cv calculation)
    print("\n3. Testing Air control valve sizing (check Cv calculation):")
    try:
        result = calculate_gas_control_valve(
            flow_rate_norm_m3_hr=1000,
            fluid_name="Air",
            inlet_pressure_psi=100,
            outlet_pressure_psi=50,
            inlet_temperature_c=20,
            valve_type="globe"
        )
        result_dict = json.loads(result)
        if "errors" in result_dict:
            print(f"  WARNING: {result_dict['errors']}")
        else:
            print("  SUCCESS: No errors")
        
        if "required_cv_us" in result_dict:
            cv = result_dict["required_cv_us"]
            print(f"  Required Cv: {cv:.1f} (should be reasonable, not 41,000+)")
            if cv < 1000:
                print("  SUCCESS: Cv value is reasonable now!")
            else:
                print("  ERROR: Cv value still seems too high")
        
        if "estimated_valve_size" in result_dict:
            print(f"  Estimated size: {result_dict['estimated_valve_size']}")
    except Exception as e:
        print(f"  ERROR: {e}")
    
    print("\n" + "="*60)
    print("TEST COMPLETE")
    print("="*60)

if __name__ == "__main__":
    test_gas_tools_fixed()