#!/usr/bin/env python3
"""
Integration tests for fluids MCP server.
Tests the tools work correctly for typical engineering scenarios.
"""

import json
import pytest

# Import the tools (package installed via pip install -e .)
from tools.gas_pipe_pressure_drop import calculate_gas_pipe_pressure_drop
from tools.gas_control_valve import calculate_gas_control_valve
from tools.blower_compressor import calculate_blower_compressor_requirements
from tools.pipe_pressure_drop import calculate_pipe_pressure_drop
from tools.pump_requirements import calculate_pump_requirements
from tools.fluid_properties import get_fluid_properties


class TestWastewaterTreatmentScenarios:
    """Test realistic wastewater treatment scenarios."""
    
    def test_complete_pump_system_design(self):
        """Test complete pump system design for water transfer."""
        result = calculate_pump_requirements(
            fluid_name="Water",
            flow_rate_gpm=500,
            temperature_c=20,
            pipe_diameter_in=6,
            material="Steel",
            static_suction_head_ft=-10,
            static_discharge_head_ft=100,
            suction_pipe_length_ft=50,
            discharge_pipe_length_ft=1000,
            suction_fittings=[
                {"type": "90_elbow", "quantity": 2},
                {"type": "entrance_sharp", "quantity": 1}
            ],
            discharge_fittings=[
                {"type": "90_elbow", "quantity": 4},
                {"type": "check_valve_swing", "quantity": 1},
                {"type": "exit_normal", "quantity": 1}
            ]
        )
        result_dict = json.loads(result)
        
        # Should have all required pump sizing outputs
        assert "total_dynamic_head_ft" in result_dict
        assert "npsha_ft" in result_dict
        assert "estimated_power_kw" in result_dict
        
        # Values should be reasonable for this scenario
        tdh = result_dict["total_dynamic_head_ft"]
        npsha = result_dict["npsha_ft"]
        power = result_dict["estimated_power_kw"]
        
        assert 120 < tdh < 140, f"TDH {tdh} ft seems unrealistic"
        assert npsha > 15, f"NPSHa {npsha} ft too low"
        assert 15 < power < 25, f"Power {power} kW seems unrealistic"
    
    def test_aeration_system_design(self):
        """Test complete aeration system for wastewater treatment."""
        # Calculate blower requirements
        blower_result = calculate_blower_compressor_requirements(
            flow_rate_norm_m3_hr=5000,  # Large aeration system
            fluid_name="Air",
            inlet_pressure_psi=14.7,    # Atmospheric
            outlet_pressure_psi=22.0,   # For 5m water depth + losses
            inlet_temperature_c=20,
            efficiency=0.75
        )
        blower_dict = json.loads(blower_result)
        
        assert "required_power_kw" in blower_dict
        assert "actual_discharge_temp_c" in blower_dict
        
        power = blower_dict["required_power_kw"]
        temp = blower_dict["actual_discharge_temp_c"]
        
        # Typical aeration power: 20-30 kW per 1000 m³/hr
        assert 80 < power < 120, f"Blower power {power} kW unrealistic for aeration"
        assert 50 < temp < 80, f"Discharge temperature {temp}°C seems unrealistic"
    
    def test_biogas_system_design(self):
        """Test biogas piping and control system."""
        # Test biogas pipe pressure drop
        pipe_result = calculate_gas_pipe_pressure_drop(
            flow_rate_norm_m3_hr=500,   # Typical biogas flow
            fluid_name="Methane",       # Biogas approximation
            inlet_pressure=103000,      # ~1.5 psig
            pipe_length=200,            # 200m pipe run
            nominal_size_in=6,          # 6" pipe
            schedule="40",
            material="Steel",
            temperature_c=35            # Digester temperature
        )
        pipe_dict = json.loads(pipe_result)
        
        # Should calculate successfully without property lookup failures
        assert "pressure_drop_total_pa" in pipe_dict or "outlet_pressure" in pipe_dict
        
        # Test biogas control valve
        valve_result = calculate_gas_control_valve(
            flow_rate_norm_m3_hr=200,   # Smaller control flow
            fluid_name="Methane",
            inlet_pressure_psi=2.0,     # 2 psig
            outlet_pressure_psi=0.5,    # 0.5 psig
            inlet_temperature_c=35,
            valve_type="butterfly"
        )
        valve_dict = json.loads(valve_result)
        
        if "required_cv_us" in valve_dict:
            cv = valve_dict["required_cv_us"]
            # Should be reasonable for biogas control
            assert 1 < cv < 50, f"Cv {cv} unrealistic for biogas control"


class TestIndustrialScenarios:
    """Test industrial fluid system scenarios."""
    
    def test_water_cooling_system(self):
        """Test water cooling system pipe sizing."""
        result = calculate_pipe_pressure_drop(
            fluid_name="Water",
            flow_rate_gpm=1000,  # Large cooling system
            temperature_c=30,    # Warm return water
            pipe_length_ft=500,
            pipe_diameter_in=8,
            material="Steel",
            fittings=[
                {"type": "90_elbow", "quantity": 6},
                {"type": "gate_valve", "quantity": 3},
                {"type": "tee_run_through", "quantity": 2}
            ]
        )
        result_dict = json.loads(result)
        
        assert "pressure_drop_total_psi" in result_dict
        assert "flow_velocity_m_s" in result_dict
        
        velocity = result_dict["flow_velocity_m_s"]
        pressure_drop = result_dict["pressure_drop_total_psi"]
        
        # Typical cooling water velocities: 1.5-3 m/s
        assert 1.5 < velocity < 3.5, f"Velocity {velocity} m/s unrealistic"
        assert pressure_drop > 0
    
    def test_gas_distribution_system(self):
        """Test natural gas distribution piping."""
        result = calculate_gas_pipe_pressure_drop(
            flow_rate_norm_m3_hr=2000,  # Industrial gas flow
            fluid_name="Methane",
            inlet_pressure=300000,       # 3 bar
            pipe_length=1000,           # 1 km distribution
            nominal_size_in=8,          # Large pipe
            schedule="40",
            material="Steel",
            temperature_c=15,           # Ground temperature
            method="Weymouth"
        )
        result_dict = json.loads(result)
        
        # Should have valid pressure drop calculation
        assert "pressure_drop_total_pa" in result_dict or "outlet_pressure" in result_dict


class TestFluidPropertyAccuracy:
    """Test fluid property lookup accuracy."""
    
    def test_common_fluid_properties(self):
        """Test property lookup for common engineering fluids."""
        fluids_to_test = [
            ("Water", 20, 1.01325),
            ("Air", 20, 1.01325), 
            ("Methane", 25, 1.01325),
            ("CarbonDioxide", 20, 1.01325)
        ]
        
        for fluid_name, temp, pressure in fluids_to_test:
            result = get_fluid_properties(
                fluid_name=fluid_name,
                temperature_c=temp,
                pressure_bar=pressure
            )
            result_dict = json.loads(result)
            
            # Should not have lookup errors
            assert "errors" not in result_dict or not result_dict["errors"]
            
            # Should have essential properties
            assert "density_kg_m3" in result_dict
            assert "dynamic_viscosity_pa_s" in result_dict
            assert "molecular_weight_kg_mol" in result_dict
            
            # Values should be reasonable
            assert result_dict["density_kg_m3"] > 0
            assert result_dict["dynamic_viscosity_pa_s"] > 0
            assert result_dict["molecular_weight_kg_mol"] > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])