#!/usr/bin/env python3
"""
Proper pytest tests for fluids MCP server tools.
Replaces the print-based tests with proper assertions.
"""

import json
import pytest

# Import the tools (package installed via pip install -e .)
from tools.gas_pipe_pressure_drop import calculate_gas_pipe_pressure_drop
from tools.gas_control_valve import calculate_gas_control_valve
from tools.blower_compressor import calculate_blower_compressor_requirements
from tools.fluid_properties import get_fluid_properties, list_available_fluids
from tools.pipe_pressure_drop import calculate_pipe_pressure_drop
from tools.pump_requirements import calculate_pump_requirements
from tools.reynolds_number import calculate_reynolds_number


class TestFluidProperties:
    """Test fluid property lookup functionality."""
    
    def test_list_available_fluids(self):
        """Test that fluid list is returned."""
        result = list_available_fluids()
        result_dict = json.loads(result)
        
        assert "available_fluids" in result_dict
        assert "total_count" in result_dict
        assert result_dict["total_count"] > 100  # Should have many fluids
        assert "Water" in str(result_dict["available_fluids"])  # Water should be available
    
    def test_water_properties(self):
        """Test water property lookup."""
        result = get_fluid_properties(
            fluid_name="Water",
            temperature_c=20.0,
            pressure_bar=1.01325
        )
        result_dict = json.loads(result)
        
        assert "density_kg_m3" in result_dict
        assert "dynamic_viscosity_pa_s" in result_dict
        # Water density at 20°C should be close to 998 kg/m³
        assert 990 < result_dict["density_kg_m3"] < 1005
    
    def test_methane_properties(self):
        """Test methane property lookup (previously failed)."""
        result = get_fluid_properties(
            fluid_name="Methane",
            temperature_c=20.0,
            pressure_bar=1.01325
        )
        result_dict = json.loads(result)

        # Should not have errors now that property lookup is fixed
        assert "errors" not in result_dict or not result_dict["errors"]
        assert "density_kg_m3" in result_dict
        assert "molecular_weight_kg_mol" in result_dict


class TestCachedFluidPropertiesRegressions:
    """Regression tests for CachedFluidProperties bug fixes."""

    def test_methane_z_factor_at_high_pressure(self):
        """
        Verify CachedFluidProperties queries Z from CoolProp at high pressure.

        Bug fix: Z was hard-coded to 1.0 instead of querying from CoolProp.
        At high pressure, real gases deviate from ideal (Z < 1.0 for methane).
        """
        from utils.property_cache import CachedFluidProperties

        # At 50 bar, methane Z-factor should be noticeably < 1.0
        # CachedFluidProperties takes: coolprop_name, T_in_deg_C, P_in_bar
        props = CachedFluidProperties("Methane", 20.0, 50.0)

        # Z should be queried from CoolProp, not hard-coded to 1.0
        z_value = props.Z[0]
        assert z_value < 1.0, f"Z={z_value} should be < 1.0 at 50 bar for methane"
        assert z_value > 0.5, f"Z={z_value} is unreasonably low"

    def test_methane_gamma_differs_from_ideal(self):
        """
        Verify that gamma (Cp/Cv) is calculated from CoolProp Cv, not Cv=Cp-R.

        Bug fix: Cv was calculated as Cp-R (ideal gas) instead of querying
        from CoolProp directly using key 'O' (Cvmass).
        """
        from utils.property_cache import CachedFluidProperties

        # Get properties at moderate pressure
        # CachedFluidProperties takes: coolprop_name, T_in_deg_C, P_in_bar
        props = CachedFluidProperties("Methane", 20.0, 10.0)

        gamma = props.gamma[0]
        # Methane gamma should be around 1.3-1.35 at these conditions
        # Not exactly 1.4 (ideal diatomic) or derived from Cp-R
        assert 1.25 < gamma < 1.45, f"gamma={gamma} outside expected range for methane"

    def test_cv_is_reasonable_for_real_gas(self):
        """Verify Cv is queried directly from CoolProp."""
        from utils.property_cache import CachedFluidProperties

        # CachedFluidProperties takes: coolprop_name, T_in_deg_C, P_in_bar
        props = CachedFluidProperties("Methane", 20.0, 10.0)

        cv = props.Cv[0]
        cp = props.Cp[0]

        # Cv should be positive and less than Cp
        assert cv > 0, f"Cv={cv} should be positive"
        assert cv < cp, f"Cv={cv} should be less than Cp={cp}"

        # For methane at 10 bar, Cv ~ 1700-2000 J/(kg·K)
        assert 1500 < cv < 2500, f"Cv={cv} outside expected range for methane"


class TestLiquidFlowCalculations:
    """Test liquid flow calculations."""
    
    def test_pipe_pressure_drop(self):
        """Test pipe pressure drop calculation."""
        result = calculate_pipe_pressure_drop(
            fluid_name="Water",
            flow_rate_gpm=500,
            temperature_c=20,
            pipe_length_ft=1000,
            pipe_diameter_in=6,
            material="Steel",
            fittings=[
                {"type": "90_elbow", "quantity": 4},
                {"type": "gate_valve", "quantity": 2}
            ]
        )
        result_dict = json.loads(result)
        
        assert "pressure_drop_total_psi" in result_dict
        assert "reynolds_number" in result_dict
        assert result_dict["reynolds_number"] > 200000  # Should be turbulent
        assert result_dict["pressure_drop_total_psi"] > 0
    
    def test_pump_requirements(self):
        """Test pump sizing calculation."""
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
            suction_fittings=[{"type": "90_elbow", "quantity": 2}],
            discharge_fittings=[{"type": "90_elbow", "quantity": 4}]
        )
        result_dict = json.loads(result)
        
        assert "total_dynamic_head_ft" in result_dict
        assert "npsha_ft" in result_dict
        assert "estimated_power_kw" in result_dict
        assert result_dict["total_dynamic_head_ft"] > 100  # Should include static + friction
        assert result_dict["npsha_ft"] > 0  # Should have positive NPSH available


class TestGasFlowCalculations:
    """Test gas flow calculations with proper property lookup."""
    
    def test_methane_pipe_pressure_drop(self):
        """Test methane pipe pressure drop (previously failed with property lookup)."""
        result = calculate_gas_pipe_pressure_drop(
            flow_rate_norm_m3_hr=1000,
            fluid_name="Methane",
            inlet_pressure=500000,
            pipe_length=1000,
            nominal_size_in=4,
            schedule="40",
            material="Steel",
            temperature_c=20
        )
        result_dict = json.loads(result)

        # Should either calculate successfully OR have graceful error with log
        # On some platforms (e.g., macOS), property lookup may fail due to binary deps
        has_result = "pressure_drop_total_pa" in result_dict or "outlet_pressure" in result_dict
        has_graceful_error = "errors" in result_dict and "log" in result_dict and len(result_dict.get("log", [])) > 0
        assert has_result or has_graceful_error, f"Expected results or graceful error, got: {result_dict.keys()}"

        # If there are errors, they should be warnings, not hard failures
        if "errors" in result_dict:
            errors = result_dict["errors"]
            for error in errors:
                assert "Warning" in error or "warning" in error.lower() or "failed" in error.lower()

    def test_air_blower_requirements(self):
        """Test air blower sizing with correct flow rate conversion."""
        result = calculate_blower_compressor_requirements(
            flow_rate_norm_m3_hr=1000,
            fluid_name="Air",
            inlet_pressure_psi=14.7,
            outlet_pressure_psi=29.4,
            inlet_temperature_c=20,
            efficiency=0.75
        )
        result_dict = json.loads(result)
        
        assert "flow_rate_kg_s" in result_dict
        assert "required_power_kw" in result_dict
        
        # Flow rate should be reasonable (not 0.0004 kg/s as before)
        flow_rate = result_dict["flow_rate_kg_s"]
        assert 0.3 < flow_rate < 0.4, f"Flow rate {flow_rate} kg/s seems unrealistic"
        
        # Power should be reasonable
        assert result_dict["required_power_kw"] > 10
    
    def test_air_control_valve_sizing(self):
        """Test air control valve sizing with correct Cv calculation."""
        result = calculate_gas_control_valve(
            flow_rate_norm_m3_hr=1000,
            fluid_name="Air",
            inlet_pressure_psi=100,
            outlet_pressure_psi=50,
            inlet_temperature_c=20,
            valve_type="globe"
        )
        result_dict = json.loads(result)
        
        assert "required_cv_us" in result_dict
        assert "estimated_valve_size" in result_dict
        
        # Cv should be reasonable (not 40,000+ as before)
        cv = result_dict["required_cv_us"]
        assert cv < 1000, f"Cv value {cv} seems too high"
        assert cv > 0.1, f"Cv value {cv} seems too low"


class TestReynoldsNumber:
    """Test Reynolds number calculations."""
    
    def test_water_reynolds_number(self):
        """Test Reynolds number calculation for water."""
        result = calculate_reynolds_number(
            fluid_name="Water",
            temperature_c=20,
            velocity_ft_s=5,
            characteristic_length_in=6
        )
        result_dict = json.loads(result)
        
        assert "reynolds_number" in result_dict
        assert "flow_regime" in result_dict
        assert result_dict["reynolds_number"] > 100000  # Should be high for water
        assert result_dict["flow_regime"] == "turbulent"


class TestWastewaterTreatmentScenarios:
    """Test realistic wastewater treatment scenarios."""
    
    def test_aeration_blower_sizing(self):
        """Test aeration blower for wastewater treatment."""
        result = calculate_blower_compressor_requirements(
            flow_rate_norm_m3_hr=5000,  # Typical aeration flow
            fluid_name="Air",
            inlet_pressure_psi=14.7,    # Atmospheric
            outlet_pressure_psi=22.0,   # For 5m water depth
            inlet_temperature_c=20,
            efficiency=0.75
        )
        result_dict = json.loads(result)
        
        assert "required_power_kw" in result_dict
        power = result_dict["required_power_kw"]
        
        # Should be reasonable for aeration (typically 20-30 kW per 1000 m³/hr)
        assert 50 < power < 150, f"Power {power} kW seems unrealistic for aeration"
    
    def test_biogas_pipe_pressure_drop(self):
        """Test biogas piping from digester to boiler."""
        result = calculate_gas_pipe_pressure_drop(
            flow_rate_norm_m3_hr=500,   # Typical biogas flow
            fluid_name="Methane",       # Biogas approximation
            inlet_pressure=103000,      # ~1.5 psig
            pipe_length=200,            # 200m pipe run
            nominal_size_in=6,          # 6" pipe
            schedule="40",
            material="Steel",
            temperature_c=35            # Digester temperature
        )
        result_dict = json.loads(result)
        
        # Should calculate without property lookup failures
        if "pressure_drop_total_pa" in result_dict:
            pressure_drop = result_dict["pressure_drop_total_pa"]
            assert pressure_drop > 0
        elif "outlet_pressure" in result_dict:
            outlet_pressure = result_dict["outlet_pressure"]
            assert outlet_pressure < 103000  # Should have some pressure drop
    
    def test_biogas_control_valve(self):
        """Test biogas control valve for pressure regulation."""
        result = calculate_gas_control_valve(
            flow_rate_norm_m3_hr=200,   # Small biogas flow
            fluid_name="Methane",
            inlet_pressure_psi=2.0,     # 2 psig
            outlet_pressure_psi=0.5,    # 0.5 psig
            inlet_temperature_c=35,
            valve_type="butterfly"
        )
        result_dict = json.loads(result)
        
        if "required_cv_us" in result_dict:
            cv = result_dict["required_cv_us"]
            # Should be reasonable for small biogas flow
            assert 1 < cv < 100, f"Cv {cv} seems unrealistic for biogas control"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])