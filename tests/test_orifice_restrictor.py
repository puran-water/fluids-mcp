#!/usr/bin/env python3
"""
Tests for orifice restrictor sizing functionality.

Validates the orifice_restrictor tool against:
- Manual calculations using Q = Cd * A * sqrt(2*dP/rho)
- User's RO permeate throttling examples
- ISO 5167 validity limits
"""

import json
import math
import pytest

from tools.orifice_restrictor import (
    calculate_orifice_sizing,
    calculate_orifice_analysis_dp,
    calculate_orifice_analysis_flow,
    generate_plate_kit,
    check_iso5167_validity,
)
from omnitools.orifice_restrictor import orifice_restrictor


class TestOrificeSizing:
    """Test orifice sizing calculations (Q + dP -> d)."""

    def test_ro_stage1_sizing(self):
        """
        User's RO Stage 1 example:
        87.2 m³/h @ 2.0 bar -> ~50mm orifice (with Cd=0.62, rho=1000)

        Note: ISO 5167 Reader-Harris-Gallagher Cd varies with beta/Re,
        so result may differ slightly from simple Cd=0.62 calculation.
        """
        result = calculate_orifice_sizing(
            flow_rate_m3_hr=87.2,
            pressure_drop_bar=2.0,
            nominal_size_in=4,  # 4" pipe
            schedule="40",
            fluid_name="Water",
            temperature_c=25,
            phase="liquid",
            size_units="mm"
        )
        result_dict = json.loads(result)

        # Should not have critical errors
        assert "error" not in result_dict, f"Sizing failed: {result_dict}"

        # Check orifice diameter is close to expected ~50mm
        d_mm = result_dict.get("orifice_diameter_mm")
        assert d_mm is not None, "Missing orifice_diameter_mm"

        # Allow ±15% tolerance for correlation differences
        assert 40 < d_mm < 60, f"Expected ~50mm, got {d_mm}mm"

        # Check beta ratio is reasonable
        beta = result_dict.get("beta_ratio")
        assert beta is not None
        assert 0.1 < beta < 0.75, f"Beta {beta} outside ISO 5167 range"

    def test_ro_stage2_sizing(self):
        """
        User's RO Stage 2 example:
        30.0 m³/h @ 1.0 bar -> ~35mm orifice
        """
        result = calculate_orifice_sizing(
            flow_rate_m3_hr=30.0,
            pressure_drop_bar=1.0,
            nominal_size_in=3,  # Assuming 3" header for stage 2
            schedule="40",
            fluid_name="Water",
            temperature_c=25,
            phase="liquid",
            size_units="mm"
        )
        result_dict = json.loads(result)

        assert "error" not in result_dict
        d_mm = result_dict.get("orifice_diameter_mm")
        assert d_mm is not None

        # Allow ±15% tolerance
        assert 28 < d_mm < 45, f"Expected ~35mm, got {d_mm}mm"

    def test_sizing_with_imperial_units(self):
        """Test sizing with GPM and psi inputs."""
        result = calculate_orifice_sizing(
            flow_rate_gpm=100,
            pressure_drop_psi=10,
            nominal_size_in=2,
            schedule="40",
            fluid_name="Water",
            temperature_c=25,
            phase="liquid",
            size_units="inch"
        )
        result_dict = json.loads(result)

        assert "error" not in result_dict
        assert "orifice_diameter_in" in result_dict
        assert result_dict["orifice_diameter_in"] > 0

    def test_sizing_with_direct_properties(self):
        """Test sizing with directly provided fluid properties."""
        result = calculate_orifice_sizing(
            flow_rate_m3_hr=50,
            pressure_drop_bar=1.5,
            pipe_diameter=0.1,  # 100mm
            fluid_density=998,  # kg/m³
            fluid_viscosity=0.001,  # Pa·s
            phase="liquid",
            size_units="mm"
        )
        result_dict = json.loads(result)

        assert "error" not in result_dict
        assert "orifice_diameter_mm" in result_dict

    def test_sizing_missing_flow_rate_error(self):
        """Test that missing flow rate returns error."""
        result = calculate_orifice_sizing(
            pressure_drop_bar=2.0,
            nominal_size_in=4,
            fluid_name="Water",
            temperature_c=25,
        )
        result_dict = json.loads(result)

        assert "errors" in result_dict
        assert any("flow_rate" in str(e).lower() for e in result_dict["errors"])


class TestOrificeAnalysisDP:
    """Test pressure drop analysis (Q + d -> dP)."""

    def test_analysis_dp_roundtrip(self):
        """
        Roundtrip test: Size orifice, then verify dP calculation matches.
        """
        # First, size an orifice
        sizing_result = json.loads(calculate_orifice_sizing(
            flow_rate_m3_hr=50,
            pressure_drop_bar=1.0,
            nominal_size_in=3,
            fluid_name="Water",
            temperature_c=25,
            phase="liquid",
        ))

        assert "error" not in sizing_result
        d_mm = sizing_result["orifice_diameter_mm"]

        # Now calculate dP with the sized orifice
        analysis_result = json.loads(calculate_orifice_analysis_dp(
            flow_rate_m3_hr=50,
            orifice_diameter_mm=d_mm,
            nominal_size_in=3,
            fluid_name="Water",
            temperature_c=25,
            phase="liquid",
        ))

        assert "error" not in analysis_result
        calc_dp_bar = analysis_result.get("pressure_drop_bar")
        assert calc_dp_bar is not None

        # Should be close to original target (within 10%)
        assert 0.9 < calc_dp_bar < 1.1, f"Expected ~1.0 bar, got {calc_dp_bar} bar"

    def test_analysis_dp_with_imperial(self):
        """Test dP analysis with imperial inputs."""
        result = calculate_orifice_analysis_dp(
            flow_rate_gpm=100,
            orifice_diameter_in=0.75,
            nominal_size_in=2,
            fluid_name="Water",
            temperature_c=25,
            phase="liquid",
        )
        result_dict = json.loads(result)

        assert "error" not in result_dict
        assert "pressure_drop_psi" in result_dict
        assert result_dict["pressure_drop_psi"] > 0


class TestOrificeAnalysisFlow:
    """Test flow rate analysis (d + dP -> Q)."""

    def test_analysis_flow_roundtrip(self):
        """
        Roundtrip test: Calculate flow from known orifice and dP.
        """
        # Known configuration
        d_mm = 50
        dp_bar = 2.0

        result = json.loads(calculate_orifice_analysis_flow(
            pressure_drop_bar=dp_bar,
            orifice_diameter_mm=d_mm,
            nominal_size_in=4,
            fluid_name="Water",
            temperature_c=25,
            phase="liquid",
        ))

        assert "error" not in result
        flow_m3_hr = result.get("flow_rate_m3_hr")
        assert flow_m3_hr is not None
        assert flow_m3_hr > 0

        # Verify with analysis_dp
        dp_check = json.loads(calculate_orifice_analysis_dp(
            flow_rate_m3_hr=flow_m3_hr,
            orifice_diameter_mm=d_mm,
            nominal_size_in=4,
            fluid_name="Water",
            temperature_c=25,
            phase="liquid",
        ))

        assert "error" not in dp_check
        dp_calc = dp_check.get("pressure_drop_bar")
        assert abs(dp_calc - dp_bar) < 0.05, f"Expected {dp_bar} bar, got {dp_calc} bar"


class TestPlateKit:
    """Test plate kit generation."""

    def test_plate_kit_generation(self):
        """Test generating a 5-plate kit with ±20% diameter variation."""
        result = json.loads(generate_plate_kit(
            flow_rate_m3_hr=87.2,
            pressure_drop_bar=2.0,
            nominal_size_in=4,
            fluid_name="Water",
            temperature_c=25,
            diameter_variability_percent=20.0,
            kit_sizes=5,
            phase="liquid",
        ))

        assert "error" not in result
        assert "plate_kit" in result
        assert len(result["plate_kit"]) == 5

        # Check design orifice is reasonable
        design_d_mm = result.get("design_orifice_mm")
        assert design_d_mm is not None
        assert 40 < design_d_mm < 60

        # Check plate kit structure
        for plate in result["plate_kit"]:
            assert "plate_number" in plate
            assert "orifice_diameter_mm" in plate
            assert "pressure_drop_bar" in plate or "pressure_drop_bar" == "Error"
            assert "beta_ratio" in plate

        # Design plate should be marked
        design_plates = [p for p in result["plate_kit"] if p.get("is_design_plate")]
        assert len(design_plates) == 1

        # Verify dP increases as diameter decreases
        diameters = [p["orifice_diameter_mm"] for p in result["plate_kit"]]
        dps = [p.get("pressure_drop_bar") for p in result["plate_kit"]]

        # Filter out errors
        valid_pairs = [(d, dp) for d, dp in zip(diameters, dps) if dp != "Error"]
        if len(valid_pairs) >= 2:
            # Smaller orifice = higher dP (ΔP ∝ 1/d⁴ approximately)
            for i in range(len(valid_pairs) - 1):
                d1, dp1 = valid_pairs[i]
                d2, dp2 = valid_pairs[i + 1]
                if d1 < d2:  # d1 smaller
                    assert dp1 > dp2, "Smaller orifice should have higher dP"

    def test_plate_kit_interpretation(self):
        """Test that plate kit includes interpretation text."""
        result = json.loads(generate_plate_kit(
            flow_rate_m3_hr=50,
            pressure_drop_bar=1.0,
            nominal_size_in=3,
            fluid_name="Water",
            temperature_c=25,
            kit_sizes=3,
            phase="liquid",
        ))

        assert "interpretation" in result
        assert len(result["interpretation"]) > 20  # Should have meaningful text


class TestOmnitoolInterface:
    """Test the unified omnitool interface."""

    def test_sizing_mode(self):
        """Test omnitool in sizing mode."""
        result = json.loads(orifice_restrictor(
            mode="sizing",
            phase="liquid",
            flow_rate_m3_hr=50,
            pressure_drop_bar=1.0,
            nominal_size_in=3,
            fluid_name="Water",
            temperature_c=25,
        ))

        assert "error" not in result
        assert result.get("mode") == "sizing"
        assert "orifice_diameter_mm" in result

    def test_analysis_dp_mode(self):
        """Test omnitool in analysis_dp mode."""
        result = json.loads(orifice_restrictor(
            mode="analysis_dp",
            phase="liquid",
            flow_rate_m3_hr=50,
            orifice_diameter_mm=40,
            nominal_size_in=3,
            fluid_name="Water",
            temperature_c=25,
        ))

        assert "error" not in result
        assert result.get("mode") == "analysis_dp"
        assert "pressure_drop_bar" in result

    def test_plate_kit_mode(self):
        """Test omnitool in plate_kit mode."""
        result = json.loads(orifice_restrictor(
            mode="plate_kit",
            phase="liquid",
            flow_rate_m3_hr=50,
            pressure_drop_bar=1.0,
            nominal_size_in=3,
            fluid_name="Water",
            temperature_c=25,
        ))

        assert "error" not in result
        assert result.get("mode") == "plate_kit"
        assert "plate_kit" in result

    def test_invalid_mode_error(self):
        """Test that invalid mode returns error."""
        result = json.loads(orifice_restrictor(
            mode="invalid_mode",
            flow_rate_m3_hr=50,
            pressure_drop_bar=1.0,
        ))

        assert "error" in result
        assert "invalid" in result["error"].lower()


class TestISO5167Validity:
    """Test ISO 5167 validity checking."""

    def test_valid_configuration(self):
        """Test that valid configuration has no warnings."""
        warnings = check_iso5167_validity(
            D=0.1,      # 100mm pipe
            D2=0.05,    # 50mm orifice -> beta=0.5
            dP=100000,  # 100 kPa (1 bar)
            Re_D=100000,
            taps='D'
        )
        assert len(warnings) == 0

    def test_beta_too_high(self):
        """Test warning for beta > 0.75."""
        warnings = check_iso5167_validity(
            D=0.1,
            D2=0.08,  # beta = 0.8
            dP=100000,
        )
        assert any("beta" in w.lower() for w in warnings)

    def test_beta_too_low(self):
        """Test warning for beta < 0.1."""
        warnings = check_iso5167_validity(
            D=0.1,
            D2=0.005,  # beta = 0.05
            dP=100000,
        )
        assert any("beta" in w.lower() for w in warnings)

    def test_dp_too_high(self):
        """Test warning for dP > 250 kPa."""
        warnings = check_iso5167_validity(
            D=0.1,
            D2=0.05,
            dP=300000,  # 300 kPa
        )
        assert any("250" in w for w in warnings)

    def test_orifice_too_small(self):
        """Test warning for orifice < 12.5mm."""
        warnings = check_iso5167_validity(
            D=0.025,  # 25mm pipe
            D2=0.01,  # 10mm orifice
            dP=100000,
        )
        assert any("12.5" in w for w in warnings)


class TestManualCalculationValidation:
    """
    Validate against manual orifice calculations.

    Manual formula: Q = Cd * A * sqrt(2 * dP / rho)
    Solving for d: d = sqrt(4 * Q / (Cd * pi * sqrt(2 * dP / rho)))
    """

    def test_manual_calculation_match(self):
        """
        Compare tool result to manual calculation with Cd=0.62.

        Note: Tool uses ISO 5167 Cd (varies with beta/Re), so expect
        some difference, but should be within engineering tolerance.
        """
        # Given conditions
        Q_m3_hr = 87.2
        Q_m3_s = Q_m3_hr / 3600  # 0.02422 m³/s
        dP_Pa = 200000  # 2 bar
        rho = 1000  # kg/m³
        Cd_simple = 0.62

        # Manual calculation
        # Q = Cd * A * sqrt(2*dP/rho)
        # A = Q / (Cd * sqrt(2*dP/rho))
        # d = sqrt(4*A/pi)

        sqrt_term = math.sqrt(2 * dP_Pa / rho)
        A_manual = Q_m3_s / (Cd_simple * sqrt_term)
        d_manual = math.sqrt(4 * A_manual / math.pi)
        d_manual_mm = d_manual * 1000

        # Tool calculation
        result = json.loads(calculate_orifice_sizing(
            flow_rate_m3_hr=Q_m3_hr,
            pressure_drop_bar=2.0,
            pipe_diameter=0.1,  # 100mm pipe
            fluid_density=rho,
            fluid_viscosity=0.001,  # 1 cP
            phase="liquid",
            size_units="mm"
        ))

        assert "error" not in result
        d_tool_mm = result["orifice_diameter_mm"]

        # ISO 5167 Cd varies with beta (typically 0.59-0.65 for this range)
        # Allow up to 15% difference from simple Cd=0.62 calculation
        percent_diff = abs(d_tool_mm - d_manual_mm) / d_manual_mm * 100
        assert percent_diff < 15, (
            f"Tool result {d_tool_mm:.1f}mm differs from manual {d_manual_mm:.1f}mm "
            f"by {percent_diff:.1f}%"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
