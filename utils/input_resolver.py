"""
Shared input resolution system using Pydantic models.

This module eliminates code duplication across tools by providing standardized
input validation, unit conversion, and default handling.
"""

from typing import Optional, Dict, Any, Union, List
from pydantic import BaseModel, Field, validator, root_validator
import logging

from .constants import (
    GPM_to_M3S, INCH_to_M, FT_to_M, PSI_to_PA, 
    LBFT3_to_KGM3, CENTIPOISE_to_PAS, DEG_C_to_K, G_GRAVITY
)

logger = logging.getLogger("fluids-mcp.input_resolver")


class FlowRateInput(BaseModel):
    """Standardized flow rate input with automatic unit conversion."""
    
    # SI units (preferred)
    m3_s: Optional[float] = Field(None, description="Flow rate in m³/s")
    kg_s: Optional[float] = Field(None, description="Mass flow rate in kg/s")
    
    # Imperial units
    gpm: Optional[float] = Field(None, description="Flow rate in US GPM")
    
    # Standard conditions
    norm_m3_hr: Optional[float] = Field(None, description="Flow rate at Normal conditions (0°C, 1 atm) in m³/hr")
    std_m3_hr: Optional[float] = Field(None, description="Flow rate at Standard conditions (15°C, 1 atm) in m³/hr")

    def get_m3_s(self) -> Optional[float]:
        """Get flow rate in m³/s, converting from other units if needed."""
        if self.m3_s is not None:
            return self.m3_s
        elif self.gpm is not None:
            return self.gpm * GPM_to_M3S
        elif self.norm_m3_hr is not None:
            return self.norm_m3_hr / 3600.0  # Convert hr to s
        elif self.std_m3_hr is not None:
            return self.std_m3_hr / 3600.0  # Convert hr to s
        return None

    def get_kg_s(self) -> Optional[float]:
        """Get mass flow rate in kg/s."""
        return self.kg_s

    def get_source(self) -> str:
        """Get description of which input was used."""
        if self.m3_s is not None:
            return "SI (m³/s)"
        elif self.gpm is not None:
            return f"Imperial (GPM) -> {self.get_m3_s():.6f} m³/s"
        elif self.norm_m3_hr is not None:
            return f"Normal conditions (m³/hr) -> {self.get_m3_s():.6f} m³/s"
        elif self.std_m3_hr is not None:
            return f"Standard conditions (m³/hr) -> {self.get_m3_s():.6f} m³/s"
        elif self.kg_s is not None:
            return "SI (kg/s)"
        return "Not specified"


class PressureInput(BaseModel):
    """Standardized pressure input with automatic unit conversion."""
    
    # SI units (preferred)
    pa: Optional[float] = Field(None, description="Pressure in Pa (absolute)")
    
    # Imperial units  
    psi: Optional[float] = Field(None, description="Pressure in psi (absolute)")

    def get_pa(self) -> Optional[float]:
        """Get pressure in Pa, converting from other units if needed."""
        if self.pa is not None:
            return self.pa
        elif self.psi is not None:
            return self.psi * PSI_to_PA
        return None

    def get_source(self) -> str:
        """Get description of which input was used."""
        if self.pa is not None:
            return "SI (Pa)"
        elif self.psi is not None:
            return f"Imperial (psi) -> {self.get_pa():.1f} Pa"
        return "Not specified"


class DimensionInput(BaseModel):
    """Standardized dimension input with automatic unit conversion."""
    
    # SI units (preferred)
    m: Optional[float] = Field(None, description="Dimension in meters")
    
    # Imperial units
    inch: Optional[float] = Field(None, description="Dimension in inches")
    ft: Optional[float] = Field(None, description="Dimension in feet")

    def get_m(self) -> Optional[float]:
        """Get dimension in meters, converting from other units if needed."""
        if self.m is not None:
            return self.m
        elif self.inch is not None:
            return self.inch * INCH_to_M
        elif self.ft is not None:
            return self.ft * FT_to_M
        return None

    def get_source(self) -> str:
        """Get description of which input was used."""
        if self.m is not None:
            return "SI (m)"
        elif self.inch is not None:
            return f"Imperial (in) -> {self.get_m():.6f} m"
        elif self.ft is not None:
            return f"Imperial (ft) -> {self.get_m():.6f} m"
        return "Not specified"


class GasPropertiesInput(BaseModel):
    """Standardized gas properties input with validation and defaults."""
    
    # Direct property inputs
    mw: Optional[float] = Field(None, description="Molecular weight in kg/kmol", gt=0)
    gamma: Optional[float] = Field(None, description="Specific heat ratio Cp/Cv", gt=1.0, le=2.0)
    z_factor: Optional[float] = Field(None, description="Compressibility factor", gt=0, le=2.0)
    viscosity: Optional[float] = Field(None, description="Dynamic viscosity in Pa·s", gt=0)
    
    # Lookup inputs
    fluid_name: Optional[str] = Field(None, description="Fluid name for property lookup")
    composition_mol: Optional[Dict[str, float]] = Field(None, description="Molar composition")
    
    # Default handling
    allow_defaults: bool = Field(False, description="Allow default properties (may cause >10% error)")

    @validator('composition_mol')
    def validate_composition(cls, v):
        """Validate that molar composition sums to approximately 1.0."""
        if v is not None:
            total = sum(v.values())
            if not (0.95 <= total <= 1.05):
                raise ValueError(f"Molar composition must sum to ~1.0, got {total:.3f}")
        return v

    def resolve_properties(self, temperature_c: float, pressure_pa: Optional[float] = None) -> Dict[str, Any]:
        """
        Resolve gas properties with clear error handling.
        
        Returns dict with:
        - mw: Molecular weight in kg/kmol
        - gamma: Specific heat ratio
        - z_factor: Compressibility factor
        - viscosity: Dynamic viscosity in Pa·s (if available)
        - source: Description of property source
        - warnings: List of warnings
        - errors: List of errors
        """
        result = {
            "mw": self.mw,
            "gamma": self.gamma, 
            "z_factor": self.z_factor,
            "viscosity": self.viscosity,
            "source": "Direct input",
            "warnings": [],
            "errors": []
        }
        
        # Try property lookup if direct values missing
        if any(x is None for x in [self.mw, self.gamma, self.z_factor]) and self.fluid_name:
            try:
                from ..utils.import_helpers import FLUIDPROP_AVAILABLE, FluidProperties
                if FLUIDPROP_AVAILABLE and FluidProperties is not None:
                    lookup_p_bar = (pressure_pa or 101325.0) / 100000.0
                    fluid_props = FluidProperties(
                        coolprop_name=self.fluid_name,
                        T_in_deg_C=temperature_c,
                        P_in_bar=lookup_p_bar
                    )
                    
                    if self.mw is None:
                        result["mw"] = float(fluid_props.MW) * 1000.0  # Convert to kg/kmol
                    if self.gamma is None and hasattr(fluid_props, 'gamma'):
                        result["gamma"] = float(fluid_props.gamma[0])
                    if self.z_factor is None and hasattr(fluid_props, 'Z'):
                        result["z_factor"] = float(fluid_props.Z[0])
                    if self.viscosity is None and hasattr(fluid_props, 'eta'):
                        result["viscosity"] = float(fluid_props.eta[0])
                    
                    result["source"] = f"Lookup ({self.fluid_name})"
                else:
                    result["warnings"].append("FluidProp unavailable for property lookup")
            except Exception as e:
                result["warnings"].append(f"Property lookup failed: {e}")
        
        # Apply defaults if explicitly allowed
        if self.allow_defaults:
            if result["mw"] is None:
                result["mw"] = 28.96  # Air
                result["warnings"].append("Using default MW=28.96 (Air) - may cause >10% error for non-air gases")
            if result["gamma"] is None:
                result["gamma"] = 1.4  # Air
                result["warnings"].append("Using default gamma=1.4 (Air) - may cause >10% error for non-air gases")
            if result["z_factor"] is None:
                result["z_factor"] = 1.0  # Ideal gas
                result["warnings"].append("Using default Z=1.0 (Ideal) - may cause errors at high pressures")
        else:
            # Fail fast if critical properties missing
            if result["mw"] is None:
                result["errors"].append("Missing gas molecular weight. Provide mw, fluid_name, or set allow_defaults=True")
            if result["gamma"] is None:
                result["errors"].append("Missing gas specific heat ratio. Provide gamma, fluid_name, or set allow_defaults=True")
            if result["z_factor"] is None:
                result["errors"].append("Missing gas compressibility factor. Provide z_factor, fluid_name, or set allow_defaults=True")
        
        return result


class InputResolver:
    """
    Centralized input resolution with consistent logging and error handling.
    
    Eliminates code duplication across all tools.
    """
    
    def __init__(self, tool_name: str):
        self.tool_name = tool_name
        self.results_log: List[str] = []
        self.error_log: List[str] = []

    def resolve_flow_rate(self, **kwargs) -> Optional[float]:
        """Resolve flow rate to m³/s with logging."""
        flow_input = FlowRateInput(**kwargs)
        result = flow_input.get_m3_s()
        
        if result is not None:
            self.results_log.append(f"Flow rate: {flow_input.get_source()}")
        else:
            self.error_log.append("Missing flow rate input")
        
        return result

    def resolve_mass_flow_rate(self, **kwargs) -> Optional[float]:
        """Resolve mass flow rate to kg/s with logging."""
        flow_input = FlowRateInput(**kwargs)
        result = flow_input.get_kg_s()
        
        if result is not None:
            self.results_log.append(f"Mass flow rate: {flow_input.get_source()}")
        else:
            self.error_log.append("Missing mass flow rate input")
        
        return result

    def resolve_pressure(self, name: str, **kwargs) -> Optional[float]:
        """Resolve pressure to Pa with logging."""
        pressure_input = PressureInput(**kwargs)
        result = pressure_input.get_pa()
        
        if result is not None:
            self.results_log.append(f"{name}: {pressure_input.get_source()}")
        else:
            self.error_log.append(f"Missing {name.lower()} input")
        
        return result

    def resolve_dimension(self, name: str, **kwargs) -> Optional[float]:
        """Resolve dimension to meters with logging."""
        dim_input = DimensionInput(**kwargs)
        result = dim_input.get_m()
        
        if result is not None:
            self.results_log.append(f"{name}: {dim_input.get_source()}")
        else:
            self.error_log.append(f"Missing {name.lower()} input")
        
        return result

    def resolve_gas_properties(self, temperature_c: float, pressure_pa: Optional[float] = None, **kwargs) -> Dict[str, Any]:
        """Resolve gas properties with comprehensive validation."""
        gas_input = GasPropertiesInput(**kwargs)
        result = gas_input.resolve_properties(temperature_c, pressure_pa)
        
        # Add to our logs
        self.results_log.append(f"Gas properties source: {result['source']}")
        self.results_log.extend(result["warnings"])
        self.error_log.extend(result["errors"])
        
        return result

    def get_logs(self) -> Dict[str, List[str]]:
        """Get accumulated logs."""
        return {
            "log": self.results_log.copy(),
            "errors": self.error_log.copy()
        }