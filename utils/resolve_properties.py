"""
Centralized fluid property resolver.

Single shared resolver that all tools call. Resolution order:
1. CoolProp via CachedFluidProperties (primary — always available)
2. FluidProperties via fluidprop (secondary — optional plugin)
3. thermo.Chemical fallback (tertiary — for NaN/missing values)
"""

import logging
import math
from dataclasses import dataclass, field
from typing import Optional, List

logger = logging.getLogger("fluids-mcp.resolve_properties")


@dataclass
class ResolvedLiquidProperties:
    """Resolved liquid properties with provenance tracking."""
    density: Optional[float] = None
    viscosity: Optional[float] = None
    kinematic_viscosity: Optional[float] = None
    vapor_pressure: Optional[float] = None
    critical_pressure: Optional[float] = None
    source: Optional[str] = None
    log: List[str] = field(default_factory=list)


@dataclass
class ResolvedGasProperties:
    """Resolved gas properties with provenance tracking."""
    mw: Optional[float] = None
    gamma: Optional[float] = None
    z_factor: Optional[float] = None
    viscosity: Optional[float] = None
    density: Optional[float] = None
    cp: Optional[float] = None
    cv: Optional[float] = None
    source: Optional[str] = None
    log: List[str] = field(default_factory=list)


def _is_bad(x):
    """Check if a value is None, NaN, or otherwise invalid."""
    try:
        return x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x)))
    except Exception:
        return True


def resolve_liquid_properties(
    fluid_name: str,
    temperature_c: float,
    pressure_bar: float = 1.01325,
) -> Optional[ResolvedLiquidProperties]:
    """
    Resolve liquid fluid properties using CoolProp (primary) with fluidprop fallback.

    Args:
        fluid_name: Fluid name (will be mapped through aliases)
        temperature_c: Temperature in Celsius
        pressure_bar: Pressure in bar

    Returns:
        ResolvedLiquidProperties or None if all lookups fail
    """
    from .fluid_aliases import map_fluid_name

    result = ResolvedLiquidProperties()
    mapped_name = map_fluid_name(fluid_name)

    # Primary: CoolProp via CachedFluidProperties
    try:
        from .property_cache import CachedFluidProperties
        props = CachedFluidProperties(
            coolprop_name=mapped_name,
            T_in_deg_C=temperature_c,
            P_in_bar=pressure_bar,
        )
        result.density = float(props.rho[0])
        result.viscosity = float(props.eta[0])
        result.kinematic_viscosity = float(props.nu[0])
        result.source = "CoolProp"
        result.log.append(f"Properties via CoolProp for {mapped_name} at {temperature_c}°C, {pressure_bar} bar")

        # Vapor pressure and critical pressure via CoolProp
        try:
            from .property_cache import cached_saturation_pressure, cached_critical_pressure
            result.vapor_pressure = cached_saturation_pressure(mapped_name, temperature_c)
            result.log.append(f"Vapor pressure via CoolProp: {result.vapor_pressure:.2f} Pa")
        except Exception:
            pass
        try:
            from .property_cache import cached_critical_pressure
            result.critical_pressure = cached_critical_pressure(mapped_name)
            result.log.append(f"Critical pressure via CoolProp: {result.critical_pressure:.2f} Pa")
        except Exception:
            pass

        return result
    except Exception as e:
        result.log.append(f"CoolProp lookup failed for {mapped_name}: {e}")

    # Secondary: FluidProperties via fluidprop
    try:
        from .import_helpers import FLUIDPROP_AVAILABLE, FluidProperties, FLUID_SELECTION
        if FLUIDPROP_AVAILABLE and FluidProperties is not None:
            actual_name = mapped_name
            if FLUID_SELECTION is not None:
                try:
                    valid_fluids = [f[0] for f in FLUID_SELECTION if f is not None and hasattr(f, '__getitem__')]
                except (TypeError, IndexError):
                    valid_fluids = []
                if valid_fluids and actual_name not in valid_fluids:
                    match = next((f for f in valid_fluids if f.lower() == actual_name.lower()), None)
                    if match:
                        actual_name = match

            fp = FluidProperties(
                coolprop_name=actual_name,
                T_in_deg_C=temperature_c,
                P_in_bar=pressure_bar,
            )
            result.density = float(fp.rho[0])
            result.viscosity = float(fp.eta[0])
            result.kinematic_viscosity = float(fp.nu[0])
            result.source = "fluidprop"
            result.log.append(f"Properties via fluidprop for {actual_name}")
            return result
    except Exception as e:
        result.log.append(f"fluidprop lookup failed: {e}")

    return None


def resolve_gas_properties(
    fluid_name: str,
    temperature_c: float,
    pressure_bar: float = 1.01325,
) -> Optional[ResolvedGasProperties]:
    """
    Resolve gas fluid properties using CoolProp (primary) with thermo.Chemical fallback.

    Args:
        fluid_name: Fluid name (will be mapped through aliases)
        temperature_c: Temperature in Celsius
        pressure_bar: Pressure in bar

    Returns:
        ResolvedGasProperties or None if all lookups fail
    """
    from .fluid_aliases import map_fluid_name

    result = ResolvedGasProperties()
    mapped_name = map_fluid_name(fluid_name)
    T_K = temperature_c + 273.15

    # Primary: CoolProp via CachedFluidProperties
    try:
        from .property_cache import CachedFluidProperties
        props = CachedFluidProperties(
            coolprop_name=mapped_name,
            T_in_deg_C=temperature_c,
            P_in_bar=pressure_bar,
        )
        result.mw = float(props.MW)
        result.gamma = float(props.gamma[0])
        result.z_factor = float(props.Z[0])
        result.viscosity = float(props.eta[0])
        result.density = float(props.rho[0])
        result.cp = float(props.Cp[0])
        result.cv = float(props.Cv[0])
        result.source = "CoolProp"
        result.log.append(f"Gas properties via CoolProp for {mapped_name} at {temperature_c}°C, {pressure_bar} bar")
        return result
    except Exception as e:
        result.log.append(f"CoolProp gas lookup failed for {mapped_name}: {e}")

    # Secondary: FluidProperties via fluidprop
    try:
        from .import_helpers import FLUIDPROP_AVAILABLE, FluidProperties, FLUID_SELECTION
        if FLUIDPROP_AVAILABLE and FluidProperties is not None:
            actual_name = mapped_name
            if FLUID_SELECTION is not None:
                try:
                    valid_fluids = [f[0] for f in FLUID_SELECTION if f is not None and hasattr(f, '__getitem__')]
                except (TypeError, IndexError):
                    valid_fluids = []
                if valid_fluids and actual_name not in valid_fluids:
                    match = next((f for f in valid_fluids if f.lower() == actual_name.lower()), None)
                    if match:
                        actual_name = match

            fp = FluidProperties(
                coolprop_name=actual_name,
                T_in_deg_C=temperature_c,
                P_in_bar=pressure_bar,
            )
            result.mw = float(fp.MW)
            result.viscosity = float(fp.eta[0])
            result.density = float(fp.rho[0])
            if hasattr(fp, 'gamma'):
                result.gamma = float(fp.gamma[0])
            elif hasattr(fp, 'Cp') and hasattr(fp, 'Cv'):
                cv = float(fp.Cv[0])
                if cv > 0:
                    result.gamma = float(fp.Cp[0]) / cv
            if hasattr(fp, 'Z'):
                result.z_factor = float(fp.Z[0])
            if hasattr(fp, 'Cp'):
                result.cp = float(fp.Cp[0])
            if hasattr(fp, 'Cv'):
                result.cv = float(fp.Cv[0])
            result.source = "fluidprop"
            result.log.append(f"Gas properties via fluidprop for {actual_name}")

            # Check if we still have gaps — fall through to thermo if so
            if not any(_is_bad(v) for v in [result.mw, result.gamma, result.z_factor, result.viscosity]):
                return result
    except Exception as e:
        result.log.append(f"fluidprop gas lookup failed: {e}")

    # Tertiary: thermo.Chemical fallback for missing values
    try:
        from thermo.chemical import Chemical
        import re as _re

        cand = [mapped_name]
        try:
            spaced = _re.sub(r'(?<!^)(?=[A-Z])', ' ', mapped_name)
            cand.extend([spaced, spaced.lower(), mapped_name.lower()])
        except Exception:
            pass
        formula_map = {
            'sulfurdioxide': 'SO2', 'sulfur dioxide': 'SO2', 'sulphur dioxide': 'SO2',
            'carbondioxide': 'CO2', 'carbon dioxide': 'CO2',
            'hydrogensulfide': 'H2S', 'hydrogen sulfide': 'H2S',
        }
        cand.extend([formula_map.get(c, c) for c in list(cand)])

        chem = None
        for ident in cand:
            try:
                chem = Chemical(ident, T=T_K, P=pressure_bar * 1e5)
                if chem and chem.P is not None:
                    break
            except Exception:
                chem = None

        if chem is not None:
            if _is_bad(result.mw) and getattr(chem, 'MW', None):
                result.mw = float(chem.MW)
            if _is_bad(result.gamma):
                try:
                    cpg = float(chem.Cpg)
                    cvg = float(chem.Cvg)
                    if cvg > 0:
                        result.gamma = cpg / cvg
                except Exception:
                    pass
            if _is_bad(result.z_factor):
                try:
                    z_val = float(chem.Zg)
                    if z_val > 0:
                        result.z_factor = z_val
                except Exception:
                    pass
            if _is_bad(result.viscosity):
                try:
                    mu_val = float(chem.mug)
                    if mu_val > 0:
                        result.viscosity = mu_val
                except Exception:
                    pass
            if result.source is None:
                result.source = "thermo"
            else:
                result.source += "+thermo"
            result.log.append("Used thermo.Chemical fallback for missing properties")

            if not any(_is_bad(v) for v in [result.mw, result.gamma, result.z_factor, result.viscosity]):
                return result
    except Exception as e:
        result.log.append(f"thermo.Chemical fallback failed: {e}")

    # Return what we have (may have partial data)
    if result.source is not None:
        return result
    return None
