from __future__ import annotations
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional
VALID_TEMPERATURE_RANGE_C = (0.0, 150.0)
_K0_LUBE = 0.0
_K1_LUBE = 0.34878

class InsufficientOilData(ValueError):
    """Raised when oil properties cannot be established from what was supplied."""

    def __init__(self, missing: List[str], context: str=''):
        self.missing = list(missing)
        msg = 'Lubricant properties cannot be established; missing: ' + ', '.join(self.missing)
        if context:
            msg = f'{msg}. {context}'
        super().__init__(msg)

def _walther_z(nu_cSt: float) -> float:
    return math.log10(math.log10(nu_cSt + 0.7))

@dataclass
class LubeOil:
    """A specific lubricant, defined by its datasheet anchors.

    Parameters
    ----------
    nu_40_cSt, nu_100_cSt:
        Kinematic viscosity at 40 °C and 100 °C. Both are required — see the module
        docstring for why a grade name alone will not do.
    density_15_kg_m3:
        Density at 15 °C.
    designation, provenance:
        Free text carried through to the rating output so a reviewer can see which oil and
        which document the numbers came from.
    """
    nu_40_cSt: float
    nu_100_cSt: float
    density_15_kg_m3: float
    designation: str = ''
    provenance: str = ''
    warnings: List[str] = field(default_factory=list)
    _A: float = field(init=False, repr=False, default=0.0)
    _B: float = field(init=False, repr=False, default=0.0)

    def __post_init__(self) -> None:
        missing = [name for name, val in (('nu_40_cSt', self.nu_40_cSt), ('nu_100_cSt', self.nu_100_cSt), ('density_15_kg_m3', self.density_15_kg_m3)) if not val or val <= 0]
        if missing:
            raise InsufficientOilData(missing, 'All three are on the lubricant datasheet.')
        if self.nu_100_cSt >= self.nu_40_cSt:
            raise InsufficientOilData(['nu_100_cSt < nu_40_cSt'], 'Viscosity must fall with temperature; check the two anchors are not swapped.')
        t40, t100 = (40.0 + 273.15, 100.0 + 273.15)
        z40, z100 = (_walther_z(self.nu_40_cSt), _walther_z(self.nu_100_cSt))
        self._B = (z40 - z100) / (math.log10(t100) - math.log10(t40))
        self._A = z40 + self._B * math.log10(t40)
        if not self.provenance:
            self.warnings.append('No provenance recorded for these oil properties; cite the lubricant datasheet so the basis is auditable.')

    def _range_check(self, T_C: float) -> None:
        lo, hi = VALID_TEMPERATURE_RANGE_C
        if not lo <= T_C <= hi:
            self.warnings.append(f'{T_C:.1f} °C is outside the {lo:.0f}-{hi:.0f} °C range over which these petroleum correlations are established; treat the result as an extrapolation.')

    @property
    def specific_gravity_15(self) -> float:
        return self.density_15_kg_m3 / 999.016

    @property
    def viscosity_index(self) -> Optional[float]:
        """ASTM D2270 viscosity index, when ``chemicals`` is available."""
        try:
            from chemicals import viscosity_index
            return float(viscosity_index(nu_40=self.nu_40_cSt * 1e-06, nu_100=self.nu_100_cSt * 1e-06))
        except Exception:
            return None

    def kinematic_viscosity_cSt(self, T_C: float) -> float:
        """Kinematic viscosity [cSt] by ASTM D341, exact at 40 and 100 °C."""
        self._range_check(T_C)
        log_nu = 10.0 ** (self._A - self._B * math.log10(T_C + 273.15))
        return 10.0 ** log_nu - 0.7

    def density_kg_m3(self, T_C: float) -> float:
        """Density [kg/m³] by the ASTM D1250 lubricating-oil expansion form."""
        self._range_check(T_C)
        rho15 = self.density_15_kg_m3
        alpha = _K0_LUBE / rho15 ** 2 + _K1_LUBE / rho15
        dT = T_C - 15.0
        return rho15 * math.exp(-alpha * dT * (1.0 + 0.8 * alpha * dT))

    def dynamic_viscosity_Pa_s(self, T_C: float) -> float:
        """Dynamic viscosity [Pa·s] = ν · ρ."""
        return self.kinematic_viscosity_cSt(T_C) * 1e-06 * self.density_kg_m3(T_C)

    def specific_heat_J_kgK(self, T_C: float) -> float:
        """Specific heat [J/kg·K] by Cragoe (1929)."""
        self._range_check(T_C)
        T_F = T_C * 9.0 / 5.0 + 32.0
        cp_btu = (0.388 + 0.00045 * T_F) / math.sqrt(self.specific_gravity_15)
        return cp_btu * 4186.8

    def thermal_conductivity_W_mK(self, T_C: float) -> float:
        """Thermal conductivity [W/m·K] by the Cragoe petroleum-fraction form."""
        self._range_check(T_C)
        return 0.117 / self.specific_gravity_15 * (1.0 - 0.00054 * T_C)

    def prandtl(self, T_C: float) -> float:
        return self.specific_heat_J_kgK(T_C) * self.dynamic_viscosity_Pa_s(T_C) / self.thermal_conductivity_W_mK(T_C)

    def properties_at(self, T_C: float) -> Dict[str, object]:
        """Full property set at ``T_C``, with basis and provenance attached."""
        before = len(self.warnings)
        out = {'temperature_C': T_C, 'kinematic_viscosity_cSt': self.kinematic_viscosity_cSt(T_C), 'dynamic_viscosity_Pa_s': self.dynamic_viscosity_Pa_s(T_C), 'density_kg_m3': self.density_kg_m3(T_C), 'specific_heat_J_kgK': self.specific_heat_J_kgK(T_C), 'thermal_conductivity_W_mK': self.thermal_conductivity_W_mK(T_C), 'prandtl': self.prandtl(T_C), 'viscosity_index': self.viscosity_index, 'designation': self.designation, 'provenance': self.provenance, 'correlations': {'viscosity': 'ASTM D341 (Walther), fitted to nu40 + nu100', 'density': 'ASTM D1250 lubricating-oil expansion', 'specific_heat': 'Cragoe (1929)', 'thermal_conductivity': 'Cragoe (1929) petroleum fraction'}}
        new = self.warnings[before:]
        del self.warnings[before:]
        for w in new:
            if w not in self.warnings:
                self.warnings.append(w)
        out['warnings'] = list(self.warnings)
        return out
