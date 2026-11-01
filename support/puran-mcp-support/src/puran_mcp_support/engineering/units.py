from __future__ import annotations
from dataclasses import dataclass
from typing import Final
FOOT_M: Final[float] = 0.3048
POUND_KG: Final[float] = 0.45359237
PSI_PA: Final[float] = 6894.757293168361
ATM_PA: Final[float] = 101325.0
R_UNIVERSAL: Final[float] = 8.314462618
M_DRY_AIR: Final[float] = 0.0289647
M_WATER: Final[float] = 0.01801528

def moist_air_density(temperature_K: float, pressure_Pa: float, relative_humidity: float=0.0) -> float:
    """Density of moist air (kg/m3) from the ideal gas law.

    ``relative_humidity`` is a fraction in [0, 1]. Saturation vapour pressure
    uses the Magnus form, which is well inside measurement error over the 0-50 C
    range this repo operates in. Moist air is *less* dense than dry air at the
    same temperature and pressure, because water (18 g/mol) is lighter than the
    air it displaces (29 g/mol).
    """
    T = float(temperature_K)
    P = float(pressure_Pa)
    rh = min(max(float(relative_humidity), 0.0), 1.0)
    if rh > 0.0:
        t_c = T - 273.15
        p_sat = 610.94 * pow(2.718281828459045, 17.625 * t_c / (t_c + 243.04))
        p_v = rh * p_sat
    else:
        p_v = 0.0
    p_d = P - p_v
    return (p_d * M_DRY_AIR + p_v * M_WATER) / (R_UNIVERSAL * T)

@dataclass(frozen=True)
class GasVolumeBasis:
    """A named reference condition for a *standard* gas volume.

    ``density_kg_m3`` is carried explicitly rather than always derived, because
    vendor and standards bodies quote rounded conventional densities (notably
    0.0750 lb/ft3 for US blower practice) and a design must reproduce the
    vendor's arithmetic exactly, not a marginally different ideal-gas value.
    Use :meth:`derived_density` to see what the ideal gas law would give.
    """
    name: str
    temperature_K: float
    pressure_Pa: float
    relative_humidity: float
    density_kg_m3: float
    note: str = ''

    def derived_density(self) -> float:
        """Ideal-gas density at this basis, for cross-checking the quoted value."""
        return moist_air_density(self.temperature_K, self.pressure_Pa, self.relative_humidity)
SCFM_BLOWER: Final[GasVolumeBasis] = GasVolumeBasis(name='scfm_blower_ASME_CAGI', temperature_K=293.15, pressure_Pa=14.696 * PSI_PA, relative_humidity=0.36, density_kg_m3=0.075 * POUND_KG / FOOT_M ** 3, note='68 F, 14.696 psia, 36% RH. US wastewater blower and ASCE convention.')
SCF_GAS_INDUSTRY: Final[GasVolumeBasis] = GasVolumeBasis(name='scf_gas_industry', temperature_K=288.706, pressure_Pa=14.696 * PSI_PA, relative_humidity=0.0, density_kg_m3=moist_air_density(288.706, 14.696 * PSI_PA, 0.0), note='60 F, 14.696 psia, dry. NOT the aeration convention -- ~1.6% denser.')
NM3_ISO: Final[GasVolumeBasis] = GasVolumeBasis(name='Nm3_ISO_13443', temperature_K=273.15, pressure_Pa=ATM_PA, relative_humidity=0.0, density_kg_m3=moist_air_density(273.15, ATM_PA, 0.0), note='0 C, 101.325 kPa, dry.')
NM3_20C: Final[GasVolumeBasis] = GasVolumeBasis(name='Nm3_20C', temperature_K=293.15, pressure_Pa=ATM_PA, relative_humidity=0.0, density_kg_m3=moist_air_density(293.15, ATM_PA, 0.0), note='20 C, 101.325 kPa, dry. Common in European aeration datasheets.')
GAS_VOLUME_BASES: Final[dict[str, GasVolumeBasis]] = {b.name: b for b in (SCFM_BLOWER, SCF_GAS_INDUSTRY, NM3_ISO, NM3_20C)}

def scfm_to_std_m3_per_day(scfm: float, basis: GasVolumeBasis) -> float:
    """Standard cubic feet per minute -> standard cubic metres per day.

    Purely geometric: ft3 -> m3 and per-minute -> per-day. The reference
    condition does not enter, because both sides are *at the same basis*. The
    ``basis`` argument is required anyway so the caller has to state which
    standard the number came from -- it is carried into the mass conversions
    below, where it does matter.
    """
    return float(scfm) * FOOT_M ** 3 * 60.0 * 24.0

def std_to_actual_m3_per_day(std_m3_per_day: float, basis: GasVolumeBasis, field_temperature_K: float, field_pressure_Pa: float) -> float:
    """Standard volumetric flow -> actual (field-condition) volumetric flow.

    Required whenever a model wants volume *at the blower discharge or the
    diffuser*, rather than at the reference condition. QSDsan's
    ``DiffusedAeration.Q_air`` is documented as a field-condition flow and
    internally divides by an air-density factor to recover standard conditions,
    so handing it a standard flow under-reports SOTR.
    """
    return float(std_m3_per_day) * (float(field_temperature_K) / basis.temperature_K) * (basis.pressure_Pa / float(field_pressure_Pa))

def site_pressure_pa(elevation_m: float, sea_level_pressure_Pa: float=ATM_PA) -> float:
    """Barometric pressure at elevation, ISA troposphere model.

    Used for the omega (pressure) correction in oxygen-transfer calculations.
    """
    z = float(elevation_m)
    return float(sea_level_pressure_Pa) * pow(1.0 - 2.25577e-05 * z, 5.25588)
