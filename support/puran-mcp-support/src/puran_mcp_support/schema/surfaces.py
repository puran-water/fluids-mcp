"""Immutable surface-capacity datasets and basis-safe conversion helpers.

Surface capacities are deliberately kept separate from mineral equilibrium
constants.  In particular, the Mao--Yue in-situ HFO values are effective
parameters inferred with a kinetic model; they are not first-principles
Dzombak--Morel diffuse-layer site densities.
"""
from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Literal, Optional
SurfaceParameterInterpretation = Literal['thermodynamic_surface_dataset', 'kinetic_calibration_surrogate']

@dataclass(frozen=True)
class SurfaceCapacityDataset:
    """Versioned site-capacity definition on a per-mole-host basis."""
    dataset_id: str
    host_phase_id: str
    weak_sites_mol_per_mol_host: float
    strong_sites_mol_per_mol_host: float
    host_molar_mass_g_mol: float
    specific_area_m2_mol_host: float
    source: str
    parameter_interpretation: SurfaceParameterInterpretation
    oxide_history_class: str
    formation_ph: Optional[float] = None
    formation_ph_basis: Optional[str] = None

    def __post_init__(self) -> None:
        positive = (self.weak_sites_mol_per_mol_host, self.strong_sites_mol_per_mol_host, self.host_molar_mass_g_mol, self.specific_area_m2_mol_host)
        if any((not math.isfinite(value) or value <= 0 for value in positive)):
            raise ValueError('surface site, mass, and area values must be finite and positive')
        if self.formation_ph is not None and (not math.isfinite(self.formation_ph)):
            raise ValueError('formation pH must be finite')

    @property
    def total_sites_mol_per_mol_host(self) -> float:
        return self.weak_sites_mol_per_mol_host + self.strong_sites_mol_per_mol_host

    @property
    def weak_to_strong_ratio(self) -> float:
        return self.weak_sites_mol_per_mol_host / self.strong_sites_mol_per_mol_host

    def site_moles_for_host_mass(self, host_mass_g: float) -> dict[str, float]:
        """Convert a host-solid mass into weak/strong site inventories."""
        if not math.isfinite(float(host_mass_g)) or host_mass_g < 0:
            raise ValueError('host mass must be finite and nonnegative')
        host_moles = float(host_mass_g) / self.host_molar_mass_g_mol
        return {'weak': host_moles * self.weak_sites_mol_per_mol_host, 'strong': host_moles * self.strong_sites_mol_per_mol_host}

    def site_density_mol_m2(self, specific_area_m2_g: float) -> dict[str, float]:
        """Return equivalent site densities for a caller-specified solid area."""
        if not math.isfinite(float(specific_area_m2_g)) or specific_area_m2_g <= 0:
            raise ValueError('specific area must be finite and positive')
        area_per_mol = float(specific_area_m2_g) * self.host_molar_mass_g_mol
        return {'weak': self.weak_sites_mol_per_mol_host / area_per_mol, 'strong': self.strong_sites_mol_per_mol_host / area_per_mol}

def split_total_sites_preserving_ratio(total_sites_mol_per_mol_host: float, *, weak_to_strong_ratio: float=40.0) -> tuple[float, float]:
    """Return ``(weak, strong)`` while preserving an explicit site ratio."""
    total = float(total_sites_mol_per_mol_host)
    ratio = float(weak_to_strong_ratio)
    if not math.isfinite(total) or total <= 0:
        raise ValueError('total site inventory must be finite and positive')
    if not math.isfinite(ratio) or ratio <= 0:
        raise ValueError('weak-to-strong site ratio must be finite and positive')
    strong = total / (ratio + 1.0)
    return (ratio * strong, strong)
DZOMBAK_MOREL_HFO_DATASET_ID = 'hfo_dzombak_morel_preformed_v1'
DZOMBAK_MOREL_HFO = SurfaceCapacityDataset(dataset_id=DZOMBAK_MOREL_HFO_DATASET_ID, host_phase_id='HFO', weak_sites_mol_per_mol_host=0.2, strong_sites_mol_per_mol_host=0.005, host_molar_mass_g_mol=89.0, specific_area_m2_mol_host=53300.0, source='Dzombak and Morel HFO parameterization as documented by USGS PHREEQC', parameter_interpretation='thermodynamic_surface_dataset', oxide_history_class='preformed_aged')
_MAO_YUE_EFFECTIVE_SITES_BY_PH = ((6.0, 0.955), (7.0, 1.46), (8.0, 2.98))

def mao_yue_effective_total_sites(formation_ph: float) -> float:
    """Interpolate the reported effective in-situ HFO capacity over pH 6--8.

    Values outside the experimental interval are clamped.  The return value is
    intentionally named *effective total sites*: using it in a diffuse-layer
    equilibrium model is a declared calibration surrogate.
    """
    ph = float(formation_ph)
    if not math.isfinite(ph):
        raise ValueError('formation pH must be finite')
    points = _MAO_YUE_EFFECTIVE_SITES_BY_PH
    if ph <= points[0][0]:
        return points[0][1]
    if ph >= points[-1][0]:
        return points[-1][1]
    for (ph0, sites0), (ph1, sites1) in zip(points, points[1:]):
        if ph0 <= ph <= ph1:
            return sites0 + (sites1 - sites0) * (ph - ph0) / (ph1 - ph0)
    raise AssertionError('unreachable interpolation interval')

def mao_yue_fresh_hfo_dataset(formation_ph: float, *, formation_ph_basis: str) -> SurfaceCapacityDataset:
    """Build a fully identified fresh-HFO surrogate at a declared local pH."""
    if not formation_ph_basis:
        raise ValueError('formation_ph_basis is required for fresh in-situ HFO')
    total = mao_yue_effective_total_sites(formation_ph)
    weak, strong = split_total_sites_preserving_ratio(total)
    return SurfaceCapacityDataset(dataset_id='hfo_mao_yue_in_situ_effective_v1', host_phase_id='HFO', weak_sites_mol_per_mol_host=weak, strong_sites_mol_per_mol_host=strong, host_molar_mass_g_mol=89.0, specific_area_m2_mol_host=53300.0, source='Mao and Yue (2016), DOI 10.1038/srep35292', parameter_interpretation='kinetic_calibration_surrogate', oxide_history_class='fresh_in_situ', formation_ph=float(formation_ph), formation_ph_basis=formation_ph_basis)
__all__ = ['DZOMBAK_MOREL_HFO', 'DZOMBAK_MOREL_HFO_DATASET_ID', 'SurfaceCapacityDataset', 'SurfaceParameterInterpretation', 'mao_yue_effective_total_sites', 'mao_yue_fresh_hfo_dataset', 'split_total_sites_preserving_ratio']
