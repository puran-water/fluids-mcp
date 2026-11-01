"""Canonical ontology with explicit native measurement bases."""
from __future__ import annotations
import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Mapping, Optional

class BasisKind(str, Enum):
    MOL_SPECIES = 'mol_species'
    KG_COD = 'kg_COD'
    KG_ELEMENT = 'kg_element'
    EQUIVALENT_CHARGE = 'equivalent_charge'
    KG_TSS = 'kg_TSS'
    DIMENSIONLESS_INVENTORY = 'dimensionless_inventory'

class Phase(str, Enum):
    AQUEOUS = 'aqueous'
    GAS = 'gas'
    SOLID = 'solid'
    BIOMASS = 'biomass'
    SURFACE = 'surface'

class ParticleClass(str, Enum):
    TRULY_DISSOLVED = 'truly_dissolved'
    COLLOIDAL = 'colloidal'
    PARTICULATE = 'particulate'
    GAS = 'gas'

class PermeabilityClass(str, Enum):
    PERMEATING = 'permeating'
    PARTIALLY_REJECTED = 'partially_rejected'
    RETAINED = 'retained'
    NOT_APPLICABLE = 'not_applicable'
ELEMENTS = ('H', 'O', 'C', 'N', 'P', 'S', 'Fe', 'Al', 'Ca', 'Mg', 'Na', 'K', 'Cl')

@dataclass(frozen=True)
class ComponentDefinition:
    component_id: str
    basis_kind: BasisKind
    measured_as: str
    native_mass_basis: str
    elemental_coefficients: Mapping[str, float] = field(default_factory=dict)
    charge: float = 0.0
    electron_equivalent_mol_e: float = 0.0
    phase: Phase = Phase.AQUEOUS
    particle_class: ParticleClass = ParticleClass.TRULY_DISSOLVED
    permeability_class: PermeabilityClass = PermeabilityClass.PERMEATING
    molecular_or_empirical_formula: Optional[str] = None

    def __post_init__(self) -> None:
        if not self.component_id or not self.component_id.strip():
            raise ValueError('component_id must be non-empty')
        unknown = sorted(set(self.elemental_coefficients) - set(ELEMENTS))
        if unknown:
            raise ValueError(f'{self.component_id}: unsupported elements {unknown}')
        values = [*self.elemental_coefficients.values(), self.charge, self.electron_equivalent_mol_e]
        if not all((math.isfinite(float(value)) for value in values)):
            raise ValueError(f'{self.component_id}: coefficients must be finite')

    @property
    def cod_kg_o2(self) -> float:
        """COD observable derived from the one redox-equivalent convention.

        One mole O2 accepts four mole electrons, hence one mol e- corresponds
        to 8 g O2 (0.008 kg O2).
        """
        return self.electron_equivalent_mol_e * 0.008

@dataclass(frozen=True)
class BasisMapEntry:
    native_component_id: str
    chemical_total_id: str
    factor: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.factor):
            raise ValueError('basis-map factor must be finite')

@dataclass(frozen=True)
class BasisMap:
    map_id: str
    entries: tuple[BasisMapEntry, ...]

    def to_chemical_totals(self, native_state: Mapping[str, float]) -> dict[str, float]:
        totals: dict[str, float] = {}
        for entry in self.entries:
            value = float(native_state.get(entry.native_component_id, 0.0))
            totals[entry.chemical_total_id] = totals.get(entry.chemical_total_id, 0.0) + entry.factor * value
        return totals

    def to_native_extents(self, chemical_extents: Mapping[str, float]) -> dict[str, float]:
        """Map chemical-total extents back when the map is one-to-one.

        Many-to-one biological aggregation is intentionally non-invertible and
        must be handled by an explicit engine adapter.
        """
        by_total: dict[str, list[BasisMapEntry]] = {}
        for entry in self.entries:
            by_total.setdefault(entry.chemical_total_id, []).append(entry)
        ambiguous = [key for key, rows in by_total.items() if len(rows) != 1]
        if ambiguous:
            raise ValueError(f'basis map {self.map_id!r} is not invertible for totals {ambiguous}')
        return {rows[0].native_component_id: float(chemical_extents.get(total, 0.0)) / rows[0].factor for total, rows in by_total.items() if rows[0].factor != 0}
__all__ = ['BasisKind', 'BasisMap', 'BasisMapEntry', 'ComponentDefinition', 'ELEMENTS', 'ParticleClass', 'PermeabilityClass', 'Phase']
