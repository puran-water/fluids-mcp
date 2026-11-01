"""Versioned thermodynamic, phase, surface, kinetics, and reagent definitions."""
from __future__ import annotations
import hashlib
import json
import math
from dataclasses import asdict, dataclass, field
from typing import Literal, Mapping, Optional
from .expression import Expression
from .ontology import BasisMap, ComponentDefinition
from .reaction import ReactionDefinition, validate_reaction_conservation
PHREEQC_PHASE_GAS_CONSTANT_J_MOL_K = 8.3147

@dataclass(frozen=True)
class PhaseDefinition:
    phase_id: str
    formula: str
    dissolution_reaction: Mapping[str, float]
    log_k_25c: float
    source: str
    metastability_class: str
    temperature_correction_id: Optional[str] = None
    delta_h_j_mol: Optional[float] = None
    analytical_expression: Optional[tuple[float, ...]] = None

    def __post_init__(self) -> None:
        values = (self.log_k_25c, self.delta_h_j_mol)
        if any((value is not None and (not math.isfinite(float(value))) for value in values)):
            raise ValueError(f'phase {self.phase_id!r} constants must be finite')
        if self.delta_h_j_mol is not None and (not self.temperature_correction_id):
            raise ValueError(f'phase {self.phase_id!r} has an enthalpy but no temperature_correction_id')
        if self.analytical_expression is not None:
            coefficients = tuple((float(value) for value in self.analytical_expression))
            if not 1 <= len(coefficients) <= 6:
                raise ValueError(f'phase {self.phase_id!r} analytical_expression must carry 1-6 coefficients')
            if any((not math.isfinite(value) for value in coefficients)):
                raise ValueError(f'phase {self.phase_id!r} analytical_expression must be finite')
            if self.delta_h_j_mol is not None:
                raise ValueError(f'phase {self.phase_id!r} declares BOTH delta_h and an analytical expression — one temperature law per phase')
            if not self.temperature_correction_id:
                raise ValueError(f'phase {self.phase_id!r} has an analytical expression but no temperature_correction_id')
            anchor = self.log10_k(298.15)
            if abs(anchor - self.log_k_25c) > 0.005:
                raise ValueError(f'phase {self.phase_id!r} analytical expression evaluates to {anchor:.4f} at 25 C but log_k_25c declares {self.log_k_25c:.4f}')

    def log10_k(self, temperature_k: float) -> float:
        """log10 K at ``temperature_k`` under this phase's declared law."""
        temperature = float(temperature_k)
        if not math.isfinite(temperature) or temperature <= 0.0:
            raise ValueError('temperature_k must be finite and positive')
        if self.analytical_expression is not None:
            a = tuple((float(value) for value in self.analytical_expression)) + (0.0,) * 6
            return a[0] + a[1] * temperature + a[2] / temperature + a[3] * math.log10(temperature) + a[4] / temperature ** 2 + a[5] * temperature ** 2
        if self.delta_h_j_mol is None:
            return self.log_k_25c
        return self.log_k_25c + self.delta_h_j_mol / PHREEQC_PHASE_GAS_CONSTANT_J_MOL_K * (1.0 / 298.15 - 1.0 / temperature) / math.log(10.0)

    def equilibrium_constant(self, temperature_k: float) -> float:
        """Return K for the declared dissolution reaction at ``temperature_k``."""
        return 10.0 ** self.log10_k(temperature_k)

@dataclass(frozen=True)
class SurfaceDefinition:
    surface_id: str
    host_phase_id: str
    site_density_mol_per_mol_host: Mapping[str, float]
    protonation_constants: Mapping[str, float]
    binding_constants: Mapping[str, float]
    electrostatic_model: str
    destroy_sites_with_host: bool = True
    surface_dataset_id: str = 'unassigned'
    source: str = 'unassigned'
    parameter_interpretation: str = 'unassigned'
    oxide_history_class: str = 'unassigned'
    site_ratio_policy: str = 'independent'
    specific_area_m2_mol_host: Optional[float] = None
    host_inventory_component_id: Optional[str] = None

@dataclass(frozen=True)
class ConditionalLigandDefinition:
    """Conditional dissolved-organic complexation and site competition.

    This is deliberately a reduced, differentiable model.  ``Ligand_sites``
    is a transported analytical binding-site inventory; the conditional beta
    values partition total dissolved Fe without pretending that an unresolved
    wastewater ligand has a unique molecular formula.  A calibrated or
    bench-tested definition is still required before recommendation use.
    """
    ligand_model_id: str
    ligand_component_id: str
    total_basis: str
    beta_fe_ii_l_mol: float
    beta_fe_iii_l_mol: float
    surface_inhibition_constant_mol_l: float
    source: str
    parameter_interpretation: Literal['screening_prior', 'bench_calibrated_conditional_model', 'plant_calibrated_conditional_model'] = 'screening_prior'
    calibrated: bool = False
    minimum_ph: float = 0.0
    maximum_ph: float = 14.0
    maximum_ligand_sites_mol_l: float = 1.0

    def __post_init__(self) -> None:
        if not all((self.ligand_model_id, self.ligand_component_id, self.total_basis, self.source)):
            raise ValueError('conditional ligand identifiers, basis, and source are required')
        nonnegative = (self.beta_fe_ii_l_mol, self.beta_fe_iii_l_mol, self.surface_inhibition_constant_mol_l, self.maximum_ligand_sites_mol_l)
        if any((not math.isfinite(float(value)) or value < 0.0 for value in nonnegative)):
            raise ValueError('conditional ligand parameters must be finite and nonnegative')
        if self.surface_inhibition_constant_mol_l <= 0.0:
            raise ValueError('surface inhibition constant must be positive')
        if not (math.isfinite(float(self.minimum_ph)) and math.isfinite(float(self.maximum_ph)) and (0.0 <= self.minimum_ph < self.maximum_ph <= 14.0)):
            raise ValueError('conditional ligand pH validity interval is invalid')

@dataclass(frozen=True)
class ReagentDefinition:
    reagent_id: str
    formula: str
    component_stoichiometry: Mapping[str, float]
    molar_mass_kg_mol: Optional[float] = None
    assay_mass_fraction_mean: Optional[float] = None
    assay_mass_fraction_standard_uncertainty: Optional[float] = None
    density_kg_l_mean: Optional[float] = None
    density_kg_l_standard_uncertainty: Optional[float] = None
    free_acidity_eq_per_kg: Optional[float] = None
    free_acidity_status: str = 'not_measured'

    def __post_init__(self) -> None:
        if not self.reagent_id or not self.formula or (not self.component_stoichiometry):
            raise ValueError('reagent ID, formula, and component stoichiometry are required')
        if self.molar_mass_kg_mol is not None and (not math.isfinite(float(self.molar_mass_kg_mol)) or self.molar_mass_kg_mol <= 0):
            raise ValueError('reagent molar mass must be finite and positive')
        if any((not math.isfinite(float(coefficient)) or coefficient < 0 for coefficient in self.component_stoichiometry.values())):
            raise ValueError('reagent source coefficients must be finite and nonnegative')

@dataclass(frozen=True)
class KineticDefinition:
    kinetic_id: str
    expression: Expression
    parameter_provenance: Mapping[str, str] = field(default_factory=dict)
    reaction_id: Optional[str] = None
    rate_role: str = 'unassigned'

@dataclass(frozen=True)
class ChemistryKernelDefinition:
    schema_version: str
    kernel_id: str
    components: tuple[ComponentDefinition, ...]
    reactions: tuple[ReactionDefinition, ...]
    phases: tuple[PhaseDefinition, ...]
    surfaces: tuple[SurfaceDefinition, ...] = ()
    ligand_models: tuple[ConditionalLigandDefinition, ...] = ()
    kinetics: tuple[KineticDefinition, ...] = ()
    reagents: tuple[ReagentDefinition, ...] = ()
    basis_maps: tuple[BasisMap, ...] = ()
    thermodynamic_dataset_id: str = 'unassigned'
    aqueous_species_set_id: str = 'unassigned'
    phase_definition_set_id: str = 'unassigned'
    surface_model_id: str = 'none'
    phase_selection_profile_id: str = 'unassigned'
    kinetics_profile_id: str = 'unassigned'
    activity_model_id: str = 'ideal'
    ligand_model_id: str = 'none'
    redox_reaction_set_id: str = 'unassigned'
    commissioning_status: str = 'uncommissioned'
    recommendation_capable: bool = False

    @property
    def sha256(self) -> str:
        payload = asdict(self)
        keyed_collections = {'components': 'component_id', 'reactions': 'reaction_id', 'phases': 'phase_id', 'surfaces': 'surface_id', 'ligand_models': 'ligand_model_id', 'kinetics': 'kinetic_id', 'reagents': 'reagent_id', 'basis_maps': 'map_id'}
        for collection, key in keyed_collections.items():
            payload[collection] = sorted(payload.get(collection, []), key=lambda row: row[key])
        for row in payload.get('kinetics', []):
            kinetic = next((k for k in self.kinetics if k.kinetic_id == row['kinetic_id']))
            row['expression'] = kinetic.expression.to_dict()
        canonical = json.dumps(payload, sort_keys=True, separators=(',', ':'), allow_nan=False, default=str)
        return hashlib.sha256(canonical.encode('utf-8')).hexdigest()

class KernelDefinitionError(ValueError):
    """Raised when a shared chemistry kernel has ambiguous identifiers or links."""

def _require_unique(label: str, identifiers: list[str]) -> None:
    duplicates = sorted((identifier for identifier in set(identifiers) if identifiers.count(identifier) > 1))
    if duplicates:
        raise KernelDefinitionError(f'duplicate {label} identifiers: {duplicates}')

def validate_kernel_definition(kernel: ChemistryKernelDefinition) -> None:
    """Validate registry identity, referential integrity, and conservation."""
    collections = (('component', [item.component_id for item in kernel.components]), ('reaction', [item.reaction_id for item in kernel.reactions]), ('phase', [item.phase_id for item in kernel.phases]), ('surface', [item.surface_id for item in kernel.surfaces]), ('ligand-model', [item.ligand_model_id for item in kernel.ligand_models]), ('kinetic', [item.kinetic_id for item in kernel.kinetics]), ('reagent', [item.reagent_id for item in kernel.reagents]), ('basis-map', [item.map_id for item in kernel.basis_maps]))
    for label, identifiers in collections:
        _require_unique(label, identifiers)
    components = {component.component_id: component for component in kernel.components}
    phases = {phase.phase_id: phase for phase in kernel.phases}
    for reaction in kernel.reactions:
        validate_reaction_conservation(reaction, components).assert_passed()
    for phase in kernel.phases:
        if phase.phase_id not in components:
            raise KernelDefinitionError(f'phase {phase.phase_id!r} has no same-ID inventory component')
        unknown = sorted(set(phase.dissolution_reaction) - set(components))
        if unknown:
            raise KernelDefinitionError(f'phase {phase.phase_id!r} references unknown components {unknown}')
        validate_reaction_conservation(ReactionDefinition(reaction_id=f'{phase.phase_id}_dissolution', stoichiometry=phase.dissolution_reaction), components).assert_passed()
    for surface in kernel.surfaces:
        if surface.host_phase_id not in phases:
            raise KernelDefinitionError(f'surface {surface.surface_id!r} references unknown host phase {surface.host_phase_id!r}')
        inventory_id = surface.host_inventory_component_id or surface.host_phase_id
        if inventory_id not in components:
            raise KernelDefinitionError(f'surface {surface.surface_id!r} references unknown host inventory {inventory_id!r}')
    for ligand in kernel.ligand_models:
        if ligand.ligand_component_id not in components:
            raise KernelDefinitionError(f'ligand model {ligand.ligand_model_id!r} references unknown component {ligand.ligand_component_id!r}')
    reactions = {reaction.reaction_id for reaction in kernel.reactions}
    for kinetic in kernel.kinetics:
        if kinetic.reaction_id is not None and kinetic.reaction_id not in reactions:
            raise KernelDefinitionError(f'kinetic {kinetic.kinetic_id!r} references unknown reaction {kinetic.reaction_id!r}')
    for reagent in kernel.reagents:
        unknown = sorted(set(reagent.component_stoichiometry) - set(components))
        if unknown:
            raise KernelDefinitionError(f'reagent {reagent.reagent_id!r} references unknown components {unknown}')
__all__ = ['ChemistryKernelDefinition', 'ConditionalLigandDefinition', 'KernelDefinitionError', 'KineticDefinition', 'PHREEQC_PHASE_GAS_CONSTANT_J_MOL_K', 'PhaseDefinition', 'ReagentDefinition', 'SurfaceDefinition', 'validate_kernel_definition']
