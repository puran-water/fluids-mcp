"""Reaction definitions and exact elemental/charge/redox validation."""
from __future__ import annotations
import hashlib
import json
from dataclasses import dataclass, field
from typing import Mapping
from .ontology import ELEMENTS, ComponentDefinition

@dataclass(frozen=True)
class ReactionDefinition:
    reaction_id: str
    stoichiometry: Mapping[str, float]
    parameter_ids: tuple[str, ...] = ()
    provenance: Mapping[str, str] = field(default_factory=dict)

    @property
    def sha256(self) -> str:
        canonical = json.dumps({'reaction_id': self.reaction_id, 'stoichiometry': dict(sorted(self.stoichiometry.items())), 'parameter_ids': self.parameter_ids, 'provenance': dict(sorted(self.provenance.items()))}, sort_keys=True, separators=(',', ':'), allow_nan=False)
        return hashlib.sha256(canonical.encode('utf-8')).hexdigest()

@dataclass(frozen=True)
class ConservationReport:
    reaction_id: str
    residuals: Mapping[str, float]
    tolerance: float

    @property
    def passed(self) -> bool:
        return all((abs(value) <= self.tolerance for value in self.residuals.values()))

    def assert_passed(self) -> None:
        if not self.passed:
            failed = {key: value for key, value in self.residuals.items() if abs(value) > self.tolerance}
            raise ValueError(f'reaction {self.reaction_id!r} is not conserved: {failed}')

def validate_reaction_conservation(reaction: ReactionDefinition, components: Mapping[str, ComponentDefinition], *, tolerance: float=1e-12) -> ConservationReport:
    unknown = sorted(set(reaction.stoichiometry) - set(components))
    if unknown:
        raise KeyError(f'reaction {reaction.reaction_id!r} references unknown components {unknown}')
    residuals = {element: 0.0 for element in ELEMENTS}
    residuals['charge'] = 0.0
    residuals['electron_equivalent'] = 0.0
    for component_id, coefficient in reaction.stoichiometry.items():
        component = components[component_id]
        for element in ELEMENTS:
            residuals[element] += coefficient * float(component.elemental_coefficients.get(element, 0.0))
        residuals['charge'] += coefficient * component.charge
        residuals['electron_equivalent'] += coefficient * component.electron_equivalent_mol_e
    return ConservationReport(reaction.reaction_id, residuals, tolerance)
__all__ = ['ConservationReport', 'ReactionDefinition', 'validate_reaction_conservation']
