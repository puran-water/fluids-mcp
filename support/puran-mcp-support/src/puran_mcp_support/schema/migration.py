"""Versioned native-state migration contracts for shared chemistry kernels."""
from __future__ import annotations
from dataclasses import asdict, dataclass
from typing import Literal
NativeStateRole = Literal['differential_state', 'algebraic_state', 'implicit_reservoir']
NativeStateImplementationStatus = Literal['governing', 'partial_native', 'sidecar_only', 'missing', 'implicit_reservoir']

@dataclass(frozen=True)
class NativeStateRequirement:
    """One canonical inventory or closure required in an engine-native model."""
    requirement_id: str
    canonical_state_ids: tuple[str, ...]
    required_role: NativeStateRole
    implementation_status: NativeStateImplementationStatus
    native_state_ids: tuple[str, ...] = ()
    basis_or_semantics: str = ''
    blocker: str | None = None

    def __post_init__(self) -> None:
        if not self.requirement_id.strip() or not self.canonical_state_ids:
            raise ValueError('native-state requirement ID and canonical states are required')
        if len(set(self.canonical_state_ids)) != len(self.canonical_state_ids):
            raise ValueError('canonical state IDs must be unique within a requirement')
        if len(set(self.native_state_ids)) != len(self.native_state_ids):
            raise ValueError('native state IDs must be unique within a requirement')
        if self.implementation_status in {'governing', 'partial_native'} and (not self.native_state_ids):
            raise ValueError('governing or partial native requirements must identify native states')
        if self.implementation_status == 'implicit_reservoir' and self.required_role != 'implicit_reservoir':
            raise ValueError('implicit-reservoir status is valid only for an implicit-reservoir role')
        if not self.satisfied and (not self.blocker):
            raise ValueError('an unsatisfied native-state requirement needs a blocker')

    @property
    def satisfied(self) -> bool:
        return self.implementation_status in {'governing', 'implicit_reservoir'}

    def public_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload['satisfied'] = self.satisfied
        return payload

@dataclass(frozen=True)
class NativeStateMigrationContract:
    """Fail-closed compatibility record for a checkpoint-breaking integration."""
    migration_id: str
    engine_id: str
    model_variant: str
    kernel_id: str
    kernel_hash: str
    source_layout_version: str
    target_layout_version: str
    checkpoint_policy: Literal['cold_reseed', 'explicit_versioned_migration']
    requirements: tuple[NativeStateRequirement, ...]
    shared_material_source_terms_coupled: bool
    shared_charge_acid_base_coupled: bool
    serial_zone_topology_coupled: bool
    initialization_and_scaling_complete: bool
    restart_migration_implemented: bool

    def __post_init__(self) -> None:
        required_text = (self.migration_id, self.engine_id, self.model_variant, self.kernel_id, self.source_layout_version, self.target_layout_version)
        if any((not value.strip() for value in required_text)):
            raise ValueError('native-state migration identifiers are required')
        if len(self.kernel_hash) != 64 or any((character not in '0123456789abcdef' for character in self.kernel_hash.lower())):
            raise ValueError('kernel_hash must be a 64-character hexadecimal SHA-256')
        if not self.requirements:
            raise ValueError('native-state migration requires at least one state requirement')
        requirement_ids = [row.requirement_id for row in self.requirements]
        if len(set(requirement_ids)) != len(requirement_ids):
            raise ValueError('native-state migration requirement IDs must be unique')

    @property
    def unsatisfied_requirements(self) -> tuple[str, ...]:
        return tuple((row.requirement_id for row in self.requirements if not row.satisfied))

    @property
    def native_coupling_complete(self) -> bool:
        return bool(not self.unsatisfied_requirements and self.shared_material_source_terms_coupled and self.shared_charge_acid_base_coupled and self.serial_zone_topology_coupled and self.initialization_and_scaling_complete and self.restart_migration_implemented)

    @property
    def blockers(self) -> tuple[str, ...]:
        blockers = [row.blocker for row in self.requirements if not row.satisfied and row.blocker is not None]
        flag_blockers = ((self.shared_material_source_terms_coupled, 'shared_chemistry_source_terms_not_governing_native_material_balances'), (self.shared_charge_acid_base_coupled, 'shared_chemistry_not_coupled_to_native_charge_and_acid_base_closure'), (self.serial_zone_topology_coupled, 'serial_chemistry_zones_not_governing_native_topology'), (self.initialization_and_scaling_complete, 'native_chemistry_initialization_and_scaling_incomplete'), (self.restart_migration_implemented, 'native_state_restart_migration_not_implemented'))
        blockers.extend((blocker for passed, blocker in flag_blockers if not passed))
        return tuple(dict.fromkeys(blockers))

    def public_dict(self) -> dict[str, object]:
        return {'schema_version': 'native_state_migration_contract/v1', 'migration_id': self.migration_id, 'engine_id': self.engine_id, 'model_variant': self.model_variant, 'kernel_id': self.kernel_id, 'kernel_hash': self.kernel_hash, 'source_layout_version': self.source_layout_version, 'target_layout_version': self.target_layout_version, 'checkpoint_policy': self.checkpoint_policy, 'requirements': [row.public_dict() for row in self.requirements], 'unsatisfied_requirements': list(self.unsatisfied_requirements), 'shared_material_source_terms_coupled': self.shared_material_source_terms_coupled, 'shared_charge_acid_base_coupled': self.shared_charge_acid_base_coupled, 'serial_zone_topology_coupled': self.serial_zone_topology_coupled, 'initialization_and_scaling_complete': self.initialization_and_scaling_complete, 'restart_migration_implemented': self.restart_migration_implemented, 'native_coupling_complete': self.native_coupling_complete, 'recommendation_eligible': False, 'blockers': list(self.blockers)}
__all__ = ['NativeStateImplementationStatus', 'NativeStateMigrationContract', 'NativeStateRequirement', 'NativeStateRole']
