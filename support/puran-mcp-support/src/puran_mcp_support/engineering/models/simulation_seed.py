from __future__ import annotations
import math
from typing import Any, Literal
from uuid import UUID
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator
from puran_mcp_support.engineering.contracts.proposal_inputs import CanonicalCredibilityStamp
EngineKey = Literal['ASM2d', 'ADM1', 'WaterTAP-ASM1', 'WaterTAP-ASM2D', 'WaterTAP-ModifiedASM2D', 'WaterTAP-ADM1', 'WaterTAP-ModifiedADM1', 'qsdsan-ASM2D', 'qsdsan-ModifiedASM2D', 'qsdsan-ADM1', 'qsdsan-ModifiedADM1', 'WaterTAP-RO', 'WaterTAP-NF', 'WaterTAP-evaporator', 'WaterTAP-IX', 'WaterTAP-membrane.gas_separation', 'phys_chem.coag_floc', 'ro.conventional', 'mvr.evaporator', 'atfd.thin_film_dryer', 'physico-chem', 'sludge', 'odor-control', 'none']
PortRole = Literal['feed', 'effluent', 'forward', 'permeate', 'concentrate', 'retentate', 'distillate', 'brine', 'filtrate', 'sludge', 'biogas', 'regen_waste', 'makeup', 'return', 'recirculation', 'vent']
_QSDSAN_ENGINE_KEYS = frozenset({'ASM2d', 'ADM1', 'physico-chem', 'sludge'})
_ASM2D_INFLUENT_REQUIRED_ENGINE_KEYS = frozenset({'ASM2d', 'WaterTAP-ASM2D', 'WaterTAP-ModifiedASM2D', 'qsdsan-ASM2D', 'qsdsan-ModifiedASM2D'})
_ADM1_INFLUENT_REQUIRED_ENGINE_KEYS = frozenset({'ADM1', 'WaterTAP-ADM1', 'WaterTAP-ModifiedADM1', 'qsdsan-ADM1', 'qsdsan-ModifiedADM1'})

def _require_non_empty(value: str, *, field_name: str) -> str:
    """Return a stripped string and reject empty values."""
    text = value.strip()
    if not text:
        raise ValueError(f'`{field_name}` must be non-empty.')
    return text

class SimulationSeedConnection(BaseModel):
    """Downstream reference for an outlet."""
    model_config = ConfigDict(extra='forbid', str_strip_whitespace=True)
    downstream_seed_id: str
    downstream_inlet_role: PortRole
    tear_stream: bool = False
    tear_rationale: str | None = None

    @field_validator('downstream_seed_id')
    @classmethod
    def _validate_downstream_seed_id(cls, value: str) -> str:
        return _require_non_empty(value, field_name='downstream_seed_id')

    @field_validator('tear_rationale')
    @classmethod
    def _validate_tear_rationale(cls, value: str | None) -> str | None:
        if value is None:
            return None
        return _require_non_empty(value, field_name='tear_rationale')

    @model_validator(mode='after')
    def _require_tear_rationale(self) -> 'SimulationSeedConnection':
        if self.tear_stream and self.tear_rationale is None:
            raise ValueError('tear_stream=True requires `tear_rationale`.')
        return self

class SimulationSeedPort(BaseModel):
    """A single inlet or outlet on a seed, with canonical stream_id.

    `stream_id` must be the same canonical id that the downstream engine's
    seed references on its matching inlet and that the Stage 7b refined
    stream_snapshot writer persists.
    """
    model_config = ConfigDict(extra='forbid', str_strip_whitespace=True)
    role: PortRole
    stream_id: str
    connects_to: SimulationSeedConnection | None = None

    @field_validator('stream_id')
    @classmethod
    def _validate_stream_id(cls, value: str) -> str:
        return _require_non_empty(value, field_name='stream_id')

class SimulationSeedTopology(BaseModel):
    """v2 topology block with typed inlet/outlet ports."""
    model_config = ConfigDict(extra='forbid', str_strip_whitespace=True)
    inlets: list[SimulationSeedPort]
    outlets: list[SimulationSeedPort]

    @model_validator(mode='after')
    def _require_unique_ports(self) -> 'SimulationSeedTopology':
        all_ports = self.inlets + self.outlets
        ids = [p.stream_id for p in all_ports]
        if len(ids) != len(set(ids)):
            raise ValueError(f'SimulationSeedTopology: duplicate stream_id across inlets+outlets: {sorted(ids)}')
        return self

class SimulationSeedInfluentComposition(BaseModel):
    """Fractionated plant-boundary influent carried by biological seeds."""
    model_config = ConfigDict(extra='forbid', str_strip_whitespace=True)
    species_kg_m3: dict[str, float]
    flow_m3_d: float = Field(gt=0)
    temperature_K: float = Field(gt=0)
    alkalinity_mol_per_m3: float = Field(ge=0)
    fractionation_basis: Literal['asm2d', 'adm1']
    source_influent_artifact_id: UUID

    @field_validator('species_kg_m3')
    @classmethod
    def _validate_species(cls, value: dict[str, float]) -> dict[str, float]:
        if not value:
            raise ValueError('influent_composition.species_kg_m3 must be non-empty.')
        normalized: dict[str, float] = {}
        for key, raw in value.items():
            component = _require_non_empty(str(key), field_name='species_kg_m3 component')
            try:
                numeric = float(raw)
            except (TypeError, ValueError) as exc:
                raise ValueError(f'influent_composition.species_kg_m3[{component!r}] must be numeric.') from exc
            if not math.isfinite(numeric) or numeric < 0:
                raise ValueError(f'influent_composition.species_kg_m3[{component!r}] must be finite and >= 0.')
            normalized[component] = numeric
        return normalized

    @model_validator(mode='after')
    def _validate_aggregate_mass(self) -> 'SimulationSeedInfluentComposition':
        from puran_mcp_support.engineering.heuristics.asm2d_influent_fractionation import adm1_cod_kg_m3, adm1_tn_kg_m3, adm1_tp_kg_m3, asm2d_cod_kg_m3, asm2d_tn_kg_m3, asm2d_tp_kg_m3
        if self.fractionation_basis == 'asm2d':
            cod = asm2d_cod_kg_m3(self.species_kg_m3)
            tn = asm2d_tn_kg_m3(self.species_kg_m3)
            tp = asm2d_tp_kg_m3(self.species_kg_m3)
        else:
            cod = adm1_cod_kg_m3(self.species_kg_m3)
            tn = adm1_tn_kg_m3(self.species_kg_m3)
            tp = adm1_tp_kg_m3(self.species_kg_m3)
        if cod <= 0 or tn <= 0 or tp <= 0:
            raise ValueError('influent_composition species must imply positive COD, TN, and TP.')
        return self

class SimulationSeed(BaseModel):
    """Validated Stage-5 simulation handoff payload."""
    model_config = ConfigDict(extra='forbid', str_strip_whitespace=True)
    schema_version: str = 'simulation_seed/v2'
    seed_id: str
    engine_key: EngineKey
    unit_type: str | None = None
    qsdsan_class: str | None = None
    watertap_class: str | None = None
    unit_id_hint: str | None = None
    kinetic_model: str | None = None
    engineering_inputs: dict[str, Any] = Field(default_factory=dict)
    constructor_kwargs: dict[str, Any] = Field(default_factory=dict)
    sizing_reference: dict[str, Any] | None = None
    sizing_outputs: dict[str, dict[str, Any]] = Field(default_factory=dict)
    seed_bindings: dict[str, str] = Field(default_factory=dict)
    stream_topology: SimulationSeedTopology
    influent_component_hints: dict[str, Any] = Field(default_factory=dict)
    influent_composition: SimulationSeedInfluentComposition | None = None
    credibility: CanonicalCredibilityStamp
    notes: list[str] = Field(default_factory=list)

    @field_validator('schema_version', mode='before')
    @classmethod
    def _reject_null_schema_version(cls, value: Any) -> Any:
        if value is None:
            raise ValueError('`schema_version` may not be null.')
        return value

    @field_validator('seed_id')
    @classmethod
    def _validate_seed_id(cls, value: str) -> str:
        return _require_non_empty(value, field_name='seed_id')

    @field_validator('kinetic_model')
    @classmethod
    def _validate_kinetic_model(cls, value: str | None) -> str | None:
        if value is None:
            return None
        return _require_non_empty(value, field_name='kinetic_model')

    @field_validator('unit_type', 'qsdsan_class', 'watertap_class', 'unit_id_hint')
    @classmethod
    def _validate_optional_strings(cls, value: str | None, info) -> str | None:
        if value is None:
            return None
        return _require_non_empty(value, field_name=info.field_name)

    @model_validator(mode='after')
    def _validate_seed(self) -> 'SimulationSeed':
        if self.schema_version != 'simulation_seed/v2':
            raise ValueError(f"Unknown simulation_seed schema_version: {self.schema_version}. Only 'simulation_seed/v2' is supported.")
        port_count = len(self.stream_topology.inlets) + len(self.stream_topology.outlets)
        if self.engine_key != 'none' and port_count < 1:
            raise ValueError("Simulation seeds with engine_key!='none' require at least one topology port.")
        if self.engine_key == 'none':
            if self.unit_type is not None:
                raise ValueError("Simulation seeds with engine_key='none' must set unit_type=None.")
            if self.qsdsan_class is not None:
                raise ValueError("Simulation seeds with engine_key='none' must set qsdsan_class=None.")
            if self.watertap_class is not None:
                raise ValueError("Simulation seeds with engine_key='none' must set watertap_class=None.")
            return self
        if self.engine_key.startswith('WaterTAP-'):
            if self.watertap_class is None:
                raise ValueError('WaterTAP simulation seeds require `watertap_class`.')
            if self.unit_type is None:
                raise ValueError('WaterTAP simulation seeds require `unit_type`.')
            self._validate_boundary_influent_composition()
            return self
        if _is_qsdsan_engine_key(self.engine_key):
            if self.qsdsan_class is None:
                raise ValueError('QSDsan simulation seeds require `qsdsan_class`.')
            if self.unit_type is None:
                raise ValueError('QSDsan simulation seeds require `unit_type`.')
        self._validate_boundary_influent_composition()
        return self

    def _validate_boundary_influent_composition(self) -> None:
        basis = _required_fractionation_basis(self.engine_key)
        if basis is None:
            return
        requires = any((inlet.role == 'feed' and _is_plant_boundary_stream(inlet.stream_id) for inlet in self.stream_topology.inlets))
        if not requires:
            return
        if self.influent_composition is None:
            raise ValueError('E951 stage7b_seed_influent_missing: boundary ASM2D/ADM1 seeds must carry influent_composition populated by Stage 5 emit.')
        if self.influent_composition.fractionation_basis != basis:
            raise ValueError(f'influent_composition.fractionation_basis must match biological engine family {basis!r}.')
        raw_basis = self.engineering_inputs.get('influent_basis')
        if not isinstance(raw_basis, dict):
            raw_basis = self.constructor_kwargs.get('influent_basis')
        if not isinstance(raw_basis, dict):
            return
        from puran_mcp_support.engineering.heuristics.asm2d_influent_fractionation import TenderInfluentBasis, validate_adm1_mass_closure, validate_asm2d_mass_closure
        influent = TenderInfluentBasis.model_validate(raw_basis)
        if basis == 'asm2d':
            validate_asm2d_mass_closure(self.influent_composition.species_kg_m3, influent)
        else:
            validate_adm1_mass_closure(self.influent_composition.species_kg_m3, influent)

def _is_qsdsan_engine_key(engine_key: str) -> bool:
    """Return True for engine families compiled by qsdsan-engine-mcp."""
    normalized = engine_key.lower()
    return engine_key in _QSDSAN_ENGINE_KEYS or normalized.startswith('qsdsan.') or normalized.startswith('qsdsan-')

def _is_plant_boundary_stream(stream_id: str) -> bool:
    normalized = stream_id.strip().lower().replace('_', '-')
    if not normalized.startswith('plant:'):
        return False
    return any((token in normalized for token in ('influent', 'feed')))

def _required_fractionation_basis(engine_key: str) -> Literal['asm2d', 'adm1'] | None:
    if engine_key in _ASM2D_INFLUENT_REQUIRED_ENGINE_KEYS:
        return 'asm2d'
    if engine_key in _ADM1_INFLUENT_REQUIRED_ENGINE_KEYS:
        return 'adm1'
    return None
