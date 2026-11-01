from __future__ import annotations
from datetime import datetime
from enum import Enum
from typing import Any, Literal
from uuid import UUID
from pydantic import BaseModel, Field, model_validator

class RiskSeverity(str, Enum):
    """Severity ladder for simulation-layer risk flags.

    - ``info``: notable but does not require mitigation; never blocks compile.
    - ``warning``: design attention item; may be proposal-visible if
      ``RiskFlag.client_safe``; never blocks compile.
    - ``critical``: design-risk hard stop. Blocks compile unconditionally via
      stage_simulation_risk_gate (E1901). PR5e intentionally does not ship a
      waiver workflow in-band; if reviewer disposition is needed later, add a
      separate ``simulation_risk_review`` artifact rather than mutating
      engine-produced diagnostics.
    """
    info = 'info'
    warning = 'warning'
    critical = 'critical'

class EvidenceRef(BaseModel):
    """Back-reference to a diagnostic field or timeseries segment.

    ``path`` is a JSONPath-like pointer into ``EngineDiagnostics`` or the
    inline timeseries. ``artifact_uri`` points at the full-resolution
    ``simulation_timeseries`` artifact when the evidence is a trajectory.
    """
    path: str
    summary: str | None = None
    artifact_uri: str | None = None

class RiskFlag(BaseModel):
    """Normalized simulation risk flag surfaced at the engine layer.

    Codes follow the pattern ``RISK-<FAMILY>-<KIND>-<NN>`` and are defined in
    ``engineering_utils.diagnostics.risk_rules``. Example codes:
    ``RISK-ASM-INHIB-SNH-01``, ``RISK-RO-SCALE-CaSO4-02``,
    ``RISK-THERMAL-SUPERSAT-01``, ``RISK-IX-LEAK-01``,
    ``RISK-PHYSCHEM-DISINFECT-CT-01``, ``RISK-SLUDGE-POLYMER-01``,
    ``RISK-ODOR-BREAKTHROUGH-01``.
    """
    code: str
    severity: RiskSeverity
    title: str
    summary: str
    evidence: list[EvidenceRef] = Field(default_factory=list)
    proposal_narrative_hint: str
    mitigation_hint: str | None = None
    client_safe: bool = True

class SolverMetadata(BaseModel):
    """Solver hygiene metadata. Never proposal-facing by default."""
    solver_name: str | None = None
    termination_condition: str | None = None
    solve_time_s: float | None = None
    iterations: int | None = None
    dof: int | None = None
    internal_only: bool = True

class ConvergenceSummary(BaseModel):
    status: Literal['converged', 'converged_with_warnings', 'failed', 'not_run']
    message: str | None = None
    residual_norm: float | None = None

class TimeseriesSeries(BaseModel):
    """A single downsampled trajectory for inline rendering.

    The full-resolution counterpart lives in the parent
    ``EngineTimeseries.timeseries_full_uri`` / ``timeseries_full_artifact_id``.
    """
    series_id: str
    unit_op: str | None = None
    variable: str
    units: str | None = None
    values: list[float]

class EngineTimeseries(BaseModel):
    """Downsampled per-engine trajectory payload embedded in the report.

    The inline series are client-safe and capped for proposal rendering.
    ``timeseries_full_uri`` / ``timeseries_full_artifact_id`` reference the
    comprehensive full-resolution ``simulation_timeseries`` artifact that
    carries every state variable the engine exposed (see PR5e "comprehensive
    export" rule — temporal downsampling is the only compression applied).

    Invariants (enforced):
    - Every ``series[*].values`` has the same length as ``time``.
    - ``series_id`` is unique across ``series``.
    - ``timeseries_full_uri`` and ``timeseries_full_artifact_id`` are either
      both present or both absent (paired).
    """
    time: list[float]
    time_unit: Literal['s', 'min', 'h', 'd']
    downsample_method: Literal['lttb', 'uniform', 'event_preserving']
    series: list[TimeseriesSeries] = Field(default_factory=list)
    timeseries_full_uri: str | None = None
    timeseries_full_artifact_id: str | None = None
    timeseries_full_content_hash: str | None = None

    @model_validator(mode='after')
    def _check_invariants(self) -> 'EngineTimeseries':
        expected_len = len(self.time)
        for s in self.series:
            if len(s.values) != expected_len:
                raise ValueError(f"EngineTimeseries: series '{s.series_id}' has {len(s.values)} values but time axis has {expected_len} points")
        ids = [s.series_id for s in self.series]
        dup = {sid for sid in ids if ids.count(sid) > 1}
        if dup:
            raise ValueError(f'EngineTimeseries: duplicate series_id(s): {sorted(dup)}')
        has_uri = self.timeseries_full_uri is not None
        has_id = self.timeseries_full_artifact_id is not None
        has_hash = self.timeseries_full_content_hash is not None
        if has_uri != has_id:
            raise ValueError('EngineTimeseries: timeseries_full_uri and timeseries_full_artifact_id must be either both set or both None')
        if has_hash and (not (has_uri and has_id)):
            raise ValueError('EngineTimeseries: timeseries_full_content_hash requires timeseries_full_uri and timeseries_full_artifact_id')
        return self

class UniversalDiagnostics(BaseModel):
    """Engine-family-agnostic diagnostic block."""
    warnings: list[str] = Field(default_factory=list)
    convergence: ConvergenceSummary
    solver_metadata: SolverMetadata | None = None
    coverage_status: Literal['full', 'partial', 'proxy', 'unavailable'] = 'partial'
    coverage_notes: list[str] = Field(default_factory=list)

class EngineDiagnostics(BaseModel):
    """Engine-family-specific diagnostic envelope.

    ``universal`` is always present. Exactly one of the family-specific
    blocks should be populated per ``EngineReport``, matching
    ``EngineReport.engine_family``.

    Family blocks remain ``dict[str, Any]`` while producer contracts stabilize,
    but biological engine reports must now carry a populated ``asm_bio`` block.
    """
    universal: UniversalDiagnostics
    topology_source: str | None = None
    asm_bio: dict[str, Any] | None = None
    ro: dict[str, Any] | None = None
    thermal: dict[str, Any] | None = None
    ix: dict[str, Any] | None = None
    physchem: dict[str, Any] | None = None
    sludge: dict[str, Any] | None = None
    odor: dict[str, Any] | None = None
    process_inhibition: dict[str, Any] | None = None

class EngineReport(BaseModel):
    """One engine's contribution to the consolidated simulation_report.

    ``predicted_streams`` carries multi-outlet streams (permeate, retentate,
    distillate, brine, regen effluent, etc.) — a single ``predicted_effluent``
    is structurally insufficient for RO/evap/IX. ``predicted_effluent`` is
    the train-terminal outlet from this engine's perspective and is what
    feeds the plant-level merge.
    """
    engine: str
    engine_family: Literal['bio', 'ro', 'thermal', 'ix', 'physchem', 'sludge', 'odor']
    model_type: str
    unit_ops: list[str] = Field(default_factory=list)
    flowsheet_svg: str | None = None
    predicted_streams: dict[str, dict[str, Any]] = Field(default_factory=dict)
    predicted_effluent: dict[str, Any] | None = None
    diagnostics: EngineDiagnostics
    risk_flags: list[RiskFlag] = Field(default_factory=list)
    timeseries: EngineTimeseries | None = None

    @model_validator(mode='after')
    def _require_bio_diagnostics(self) -> 'EngineReport':
        if self.engine_family != 'bio':
            return self
        if not isinstance(self.flowsheet_svg, str) or not self.flowsheet_svg.strip():
            raise ValueError('EngineReport bio rows require non-empty flowsheet_svg')
        if not isinstance(self.diagnostics.asm_bio, dict) or not self.diagnostics.asm_bio:
            raise ValueError('EngineReport bio rows require non-empty diagnostics.asm_bio')
        return self

class SimulationReport(BaseModel):
    """Stage-(c) consolidated simulation report for a plant.

    ``predicted_effluent`` is the reconciled plant-terminal outlet, merged
    from each engine's ``predicted_effluent`` (downstream engine wins for
    overlapping parameters). The compliance gate reads from this field.

    Stage is fixed to ``"c"`` — this model represents the consolidated
    multi-engine contract only.
    """
    project_ref: str
    stage: Literal['c'] = 'c'
    plant_slug: str | None = None
    workflow_run_id: UUID | None = None
    process_balance_run_id: UUID | None = None
    process_balance_result_hash: str | None = Field(default=None, pattern='^[0-9a-f]{64}$')
    generated_at: datetime
    unit_ops: list[str] = Field(default_factory=list)
    predicted_effluent: dict[str, Any] = Field(default_factory=dict)
    engine_reports: list[EngineReport] = Field(default_factory=list)
    metadata: dict[str, Any] = Field(default_factory=dict)
