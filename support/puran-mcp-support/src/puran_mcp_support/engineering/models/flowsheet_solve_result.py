from __future__ import annotations
from typing import Any, Literal
from pydantic import BaseModel, ConfigDict, Field
from puran_mcp_support.engineering.models.solver_infeasibility_report import SolverInfeasibilityReport

class ResidualBasis(BaseModel):
    """Structured breakdown of the residuals behind a near-feasible acceptance.

    Populated ONLY when a solve is accepted ``near_feasible`` (the square ASM2D
    recycle case where IPOPT terminates non-optimal but every equation and bound
    is satisfied within the proposal tolerance). It records *why* the iterate was
    deemed acceptable so a downstream consumer can audit the trade rather than
    reading a non-optimal solve as a strict optimum. Every field is optional so a
    caller can backfill whatever the diagnostics surfaced.
    """
    model_config = ConfigDict(extra='forbid')
    max_equality_residual: float | None = None
    worst_constraint_name: str | None = None
    constraint_residual_tolerance: float | None = None
    constraint_residual_count: int | None = None
    max_bound_violation: float | None = None
    worst_bound_variable: str | None = None
    bound_violation_tolerance: float | None = None
    bound_violation_count: int | None = None

class SolveResult(BaseModel):
    """Structured response from watertap-engine-mcp.solve_flowsheet."""
    model_config = ConfigDict(extra='forbid')
    status: Literal['optimal', 'non_physical', 'infeasible_structured', 'invalid_property_range', 'convergence_budget_exceeded', 'topology_template_not_found', 'feedstock_class_not_found']
    solution_quality: Literal['strict_optimal', 'near_feasible', 'non_physical', 'infeasible_structured'] | None = None
    residual_basis: ResidualBasis | dict[str, Any] | None = None
    physicality: dict[str, Any] | None = None
    solve_seconds: float = 0.0
    effluent: dict[str, Any] = Field(default_factory=dict)
    design_parameters_settled: dict[str, Any] = Field(default_factory=dict)
    stream_closure: dict[str, Any] = Field(default_factory=dict)
    diagnostics: dict[str, Any] = Field(default_factory=dict)
    solver_infeasibility_report: SolverInfeasibilityReport | dict[str, Any] | None = None
    flowsheet_svg: str = ''
    artifact_id: str = ''
    artifact_path: str = ''
    warnings: list[str] = Field(default_factory=list)
