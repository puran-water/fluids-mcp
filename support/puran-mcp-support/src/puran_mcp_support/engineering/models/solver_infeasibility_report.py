from __future__ import annotations
from datetime import datetime, timezone
UTC = timezone.utc
from typing import Any, Literal
from pydantic import BaseModel, ConfigDict, Field, model_validator

class SpeciesMassBalanceImbalance(BaseModel):
    """One CSTR/species mass-balance residual at the best iterate."""
    model_config = ConfigDict(extra='allow')
    unit_id: str = Field(..., min_length=1)
    species: str = Field(..., min_length=1)
    constraint_name: str = Field(..., min_length=1)
    residual: float
    percent_imbalance: float | None = None
    body_value: float | None = None
    bound: float | None = None

class VariableAtBound(BaseModel):
    """One variable sitting on or outside a numerical/physical bound."""
    model_config = ConfigDict(extra='allow')
    variable_name: str = Field(..., min_length=1)
    value: float | None = None
    lower_bound: float | None = None
    upper_bound: float | None = None
    bound_type: Literal['lower', 'upper', 'below_lower', 'above_upper', 'fixed']
    distance_to_bound: float | None = None

class ConstraintResidualSummary(BaseModel):
    """Largest active constraint residual at the best iterate."""
    model_config = ConfigDict(extra='allow')
    constraint_name: str = Field(..., min_length=1)
    residual: float
    body_value: float | None = None
    bound: float | None = None
    is_equality: bool = True

class SolverInfeasibilityReport(BaseModel):
    """Machine-readable failure artifact for non-converged WaterTAP solves."""
    model_config = ConfigDict(extra='allow')
    schema_version: Literal['solver_infeasibility_report/v1'] = 'solver_infeasibility_report/v1'
    generated_at: datetime = Field(default_factory=lambda: datetime.now(UTC))
    engine: str = 'watertap-engine-mcp'
    model_family: str = 'unknown'
    convergence_status: Literal['solver_infeasible_with_diagnostic'] = 'solver_infeasible_with_diagnostic'
    solver_status: str | None = None
    termination_condition: str | None = None
    message: str | None = None
    species_mass_balance_imbalances: list[SpeciesMassBalanceImbalance] = Field(default_factory=list)
    variables_at_bounds: list[VariableAtBound] = Field(default_factory=list)
    max_constraint_residual: ConstraintResidualSummary | None = None
    most_likely_cause: str = Field(..., min_length=1)
    suggested_remediation: list[str] = Field(default_factory=list)
    diagnostic_basis: dict[str, Any] = Field(default_factory=dict)

    @model_validator(mode='after')
    def _require_remediation(self) -> 'SolverInfeasibilityReport':
        if not any((str(item).strip() for item in self.suggested_remediation)):
            raise ValueError('suggested_remediation must contain at least one item')
        return self
