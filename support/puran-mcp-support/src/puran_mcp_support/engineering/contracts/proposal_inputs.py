from __future__ import annotations
from typing import Any, Literal, TypeVar
from pydantic import AliasChoices, BaseModel, ConfigDict, Field, ValidationError, field_validator, model_validator

class ProposalContractModel(BaseModel):
    """Base config for all proposal input contracts."""
    model_config = ConfigDict(extra='allow', populate_by_name=True)

class CanonicalCredibilityStamp(ProposalContractModel):
    status: str
    grade: str
    basis: str
    model_source: str | None = None
    calibration_dataset_ref: str | None = None
    plant_reference: str | None = None
    validity_range: str | None = None
    reviewer: str | None = None
    evidence_artifact_ref: str | None = None
    vendor_artifact_ref: str | None = None
    evidence_gate_reason: str | None = None

    @field_validator('grade', mode='before')
    @classmethod
    def _normalize_grade(cls, value: Any) -> Any:
        """Use the JSON-schema spelling while accepting legacy DB spelling."""
        if value == 'order_of_magnitude':
            return 'order-of-magnitude'
        return value
