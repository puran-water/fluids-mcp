from __future__ import annotations
from typing import Literal
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator
MODEL_KEY = 'watertap.asm2d_mbr/v1'
COMPONENT_NAMESPACE = 'watertap.modified_asm2d/v1'

class NativeBiologyParameterField(BaseModel):
    """One source-prescribed native parameter, with its native unit."""
    model_config = ConfigDict(extra='forbid', frozen=True, allow_inf_nan=False)
    value: float
    unit: str = Field(min_length=1)
    source_clause_ref: str = Field(min_length=1)

    @field_validator('value', mode='before')
    @classmethod
    def _numeric_value(cls, value):
        if isinstance(value, bool) or not isinstance(value, (float, int)):
            raise ValueError('native parameter values must be explicit numbers')
        return value

    @field_validator('unit', 'source_clause_ref')
    @classmethod
    def _trimmed_text(cls, value):
        if value != value.strip() or not value:
            raise ValueError('native parameter units and source clauses must be nonempty trimmed text')
        return value

class NativeBiologyParameterBasis(BaseModel):
    """Immutable evidence for the complete temperature/parameter prescription."""
    model_config = ConfigDict(extra='forbid', frozen=True, allow_inf_nan=False)
    schema_version: Literal['native_biology_parameter_basis/v1']
    model_key: Literal['watertap.asm2d_mbr/v1']
    component_namespace: Literal['watertap.modified_asm2d/v1']
    fields: dict[str, NativeBiologyParameterField]
