from __future__ import annotations
import hashlib
import json
import math
from datetime import datetime
from enum import Enum
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator
from puran_mcp_support.engineering.heuristics.ro.element_design_basis_8inch import DeviationAuthority

class ParameterTier(str, Enum):
    CATALOG = 'catalog'
    VENDOR_PROJECTION_FIT = 'vendor_projection_fit'
    HOUSE_DEFAULT = 'house_default'
    SITE_CALIBRATED = 'site_calibrated'
    FAMILY_PRIOR_CANDIDATE = 'family_prior_candidate'
    HOUSE_RULED = 'house_ruled'

class ParameterValue(BaseModel):
    model_config = ConfigDict(extra='forbid', frozen=True)
    value: float
    unit: str = Field(min_length=1, description="'1' for dimensionless")
    uncertainty: float | None = None
    note: str = ''

    @field_validator('value', 'uncertainty')
    @classmethod
    def _finite(cls, value: float | None) -> float | None:
        if value is not None and (not math.isfinite(value)):
            raise ValueError('parameter values must be finite')
        return value

class ValidityWindow(BaseModel):
    model_config = ConfigDict(extra='forbid', frozen=True)
    start: datetime
    end: datetime | None = None

    @model_validator(mode='after')
    def _ordered(self) -> 'ValidityWindow':
        if self.start.tzinfo is None or (self.end is not None and self.end.tzinfo is None):
            raise ValueError('validity windows are timezone-aware')
        if self.end is not None and self.end <= self.start:
            raise ValueError('a validity window must end after it starts')
        return self

    def covers(self, at: datetime) -> bool:
        if at.tzinfo is None:
            raise ValueError('compare a validity window with a timezone-aware instant')
        return self.start <= at and (self.end is None or at < self.end)

class ParameterProvenance(BaseModel):
    model_config = ConfigDict(extra='forbid', frozen=True)
    method: str = Field(min_length=1, description='how the values were obtained (fit, datasheet, standard ...)')
    produced_by: str = Field(min_length=1, description='the run, script or document that produced them')
    evidence_refs: tuple[str, ...] = ()
    ruling_ref: str | None = None
    authority: DeviationAuthority | None = None
    promoted_from: str | None = None
    source_object: str | None = None
    source_sha256: str | None = None

class PriorMember(BaseModel):
    """One site calibration a family prior aggregated."""
    model_config = ConfigDict(extra='forbid', frozen=True)
    asset_uid: str
    site: str
    window: ValidityWindow
    set_hash: str = Field(pattern='^[0-9a-f]{64}$')
    produced_by: str

class ParameterSpread(BaseModel):
    """The spread of one parameter over a prior's members, in the parameter's own unit."""
    model_config = ConfigDict(extra='forbid', frozen=True)
    n: int = Field(ge=2)
    mean: float
    median: float
    minimum: float
    maximum: float
    std: float = Field(ge=0.0)
    cv: float | None = None
    unit: str

class PriorAggregate(BaseModel):
    """How a ``family_prior_candidate`` was formed: its members (one site calibration per asset), the spread of each
    parameter, and the parameters left out because not every member carries them."""
    model_config = ConfigDict(extra='forbid', frozen=True)
    method: str = Field(min_length=1)
    members: tuple[PriorMember, ...] = Field(min_length=2)
    spread: dict[str, ParameterSpread]
    sites: tuple[str, ...]
    excluded_parameters: tuple[str, ...] = ()

    @model_validator(mode='after')
    def _one_set_per_asset(self) -> 'PriorAggregate':
        assets = [m.asset_uid for m in self.members]
        if len(set(assets)) != len(assets):
            raise ValueError('a family prior takes one site calibration per asset')
        if tuple(sorted({m.site for m in self.members})) != self.sites:
            raise ValueError("a family prior's sites are its members' sites, sorted")
        return self

class ParameterSet(BaseModel):
    model_config = ConfigDict(extra='forbid', frozen=True)
    schema_version: str = 'process_model_parameter_set/v1'
    model_key: str
    tier: ParameterTier
    values: dict[str, ParameterValue] = Field(min_length=1)
    provenance: ParameterProvenance
    subject: str | None = None
    settings: dict[str, str | bool] = Field(default_factory=dict)
    asset_uid: str | None = None
    site: str | None = None
    window: ValidityWindow | None = None
    window_basis: str | None = None
    aggregate: PriorAggregate | None = None

    @model_validator(mode='after')
    def _tier_carries_its_scope(self) -> 'ParameterSet':
        if self.tier is ParameterTier.SITE_CALIBRATED:
            if not (self.asset_uid and self.site and self.window):
                raise ValueError('a site_calibrated set names its site, asset_uid and validity window')
            if not self.asset_uid.startswith(f'{self.site}/'):
                raise ValueError(f'asset {self.asset_uid!r} is not an asset of site {self.site!r}')
        elif self.asset_uid or self.site:
            raise ValueError(f'a {self.tier.value} set is not bound to one asset or site')
        if self.tier is ParameterTier.HOUSE_RULED:
            if not (self.provenance.ruling_ref or '').strip():
                raise ValueError('a house_ruled set cites its ruling (provenance.ruling_ref)')
            authority = self.provenance.authority
            if authority is None or authority.status != 'ruled':
                raise ValueError("a house_ruled set carries the ruled authority that promoted it (provenance.authority, status 'ruled'); a proposal pending a ruling is not one")
            if authority.ruling_ref != self.provenance.ruling_ref:
                raise ValueError('provenance.ruling_ref and provenance.authority.ruling_ref cite the same ruling')
            if not self.provenance.promoted_from:
                raise ValueError('a house_ruled set names the family-prior candidate it promoted (provenance.promoted_from)')
        else:
            if self.provenance.ruling_ref:
                raise ValueError(f'a {self.tier.value} set cites no ruling; only house_ruled does')
            if self.provenance.authority is not None or self.provenance.promoted_from:
                raise ValueError(f'a {self.tier.value} set carries no promotion authority; only house_ruled does')
        if self.tier is ParameterTier.FAMILY_PRIOR_CANDIDATE:
            if self.aggregate is None:
                raise ValueError('a family_prior_candidate set records how it was aggregated (aggregate)')
            if set(self.aggregate.spread) != set(self.values):
                raise ValueError('a family prior reports the spread of every value it carries, and of nothing else')
        elif self.aggregate is not None:
            raise ValueError(f'a {self.tier.value} set is not an aggregate of site calibrations')
        return self

    def set_hash(self) -> str:
        """sha256 of the set's canonical JSON: the identity a prior member and a promotion cite."""
        canonical = json.dumps(self.model_dump(mode='json'), sort_keys=True, separators=(',', ':'))
        return hashlib.sha256(canonical.encode()).hexdigest()

    def to_json(self) -> str:
        return self.model_dump_json(indent=1)

    @classmethod
    def from_json(cls, text: str) -> 'ParameterSet':
        return cls.model_validate_json(text)
