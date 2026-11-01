from __future__ import annotations
import math
from enum import Enum
from pydantic import field_validator, model_validator
from puran_mcp_support.engineering.pydantic_base import StrictFiniteBaseModel

class Distribution(str, Enum):
    """Sampling distribution for a factor.

    ``uniform`` samples linearly on ``[low, high]``. ``loguniform`` samples
    uniformly on ``[log10(low), log10(high)]`` (requires ``low > 0``); the driver
    exponentiates back to engineering units before evaluation.
    """
    uniform = 'uniform'
    loguniform = 'loguniform'

class Factor(StrictFiniteBaseModel):
    """A single varied input dimension."""
    name: str
    low: float
    high: float
    distribution: Distribution = Distribution.uniform
    nominal: float | None = None
    description: str | None = None
    units: str | None = None

    @field_validator('name')
    @classmethod
    def _name_nonempty(cls, value: str) -> str:
        if not value or not value.strip():
            raise ValueError('Factor.name must be a non-empty string')
        return value

    @model_validator(mode='after')
    def _check_bounds(self) -> 'Factor':
        if not self.high > self.low:
            raise ValueError(f'Factor {self.name!r}: high ({self.high}) must exceed low ({self.low})')
        if self.distribution is Distribution.loguniform and self.low <= 0:
            raise ValueError(f'Factor {self.name!r}: loguniform requires low > 0 (got {self.low})')
        if self.nominal is not None and (not self.low <= self.nominal <= self.high):
            raise ValueError(f'Factor {self.name!r}: nominal ({self.nominal}) outside [low, high]')
        return self

    @property
    def is_log(self) -> bool:
        return self.distribution is Distribution.loguniform

    @property
    def nominal_value(self) -> float:
        """Nominal value; defaults to the (geometric for log) midpoint."""
        if self.nominal is not None:
            return self.nominal
        if self.is_log:
            return float(math.sqrt(self.low * self.high))
        return float(0.5 * (self.low + self.high))

    def salib_bounds(self) -> list[float]:
        """Bounds in the space SALib samples (log10 space for loguniform)."""
        if self.is_log:
            return [math.log10(self.low), math.log10(self.high)]
        return [self.low, self.high]

class GSAMetadata(StrictFiniteBaseModel):
    """Provenance carried through to the result schema.

    ``input_hash`` is a placeholder the per-engine adapter fills with a hash of
    the base flowsheet / influent so a study is reproducible against a known
    model state.
    """
    engine: str | None = None
    model_version: str | None = None
    input_hash: str | None = None
    label: str | None = None
    notes: str | None = None

class GSAProblem(StrictFiniteBaseModel):
    """Typed, serialisable GSA spec: factors + KPI names + metadata."""
    factors: list[Factor]
    outputs: list[str]
    metadata: GSAMetadata = GSAMetadata()

    @model_validator(mode='after')
    def _check_unique(self) -> 'GSAProblem':
        if not self.factors:
            raise ValueError('GSAProblem requires at least one factor')
        names = [f.name for f in self.factors]
        if len(set(names)) != len(names):
            raise ValueError(f'Factor names must be unique, got {names}')
        if not self.outputs:
            raise ValueError('GSAProblem requires at least one output (KPI) name')
        if len(set(self.outputs)) != len(self.outputs):
            raise ValueError(f'Output names must be unique, got {self.outputs}')
        return self

    @property
    def factor_names(self) -> list[str]:
        return [f.name for f in self.factors]

    @property
    def log_mask(self) -> list[bool]:
        return [f.is_log for f in self.factors]

    def to_salib_problem(self) -> dict:
        """The SALib problem dict: ``{num_vars, names, bounds}``.

        Bounds for loguniform factors are in log10 space; the driver maps samples
        back with ``10 ** value`` before ``transform``/``evaluate``.
        """
        return {'num_vars': len(self.factors), 'names': self.factor_names, 'bounds': [f.salib_bounds() for f in self.factors]}

    def nominal_inputs(self) -> dict[str, float]:
        """Each factor at its nominal value (engineering units, not log space)."""
        return {f.name: f.nominal_value for f in self.factors}
