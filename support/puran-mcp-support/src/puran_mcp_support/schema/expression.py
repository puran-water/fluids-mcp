"""Typed, serializable expression graph shared by numeric and Pyomo adapters."""
from __future__ import annotations
import hashlib
import json
from dataclasses import dataclass
from typing import Any, Mapping, Union
Scalar = Union[int, float]

class Expression:

    def to_dict(self) -> dict[str, Any]:
        raise NotImplementedError

    @property
    def sha256(self) -> str:
        canonical = json.dumps(self.to_dict(), sort_keys=True, separators=(',', ':'), allow_nan=False)
        return hashlib.sha256(canonical.encode('utf-8')).hexdigest()

def as_expression(value: 'Expression | Scalar') -> 'Expression':
    return value if isinstance(value, Expression) else Constant(float(value))

@dataclass(frozen=True)
class Constant(Expression):
    value: float

    def to_dict(self) -> dict[str, Any]:
        return {'op': 'constant', 'value': self.value}

@dataclass(frozen=True)
class Symbol(Expression):
    name: str

    def to_dict(self) -> dict[str, Any]:
        return {'op': 'symbol', 'name': self.name}

@dataclass(frozen=True)
class Unary(Expression):
    op: str
    value: Expression

    def to_dict(self) -> dict[str, Any]:
        return {'op': self.op, 'value': self.value.to_dict()}

@dataclass(frozen=True)
class Binary(Expression):
    op: str
    left: Expression
    right: Expression

    def to_dict(self) -> dict[str, Any]:
        return {'op': self.op, 'left': self.left.to_dict(), 'right': self.right.to_dict()}

@dataclass(frozen=True)
class Nary(Expression):
    op: str
    values: tuple[Expression, ...]

    def to_dict(self) -> dict[str, Any]:
        return {'op': self.op, 'values': [value.to_dict() for value in self.values]}

@dataclass(frozen=True)
class SaturationRatio(Expression):
    activities: Mapping[str, Expression]
    stoichiometry: Mapping[str, float]
    equilibrium_constant: Expression

    def to_dict(self) -> dict[str, Any]:
        return {'op': 'saturation_ratio', 'activities': {key: self.activities[key].to_dict() for key in sorted(self.activities)}, 'stoichiometry': {key: self.stoichiometry[key] for key in sorted(self.stoichiometry)}, 'equilibrium_constant': self.equilibrium_constant.to_dict()}

def add(*values: Expression | Scalar) -> Expression:
    return Nary('add', tuple((as_expression(value) for value in values)))

def multiply(*values: Expression | Scalar) -> Expression:
    return Nary('multiply', tuple((as_expression(value) for value in values)))

def power(left: Expression | Scalar, right: Expression | Scalar) -> Expression:
    return Binary('power', as_expression(left), as_expression(right))

def exp(value: Expression | Scalar) -> Expression:
    return Unary('exp', as_expression(value))

def log(value: Expression | Scalar) -> Expression:
    return Unary('log', as_expression(value))

def smooth_positive(value: Expression | Scalar, epsilon: Expression | Scalar) -> Expression:
    return Binary('smooth_positive', as_expression(value), as_expression(epsilon))

def smooth_minimum(left: Expression | Scalar, right: Expression | Scalar, epsilon: Expression | Scalar) -> Expression:
    return Nary('smooth_minimum', (as_expression(left), as_expression(right), as_expression(epsilon)))

def monod(substrate: Expression | Scalar, half_saturation: Expression | Scalar) -> Expression:
    return Binary('monod', as_expression(substrate), as_expression(half_saturation))

def noncompetitive_inhibition(inhibitor: Expression | Scalar, inhibition_constant: Expression | Scalar) -> Expression:
    return Binary('noncompetitive_inhibition', as_expression(inhibitor), as_expression(inhibition_constant))

def nucleation_rate(k_n: Expression | Scalar, volume: Expression | Scalar, saturation_index: Expression | Scalar, critical_si: Expression | Scalar, exponent: Expression | Scalar, epsilon: Expression | Scalar) -> Expression:
    driving = smooth_positive(add(saturation_index, multiply(-1.0, critical_si)), epsilon)
    return multiply(k_n, volume, power(driving, exponent))

def growth_rate(k_g: Expression | Scalar, area: Expression | Scalar, saturation_ratio: Expression | Scalar, exponent: Expression | Scalar, epsilon: Expression | Scalar) -> Expression:
    driving = smooth_positive(add(saturation_ratio, -1.0), epsilon)
    return multiply(k_g, area, power(driving, exponent))

def dissolution_rate(k_d: Expression | Scalar, area: Expression | Scalar, saturation_ratio: Expression | Scalar, exponent: Expression | Scalar, epsilon: Expression | Scalar) -> Expression:
    driving = smooth_positive(add(1.0, multiply(-1.0, saturation_ratio)), epsilon)
    return multiply(-1.0, k_d, area, power(driving, exponent))
__all__ = ['Binary', 'Constant', 'Expression', 'Nary', 'SaturationRatio', 'Symbol', 'Unary', 'add', 'as_expression', 'dissolution_rate', 'exp', 'growth_rate', 'log', 'monod', 'multiply', 'noncompetitive_inhibition', 'nucleation_rate', 'power', 'smooth_minimum', 'smooth_positive']
