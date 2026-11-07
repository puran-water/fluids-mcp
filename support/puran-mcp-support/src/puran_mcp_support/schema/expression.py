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
FACTOR_NAMES: dict[str, str] = {}

def named_factor(name: str, expression: Expression) -> Expression:
    """``expression`` unchanged, its name recorded in `FACTOR_NAMES`; one factor never carries two names."""
    key = expression.sha256
    if FACTOR_NAMES.setdefault(key, name) != name:
        raise ValueError(f'factor {name!r} is already named {FACTOR_NAMES[key]!r}')
    return expression

def _symbol_names(expression: Expression) -> list[str]:
    found: set[str] = set()

    def walk(node):
        if isinstance(node, dict):
            if node.get('op') == 'symbol':
                found.add(node['name'])
            for value in node.values():
                walk(value)
        elif isinstance(node, list):
            for value in node:
                walk(value)
    walk(expression.to_dict())
    return sorted(found)

def product_factors(expression: Expression) -> tuple[tuple[str, Expression], ...]:
    """``expression`` as the factors of its product, in order, each ``(name, factor)``: a named factor
    (`named_factor`) as itself; a nested ``multiply`` flattened; any other factor named by its operator and the symbols
    it reads (``monod(K_so4_hSRB, S_SO4)``), never a guessed name. The product of the factors is ``expression``."""
    name = FACTOR_NAMES.get(expression.sha256)
    if name is not None:
        return ((name, expression),)
    if isinstance(expression, Nary) and expression.op == 'multiply':
        return tuple((factor for value in expression.values for factor in product_factors(value)))
    if isinstance(expression, Symbol):
        return ((expression.name, expression),)
    if isinstance(expression, Constant):
        return ((f'constant {expression.value:g}', expression),)
    op = getattr(expression, 'op', None) or 'saturation_ratio'
    return ((f'{op}({', '.join(_symbol_names(expression))})', expression),)
__all__ = ['Binary', 'Constant', 'Expression', 'FACTOR_NAMES', 'Nary', 'SaturationRatio', 'Symbol', 'Unary', 'add', 'as_expression', 'dissolution_rate', 'exp', 'growth_rate', 'log', 'monod', 'multiply', 'named_factor', 'noncompetitive_inhibition', 'nucleation_rate', 'power', 'product_factors', 'smooth_minimum', 'smooth_positive']
