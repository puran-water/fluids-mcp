"""Pure-Python evaluator for the shared expression graph."""
from __future__ import annotations
import math
from typing import Mapping
from puran_mcp_support.schema.expression import Binary, Constant, Expression, Nary, SaturationRatio, Symbol, Unary

class ExpressionDomainError(ValueError):
    pass

def evaluate(expression: Expression, symbols: Mapping[str, float]) -> float:
    if isinstance(expression, Constant):
        return float(expression.value)
    if isinstance(expression, Symbol):
        if expression.name not in symbols:
            raise KeyError(f'missing expression symbol {expression.name!r}')
        value = float(symbols[expression.name])
        if not math.isfinite(value):
            raise ExpressionDomainError(f'symbol {expression.name!r} is non-finite')
        return value
    if isinstance(expression, Unary):
        value = evaluate(expression.value, symbols)
        if expression.op == 'exp':
            return math.exp(value)
        if expression.op == 'log':
            if value <= 0:
                raise ExpressionDomainError('log operand must be positive')
            return math.log(value)
        raise ValueError(f'unknown unary op {expression.op!r}')
    if isinstance(expression, Binary):
        left = evaluate(expression.left, symbols)
        right = evaluate(expression.right, symbols)
        if expression.op == 'power':
            return left ** right
        if expression.op == 'monod':
            denominator = left + right
            if denominator <= 0:
                raise ExpressionDomainError('Monod denominator must be positive')
            return left / denominator
        if expression.op == 'noncompetitive_inhibition':
            denominator = left + right
            if denominator <= 0:
                raise ExpressionDomainError('inhibition denominator must be positive')
            return right / denominator
        if expression.op == 'smooth_positive':
            if right <= 0:
                raise ExpressionDomainError('smoothing epsilon must be positive')
            return 0.5 * (left + math.sqrt(left * left + right * right))
        raise ValueError(f'unknown binary op {expression.op!r}')
    if isinstance(expression, Nary):
        values = [evaluate(value, symbols) for value in expression.values]
        if expression.op == 'add':
            return sum(values)
        if expression.op == 'multiply':
            return math.prod(values)
        if expression.op == 'smooth_minimum':
            if len(values) != 3 or values[2] <= 0:
                raise ExpressionDomainError('smooth_minimum requires left, right, epsilon>0')
            left, right, epsilon = values
            return 0.5 * (left + right - math.sqrt((left - right) ** 2 + epsilon ** 2))
        raise ValueError(f'unknown n-ary op {expression.op!r}')
    if isinstance(expression, SaturationRatio):
        iap = 1.0
        for species, coefficient in expression.stoichiometry.items():
            if species not in expression.activities:
                raise KeyError(f'missing activity expression for {species!r}')
            activity = evaluate(expression.activities[species], symbols)
            if activity < 0:
                raise ExpressionDomainError(f'activity {species!r} is negative')
            iap *= activity ** coefficient
        equilibrium_constant = evaluate(expression.equilibrium_constant, symbols)
        if equilibrium_constant <= 0:
            raise ExpressionDomainError('equilibrium constant must be positive')
        return iap / equilibrium_constant
    raise TypeError(f'unsupported expression type {type(expression).__name__}')

def directional_derivative(expression: Expression, symbols: Mapping[str, float], direction: Mapping[str, float], *, step: float=1e-07) -> float:
    plus = {key: float(value) + step * direction.get(key, 0.0) for key, value in symbols.items()}
    minus = {key: float(value) - step * direction.get(key, 0.0) for key, value in symbols.items()}
    return (evaluate(expression, plus) - evaluate(expression, minus)) / (2.0 * step)
__all__ = ['ExpressionDomainError', 'directional_derivative', 'evaluate']
