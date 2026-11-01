from __future__ import annotations
from enum import Enum
from typing import Iterable

class UnitOpType(str, Enum):
    """Canonical unit-operation categories for cross-engine handoffs."""
    MEMBRANE = 'membrane'
    MEMBRANE_ZERO_ORDER = 'membrane_zero_order'
    THERMAL = 'thermal'
    PUMP = 'pump'
    PUMP_ZERO_ORDER = 'pump_zero_order'
    ENERGY_RECOVERY = 'energy_recovery'
    MIXER = 'mixer'
    SPLITTER = 'splitter'
    FEED = 'feed'
    FEED_ZERO_ORDER = 'feed_zero_order'
    PRODUCT = 'product'
    CRYSTALLIZER = 'crystallizer'
    BIOLOGICAL = 'biological'
    REACTOR = 'reactor'
    SEPARATOR = 'separator'
    CLARIFIER = 'clarifier'
    SLUDGE = 'sludge'
    JUNCTION = 'junction'
    UTILITY = 'utility'
    PRETREATMENT = 'pretreatment'
    OTHER = 'other'
UNIT_OP_TYPE_VALUES: tuple[str, ...] = tuple((member.value for member in UnitOpType))

def normalize_unit_op_type(value: UnitOpType | str) -> UnitOpType:
    """Normalize enum names or values into a canonical ``UnitOpType``."""
    if isinstance(value, UnitOpType):
        return value
    token = str(value).strip()
    if not token:
        raise ValueError('Unit operation type cannot be empty')
    try:
        return UnitOpType(token)
    except ValueError:
        normalized_name = token.upper().replace('-', '_').replace(' ', '_')
        try:
            return UnitOpType[normalized_name]
        except KeyError as exc:
            valid = ', '.join(UNIT_OP_TYPE_VALUES)
            raise ValueError(f'Unknown unit operation type {value!r}. Valid values: {valid}') from exc

def unit_op_type_values(types: Iterable[UnitOpType] | None=None) -> tuple[str, ...]:
    """Return stable string values for all or selected canonical unit types."""
    members = UnitOpType if types is None else types
    return tuple((member.value for member in members))
