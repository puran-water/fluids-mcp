"""Simultaneous, positivity-preserving reaction-extent application."""
from __future__ import annotations
import math
from dataclasses import dataclass, field
from typing import Mapping, Sequence
from puran_mcp_support.schema import ReactionDefinition

@dataclass(frozen=True)
class ExtentStepResult:
    state: Mapping[str, float]
    applied_extents: Mapping[str, float]
    common_scale: float
    closure_residuals: Mapping[str, float]
    reservoir_exchanges: Mapping[str, float] = field(default_factory=dict)

def constrained_reaction_step(state: Mapping[str, float], reactions: Sequence[ReactionDefinition], requested_extents: Mapping[str, float], *, nonnegative_tolerance: float=1e-12, implicit_reservoir_component_ids: Sequence[str]=()) -> ExtentStepResult:
    """Apply coupled reaction extents with one availability scale.

    All reactions see the same pre-step aqueous state. If their combined
    consumption would exceed any available component, every extent is scaled
    by the same limiting factor. This is order-independent and neither clips
    nor creates mass after the step.
    """
    reaction_by_id = {reaction.reaction_id: reaction for reaction in reactions}
    unknown = sorted(set(requested_extents) - set(reaction_by_id))
    if unknown:
        raise KeyError(f'requested extents reference unknown reactions {unknown}')
    if any((not math.isfinite(float(value)) for value in state.values())):
        raise ValueError('state contains non-finite values')
    if any((float(value) < -nonnegative_tolerance for value in state.values())):
        raise ValueError('state contains a materially negative inventory')
    delta: dict[str, float] = {}
    for reaction_id, extent in requested_extents.items():
        extent = float(extent)
        if not math.isfinite(extent):
            raise ValueError(f'extent for {reaction_id!r} is non-finite')
        reaction = reaction_by_id[reaction_id]
        for component_id, coefficient in reaction.stoichiometry.items():
            delta[component_id] = delta.get(component_id, 0.0) + coefficient * extent
    reservoir_ids = frozenset((str(component_id) for component_id in implicit_reservoir_component_ids))
    unknown_reservoirs = sorted(reservoir_ids - (set(state) | {component_id for reaction in reactions for component_id in reaction.stoichiometry}))
    if unknown_reservoirs:
        raise KeyError(f'implicit reservoirs reference unknown components {unknown_reservoirs}')
    common_scale = 1.0
    for component_id, change in delta.items():
        if change >= 0 or component_id in reservoir_ids:
            continue
        available = max(float(state.get(component_id, 0.0)), 0.0)
        common_scale = min(common_scale, available / -change if change else 1.0)
    common_scale = max(0.0, min(1.0, common_scale))
    applied = {reaction_id: float(extent) * common_scale for reaction_id, extent in requested_extents.items()}
    new_state = dict(state)
    closure: dict[str, float] = {}
    reservoir_exchanges: dict[str, float] = {}
    for component_id, change in delta.items():
        if component_id in reservoir_ids:
            new_state[component_id] = float(state.get(component_id, 0.0))
            reservoir_exchanges[component_id] = -common_scale * change
            closure[component_id] = 0.0
            continue
        expected = float(state.get(component_id, 0.0)) + common_scale * change
        if expected < 0 and expected >= -nonnegative_tolerance:
            expected = 0.0
        if expected < -nonnegative_tolerance:
            raise RuntimeError(f'constrained extent step produced {component_id}={expected:g}')
        new_state[component_id] = expected
        closure[component_id] = new_state[component_id] - (float(state.get(component_id, 0.0)) + common_scale * change)
    return ExtentStepResult(new_state, applied, common_scale, closure, reservoir_exchanges)
__all__ = ['ExtentStepResult', 'constrained_reaction_step']
