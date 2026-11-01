"""Positive inventory stepping for dynamic CSTR/MBR chemistry adapters."""
from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Any, Mapping, Sequence
from puran_mcp_support.schema import ComponentDefinition, ReactionDefinition
from .extents import ExtentStepResult, constrained_reaction_step

@dataclass(frozen=True)
class InventoryStepResult:
    state: Mapping[str, float]
    reaction_step: ExtentStepResult
    hydraulic_state: Mapping[str, float]
    wasting_removed: Mapping[str, float]
    balance_report: Mapping[str, Any] | None = None
_BALANCE_PROPERTIES = ('H', 'O', 'C', 'N', 'P', 'S', 'Fe', 'Al', 'Ca', 'Mg', 'Na', 'K', 'Cl', 'charge', 'electron_equivalent')

def _step_balance_report(*, initial_state: Mapping[str, float], hydraulic_state: Mapping[str, float], final_state: Mapping[str, float], wasting_removed: Mapping[str, float], reservoir_exchanges: Mapping[str, float], component_definitions: Sequence[ComponentDefinition], absolute_tolerance: float=1e-10, relative_tolerance: float=1e-10) -> dict[str, Any]:
    """Verify that internal reactions neither create nor destroy properties."""
    component_by_id = {component.component_id: component for component in component_definitions}
    if len(component_by_id) != len(component_definitions):
        raise ValueError('component_definitions contain duplicate component IDs')
    state_ids = set(initial_state) | set(hydraulic_state) | set(final_state)
    unknown = sorted(state_ids - set(component_by_id))
    if unknown:
        raise KeyError(f'balance audit lacks component definitions for {unknown}')
    net_internal_change = {component_id: float(final_state.get(component_id, 0.0)) - float(initial_state.get(component_id, 0.0)) - (float(hydraulic_state.get(component_id, 0.0)) - float(initial_state.get(component_id, 0.0))) + float(wasting_removed.get(component_id, 0.0)) - float(reservoir_exchanges.get(component_id, 0.0)) for component_id in state_ids | set(reservoir_exchanges)}

    def coefficient(component_id: str, property_id: str) -> float:
        component = component_by_id[component_id]
        if property_id == 'charge':
            return float(component.charge)
        if property_id == 'electron_equivalent':
            return float(component.electron_equivalent_mol_e)
        return float(component.elemental_coefficients.get(property_id, 0.0))
    residuals: dict[str, float] = {}
    tolerances: dict[str, float] = {}
    property_passed: dict[str, bool] = {}
    for property_id in _BALANCE_PROPERTIES:
        contributions = [coefficient(component_id, property_id) * change for component_id, change in net_internal_change.items()]
        residual = float(sum(contributions))
        scale = float(sum((abs(item) for item in contributions)))
        tolerance = absolute_tolerance + relative_tolerance * scale
        residuals[property_id] = residual
        tolerances[property_id] = tolerance
        property_passed[property_id] = math.isfinite(residual) and abs(residual) <= tolerance
    return {'basis': 'moles_or_equivalents_per_inventory_volume_over_step', 'residuals': residuals, 'tolerances': tolerances, 'property_passed': property_passed, 'passed': all(property_passed.values()), 'implicit_reservoir_exchanges': dict(reservoir_exchanges)}

def advance_cstr_inventory(state: Mapping[str, float], *, reactions: Sequence[ReactionDefinition], rates_mol_per_volume_time: Mapping[str, float], inlet: Mapping[str, float], dose_source_mol_per_volume_time: Mapping[str, float], aqueous_component_ids: Sequence[str], retained_solid_ids: Sequence[str], component_definitions: Sequence[ComponentDefinition] | None=None, dt: float, hydraulic_residence_time: float, solids_retention_time: float, hydraulic_exchange_fractions: Mapping[str, float] | None=None, hydraulic_inlet_fractions: Mapping[str, float] | None=None, hydraulic_outlet_fractions: Mapping[str, float] | None=None, implicit_reservoir_component_ids: Sequence[str]=()) -> InventoryStepResult:
    """Advance a well-mixed inventory without post-step concentration clipping.

    Hydraulic exchange and nonnegative dose sources are integrated analytically,
    all chemical extents are then applied simultaneously through the shared
    availability constraint, and retained-solid wasting is integrated exactly.
    This Lie-split step is positivity preserving; timestep convergence remains a
    caller responsibility for accuracy.
    """
    values = (dt, hydraulic_residence_time, solids_retention_time)
    if any((not math.isfinite(float(value)) or value <= 0 for value in values)):
        raise ValueError('dt, HRT, and SRT must be finite and positive')
    if any((float(value) < 0 for value in dose_source_mol_per_volume_time.values())):
        raise ValueError('dose sources must be nonnegative; represent sinks as reactions')
    hydraulic = {key: float(value) for key, value in state.items()}
    exchange = {component_id: 1.0 for component_id in aqueous_component_ids}
    exchange.update({component_id: 0.0 for component_id in retained_solid_ids if component_id not in exchange})
    for component_id, fraction in (hydraulic_exchange_fractions or {}).items():
        converted = float(fraction)
        if not math.isfinite(converted) or not 0.0 <= converted <= 1.0:
            raise ValueError('hydraulic exchange fractions must be finite in [0, 1]')
        exchange[str(component_id)] = converted
    inlet_fractions = dict(exchange)
    outlet_fractions = dict(exchange)
    for supplied, destination, label in ((hydraulic_inlet_fractions or {}, inlet_fractions, 'hydraulic inlet'), (hydraulic_outlet_fractions or {}, outlet_fractions, 'hydraulic outlet')):
        for component_id, fraction in supplied.items():
            converted = float(fraction)
            if not math.isfinite(converted) or not 0.0 <= converted <= 1.0:
                raise ValueError(f'{label} fractions must be finite in [0, 1]')
            destination[str(component_id)] = converted
    hydraulic_ids = set(inlet_fractions) | set(outlet_fractions)
    unknown_exchange = sorted(hydraulic_ids - set(state) - set(inlet) - set(dose_source_mol_per_volume_time))
    if unknown_exchange:
        raise KeyError(f'hydraulic exchange references unknown components {unknown_exchange}')
    for component_id in hydraulic_ids:
        inlet_fraction = inlet_fractions.get(component_id, 0.0)
        outlet_fraction = outlet_fractions.get(component_id, 0.0)
        current = float(state.get(component_id, 0.0))
        feed = float(inlet.get(component_id, 0.0))
        source = float(dose_source_mol_per_volume_time.get(component_id, 0.0))
        if min(current, feed, source) < 0:
            raise ValueError(f'negative hydraulic inventory/source for {component_id!r}')
        inlet_rate = inlet_fraction * feed / float(hydraulic_residence_time)
        if outlet_fraction == 0.0:
            hydraulic[component_id] = current + (inlet_rate + source) * float(dt)
        else:
            effective_rate = outlet_fraction / float(hydraulic_residence_time)
            decay = math.exp(-effective_rate * float(dt))
            equilibrium = (inlet_rate + source) / effective_rate
            hydraulic[component_id] = equilibrium + (current - equilibrium) * decay
    requested_extents = {reaction_id: float(rate) * float(dt) for reaction_id, rate in rates_mol_per_volume_time.items()}
    reacted = constrained_reaction_step(hydraulic, reactions, requested_extents, implicit_reservoir_component_ids=implicit_reservoir_component_ids)
    final_state = dict(reacted.state)
    wasting_removed: dict[str, float] = {}
    wasting_decay = math.exp(-float(dt) / float(solids_retention_time))
    for component_id in retained_solid_ids:
        before = float(final_state.get(component_id, 0.0))
        after = before * wasting_decay
        final_state[component_id] = after
        wasting_removed[component_id] = before - after
    if any((value < -1e-12 or not math.isfinite(value) for value in final_state.values())):
        raise RuntimeError('inventory step produced a nonphysical state')
    balance_report = None
    if component_definitions is not None:
        balance_report = _step_balance_report(initial_state=state, hydraulic_state=hydraulic, final_state=final_state, wasting_removed=wasting_removed, reservoir_exchanges=reacted.reservoir_exchanges, component_definitions=component_definitions)
    return InventoryStepResult(final_state, reacted, hydraulic, wasting_removed, balance_report)
__all__ = ['InventoryStepResult', 'advance_cstr_inventory']
