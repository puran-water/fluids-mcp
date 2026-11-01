"""Independent whole-window inventory balance utilities.

The dynamic engines own their native state bases and prepare property
inventories/fluxes.  This module owns only the numerical ledger so every
adapter uses the same sign convention and tolerance rule.
"""
from __future__ import annotations
import math
from typing import Any, Mapping, Sequence

def _finite_series(values: Sequence[float], *, label: str, expected_length: int) -> tuple[float, ...]:
    if len(values) != expected_length:
        raise ValueError(f'{label} length {len(values)} != time length {expected_length}')
    converted = tuple((float(value) for value in values))
    if any((not math.isfinite(value) for value in converted)):
        raise ValueError(f'{label} contains nonfinite values')
    return converted

def _trapezoid(times: Sequence[float], values: Sequence[float]) -> float:
    return float(sum((0.5 * (left + right) * (t_right - t_left) for left, right, t_left, t_right in zip(values, values[1:], times, times[1:]))))

def _simpson_nonuniform(times: Sequence[float], values: Sequence[float]) -> float:
    """Integrate quadratic pairs on a nonuniform grid, with one terminal trap."""
    total = 0.0
    paired_stop = len(times) if len(times) % 2 else len(times) - 1
    for index in range(0, paired_stop - 2, 2):
        h0 = times[index + 1] - times[index]
        h1 = times[index + 2] - times[index + 1]
        span = h0 + h1
        total += span / 6.0 * ((2.0 - h1 / h0) * values[index] + span * span / (h0 * h1) * values[index + 1] + (2.0 - h0 / h1) * values[index + 2])
    if len(times) % 2 == 0:
        total += 0.5 * (values[-2] + values[-1]) * (times[-1] - times[-2])
    return float(total)

def integrate_time_series(times: Sequence[float], values: Sequence[float], *, quadrature_method: str='trapezoid') -> float:
    """Integrate one finite time series using the governed ledger quadrature.

    This small public helper lets engine adapters report disaggregated outlet
    terms (for example, permeate versus sludge wasting) with exactly the same
    quadrature used by :func:`audit_dynamic_inventory_window`.
    """
    converted_times = tuple((float(value) for value in times))
    if len(converted_times) < 2:
        raise ValueError('at least two time samples are required')
    if any((not math.isfinite(value) for value in converted_times)):
        raise ValueError('times contain nonfinite values')
    if any((right <= left for left, right in zip(converted_times, converted_times[1:]))):
        raise ValueError('times must be strictly increasing')
    converted_values = _finite_series(values, label='values', expected_length=len(converted_times))
    quadrature = {'trapezoid': _trapezoid, 'simpson_nonuniform': _simpson_nonuniform}.get(str(quadrature_method))
    if quadrature is None:
        raise ValueError("quadrature_method must be 'trapezoid' or 'simpson_nonuniform'")
    return quadrature(converted_times, converted_values)

def audit_dynamic_inventory_window(times: Sequence[float], *, inventory_by_property: Mapping[str, Sequence[float]], inlet_flux_by_property: Mapping[str, Sequence[float]], outlet_flux_by_property: Mapping[str, Sequence[float]], absolute_tolerance: float=1e-08, relative_tolerance: float=0.0001, inventory_unit: str='amount', flux_unit: str='amount_per_time', time_unit: str='time', quadrature_method: str='trapezoid') -> dict[str, Any]:
    """Audit ``accumulation = inlet - outlet`` over a sampled time window.

    The calculation is deliberately independent of reaction-rate expressions.
    A dynamic adapter supplies total retained inventories and all external
    property fluxes in one consistent native basis.  Internal reactions should
    therefore cancel. Fluxes are integrated over the supplied time samples and
    each property is checked with a combined absolute/relative tolerance.
    ``simpson_nonuniform`` uses quadratic panels on irregular samples and a
    terminal trapezoid when the sample count is even.
    """
    converted_times = tuple((float(value) for value in times))
    if len(converted_times) < 2:
        raise ValueError('at least two time samples are required')
    if any((not math.isfinite(value) for value in converted_times)):
        raise ValueError('times contain nonfinite values')
    if any((right <= left for left, right in zip(converted_times, converted_times[1:]))):
        raise ValueError('times must be strictly increasing')
    if not math.isfinite(float(absolute_tolerance)) or float(absolute_tolerance) < 0:
        raise ValueError('absolute_tolerance must be finite and nonnegative')
    if not math.isfinite(float(relative_tolerance)) or float(relative_tolerance) < 0:
        raise ValueError('relative_tolerance must be finite and nonnegative')
    integrate_time_series(converted_times, [0.0] * len(converted_times), quadrature_method=quadrature_method)
    properties = set(inventory_by_property)
    if properties != set(inlet_flux_by_property) or properties != set(outlet_flux_by_property):
        raise ValueError('inventory, inlet-flux, and outlet-flux property IDs must match')
    if not properties:
        raise ValueError('at least one property must be audited')
    rows: dict[str, dict[str, float | bool]] = {}
    for property_id in sorted(properties):
        inventory = _finite_series(inventory_by_property[property_id], label=f'inventory[{property_id}]', expected_length=len(converted_times))
        inlet = _finite_series(inlet_flux_by_property[property_id], label=f'inlet_flux[{property_id}]', expected_length=len(converted_times))
        outlet = _finite_series(outlet_flux_by_property[property_id], label=f'outlet_flux[{property_id}]', expected_length=len(converted_times))
        inlet_integral = integrate_time_series(converted_times, inlet, quadrature_method=quadrature_method)
        outlet_integral = integrate_time_series(converted_times, outlet, quadrature_method=quadrature_method)
        accumulation = inventory[-1] - inventory[0]
        expected_accumulation = inlet_integral - outlet_integral
        residual = accumulation - expected_accumulation
        scale = max(abs(inventory[0]), abs(inventory[-1]), abs(accumulation), abs(inlet_integral) + abs(outlet_integral))
        tolerance = float(absolute_tolerance) + float(relative_tolerance) * scale
        normalized_residual = residual / scale if scale > 0 else residual
        rows[property_id] = {'initial_inventory': inventory[0], 'final_inventory': inventory[-1], 'accumulation': accumulation, 'integrated_inlet': inlet_integral, 'integrated_outlet': outlet_integral, 'expected_accumulation': expected_accumulation, 'residual': residual, 'normalized_residual': normalized_residual, 'tolerance': tolerance, 'passed': abs(residual) <= tolerance}
    return {'status': 'complete', 'equation': 'final_inventory - initial_inventory = integrated_inlet - integrated_outlet', 'window_start': converted_times[0], 'window_end': converted_times[-1], 'window_duration': converted_times[-1] - converted_times[0], 'sample_count': len(converted_times), 'inventory_unit': inventory_unit, 'flux_unit': flux_unit, 'time_unit': time_unit, 'absolute_tolerance': float(absolute_tolerance), 'relative_tolerance': float(relative_tolerance), 'quadrature_method': str(quadrature_method), 'properties': rows, 'property_passed': {property_id: bool(row['passed']) for property_id, row in rows.items()}, 'passed': all((bool(row['passed']) for row in rows.values()))}
__all__ = ['audit_dynamic_inventory_window', 'integrate_time_series']
