"""Dependency-light endpoint contract for the independent PHREEQC audit."""
from __future__ import annotations
import math
import re
from collections.abc import Mapping, Sequence
from typing import Any
from puran_mcp_support.schema import MADM1_FECL3_FE_S_P_KERNEL, ChemistryKernelDefinition
_REQUIRED_VERTICAL_SLICE_TOTALS = frozenset({'P', 'S(6)', 'S(-2)', 'Fe(3)', 'Fe(2)', 'C(4)', 'Cl'})
_SUPPORTED_AUDIT_TOTALS = frozenset({'Al', 'C(4)', 'Ca', 'Cl', 'Fe(2)', 'Fe(3)', 'K', 'Mg', 'N(-3)', 'N(5)', 'Na', 'P', 'S(-2)', 'S(6)'})
_GAS_MODES = frozenset({'closed', 'fixed_partial_pressure', 'finite_headspace'})
_SUPPORTED_GASES = frozenset({'ch4', 'ch4(g)', 'co2', 'co2(g)', 'h2', 'h2(g)', 'h2s', 'h2s(g)', 'n2', 'n2(g)', 'nh3', 'nh3(g)', 'o2', 'o2(g)'})

def _finite(value: Any, *, field: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f'{field} must be numeric, not bool')
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f'{field} must be finite')
    return number

def required_phase_activity_ids(phase_ids: Sequence[str], *, kernel: ChemistryKernelDefinition=MADM1_FECL3_FE_S_P_KERNEL) -> tuple[str, ...]:
    """Return the exact aqueous/solvent activities used by selected phases."""
    phases = {phase.phase_id: phase for phase in kernel.phases}
    normalized = tuple((str(phase_id) for phase_id in phase_ids))
    if not normalized or len(set(normalized)) != len(normalized):
        raise ValueError('phase_ids must be nonempty and unique')
    missing = sorted(set(normalized) - set(phases))
    if missing:
        raise KeyError(f'kernel {kernel.kernel_id!r} omits phases {missing}')
    return tuple(sorted({component_id for phase_id in normalized for component_id in phases[phase_id].dissolution_reaction if component_id != phase_id}))

def selected_phase_activities(phase_ids: Sequence[str], activities: Mapping[str, float], *, solvent_water_activity: float=1.0, kernel: ChemistryKernelDefinition=MADM1_FECL3_FE_S_P_KERNEL) -> dict[str, float]:
    """Select all and only the activities consumed by the audited phases."""
    required = required_phase_activity_ids(phase_ids, kernel=kernel)
    available = {str(key): _finite(value, field=f'activity.{key}') for key, value in activities.items()}
    available['H2O'] = _finite(solvent_water_activity, field='solvent_water_activity')
    missing = sorted(set(required) - set(available))
    if missing:
        raise KeyError(f'endpoint omitted phase activities {missing}')
    selected = {component_id: available[component_id] for component_id in required}
    invalid = sorted((component_id for component_id, value in selected.items() if value <= 0.0))
    if invalid:
        raise ValueError(f'endpoint phase activities must be positive; invalid={invalid}')
    return selected

def build_endpoint_thermodynamic_audit_request(*, totals_mol_l: Mapping[str, float], temperature_c: float, ph: float, phase_saturation_log10_si: Mapping[str, float], activities: Mapping[str, float], source_mass_balance_passed: bool, kinetic_phase_ids: Sequence[str] | None=None, gas_boundary_mode: str='closed', gas_boundary: Mapping[str, float] | None=None, gas_volume_l: float | None=None, gas_total_pressure_atm: float | None=None, solvent_water_activity: float=1.0, species_coverage_fraction: float=0.0, reduced_model_database_hash: str | None=None, database: str='minteq.v4.dat', mechanism_interpretation: Mapping[str, Any] | None=None, kernel: ChemistryKernelDefinition=MADM1_FECL3_FE_S_P_KERNEL) -> dict[str, Any]:
    """Build the exact public request consumed by water-chemistry-mcp.

    This helper intentionally performs no PHREEQC execution.  It prevents the
    independently deployed QSDsan and WaterTAP adapters from inventing subtly
    different endpoint field names or selected-activity rules.
    """
    totals = {str(key): _finite(value, field=f'totals_mol_l.{key}') for key, value in totals_mol_l.items()}
    missing_totals = sorted(_REQUIRED_VERTICAL_SLICE_TOTALS - set(totals))
    unsupported_totals = sorted(set(totals) - _SUPPORTED_AUDIT_TOTALS)
    invalid_totals = sorted((key for key, value in totals.items() if value < 0.0))
    if missing_totals or unsupported_totals or invalid_totals:
        raise ValueError(f'endpoint analytical totals are incomplete or invalid; missing={missing_totals}, unsupported={unsupported_totals}, negative={invalid_totals}')
    if not isinstance(source_mass_balance_passed, bool):
        raise ValueError('source_mass_balance_passed must be boolean')
    temperature = _finite(temperature_c, field='temperature_c')
    ph_value = _finite(ph, field='ph')
    if temperature <= -273.15:
        raise ValueError('temperature_c must be above absolute zero')
    if not 0.0 <= ph_value <= 14.0:
        raise ValueError('ph must be in [0, 14]')
    saturation = {str(key): _finite(value, field=f'phase_saturation.{key}') for key, value in phase_saturation_log10_si.items()}
    phase_ids = tuple(saturation)
    required_phase_activity_ids(phase_ids, kernel=kernel)
    kinetic = phase_ids if kinetic_phase_ids is None else tuple((str(phase_id) for phase_id in kinetic_phase_ids))
    if set(kinetic) - set(phase_ids):
        raise ValueError('kinetic_phase_ids must be a subset of phase_ids')
    mode = str(gas_boundary_mode)
    if mode not in _GAS_MODES:
        raise ValueError(f'unsupported gas_boundary_mode {mode!r}')
    gases = {str(key): _finite(value, field=f'gas_boundary.{key}') for key, value in (gas_boundary or {}).items()}
    unsupported_gases = sorted((key for key in gases if key.strip().lower() not in _SUPPORTED_GASES))
    negative_gases = sorted((key for key, value in gases.items() if value < 0.0))
    if unsupported_gases or negative_gases:
        raise ValueError(f'gas boundary contains unsupported or negative entries; unsupported={unsupported_gases}, negative={negative_gases}')
    volume = None if gas_volume_l is None else _finite(gas_volume_l, field='gas_volume_l')
    pressure = None if gas_total_pressure_atm is None else _finite(gas_total_pressure_atm, field='gas_total_pressure_atm')
    if mode == 'closed':
        if gases or volume is not None or pressure is not None:
            raise ValueError('closed gas boundary cannot carry gas values')
    elif mode == 'fixed_partial_pressure':
        if not gases or any((value <= 0.0 for value in gases.values())):
            raise ValueError('fixed gas partial pressures must be positive')
        if volume is not None or pressure is not None:
            raise ValueError('fixed partial-pressure boundary cannot carry volume/pressure')
    else:
        if not gases or not any((value > 0.0 for value in gases.values())):
            raise ValueError('finite headspace requires a positive gas component')
        if volume is None or volume <= 0.0:
            raise ValueError('finite headspace requires positive gas_volume_l')
        if pressure is not None and pressure <= 0.0:
            raise ValueError('gas_total_pressure_atm must be positive')
    coverage = _finite(species_coverage_fraction, field='species_coverage_fraction')
    if not 0.0 <= coverage <= 1.0:
        raise ValueError('species_coverage_fraction must be in [0, 1]')
    if reduced_model_database_hash is not None and (not re.fullmatch('[0-9a-fA-F]{64}', str(reduced_model_database_hash))):
        raise ValueError('reduced_model_database_hash must be a SHA-256 hex digest')
    interpretation = dict(mechanism_interpretation or {})
    return {'totals_mol_l': totals, 'temperature_c': temperature, 'ph': ph_value, 'gas_boundary_mode': mode, 'gas_boundary': gases, 'gas_volume_l': volume, 'gas_total_pressure_atm': pressure, 'phase_ids': phase_ids, 'kinetic_phase_ids': kinetic, 'reduced_model_saturation': saturation, 'reduced_model_saturation_convention': 'log10_si', 'reduced_model_activities': selected_phase_activities(phase_ids, activities, solvent_water_activity=solvent_water_activity, kernel=kernel), 'source_mass_balance_passed': bool(source_mass_balance_passed), 'species_coverage_fraction': coverage, 'reduced_model_database_hash': None if reduced_model_database_hash is None else str(reduced_model_database_hash).lower(), 'reduced_model_kernel_hash': kernel.sha256, 'database': str(database), 'surface_dataset_id': interpretation.get('surface_dataset_id'), 'surface_parameter_interpretation': interpretation.get('surface_parameter_interpretation', 'source_not_supplied'), 'oxide_history_class': interpretation.get('oxide_history_class', 'unresolved'), 'oxide_history_id': interpretation.get('oxide_history_id', 'not_supplied'), 'formation_ph': interpretation.get('formation_ph'), 'formation_ph_basis': interpretation.get('formation_ph_basis', 'endpoint_pH_not_oxide_formation_pH'), 'active_metal_fraction_basis': interpretation.get('active_metal_fraction_basis', 'not_identified'), 'effective_active_fe_fraction': interpretation.get('effective_active_fe_fraction'), 'parameter_contact': interpretation.get('parameter_contact'), 'parameter_contact_evidence_id': interpretation.get('parameter_contact_evidence_id'), 'mechanism_identified': bool(interpretation.get('mechanism_identified', False)), 'mechanism_evidence_id': interpretation.get('mechanism_evidence_id'), 'extrapolation_status': interpretation.get('extrapolation_status', 'not_assessed'), 'extrapolation_reasons': list(interpretation.get('extrapolation_reasons') or ())}
__all__ = ['build_endpoint_thermodynamic_audit_request', 'required_phase_activity_ids', 'selected_phase_activities']
