"""Serial Fe-S-P chemistry train shared by dynamic and steady adapters."""
from __future__ import annotations
import hashlib
import json
import math
from collections.abc import Mapping
from dataclasses import asdict
from typing import Any
from puran_mcp_support.schema import MADM1_FECL3_FE_S_P_KERNEL, ChemistryKernelDefinition, SerialChemistryTrainDefinition, mechanism_interpretation_from_ledger
from .audit import build_endpoint_thermodynamic_audit_request
from .charge import infer_fixed_background_charge_eq_m3, reduced_charge_density_eq_m3, solve_relative_charge_balance
from .extents import constrained_reaction_step
from .reactive import REDUCED_AQUEOUS_COMPONENT_IDS, advance_reduced_chemistry_cstr, evaluate_reduced_chemistry_rates
_REQUIRED_FEED_COMPONENTS = frozenset({'Fe_III_aq', 'Fe_II_aq', 'H2_aq', 'HS_minus', 'SO4', 'PO4', 'CO3', 'Cl_minus'})
SERIAL_CHEMISTRY_CHECKPOINT_SCHEMA = 'process-chemistry-serial-checkpoint/v1'

def _finite_nonnegative(value: Any, *, field: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f'{field} must be numeric, not bool')
    number = float(value)
    if not math.isfinite(number) or number < 0.0:
        raise ValueError(f'{field} must be finite and nonnegative')
    return number

def _complete_state(supplied: Mapping[str, float], *, component_ids: frozenset[str], inlet_ph: float) -> dict[str, float]:
    unknown = sorted(set(supplied) - component_ids)
    if unknown:
        raise KeyError(f'serial chemistry state contains unknown components {unknown}')
    state = {component_id: 0.0 for component_id in component_ids}
    state.update({str(component_id): _finite_nonnegative(value, field=str(component_id)) for component_id, value in supplied.items()})
    if 'H_plus' not in supplied:
        state['H_plus'] = 10.0 ** (-inlet_ph) * 1000.0
    if 'H2O' not in supplied:
        state['H2O'] = 55500.0
    return state

def _reactive_p_inventory(state: Mapping[str, float]) -> float:
    return float(state.get('PO4', 0.0)) + 2.0 * float(state.get('Vivianite', 0.0)) + float(state.get('HFO_fresh_surface_P', 0.0)) + float(state.get('HFO_reoxidized_surface_P', 0.0)) + float(state.get('HFO_aged_surface_P', 0.0))

def _state_hash(states: Mapping[str, Mapping[str, float]]) -> str:
    canonical = json.dumps(states, sort_keys=True, separators=(',', ':'), allow_nan=False)
    return hashlib.sha256(canonical.encode('utf-8')).hexdigest()

def _stable_hash(value: Any) -> str:
    canonical = json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)
    return hashlib.sha256(canonical.encode('utf-8')).hexdigest()

def serial_chemistry_train_fingerprint(train: SerialChemistryTrainDefinition, *, kernel: ChemistryKernelDefinition=MADM1_FECL3_FE_S_P_KERNEL, phase_parameters: Mapping[str, Mapping[str, float]] | None=None, redox_parameters: Mapping[str, Mapping[str, float]] | None=None, oxidation_parameters: Mapping[str, Mapping[str, float]] | None=None, aging_parameters: Mapping[str, Mapping[str, float]] | None=None, surface_parameters: Mapping[str, Mapping[str, float]] | None=None, ligand_parameters: Mapping[str, float] | None=None, zone_parameter_overrides: Mapping[str, Mapping[str, Any]] | None=None) -> dict[str, Any]:
    """Fingerprint structural serial chemistry state semantics.

    Temperature, feed, dose, activity coefficients, and imposed pH are
    operating forcings and may change between windows.  The reaction/phase/
    surface kinetic profile and hydraulic train define the meaning of retained
    checkpoint states and therefore require a reseed when changed.
    """
    train_definition = asdict(train)
    configuration = {'phase_parameters': phase_parameters or {}, 'redox_parameters': redox_parameters or {}, 'oxidation_parameters': oxidation_parameters or {}, 'aging_parameters': aging_parameters or {}, 'surface_parameters': surface_parameters or {}, 'ligand_parameters': ligand_parameters or {}, 'zone_parameter_overrides': zone_parameter_overrides or {}}
    cell_ids = tuple((f'{zone.zone_id}:{cell_index + 1}' for zone in train.zones for cell_index in range(zone.tanks_in_series)))
    return {'schema_version': SERIAL_CHEMISTRY_CHECKPOINT_SCHEMA, 'kernel_id': kernel.kernel_id, 'kernel_hash': kernel.sha256, 'train_id': train.train_id, 'train_definition_hash': _stable_hash(train_definition), 'model_configuration_hash': _stable_hash(configuration), 'cell_ids': list(cell_ids)}

def _restore_serial_checkpoint(checkpoint: Mapping[str, Any], *, expected_fingerprint: Mapping[str, Any], component_ids: frozenset[str], inlet_ph: float) -> tuple[dict[str, dict[str, float]], float, float]:
    if not isinstance(checkpoint, Mapping):
        raise ValueError('serial chemistry checkpoint must be a mapping')
    mismatches: list[str] = []
    for field in ('schema_version', 'kernel_id', 'kernel_hash', 'train_id', 'train_definition_hash', 'model_configuration_hash', 'cell_ids'):
        if checkpoint.get(field) != expected_fingerprint.get(field):
            mismatches.append(f'{field} mismatch ({checkpoint.get(field)!r} != {expected_fingerprint.get(field)!r})')
    if checkpoint.get('persistable') is not True:
        mismatches.append('checkpoint is not marked persistable')
    raw_states = checkpoint.get('cell_states_mol_m3')
    if not isinstance(raw_states, Mapping):
        mismatches.append('cell_states_mol_m3 is missing or not a mapping')
        raw_states = {}
    expected_cells = set(expected_fingerprint['cell_ids'])
    supplied_cells = set(raw_states)
    if supplied_cells != expected_cells:
        mismatches.append(f'cell-state IDs differ (missing={sorted(expected_cells - supplied_cells)}, extra={sorted(supplied_cells - expected_cells)})')
    if mismatches:
        raise ValueError('serial chemistry checkpoint is incompatible; reseed required: ' + '; '.join(mismatches))
    states = {str(cell_id): _complete_state(state, component_ids=component_ids, inlet_ph=inlet_ph) for cell_id, state in raw_states.items()}
    actual_inventory_id = _state_hash(states)
    if checkpoint.get('final_inventory_id') != actual_inventory_id:
        raise ValueError('serial chemistry checkpoint inventory hash mismatch; reseed required')
    time_days = _finite_nonnegative(checkpoint.get('time_days', 0.0), field='serial_checkpoint.time_days')
    try:
        background = float(checkpoint.get('fixed_background_charge_eq_m3'))
    except (TypeError, ValueError) as exc:
        raise ValueError('serial chemistry checkpoint background charge is invalid; reseed required') from exc
    if not math.isfinite(background):
        raise ValueError('serial chemistry checkpoint background charge is invalid; reseed required')
    return (states, time_days, background)
_LEDGER_PROPERTIES = ('H', 'O', 'C', 'N', 'P', 'S', 'Fe', 'Al', 'Ca', 'Mg', 'Na', 'K', 'Cl', 'charge', 'electron_equivalent')

def _whole_train_balance(*, initial_states: Mapping[str, Mapping[str, float]], final_states: Mapping[str, Mapping[str, float]], cell_volumes_m3: Mapping[str, float], external_input_mol: Mapping[str, float], external_output_mol: Mapping[str, float], wasting_mol: Mapping[str, float], reservoir_exchange_mol: Mapping[str, float], kernel: ChemistryKernelDefinition, absolute_tolerance: float=1e-06, relative_tolerance: float=0.002) -> dict[str, Any]:
    component_by_id = {component.component_id: component for component in kernel.components}

    def inventory(states: Mapping[str, Mapping[str, float]], component_id: str) -> float:
        return sum((float(states[cell_id].get(component_id, 0.0)) * cell_volumes_m3[cell_id] for cell_id in states))
    component_rows: dict[str, dict[str, float]] = {}
    for component_id in component_by_id:
        initial = inventory(initial_states, component_id)
        final = inventory(final_states, component_id)
        inlet = float(external_input_mol.get(component_id, 0.0))
        outlet = float(external_output_mol.get(component_id, 0.0))
        wasted = float(wasting_mol.get(component_id, 0.0))
        reservoir = float(reservoir_exchange_mol.get(component_id, 0.0))
        residual = final - initial - inlet + outlet + wasted - reservoir
        component_rows[component_id] = {'initial_inventory_mol': initial, 'final_inventory_mol': final, 'external_input_mol': inlet, 'external_output_mol': outlet, 'wasting_output_mol': wasted, 'implicit_reservoir_input_mol': reservoir, 'residual_mol': residual}
    properties: dict[str, dict[str, float | bool]] = {}
    for property_id in _LEDGER_PROPERTIES:

        def coefficient(component_id: str) -> float:
            component = component_by_id[component_id]
            if property_id == 'charge':
                return float(component.charge)
            if property_id == 'electron_equivalent':
                return float(component.electron_equivalent_mol_e)
            return float(component.elemental_coefficients.get(property_id, 0.0))
        initial = sum((coefficient(component_id) * row['initial_inventory_mol'] for component_id, row in component_rows.items()))
        final = sum((coefficient(component_id) * row['final_inventory_mol'] for component_id, row in component_rows.items()))
        inlet = sum((coefficient(component_id) * row['external_input_mol'] for component_id, row in component_rows.items()))
        outlet = sum((coefficient(component_id) * (row['external_output_mol'] + row['wasting_output_mol']) for component_id, row in component_rows.items()))
        reservoir = sum((coefficient(component_id) * row['implicit_reservoir_input_mol'] for component_id, row in component_rows.items()))
        residual = final - initial - inlet + outlet - reservoir
        scale = max(abs(initial), abs(final), abs(inlet) + abs(outlet) + abs(reservoir), 1.0)
        tolerance = absolute_tolerance + relative_tolerance * scale
        properties[property_id] = {'initial_inventory': initial, 'final_inventory': final, 'external_input': inlet, 'external_output': outlet, 'implicit_reservoir_input': reservoir, 'residual': residual, 'tolerance': tolerance, 'passed': abs(residual) <= tolerance}
    governed = ('C', 'N', 'P', 'S', 'Fe', 'Al', 'Ca', 'Mg', 'Na', 'K', 'Cl', 'electron_equivalent')
    return {'equation': 'final - initial - external_input + permeate_output + wasting_output - reservoir_input = 0', 'integration_rule': 'trapezoidal terminal-port flux; exact constant feed/dose/source and wasting ledgers', 'absolute_tolerance': absolute_tolerance, 'relative_tolerance': relative_tolerance, 'component_rows': component_rows, 'properties': properties, 'governed_properties': governed, 'passed': all((bool(properties[property_id]['passed']) for property_id in governed)), 'acid_base_projection_properties': {property_id: bool(properties[property_id]['passed']) for property_id in ('H', 'O', 'charge')}}

def simulate_serial_chemistry_train(train: SerialChemistryTrainDefinition, *, feed_mol_m3: Mapping[str, float], inlet_ph: float, horizon_days: float, time_step_days: float, dose_mol_d: Mapping[str, float] | None=None, initial_state_by_zone: Mapping[str, Mapping[str, float]] | None=None, initial_cell_states_mol_m3: Mapping[str, Mapping[str, float]] | None=None, resume_checkpoint: Mapping[str, Any] | None=None, temperature_k: float=298.15, kernel: ChemistryKernelDefinition=MADM1_FECL3_FE_S_P_KERNEL, phase_parameters: Mapping[str, Mapping[str, float]] | None=None, redox_parameters: Mapping[str, Mapping[str, float]] | None=None, oxidation_parameters: Mapping[str, Mapping[str, float]] | None=None, aging_parameters: Mapping[str, Mapping[str, float]] | None=None, surface_parameters: Mapping[str, Mapping[str, float]] | None=None, ligand_parameters: Mapping[str, float] | None=None, activity_coefficients: Mapping[str, float] | None=None, solvent_water_activity: float=1.0, species_coverage_fraction: float=0.0, zone_parameter_overrides: Mapping[str, Mapping[str, Any]] | None=None, fixed_background_charge_eq_m3: float | None=None, imposed_ph_by_zone: Mapping[str, float] | None=None, sample_every_steps: int=1) -> dict[str, Any]:
    """Advance an ordered injection-to-MBR chemistry train.

    Biological transformations and native headspace gas transfer remain owned
    by the host mADM1 engine.  This function owns the transported mineral,
    oxide-history, surface-P, redox-donor, and counterion inventories.  Every
    nonclosed gas boundary is therefore rejected here unless the host engine
    has already supplied its liquid endpoint as the feed boundary.
    """
    ph = float(inlet_ph)
    horizon = float(horizon_days)
    dt_requested = float(time_step_days)
    temperature = float(temperature_k)
    water_activity = _finite_nonnegative(solvent_water_activity, field='solvent_water_activity')
    coverage = _finite_nonnegative(species_coverage_fraction, field='species_coverage_fraction')
    if not math.isfinite(ph) or not 0.0 <= ph <= 14.0:
        raise ValueError('inlet_ph must be finite in [0, 14]')
    if not math.isfinite(horizon) or horizon <= 0.0:
        raise ValueError('horizon_days must be finite and positive')
    if not math.isfinite(dt_requested) or dt_requested <= 0.0:
        raise ValueError('time_step_days must be finite and positive')
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise ValueError('temperature_k must be finite and positive')
    if water_activity <= 0.0:
        raise ValueError('solvent_water_activity must be positive')
    if coverage > 1.0:
        raise ValueError('species_coverage_fraction must be in [0, 1]')
    if isinstance(sample_every_steps, bool) or int(sample_every_steps) != sample_every_steps or sample_every_steps < 1:
        raise ValueError('sample_every_steps must be a positive integer')
    if any((zone.gas_boundary_mode != 'closed' for zone in train.zones)):
        raise ValueError('the shared serial mineral kernel accepts closed liquid zones only; couple finite-headspace/fixed-pressure gas transfer in native mADM1 and pass its liquid endpoint explicitly')
    component_ids = frozenset((component.component_id for component in kernel.components))
    missing_feed = sorted(_REQUIRED_FEED_COMPONENTS - set(feed_mol_m3))
    if missing_feed:
        raise KeyError(f'serial chemistry feed omitted analytical totals {missing_feed}')
    feed = _complete_state(feed_mol_m3, component_ids=component_ids, inlet_ph=ph)
    dose = {str(component_id): _finite_nonnegative(value, field=f'dose.{component_id}') for component_id, value in (dose_mol_d or {}).items()}
    unknown_dose = sorted(set(dose) - set(REDUCED_AQUEOUS_COMPONENT_IDS))
    if unknown_dose:
        raise KeyError(f'serial chemistry dose contains unknown aqueous components {unknown_dose}')
    zone_overrides = {str(key): dict(value) for key, value in (zone_parameter_overrides or {}).items()}
    unknown_zone_overrides = sorted(set(zone_overrides) - {zone.zone_id for zone in train.zones})
    if unknown_zone_overrides:
        raise KeyError(f'zone parameter overrides reference unknown zones {unknown_zone_overrides}')
    imposed = {str(key): float(value) for key, value in (imposed_ph_by_zone or {}).items()}
    unknown_imposed = sorted(set(imposed) - {zone.zone_id for zone in train.zones})
    if unknown_imposed:
        raise KeyError(f'imposed pH references unknown zones {unknown_imposed}')
    for zone in train.zones:
        if zone.local_ph_model in {'fixed', 'external_property_package'}:
            if zone.zone_id not in imposed:
                raise ValueError(f'zone {zone.zone_id!r} requires an explicit imposed pH')
            if not math.isfinite(imposed[zone.zone_id]) or not 0.0 <= imposed[zone.zone_id] <= 14.0:
                raise ValueError(f'zone {zone.zone_id!r} imposed pH must be in [0, 14]')
    retained_ids = tuple((component.component_id for component in kernel.components if component.component_id not in REDUCED_AQUEOUS_COMPONENT_IDS and component.component_id != 'H2O' and (component.permeability_class.value == 'retained')))
    initial_by_zone = {str(key): value for key, value in (initial_state_by_zone or {}).items()}
    unknown_initial = sorted(set(initial_by_zone) - {zone.zone_id for zone in train.zones})
    if unknown_initial:
        raise KeyError(f'initial state references unknown zones {unknown_initial}')
    if sum((value is not None for value in (initial_state_by_zone, initial_cell_states_mol_m3, resume_checkpoint))) > 1:
        raise ValueError('supply only one of initial_state_by_zone, initial_cell_states_mol_m3, or resume_checkpoint')
    cell_states: dict[str, dict[str, float]] = {}
    cell_volumes: dict[str, float] = {}
    cell_metadata: list[tuple[Any, int, str]] = []
    for zone in train.zones:
        for cell_index in range(zone.tanks_in_series):
            cell_id = f'{zone.zone_id}:{cell_index + 1}'
            cell_volumes[cell_id] = zone.volume_m3 / zone.tanks_in_series
            cell_metadata.append((zone, cell_index, cell_id))
    fingerprint = serial_chemistry_train_fingerprint(train, kernel=kernel, phase_parameters=phase_parameters, redox_parameters=redox_parameters, oxidation_parameters=oxidation_parameters, aging_parameters=aging_parameters, surface_parameters=surface_parameters, ligand_parameters=ligand_parameters, zone_parameter_overrides=zone_overrides)
    starting_checkpoint_time = 0.0
    checkpoint_background: float | None = None
    if resume_checkpoint is not None:
        cell_states, starting_checkpoint_time, checkpoint_background = _restore_serial_checkpoint(resume_checkpoint, expected_fingerprint=fingerprint, component_ids=component_ids, inlet_ph=ph)
    elif initial_cell_states_mol_m3 is not None:
        raw_cell_states = {str(cell_id): state for cell_id, state in initial_cell_states_mol_m3.items()}
        expected_cells = {cell_id for _, _, cell_id in cell_metadata}
        supplied_cells = set(raw_cell_states)
        if supplied_cells != expected_cells:
            raise KeyError(f'initial cell states must contain every train cell exactly; missing={sorted(expected_cells - supplied_cells)}, extra={sorted(supplied_cells - expected_cells)}')
        cell_states = {cell_id: _complete_state(raw_cell_states[cell_id], component_ids=component_ids, inlet_ph=imposed.get(zone.zone_id, ph)) for zone, _, cell_id in cell_metadata}
    else:
        for zone, _, cell_id in cell_metadata:
            seed = _complete_state(initial_by_zone.get(zone.zone_id, feed), component_ids=component_ids, inlet_ph=imposed.get(zone.zone_id, ph))
            cell_states[cell_id] = dict(seed)
    initial_cell_states = {cell_id: dict(state) for cell_id, state in cell_states.items()}
    background = checkpoint_background if checkpoint_background is not None and fixed_background_charge_eq_m3 is None else infer_fixed_background_charge_eq_m3(feed, measured_ph=ph, temperature_k=temperature, hydrogen_activity_coefficient=float((activity_coefficients or {}).get('H_plus', 1.0))) if fixed_background_charge_eq_m3 is None else float(fixed_background_charge_eq_m3)
    if not math.isfinite(background):
        raise ValueError('fixed background charge must be finite')
    trajectory: list[dict[str, Any]] = []
    local_balance_passed = True
    charge_passed = True
    extent_limiter_count = 0
    forced_surface_release_count = 0
    reservoir_totals: dict[str, float] = {}
    wasting_totals: dict[str, float] = {}
    external_input_totals: dict[str, float] = {}
    external_output_totals: dict[str, float] = {}
    elapsed = 0.0
    step_index = 0
    prior_terminal = dict(cell_states[cell_metadata[-1][2]])
    terminal_change_rate = math.inf
    last_charge_by_cell: dict[str, dict[str, object]] = {}
    last_rate_by_cell: dict[str, dict[str, Any]] = {}
    while elapsed < horizon - 1e-15:
        dt = min(dt_requested, horizon - elapsed)
        terminal_before_step = dict(cell_states[cell_metadata[-1][2]])
        flow = train.hydraulic_flow_m3_d
        for component_id in component_ids:
            if component_id == 'H2O':
                continue
            external_input_totals[component_id] = external_input_totals.get(component_id, 0.0) + flow * float(feed.get(component_id, 0.0)) * dt
        for component_id, mol_d in dose.items():
            external_input_totals[component_id] = external_input_totals.get(component_id, 0.0) + mol_d * dt
        external_input_totals['O2_aq'] = external_input_totals.get('O2_aq', 0.0) + sum((zone.oxygen_source_mol_m3_d * zone.volume_m3 for zone in train.zones)) * dt
        upstream = feed
        for zone, cell_index, cell_id in cell_metadata:
            cell_volume = zone.volume_m3 / zone.tanks_in_series
            cell_hrt = zone.hydraulic_residence_time_d / zone.tanks_in_series
            source: dict[str, float] = {}
            if zone.zone_kind == 'injection_micromixing' and cell_index == 0:
                source.update({component_id: mol_d * zone.dose_fraction / cell_volume for component_id, mol_d in dose.items()})
            if zone.oxygen_source_mol_m3_d:
                source['O2_aq'] = source.get('O2_aq', 0.0) + zone.oxygen_source_mol_m3_d
            exchange = {component_id: 1.0 for component_id in REDUCED_AQUEOUS_COMPONENT_IDS}
            solid_exchange = 1.0 - zone.particulate_rejection_fraction if zone.zone_kind == 'mbr_separation' else 1.0
            exchange.update({component_id: solid_exchange for component_id in retained_ids})
            options: dict[str, Any] = {'temperature_k': temperature, 'solvent_water_activity': water_activity, 'phase_parameters': phase_parameters, 'redox_parameters': redox_parameters, 'oxidation_parameters': oxidation_parameters, 'aging_parameters': aging_parameters, 'surface_parameters': surface_parameters, 'ligand_parameters': ligand_parameters, 'activity_coefficients': activity_coefficients, 'zone_volume': cell_volume}
            options.update(zone_overrides.get(zone.zone_id, {}))
            step, rate_report = advance_reduced_chemistry_cstr(cell_states[cell_id], inlet_mol_m3=upstream, dose_source_mol_m3_d=source, dt_days=dt, hydraulic_residence_time_days=cell_hrt, solids_retention_time_days=zone.solids_retention_time_d, kernel=kernel, hydraulic_exchange_fractions=exchange, hydraulic_inlet_fractions={component_id: 1.0 for component_id in retained_ids}, hydraulic_outlet_fractions={component_id: solid_exchange for component_id in retained_ids}, **options)
            updated = dict(step.state)
            if step.reaction_step.common_scale < 1.0 - 1e-12:
                extent_limiter_count += 1
            if step.balance_report is not None:
                local_balance_passed = local_balance_passed and bool(step.balance_report['passed'])
            for component_id, amount in step.reaction_step.reservoir_exchanges.items():
                reservoir_totals[component_id] = reservoir_totals.get(component_id, 0.0) + amount * cell_volume
            for component_id, amount in step.wasting_removed.items():
                wasting_totals[component_id] = wasting_totals.get(component_id, 0.0) + amount * cell_volume
            active_surface_parameters = options.get('surface_parameters') or {}
            release_extents: dict[str, float] = {}
            for history, host_id, surface_p_id in (('fresh', 'HFO', 'HFO_fresh_surface_P'), ('reoxidized', 'HFO_reoxidized', 'HFO_reoxidized_surface_P'), ('aged', 'HFO_aged', 'HFO_aged_surface_P')):
                if history not in active_surface_parameters:
                    continue
                capacity = _finite_nonnegative(active_surface_parameters[history].get('site_capacity', 0.0), field=f'{zone.zone_id}.surface.{history}.site_capacity')
                excess = max(float(updated.get(surface_p_id, 0.0)) - capacity * float(updated.get(host_id, 0.0)), 0.0)
                if excess > 1e-14:
                    release_extents[f'release_p_from_{history}_hfo'] = excess
            if release_extents:
                released = constrained_reaction_step(updated, kernel.reactions, release_extents, implicit_reservoir_component_ids=('H_plus',))
                updated = dict(released.state)
                forced_surface_release_count += len(release_extents)
                for component_id, amount in released.reservoir_exchanges.items():
                    reservoir_totals[component_id] = reservoir_totals.get(component_id, 0.0) + amount * cell_volume
            if zone.local_ph_model == 'relative_charge_balance':
                charge = solve_relative_charge_balance(updated, fixed_background_charge_eq_m3=background, temperature_k=temperature, hydrogen_activity_coefficient=float((activity_coefficients or {}).get('H_plus', 1.0)))
                last_charge_by_cell[cell_id] = charge.public_dict()
                if not charge.converged or charge.hydrogen_concentration_mol_m3 is None:
                    charge_passed = False
                    cell_states[cell_id] = updated
                    break
                updated['H_plus'] = charge.hydrogen_concentration_mol_m3
            else:
                imposed_zone_ph = imposed[zone.zone_id]
                h_activity = 10.0 ** (-imposed_zone_ph)
                updated['H_plus'] = h_activity / float((activity_coefficients or {}).get('H_plus', 1.0)) * 1000.0
                residual = reduced_charge_density_eq_m3(updated, hydrogen_activity_mol_l=h_activity, temperature_k=temperature, fixed_background_charge_eq_m3=background, hydrogen_activity_coefficient=float((activity_coefficients or {}).get('H_plus', 1.0)))
                last_charge_by_cell[cell_id] = {'status': f'{zone.local_ph_model}_imposed', 'converged': True, 'bracketed': True, 'ph': imposed_zone_ph, 'charge_residual_eq_m3': residual, 'external_titrant_or_property_package_reservoir_required': True}
            ph_projection_exchange = updated['H_plus'] - float(step.state.get('H_plus', 0.0))
            reservoir_totals['H_plus'] = reservoir_totals.get('H_plus', 0.0) + ph_projection_exchange * cell_volume
            cell_states[cell_id] = updated
            upstream = updated
            last_rate_by_cell[cell_id] = rate_report
        if not charge_passed:
            break
        elapsed += dt
        step_index += 1
        terminal = cell_states[cell_metadata[-1][2]]
        for component_id in component_ids:
            if component_id == 'H2O':
                continue
            leakage = 1.0 if component_id in REDUCED_AQUEOUS_COMPONENT_IDS else 1.0 - train.zones[-1].particulate_rejection_fraction
            average_concentration = 0.5 * (float(terminal_before_step.get(component_id, 0.0)) + float(terminal.get(component_id, 0.0)))
            external_output_totals[component_id] = external_output_totals.get(component_id, 0.0) + flow * leakage * average_concentration * dt
        terminal_change_rate = max((abs(float(terminal.get(component_id, 0.0)) - float(prior_terminal.get(component_id, 0.0))) / dt for component_id in component_ids))
        prior_terminal = dict(terminal)
        if step_index % int(sample_every_steps) == 0 or elapsed >= horizon - 1e-15:
            trajectory.append({'time_d': starting_checkpoint_time + elapsed, 'window_time_d': elapsed, 'terminal_po4_mol_m3': terminal['PO4'], 'terminal_reactive_p_mol_m3': _reactive_p_inventory(terminal), 'zone_ph': {cell_id: last_charge_by_cell[cell_id].get('ph') for _, _, cell_id in cell_metadata if cell_id in last_charge_by_cell}})
    terminal_cell_id = cell_metadata[-1][2]
    terminal = dict(cell_states[terminal_cell_id])
    mbr_zone = train.zones[-1]
    terminal_options: dict[str, Any] = {'temperature_k': temperature, 'solvent_water_activity': water_activity, 'phase_parameters': phase_parameters, 'redox_parameters': redox_parameters, 'oxidation_parameters': oxidation_parameters, 'aging_parameters': aging_parameters, 'surface_parameters': surface_parameters, 'ligand_parameters': ligand_parameters, 'activity_coefficients': activity_coefficients, 'zone_volume': mbr_zone.volume_m3 / mbr_zone.tanks_in_series}
    terminal_options.update(zone_overrides.get(mbr_zone.zone_id, {}))
    terminal_rate_report = evaluate_reduced_chemistry_rates(terminal, kernel=kernel, **terminal_options)
    rejection = mbr_zone.particulate_rejection_fraction
    permeate = {component_id: terminal.get(component_id, 0.0) if component_id in REDUCED_AQUEOUS_COMPONENT_IDS else (1.0 - rejection) * terminal.get(component_id, 0.0) for component_id in component_ids if component_id != 'H2O'}
    permeate_reactive_p = _reactive_p_inventory(permeate)
    balance = _whole_train_balance(initial_states=initial_cell_states, final_states=cell_states, cell_volumes_m3=cell_volumes, external_input_mol=external_input_totals, external_output_mol=external_output_totals, wasting_mol=wasting_totals, reservoir_exchange_mol=reservoir_totals, kernel=kernel)
    status = 'completed' if charge_passed else 'charge_balance_failed'
    blockers = ['uncommissioned kinetic parameters', 'serial mineral kernel does not by itself represent organic/biomass/polyphosphate TP', 'native mADM1 owns biological sulfur generation and finite-headspace gas transfer', 'plant calibration, independent holdout validation, and PHREEQC audit are required']
    if extent_limiter_count:
        blockers.append('reaction extent limiter activated; demonstrate timestep convergence')
    if not charge_passed:
        blockers.append('relative charge solve failed')
    if not local_balance_passed:
        blockers.append('local reaction balance failed')
    if not balance['passed']:
        blockers.append('whole-train transient balance failed; reduce the timestep and demonstrate convergence')
    final_inventory_id = _state_hash(cell_states)
    checkpoint_time = starting_checkpoint_time + elapsed
    serial_checkpoint = {**fingerprint, 'time_days': checkpoint_time, 'cell_states_mol_m3': {cell_id: dict(cell_states[cell_id]) for _, _, cell_id in cell_metadata}, 'fixed_background_charge_eq_m3': background, 'final_inventory_id': final_inventory_id, 'persistable': bool(status == 'completed' and charge_passed and local_balance_passed and balance['passed'])}
    terminal_activities = dict(terminal_rate_report['activities_mol_l'])
    terminal_phase_si = dict(terminal_rate_report['saturation_index'])
    terminal_temperature_k = float(terminal_options['temperature_k'])
    terminal_water_activity = float(terminal_options.get('solvent_water_activity', water_activity))
    terminal_ph = -math.log10(terminal_activities['H_plus'])
    mechanism_interpretation = mechanism_interpretation_from_ledger({'surface_dataset_id': kernel.surfaces[0].surface_dataset_id if kernel.surfaces else None, 'surface_parameter_interpretation': 'history_specific_uncommissioned_serial_profiles', 'oxide_history_class': 'serial_mixed_history', 'oxide_history_id': 'fresh_fes_reoxidized_and_aged_inventories_explicit_uncommissioned', 'oxide_formation_ph_basis': 'zone_charge_or_imposed_pH', 'active_metal_fraction_basis': 'not_identified', 'effective_active_fe_fraction': None, 'parameter_contact': None, 'parameter_contact_evidence_id': None, 'mechanism_identified': False, 'mechanism_evidence_id': None, 'extrapolation_status': 'not_assessed', 'extrapolation_reasons': ()})
    endpoint_translation = {'kernel_id': kernel.kernel_id, 'kernel_hash': kernel.sha256, 'totals_mol_l': {'P': terminal['PO4'] * 0.001, 'S(6)': terminal['SO4'] * 0.001, 'S(-2)': terminal['HS_minus'] * 0.001, 'Fe(3)': terminal['Fe_III_aq'] * 0.001, 'Fe(2)': terminal['Fe_II_aq'] * 0.001, 'C(4)': terminal['CO3'] * 0.001, 'Cl': terminal['Cl_minus'] * 0.001}, 'kinetic_phase_inventory_mol_l': {phase.phase_id: terminal.get(phase.phase_id, 0.0) * 0.001 for phase in kernel.phases if phase.phase_id != 'HFO'}, 'hfo_history_inventory_mol_fe_l': {'fresh': terminal.get('HFO', 0.0) * 0.001, 'fes_reoxidized': terminal.get('HFO_reoxidized', 0.0) * 0.001, 'aged': terminal.get('HFO_aged', 0.0) * 0.001}, 'p_total_basis': 'serial_reduced_PO4_only_not_plant_TP', 'temperature_c': terminal_temperature_k - 273.15, 'ph': terminal_ph, 'gas_boundary_mode': 'closed', 'gas_boundary': {}, 'analytical_totals_ready': True, 'translation_ready': True, 'audit_ready': True, 'audit_readiness_blockers': [], 'mechanism_interpretation': mechanism_interpretation, 'mapping_limitations': ['serial PO4 is modeled dissolved reactive P, not plant TP', 'conditional ligand chemistry is reduced rather than molecule-specific', 'native finite-headspace transfer precedes this closed-liquid endpoint audit', 'kinetic retained phases are not equilibrium-phase assertions']}
    solver_converged = bool(charge_passed and elapsed >= horizon - 1e-15)
    mass_balance_passed = local_balance_passed and bool(balance['passed']) and solver_converged
    audit_request = build_endpoint_thermodynamic_audit_request(totals_mol_l=endpoint_translation['totals_mol_l'], temperature_c=endpoint_translation['temperature_c'], ph=terminal_ph, phase_saturation_log10_si=terminal_phase_si, activities=terminal_activities, solvent_water_activity=terminal_water_activity, source_mass_balance_passed=mass_balance_passed, species_coverage_fraction=coverage, mechanism_interpretation=mechanism_interpretation, kernel=kernel)
    return {'status': status, 'solver_converged': solver_converged, 'chemistry_output_available': bool(solver_converged and mass_balance_passed), 'kernel_id': kernel.kernel_id, 'kernel_hash': kernel.sha256, 'train_id': train.train_id, 'hydraulic_flow_m3_d': train.hydraulic_flow_m3_d, 'elapsed_days': elapsed, 'checkpoint_time_days': checkpoint_time, 'resumed_from_checkpoint': resume_checkpoint is not None, 'starting_checkpoint_time_days': starting_checkpoint_time, 'requested_horizon_days': horizon, 'time_step_days': dt_requested, 'cell_count': len(cell_metadata), 'zone_states_mol_m3': {cell_id: dict(cell_states[cell_id]) for _, _, cell_id in cell_metadata}, 'terminal_state_mol_m3': terminal, 'mbr_permeate_state_mol_m3': permeate, 'p_metrics': {'terminal_orthophosphate_mol_p_m3': terminal.get('PO4', 0.0), 'terminal_reactive_p_mol_p_m3': _reactive_p_inventory(terminal), 'permeate_reactive_p_mol_p_m3': permeate_reactive_p, 'permeate_orthophosphate_mg_p_l': terminal.get('PO4', 0.0) * 30.973761998, 'permeate_reactive_p_mg_p_l': permeate_reactive_p * 30.973761998, 'measurement_basis': 'modeled_inorganic_and_surface_reactive_P_not_plant_TP'}, 'oxide_history_inventory_mol_fe_m3': {'fresh': terminal.get('HFO', 0.0), 'fes_reoxidized': terminal.get('HFO_reoxidized', 0.0), 'aged': terminal.get('HFO_aged', 0.0)}, 'surface_p_inventory_mol_p_m3': {'fresh': terminal.get('HFO_fresh_surface_P', 0.0), 'fes_reoxidized': terminal.get('HFO_reoxidized_surface_P', 0.0), 'aged': terminal.get('HFO_aged_surface_P', 0.0)}, 'iron_sulfide_inventory_mol_m3': {'amorphous': terminal.get('FeS_amorphous', 0.0), 'mackinawite': terminal.get('Mackinawite', 0.0)}, 'fixed_background_charge_eq_m3': background, 'charge_diagnostics': last_charge_by_cell, 'ph_charge_solve_passed': charge_passed, 'local_reaction_balance_passed': local_balance_passed, 'whole_train_balance': balance, 'mass_balance_passed': mass_balance_passed, 'implicit_reservoir_exchange_mol': reservoir_totals, 'cumulative_wasting_mol': wasting_totals, 'extent_limiter_activation_count': extent_limiter_count, 'forced_host_capacity_release_count': forced_surface_release_count, 'terminal_change_rate_mol_m3_d': terminal_change_rate, 'trajectory': trajectory, 'last_rate_reports': last_rate_by_cell, 'initial_inventory_id': _state_hash(initial_cell_states), 'final_inventory_id': final_inventory_id, 'serial_checkpoint': serial_checkpoint, 'outlet_activities': terminal_activities, 'phase_saturation_index': terminal_phase_si, 'phreeqc_endpoint_translation': endpoint_translation, 'phreeqc_audit_request': audit_request, 'phreeqc_audit_dispatch': {'tool': 'water-chemistry-mcp.audit_process_chemistry_endpoint', 'request_ready': True, 'executed_by_qsdsan': False, 'request_kernel_hash': kernel.sha256}, 'mechanism_interpretation': mechanism_interpretation, 'gas_coupling_class': 'delegated_to_native_madm1_explicit_boundary', 'mechanism_identified': False, 'recommendation_eligible': False, 'recommendation_blockers': blockers}
__all__ = ['SERIAL_CHEMISTRY_CHECKPOINT_SCHEMA', 'serial_chemistry_train_fingerprint', 'simulate_serial_chemistry_train']
