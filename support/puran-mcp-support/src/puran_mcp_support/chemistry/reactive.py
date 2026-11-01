"""Numeric execution of the shared reduced Fe-S-P chemistry kernel.

This module is the dependency-light dynamic counterpart to WaterTAP's Pyomo
compiler. It evaluates the same serialized phase, reaction, speciation, and
kinetic definitions and then advances a positivity-preserving CSTR/MBR
inventory step. QSDsan adapters can call it without importing Pyomo.
"""
from __future__ import annotations
import math
from collections.abc import Mapping
from typing import Any
from puran_mcp_support.schema import MADM1_FECL3_FE_S_P_KERNEL, ChemistryKernelDefinition, reduced_free_ion_fractions
from .evaluator import evaluate
from .inventory import InventoryStepResult, advance_cstr_inventory
REDUCED_AQUEOUS_COMPONENT_IDS = ('Fe_III_aq', 'Fe_II_aq', 'H2_aq', 'O2_aq', 'HS_minus', 'SO4', 'PO4', 'CO3', 'Cl_minus', 'Ligand_sites', 'H_plus')
_PHASE_PARAMETER_DEFAULTS = {'k_n': 0.0, 'k_g': 0.0, 'k_d': 0.0, 'si_critical': 0.0, 'n_nucleation': 1.0, 'n_growth': 1.0, 'n_dissolution': 1.0, 'area_per_inventory': 1.0, 'seed_surface_area': 0.0}
_REDOX_PARAMETER_DEFAULTS = {'k_red': 0.0, 'K_donor': 1e-06}
_OXIDATION_PARAMETER_DEFAULTS = {'k_oxidation': 0.0, 'K_oxidant': 1e-06}
_AGING_PARAMETER_DEFAULTS = {'k_aging': 0.0}
_SURFACE_PARAMETER_DEFAULTS = {'site_capacity': 0.0, 'k_adsorption': 0.0, 'K_P': 1e-06, 'K_ligand': 1.0, 'k_desorption': 0.0}
_REDOX_ACCEPTOR = {'ferric_reduction_by_hydrogen': 'Fe_III_aq', 'ferric_reduction_by_sulfide': 'Fe_III_aq', 'hfo_reduction_by_hydrogen': 'HFO', 'hfo_reduction_by_sulfide': 'HFO', 'reoxidized_hfo_reduction_by_hydrogen': 'HFO_reoxidized', 'reoxidized_hfo_reduction_by_sulfide': 'HFO_reoxidized', 'aged_hfo_reduction_by_hydrogen': 'HFO_aged', 'aged_hfo_reduction_by_sulfide': 'HFO_aged'}
_REDOX_DONOR = {reaction_id: 'H2_aq' if reaction_id.endswith('hydrogen') else 'HS_minus' for reaction_id in _REDOX_ACCEPTOR}
_OXIDATION_SOLID = {'oxidize_amorphous_fes_to_reoxidized_hfo_s0': 'FeS_amorphous', 'oxidize_amorphous_fes_to_reoxidized_hfo_sulfate': 'FeS_amorphous', 'oxidize_mackinawite_to_reoxidized_hfo_s0': 'Mackinawite', 'oxidize_mackinawite_to_reoxidized_hfo_sulfate': 'Mackinawite'}
_AGING_SOURCE = {'age_amorphous_fes_to_mackinawite': 'FeS_amorphous', 'age_fresh_hfo': 'HFO', 'age_reoxidized_hfo': 'HFO_reoxidized', 'age_fresh_hfo_surface_p': 'HFO_fresh_surface_P', 'age_reoxidized_hfo_surface_p': 'HFO_reoxidized_surface_P'}
_SURFACE_HISTORY = {'fresh': ('HFO', 'HFO_fresh_surface_P'), 'reoxidized': ('HFO_reoxidized', 'HFO_reoxidized_surface_P'), 'aged': ('HFO_aged', 'HFO_aged_surface_P')}

def _finite_nonnegative(value: Any, *, field: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0.0:
        raise ValueError(f'{field} must be finite and nonnegative')
    return number

def _normalized_profiles(kernel: ChemistryKernelDefinition, *, phase_parameters: Mapping[str, Mapping[str, float]] | None, redox_parameters: Mapping[str, Mapping[str, float]] | None, oxidation_parameters: Mapping[str, Mapping[str, float]] | None, aging_parameters: Mapping[str, Mapping[str, float]] | None, surface_parameters: Mapping[str, Mapping[str, float]] | None) -> tuple[dict[str, dict[str, float]], dict[str, dict[str, float]], dict[str, dict[str, float]], dict[str, dict[str, float]], dict[str, dict[str, float]]]:
    phase_ids = {phase.phase_id for phase in kernel.phases}
    redox_ids = set(_REDOX_ACCEPTOR) & {reaction.reaction_id for reaction in kernel.reactions}
    oxidation_ids = set(_OXIDATION_SOLID) & {reaction.reaction_id for reaction in kernel.reactions}
    aging_ids = set(_AGING_SOURCE) & {reaction.reaction_id for reaction in kernel.reactions}
    surface_ids = {history for history in _SURFACE_HISTORY if f'adsorb_p_on_{history}_hfo' in {reaction.reaction_id for reaction in kernel.reactions}}
    raw_phase = dict(phase_parameters or {})
    raw_redox = dict(redox_parameters or {})
    raw_oxidation = dict(oxidation_parameters or {})
    raw_aging = dict(aging_parameters or {})
    raw_surface = dict(surface_parameters or {})
    unknown_phases = sorted(set(raw_phase) - phase_ids)
    unknown_redox = sorted(set(raw_redox) - redox_ids)
    unknown_oxidation = sorted(set(raw_oxidation) - oxidation_ids)
    unknown_aging = sorted(set(raw_aging) - aging_ids)
    unknown_surface = sorted(set(raw_surface) - surface_ids)
    if unknown_phases or unknown_redox or unknown_oxidation or unknown_aging or unknown_surface:
        raise KeyError(f'unknown kinetic identifiers: phases={unknown_phases}, redox={unknown_redox}, oxidation={unknown_oxidation}, aging={unknown_aging}, surface={unknown_surface}')
    phases: dict[str, dict[str, float]] = {}
    for phase_id in sorted(phase_ids):
        supplied = dict(raw_phase.get(phase_id, {}))
        unknown = sorted(set(supplied) - set(_PHASE_PARAMETER_DEFAULTS))
        if unknown:
            raise KeyError(f'phase {phase_id!r} has unknown kinetic fields {unknown}')
        phases[phase_id] = {name: _finite_nonnegative(supplied.get(name, default), field=f'{phase_id}.{name}') for name, default in _PHASE_PARAMETER_DEFAULTS.items()}
    redox: dict[str, dict[str, float]] = {}
    for reaction_id in sorted(redox_ids):
        supplied = dict(raw_redox.get(reaction_id, {}))
        unknown = sorted(set(supplied) - set(_REDOX_PARAMETER_DEFAULTS))
        if unknown:
            raise KeyError(f'redox reaction {reaction_id!r} has unknown kinetic fields {unknown}')
        redox[reaction_id] = {name: _finite_nonnegative(supplied.get(name, default), field=f'{reaction_id}.{name}') for name, default in _REDOX_PARAMETER_DEFAULTS.items()}

    def normalized_group(identifiers: set[str], supplied_by_id: Mapping[str, Mapping[str, float]], defaults: Mapping[str, float], *, group: str) -> dict[str, dict[str, float]]:
        normalized: dict[str, dict[str, float]] = {}
        for identifier in sorted(identifiers):
            supplied = dict(supplied_by_id.get(identifier, {}))
            unknown = sorted(set(supplied) - set(defaults))
            if unknown:
                raise KeyError(f'{group} {identifier!r} has unknown fields {unknown}')
            normalized[identifier] = {name: _finite_nonnegative(supplied.get(name, default), field=f'{identifier}.{name}') for name, default in defaults.items()}
        return normalized
    oxidation = normalized_group(oxidation_ids, raw_oxidation, _OXIDATION_PARAMETER_DEFAULTS, group='oxidation reaction')
    aging = normalized_group(aging_ids, raw_aging, _AGING_PARAMETER_DEFAULTS, group='aging reaction')
    surfaces = normalized_group(surface_ids, raw_surface, _SURFACE_PARAMETER_DEFAULTS, group='surface history')
    if any((profile['K_oxidant'] <= 0.0 for profile in oxidation.values())):
        raise ValueError('oxidation half-saturation constants must be positive')
    if any((profile['K_P'] <= 0.0 or profile['K_ligand'] <= 0.0 for profile in surfaces.values())):
        raise ValueError('surface half-saturation/inhibition constants must be positive')
    return (phases, redox, oxidation, aging, surfaces)

def conditional_ligand_free_fractions(*, ligand_total_mol_l: float, fe_ii_total_mol_l: float, fe_iii_total_mol_l: float, beta_fe_ii_l_mol: float, beta_fe_iii_l_mol: float) -> dict[str, float]:
    """Solve a finite shared 1:1 conditional-ligand balance.

    The analytical Fe totals are unchanged.  Only their free fractions are
    returned, so this reduced model cannot create or destroy Fe while still
    exposing wastewater-organic complexation uncertainty to mineral IAPs.
    """
    values = {'ligand_total_mol_l': ligand_total_mol_l, 'fe_ii_total_mol_l': fe_ii_total_mol_l, 'fe_iii_total_mol_l': fe_iii_total_mol_l, 'beta_fe_ii_l_mol': beta_fe_ii_l_mol, 'beta_fe_iii_l_mol': beta_fe_iii_l_mol}
    converted = {key: _finite_nonnegative(value, field=key) for key, value in values.items()}
    ligand_total = converted['ligand_total_mol_l']
    if ligand_total == 0.0:
        return {'free_ligand_mol_l': 0.0, 'Fe_II_aq': 1.0, 'Fe_III_aq': 1.0, 'Fe_II_complexed_mol_l': 0.0, 'Fe_III_complexed_mol_l': 0.0}
    fe_ii = converted['fe_ii_total_mol_l']
    fe_iii = converted['fe_iii_total_mol_l']
    beta_ii = converted['beta_fe_ii_l_mol']
    beta_iii = converted['beta_fe_iii_l_mol']

    def residual(free_ligand: float) -> float:
        complex_ii = fe_ii * beta_ii * free_ligand / (1.0 + beta_ii * free_ligand)
        complex_iii = fe_iii * beta_iii * free_ligand / (1.0 + beta_iii * free_ligand)
        return free_ligand + complex_ii + complex_iii - ligand_total
    lower = 0.0
    upper = ligand_total
    for _ in range(100):
        midpoint = 0.5 * (lower + upper)
        if residual(midpoint) > 0.0:
            upper = midpoint
        else:
            lower = midpoint
    free_ligand = 0.5 * (lower + upper)
    fraction_ii = 1.0 / (1.0 + beta_ii * free_ligand)
    fraction_iii = 1.0 / (1.0 + beta_iii * free_ligand)
    return {'free_ligand_mol_l': free_ligand, 'Fe_II_aq': fraction_ii, 'Fe_III_aq': fraction_iii, 'Fe_II_complexed_mol_l': fe_ii * (1.0 - fraction_ii), 'Fe_III_complexed_mol_l': fe_iii * (1.0 - fraction_iii)}

def evaluate_reduced_chemistry_rates(state_mol_m3: Mapping[str, float], *, kernel: ChemistryKernelDefinition=MADM1_FECL3_FE_S_P_KERNEL, temperature_k: float=298.15, phase_parameters: Mapping[str, Mapping[str, float]] | None=None, redox_parameters: Mapping[str, Mapping[str, float]] | None=None, oxidation_parameters: Mapping[str, Mapping[str, float]] | None=None, aging_parameters: Mapping[str, Mapping[str, float]] | None=None, surface_parameters: Mapping[str, Mapping[str, float]] | None=None, ligand_parameters: Mapping[str, float] | None=None, activity_coefficients: Mapping[str, float] | None=None, solvent_water_activity: float=1.0, zone_volume: float=1.0) -> dict[str, Any]:
    """Evaluate shared activities, phase SI, and simultaneous reaction rates."""
    temperature = float(temperature_k)
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise ValueError('temperature_k must be finite and positive')
    water_activity = float(solvent_water_activity)
    if not math.isfinite(water_activity) or water_activity <= 0.0:
        raise ValueError('solvent_water_activity must be finite and positive')
    volume = float(zone_volume)
    if not math.isfinite(volume) or volume <= 0.0:
        raise ValueError('zone_volume must be finite and positive')
    values = {str(key): float(value) for key, value in state_mol_m3.items()}
    invalid = sorted((component_id for component_id, value in values.items() if not math.isfinite(value) or value < 0.0))
    if invalid:
        raise ValueError(f'state must be finite and nonnegative; invalid={invalid}')
    missing_aqueous = sorted(set(REDUCED_AQUEOUS_COMPONENT_IDS) - set(values))
    if missing_aqueous:
        raise KeyError(f'state omitted reduced aqueous components {missing_aqueous}')
    phase_profile, redox_profile, oxidation_profile, aging_profile, surface_profile = _normalized_profiles(kernel, phase_parameters=phase_parameters, redox_parameters=redox_parameters, oxidation_parameters=oxidation_parameters, aging_parameters=aging_parameters, surface_parameters=surface_parameters)
    gammas = {component_id: 1.0 for component_id in REDUCED_AQUEOUS_COMPONENT_IDS}
    for component_id, raw in (activity_coefficients or {}).items():
        if component_id not in gammas:
            raise KeyError(f'unknown activity coefficient component {component_id!r}')
        gamma = float(raw)
        if not math.isfinite(gamma) or gamma <= 0.0:
            raise ValueError(f'activity coefficient for {component_id!r} must be finite and positive')
        gammas[component_id] = gamma
    hydrogen_activity = max(values['H_plus'] * 0.001 * gammas['H_plus'], 1e-30)
    free_fractions = reduced_free_ion_fractions(hydrogen_activity=hydrogen_activity, temperature_k=temperature)
    ligand_definition = kernel.ligand_models[0] if kernel.ligand_models else None
    ligand_raw = dict(ligand_parameters or {})
    allowed_ligand_fields = {'beta_fe_ii_l_mol', 'beta_fe_iii_l_mol'}
    unknown_ligand = sorted(set(ligand_raw) - allowed_ligand_fields)
    if unknown_ligand:
        raise KeyError(f'unknown conditional ligand fields {unknown_ligand}')
    beta_fe_ii = _finite_nonnegative(ligand_raw.get('beta_fe_ii_l_mol', 0.0 if ligand_definition is None else ligand_definition.beta_fe_ii_l_mol), field='ligand.beta_fe_ii_l_mol')
    beta_fe_iii = _finite_nonnegative(ligand_raw.get('beta_fe_iii_l_mol', 0.0 if ligand_definition is None else ligand_definition.beta_fe_iii_l_mol), field='ligand.beta_fe_iii_l_mol')
    ligand_speciation = conditional_ligand_free_fractions(ligand_total_mol_l=values['Ligand_sites'] * 0.001, fe_ii_total_mol_l=values['Fe_II_aq'] * 0.001, fe_iii_total_mol_l=values['Fe_III_aq'] * 0.001, beta_fe_ii_l_mol=beta_fe_ii, beta_fe_iii_l_mol=beta_fe_iii)
    free_fractions.update({'Fe_II_aq': ligand_speciation['Fe_II_aq'], 'Fe_III_aq': ligand_speciation['Fe_III_aq']})
    activities = {component_id: max(values[component_id] * 0.001 * gammas[component_id] * free_fractions.get(component_id, 1.0), 1e-30) for component_id in REDUCED_AQUEOUS_COMPONENT_IDS}
    saturation_index: dict[str, float] = {}
    saturation_ratio: dict[str, float] = {}
    symbols: dict[str, float] = {'zone_volume': volume, **{f'activity.{component_id}': activity for component_id, activity in activities.items()}}
    for phase in kernel.phases:
        log_iap = 0.0
        for component_id, coefficient in phase.dissolution_reaction.items():
            if component_id == phase.phase_id:
                continue
            activity = water_activity if component_id == 'H2O' else activities[component_id]
            log_iap += float(coefficient) * math.log(activity)
        equilibrium_constant = phase.equilibrium_constant(temperature)
        si = (log_iap - math.log(equilibrium_constant)) / math.log(10.0)
        ratio = math.exp(math.log(10.0) * si)
        saturation_index[phase.phase_id] = si
        saturation_ratio[phase.phase_id] = ratio
        profile = phase_profile[phase.phase_id]
        area = profile['area_per_inventory'] * values.get(phase.phase_id, 0.0) + profile['seed_surface_area']
        for name, value in profile.items():
            symbols[f'{name}.{phase.phase_id}'] = value
        symbols[f'SI.{phase.phase_id}'] = si
        symbols[f'SI_critical.{phase.phase_id}'] = profile['si_critical']
        symbols[f'area.{phase.phase_id}'] = area
        symbols[f'saturation_ratio.{phase.phase_id}'] = ratio
    for reaction_id, profile in redox_profile.items():
        acceptor_id = _REDOX_ACCEPTOR[reaction_id]
        donor_id = _REDOX_DONOR[reaction_id]
        symbols[f'k_red.{reaction_id}'] = profile['k_red']
        symbols[f'K_donor.{reaction_id}'] = profile['K_donor']
        symbols[f'inventory.{acceptor_id}'] = values.get(acceptor_id, 0.0)
        symbols[f'activity.{donor_id}'] = activities[donor_id]
    for reaction_id, profile in oxidation_profile.items():
        solid_id = _OXIDATION_SOLID[reaction_id]
        symbols[f'k_oxidation.{reaction_id}'] = profile['k_oxidation']
        symbols[f'K_oxidant.{reaction_id}'] = profile['K_oxidant']
        symbols[f'inventory.{solid_id}'] = values.get(solid_id, 0.0)
    for reaction_id, profile in aging_profile.items():
        source_id = _AGING_SOURCE[reaction_id]
        symbols[f'k_aging.{reaction_id}'] = profile['k_aging']
        symbols[f'inventory.{source_id}'] = values.get(source_id, 0.0)
    for history, profile in surface_profile.items():
        host_id, surface_p_id = _SURFACE_HISTORY[history]
        symbols[f'site_capacity.{history}'] = profile['site_capacity']
        symbols[f'k_surface_adsorption.{history}'] = profile['k_adsorption']
        symbols[f'K_surface_P.{history}'] = profile['K_P']
        symbols[f'K_surface_ligand.{history}'] = profile['K_ligand']
        symbols[f'k_surface_desorption.{history}'] = profile['k_desorption']
        symbols[f'inventory.{host_id}'] = values.get(host_id, 0.0)
        symbols[f'inventory.{surface_p_id}'] = values.get(surface_p_id, 0.0)
    reaction_rates = {reaction.reaction_id: 0.0 for reaction in kernel.reactions}
    kinetic_rates: dict[str, float] = {}
    phase_net_rates = {phase.phase_id: 0.0 for phase in kernel.phases}
    for kinetic in kernel.kinetics:
        rate = evaluate(kinetic.expression, symbols)
        kinetic_rates[kinetic.kinetic_id] = rate
        symbols[f'rate.{kinetic.kinetic_id}'] = rate
        if kinetic.reaction_id is not None:
            reaction_rates[kinetic.reaction_id] += rate
        if kinetic.rate_role.startswith('phase_'):
            phase_id = kinetic.kinetic_id.rsplit('.', 1)[0]
            phase_net_rates[phase_id] += rate
    return {'kernel_id': kernel.kernel_id, 'kernel_hash': kernel.sha256, 'activities_mol_l': activities, 'free_ion_fraction': free_fractions, 'conditional_ligand_speciation': {**ligand_speciation, 'beta_fe_ii_l_mol': beta_fe_ii, 'beta_fe_iii_l_mol': beta_fe_iii, 'ligand_model_id': None if ligand_definition is None else ligand_definition.ligand_model_id, 'calibrated': False if ligand_definition is None else ligand_definition.calibrated}, 'saturation_index': saturation_index, 'saturation_ratio': saturation_ratio, 'kinetic_rates': kinetic_rates, 'phase_net_rates_mol_m3_d': phase_net_rates, 'reaction_rates_mol_m3_d': reaction_rates, 'activity_floor_mol_l': 1e-30}

def advance_reduced_chemistry_cstr(state_mol_m3: Mapping[str, float], *, inlet_mol_m3: Mapping[str, float], dose_source_mol_m3_d: Mapping[str, float] | None=None, dt_days: float, hydraulic_residence_time_days: float, solids_retention_time_days: float, kernel: ChemistryKernelDefinition=MADM1_FECL3_FE_S_P_KERNEL, temperature_k: float=298.15, phase_parameters: Mapping[str, Mapping[str, float]] | None=None, redox_parameters: Mapping[str, Mapping[str, float]] | None=None, oxidation_parameters: Mapping[str, Mapping[str, float]] | None=None, aging_parameters: Mapping[str, Mapping[str, float]] | None=None, surface_parameters: Mapping[str, Mapping[str, float]] | None=None, ligand_parameters: Mapping[str, float] | None=None, activity_coefficients: Mapping[str, float] | None=None, solvent_water_activity: float=1.0, zone_volume: float=1.0, solvent_water_inventory_mol_m3: float=55500.0, hydraulic_exchange_fractions: Mapping[str, float] | None=None, hydraulic_inlet_fractions: Mapping[str, float] | None=None, hydraulic_outlet_fractions: Mapping[str, float] | None=None) -> tuple[InventoryStepResult, dict[str, Any]]:
    """Advance one dynamic step with rates compiled from the shared IR."""
    working_state = {str(key): float(value) for key, value in state_mol_m3.items()}
    if 'H2O' not in working_state:
        working_state['H2O'] = _finite_nonnegative(solvent_water_inventory_mol_m3, field='solvent_water_inventory_mol_m3')
    rate_report = evaluate_reduced_chemistry_rates(working_state, kernel=kernel, temperature_k=temperature_k, phase_parameters=phase_parameters, redox_parameters=redox_parameters, oxidation_parameters=oxidation_parameters, aging_parameters=aging_parameters, surface_parameters=surface_parameters, ligand_parameters=ligand_parameters, activity_coefficients=activity_coefficients, solvent_water_activity=solvent_water_activity, zone_volume=zone_volume)
    retained_ids = tuple((component.component_id for component in kernel.components if component.component_id not in REDUCED_AQUEOUS_COMPONENT_IDS and component.component_id != 'H2O' and (component.permeability_class.value == 'retained')))
    step = advance_cstr_inventory(working_state, reactions=kernel.reactions, rates_mol_per_volume_time=rate_report['reaction_rates_mol_m3_d'], inlet=inlet_mol_m3, dose_source_mol_per_volume_time=dose_source_mol_m3_d or {}, aqueous_component_ids=REDUCED_AQUEOUS_COMPONENT_IDS, retained_solid_ids=retained_ids, component_definitions=kernel.components, dt=dt_days, hydraulic_residence_time=hydraulic_residence_time_days, solids_retention_time=solids_retention_time_days, hydraulic_exchange_fractions=hydraulic_exchange_fractions, hydraulic_inlet_fractions=hydraulic_inlet_fractions, hydraulic_outlet_fractions=hydraulic_outlet_fractions, implicit_reservoir_component_ids=('H2O', 'H_plus'))
    return (step, rate_report)

def reduced_chemistry_cstr_rhs(state_mol_m3: Mapping[str, float], *, inlet_mol_m3: Mapping[str, float], dose_source_mol_m3_d: Mapping[str, float] | None=None, hydraulic_residence_time_days: float, solids_retention_time_days: float, kernel: ChemistryKernelDefinition=MADM1_FECL3_FE_S_P_KERNEL, **chemistry_options: Any) -> tuple[dict[str, float], dict[str, Any]]:
    """Return the zero-dimensional dynamic RHS rendered from the shared IR.

    This is the continuous-time form used for QSDsan/BDF integration and is
    algebraically identical to WaterTAP's zero-accumulation constraints for
    the same state, branch, and parameter profile. Solvent water remains an
    explicitly declared implicit reservoir, matching the reduced Pyomo block.
    """
    hrt = float(hydraulic_residence_time_days)
    srt = float(solids_retention_time_days)
    if not math.isfinite(hrt) or hrt <= 0.0 or (not math.isfinite(srt)) or (srt <= 0.0):
        raise ValueError('HRT and SRT must be finite and positive')
    state = {str(key): float(value) for key, value in state_mol_m3.items()}
    source = {str(key): float(value) for key, value in (dose_source_mol_m3_d or {}).items()}
    unknown_source = sorted(set(source) - set(REDUCED_AQUEOUS_COMPONENT_IDS))
    if unknown_source:
        raise KeyError(f'dose source contains unknown aqueous components {unknown_source}')
    if any((not math.isfinite(value) or value < 0.0 for value in source.values())):
        raise ValueError('dose sources must be finite and nonnegative')
    report = evaluate_reduced_chemistry_rates(state, kernel=kernel, **chemistry_options)
    reaction_by_id = {reaction.reaction_id: reaction for reaction in kernel.reactions}
    rates = report['reaction_rates_mol_m3_d']
    retained_ids = tuple((component.component_id for component in kernel.components if component.component_id not in REDUCED_AQUEOUS_COMPONENT_IDS and component.component_id != 'H2O' and (component.permeability_class.value == 'retained')))
    derivatives: dict[str, float] = {}
    for component_id in REDUCED_AQUEOUS_COMPONENT_IDS:
        reaction_term = sum((reaction_by_id[reaction_id].stoichiometry.get(component_id, 0.0) * rate for reaction_id, rate in rates.items()))
        derivatives[component_id] = (float(inlet_mol_m3.get(component_id, 0.0)) - state[component_id]) / hrt + source.get(component_id, 0.0) + reaction_term
    for component_id in retained_ids:
        reaction_term = sum((reaction_by_id[reaction_id].stoichiometry.get(component_id, 0.0) * rate for reaction_id, rate in rates.items()))
        derivatives[component_id] = reaction_term - state.get(component_id, 0.0) / srt
    return (derivatives, report)
__all__ = ['REDUCED_AQUEOUS_COMPONENT_IDS', 'advance_reduced_chemistry_cstr', 'conditional_ligand_free_fractions', 'evaluate_reduced_chemistry_rates', 'reduced_chemistry_cstr_rhs']
