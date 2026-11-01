"""Fail-closed relative charge closure for the reduced Fe-S-P state.

The reduced chemistry state carries analytical phosphate, sulfide, and
carbonate totals rather than every protonation state.  This module projects
those totals onto their average aqueous charge at a trial hydrogen activity.
An unmodelled strong-ion background is inferred once from a measured inlet pH
and then transported unchanged through the serial chemistry train.  It is a
declared boundary term, not a hidden neutral-pH fallback.
"""
from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Mapping
from puran_mcp_support.schema import ACID_DISSOCIATION_CONSTANTS

@dataclass(frozen=True)
class ChargeBalanceResult:
    status: str
    converged: bool
    bracketed: bool
    ph: float | None
    hydrogen_activity_mol_l: float | None
    hydrogen_concentration_mol_m3: float | None
    charge_residual_eq_m3: float | None
    fixed_background_charge_eq_m3: float
    iterations: int
    lower_ph: float
    upper_ph: float

    def public_dict(self) -> dict[str, object]:
        return {'status': self.status, 'converged': self.converged, 'bracketed': self.bracketed, 'ph': self.ph, 'hydrogen_activity_mol_l': self.hydrogen_activity_mol_l, 'hydrogen_concentration_mol_m3': self.hydrogen_concentration_mol_m3, 'charge_residual_eq_m3': self.charge_residual_eq_m3, 'fixed_background_charge_eq_m3': self.fixed_background_charge_eq_m3, 'iterations': self.iterations, 'lower_ph': self.lower_ph, 'upper_ph': self.upper_ph}

def _finite_nonnegative(value: float, *, field: str) -> float:
    converted = float(value)
    if not math.isfinite(converted) or converted < 0.0:
        raise ValueError(f'{field} must be finite and nonnegative')
    return converted

def reduced_acid_base_average_charges(*, hydrogen_activity_mol_l: float, temperature_k: float) -> dict[str, float]:
    """Return charge per mole of each transported analytical total."""
    h = float(hydrogen_activity_mol_l)
    temperature = float(temperature_k)
    if not math.isfinite(h) or h <= 0.0:
        raise ValueError('hydrogen activity must be finite and positive')
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise ValueError('temperature_k must be finite and positive')
    constants = ACID_DISSOCIATION_CONSTANTS
    kp1 = constants['phosphoric_1'].ka(temperature)
    kp2 = constants['phosphoric_2'].ka(temperature)
    kp3 = constants['phosphoric_3'].ka(temperature)
    p_denominator = h ** 3 + kp1 * h ** 2 + kp1 * kp2 * h + kp1 * kp2 * kp3
    phosphate_charge = -(kp1 * h ** 2 + 2.0 * kp1 * kp2 * h + 3.0 * kp1 * kp2 * kp3) / p_denominator
    ks1 = constants['sulfide_1'].ka(temperature)
    ks2 = constants['sulfide_2'].ka(temperature)
    s_denominator = h ** 2 + ks1 * h + ks1 * ks2
    sulfide_charge = -(ks1 * h + 2.0 * ks1 * ks2) / s_denominator
    kc1 = constants['carbonic_1'].ka(temperature)
    kc2 = constants['carbonic_2'].ka(temperature)
    c_denominator = h ** 2 + kc1 * h + kc1 * kc2
    carbonate_charge = -(kc1 * h + 2.0 * kc1 * kc2) / c_denominator
    return {'PO4': phosphate_charge, 'HS_minus': sulfide_charge, 'CO3': carbonate_charge}

def reduced_charge_density_eq_m3(state_mol_m3: Mapping[str, float], *, hydrogen_activity_mol_l: float, temperature_k: float, fixed_background_charge_eq_m3: float=0.0, hydrogen_activity_coefficient: float=1.0, hydroxide_activity_coefficient: float=1.0) -> float:
    """Evaluate bulk liquid plus declared surface charge in equivalents/m3."""
    state = {str(key): _finite_nonnegative(value, field=str(key)) for key, value in state_mol_m3.items()}
    h_activity = float(hydrogen_activity_mol_l)
    background = float(fixed_background_charge_eq_m3)
    gamma_h = float(hydrogen_activity_coefficient)
    gamma_oh = float(hydroxide_activity_coefficient)
    if not math.isfinite(background):
        raise ValueError('fixed background charge must be finite')
    if not math.isfinite(gamma_h) or gamma_h <= 0.0:
        raise ValueError('hydrogen activity coefficient must be finite and positive')
    if not math.isfinite(gamma_oh) or gamma_oh <= 0.0:
        raise ValueError('hydroxide activity coefficient must be finite and positive')
    average = reduced_acid_base_average_charges(hydrogen_activity_mol_l=h_activity, temperature_k=temperature_k)
    charge = background
    charge += 3.0 * state.get('Fe_III_aq', 0.0)
    charge += 2.0 * state.get('Fe_II_aq', 0.0)
    charge -= 2.0 * state.get('SO4', 0.0)
    charge -= state.get('Cl_minus', 0.0)
    charge += average['PO4'] * state.get('PO4', 0.0)
    charge += average['HS_minus'] * state.get('HS_minus', 0.0)
    charge += average['CO3'] * state.get('CO3', 0.0)
    charge -= state.get('HFO_fresh_surface_P', 0.0)
    charge -= state.get('HFO_reoxidized_surface_P', 0.0)
    charge -= state.get('HFO_aged_surface_P', 0.0)
    h_concentration_mol_m3 = h_activity / gamma_h * 1000.0
    kw = ACID_DISSOCIATION_CONSTANTS['water'].ka(float(temperature_k))
    oh_concentration_mol_m3 = kw / h_activity / gamma_oh * 1000.0
    return charge + h_concentration_mol_m3 - oh_concentration_mol_m3

def infer_fixed_background_charge_eq_m3(state_mol_m3: Mapping[str, float], *, measured_ph: float, temperature_k: float, hydrogen_activity_coefficient: float=1.0, hydroxide_activity_coefficient: float=1.0) -> float:
    """Infer the transported unmodelled-ion charge from one measured pH."""
    ph = float(measured_ph)
    if not math.isfinite(ph) or not 0.0 <= ph <= 14.0:
        raise ValueError('measured_ph must be finite in [0, 14]')
    explicit = reduced_charge_density_eq_m3(state_mol_m3, hydrogen_activity_mol_l=10.0 ** (-ph), temperature_k=temperature_k, fixed_background_charge_eq_m3=0.0, hydrogen_activity_coefficient=hydrogen_activity_coefficient, hydroxide_activity_coefficient=hydroxide_activity_coefficient)
    return -explicit

def solve_relative_charge_balance(state_mol_m3: Mapping[str, float], *, fixed_background_charge_eq_m3: float, temperature_k: float, lower_ph: float=0.0, upper_ph: float=14.0, absolute_tolerance_eq_m3: float=1e-08, ph_tolerance: float=1e-10, max_iterations: int=200, hydrogen_activity_coefficient: float=1.0, hydroxide_activity_coefficient: float=1.0) -> ChargeBalanceResult:
    """Solve pH by bisection and return an explicit failure classification."""
    lower = float(lower_ph)
    upper = float(upper_ph)
    if not all((math.isfinite(value) for value in (lower, upper))) or not 0.0 <= lower < upper <= 14.0:
        raise ValueError('pH bounds must satisfy 0 <= lower < upper <= 14')
    if not math.isfinite(float(absolute_tolerance_eq_m3)) or absolute_tolerance_eq_m3 < 0.0:
        raise ValueError('absolute charge tolerance must be finite and nonnegative')
    if not math.isfinite(float(ph_tolerance)) or ph_tolerance <= 0.0:
        raise ValueError('pH tolerance must be finite and positive')
    if isinstance(max_iterations, bool) or int(max_iterations) != max_iterations or max_iterations < 1:
        raise ValueError('max_iterations must be a positive integer')

    def residual(ph: float) -> float:
        return reduced_charge_density_eq_m3(state_mol_m3, hydrogen_activity_mol_l=10.0 ** (-ph), temperature_k=temperature_k, fixed_background_charge_eq_m3=fixed_background_charge_eq_m3, hydrogen_activity_coefficient=hydrogen_activity_coefficient, hydroxide_activity_coefficient=hydroxide_activity_coefficient)
    low_residual = residual(lower)
    high_residual = residual(upper)
    if abs(low_residual) <= absolute_tolerance_eq_m3:
        selected_ph, selected_residual, iterations = (lower, low_residual, 0)
    elif abs(high_residual) <= absolute_tolerance_eq_m3:
        selected_ph, selected_residual, iterations = (upper, high_residual, 0)
    elif low_residual * high_residual > 0.0:
        return ChargeBalanceResult(status='charge_balance_unbracketed', converged=False, bracketed=False, ph=None, hydrogen_activity_mol_l=None, hydrogen_concentration_mol_m3=None, charge_residual_eq_m3=None, fixed_background_charge_eq_m3=float(fixed_background_charge_eq_m3), iterations=0, lower_ph=lower, upper_ph=upper)
    else:
        selected_ph = 0.5 * (lower + upper)
        selected_residual = residual(selected_ph)
        iterations = 0
        for iterations in range(1, int(max_iterations) + 1):
            selected_ph = 0.5 * (lower + upper)
            selected_residual = residual(selected_ph)
            if abs(selected_residual) <= absolute_tolerance_eq_m3 or upper - lower <= ph_tolerance:
                break
            if low_residual * selected_residual <= 0.0:
                upper = selected_ph
                high_residual = selected_residual
            else:
                lower = selected_ph
                low_residual = selected_residual
        else:
            return ChargeBalanceResult(status='charge_balance_nonconverged', converged=False, bracketed=True, ph=selected_ph, hydrogen_activity_mol_l=10.0 ** (-selected_ph), hydrogen_concentration_mol_m3=10.0 ** (-selected_ph) / hydrogen_activity_coefficient * 1000.0, charge_residual_eq_m3=selected_residual, fixed_background_charge_eq_m3=float(fixed_background_charge_eq_m3), iterations=int(max_iterations), lower_ph=float(lower_ph), upper_ph=float(upper_ph))
    h_activity = 10.0 ** (-selected_ph)
    converged = abs(selected_residual) <= absolute_tolerance_eq_m3 or upper - lower <= ph_tolerance
    return ChargeBalanceResult(status='converged' if converged else 'charge_balance_nonconverged', converged=converged, bracketed=True, ph=selected_ph, hydrogen_activity_mol_l=h_activity, hydrogen_concentration_mol_m3=h_activity / hydrogen_activity_coefficient * 1000.0, charge_residual_eq_m3=selected_residual, fixed_background_charge_eq_m3=float(fixed_background_charge_eq_m3), iterations=iterations, lower_ph=float(lower_ph), upper_ph=float(upper_ph))
__all__ = ['ChargeBalanceResult', 'infer_fixed_background_charge_eq_m3', 'reduced_acid_base_average_charges', 'reduced_charge_density_eq_m3', 'solve_relative_charge_balance']
