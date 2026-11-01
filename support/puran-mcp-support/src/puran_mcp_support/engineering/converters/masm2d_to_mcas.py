from __future__ import annotations
import logging
from typing import Any, Dict, Optional
logger = logging.getLogger(__name__)
MW_Na = 22.99
MW_K = 39.098
MW_Ca = 40.078
MW_Mg = 24.305
MW_Cl = 35.453
MW_N = 14.007
MW_P = 30.974
MW_C = 12.011
MW_HCO3 = 61.017

def masm2d_to_mcas(state: Dict[str, float], *, charge_balance: bool=True, pH: Optional[float]=None, temperature_c: float=25.0) -> Dict[str, Any]:
    """Extract MCAS ion composition from an mASM2d state dict.

    Args:
        state: mASM2d state variables in mg/L (keys: S_Na, S_K, …).
        charge_balance: If True, adjust Na+ or Cl- so that
            sum(cation equivalents) == sum(anion equivalents).
        pH: Optional pH to include in the output dict.
        temperature_c: Temperature (°C) to include in the output dict.

    Returns:
        MCAS-format dict with ``{ion}_mg_l`` keys (mg/L).
        Includes ``pH`` and ``temp`` if provided.
    """
    mcas: Dict[str, Any] = {}
    mcas['Na_mg_l'] = state.get('S_Na', 0.0)
    mcas['K_mg_l'] = state.get('S_K', 0.0)
    mcas['Ca_mg_l'] = state.get('S_Ca', 0.0)
    mcas['Mg_mg_l'] = state.get('S_Mg', 0.0)
    mcas['Cl_mg_l'] = state.get('S_Cl', 0.0)
    mcas['NH4_N_mg_l'] = state.get('S_NH4', 0.0)
    mcas['NO3_N_mg_l'] = state.get('S_NO3', 0.0)
    mcas['NO2_N_mg_l'] = 0.0
    mcas['PO4_P_mg_l'] = state.get('S_PO4', 0.0)
    s_ic = state.get('S_IC', 0.0)
    mcas['HCO3_mg_l'] = s_ic * MW_HCO3 / MW_C
    mcas['SO4_mg_l'] = 0.0
    mcas['SiO2_mg_l'] = 0.0
    mcas['Fe_mg_l'] = 0.0
    mcas['Al_mg_l'] = 0.0
    if charge_balance:
        _enforce_charge_balance(mcas)
    if pH is not None:
        mcas['pH'] = pH
    mcas['temp'] = temperature_c
    return mcas

def _enforce_charge_balance(mcas: Dict[str, float]) -> None:
    """Add Na+ or Cl- so cation meq/L == anion meq/L.

    Mutates *mcas* in place.
    """
    cations_meq = mcas['Na_mg_l'] / MW_Na + mcas['K_mg_l'] / MW_K + mcas['Ca_mg_l'] / MW_Ca * 2 + mcas['Mg_mg_l'] / MW_Mg * 2 + mcas['NH4_N_mg_l'] / MW_N
    anions_meq = mcas['Cl_mg_l'] / MW_Cl + mcas['NO3_N_mg_l'] / MW_N + mcas['NO2_N_mg_l'] / MW_N + mcas['HCO3_mg_l'] / MW_HCO3 + mcas['PO4_P_mg_l'] / MW_P * 2 + mcas['SO4_mg_l'] / 96.062 * 2
    deficit_meq = cations_meq - anions_meq
    if deficit_meq > 1e-09:
        mcas['Cl_mg_l'] += deficit_meq * MW_Cl
        logger.debug('Charge-balanced by adding %.2f meq/L Cl⁻', deficit_meq)
    elif deficit_meq < -1e-09:
        mcas['Na_mg_l'] += abs(deficit_meq) * MW_Na
        logger.debug('Charge-balanced by adding %.2f meq/L Na⁺', abs(deficit_meq))
