from __future__ import annotations
from typing import Any, Dict, List, Optional, Tuple
_MW_N = 14.007
_MW_P = 30.974
SPECIES_MAP: Dict[str, Tuple[str, int, float, float, str]] = {'Na_mg_l': ('Na_+', 1, 22.99, 1.84e-10, 'as_element'), 'K_mg_l': ('K_+', 1, 39.098, 1.25e-10, 'as_element'), 'Ca_mg_l': ('Ca_2+', 2, 40.078, 3.1e-10, 'as_element'), 'Mg_mg_l': ('Mg_2+', 2, 24.305, 3.47e-10, 'as_element'), 'Cl_mg_l': ('Cl_-', -1, 35.453, 1.21e-10, 'as_element'), 'NH4_N_mg_l': ('NH4_+', 1, 18.039, 1.25e-10, 'as_N'), 'NO3_N_mg_l': ('NO3_-', -1, 62.005, 1.29e-10, 'as_N'), 'HCO3_mg_l': ('HCO3_-', -1, 61.017, 1.18e-10, 'ion'), 'SO4_mg_l': ('SO4_2-', -2, 96.062, 2.3e-10, 'ion'), 'PO4_P_mg_l': ('HPO4_2-', -2, 95.979, 3.27e-10, 'as_P'), 'SiO2_mg_l': ('SiO2', 0, 60.084, 2.12e-10, 'ion')}

def eng_mcas_to_watertap_components(mcas: Dict[str, float], *, include_zero: bool=False) -> Dict[str, float]:
    """Convert engineering MCAS dict to canonical WaterTAP public feed components.

    Returned concentrations are in mg/L on the full-ion basis expected by the
    WaterTAP MCP `set_feed_state(..., concentration_units="mg/L")` interface.
    """
    components: Dict[str, float] = {}
    for eng_key, (wt_name, _z, mw, _r, basis) in SPECIES_MAP.items():
        conc_mg_l = float(mcas.get(eng_key, 0.0))
        if conc_mg_l <= 0 and (not include_zero):
            continue
        if basis == 'as_N':
            components[wt_name] = conc_mg_l * mw / _MW_N
        elif basis == 'as_P':
            components[wt_name] = conc_mg_l * mw / _MW_P
        else:
            components[wt_name] = conc_mg_l
    return components

def eng_mcas_to_component_charges(mcas: Dict[str, float], *, include_zero: bool=False) -> Dict[str, int]:
    """Return WaterTAP canonical component-charge map for engineering MCAS payloads."""
    charges: Dict[str, int] = {}
    for eng_key, (wt_name, z, _mw, _r, _basis) in SPECIES_MAP.items():
        conc_mg_l = float(mcas.get(eng_key, 0.0))
        if conc_mg_l > 0 or include_zero:
            charges[wt_name] = z
    return charges

def eng_mcas_to_ed_nacl_pair_state(mcas: Dict[str, float], *, flow_m3_hr: float=1.0, temperature_c: float=25.0, pressure_bar: float=1.013) -> Dict[str, Any]:
    """Collapse engineering MCAS into an approximate Na+/Cl- ED feed.

    This reducer is intentionally approximate. It preserves total dissolved
    solids mass while collapsing broader ionic chemistry onto an NaCl-equivalent
    sodium/chloride pair suitable for the narrow WaterTAP ED public profile.
    """
    total_mg_l = 0.0
    source_breakdown_mg_l: Dict[str, float] = {}
    for key, value in mcas.items():
        if not key.endswith('_mg_l'):
            continue
        if key == 'temp':
            continue
        numeric = float(value)
        if numeric <= 0:
            continue
        total_mg_l += numeric
        source_breakdown_mg_l[key] = numeric
    nacl_mw = SPECIES_MAP['Na_mg_l'][2] + SPECIES_MAP['Cl_mg_l'][2]
    na_fraction = SPECIES_MAP['Na_mg_l'][2] / nacl_mw
    cl_fraction = SPECIES_MAP['Cl_mg_l'][2] / nacl_mw
    na_mg_l = total_mg_l * na_fraction
    cl_mg_l = total_mg_l * cl_fraction
    return {'flow_vol_m3_hr': flow_m3_hr, 'temperature_C': temperature_c, 'pressure_bar': pressure_bar, 'components': {'Na_+': na_mg_l, 'Cl_-': cl_mg_l}, 'concentration_units': 'mg/L', 'concentration_basis': 'mass', 'component_charges': {'Na_+': 1, 'Cl_-': -1}, 'estimated_fields': ['Na_+', 'Cl_-'], 'available_fields': ['Na_+', 'Cl_-'], 'source_breakdown_mg_l': source_breakdown_mg_l}
