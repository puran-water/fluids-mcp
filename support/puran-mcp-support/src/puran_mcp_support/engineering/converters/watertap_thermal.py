from __future__ import annotations
from typing import Any, Dict

def _mcas_tds_mg_l(mcas: Dict[str, float]) -> float:
    total = 0.0
    for key, value in mcas.items():
        if not key.endswith('_mg_l'):
            continue
        if key in {'temp'}:
            continue
        total += float(value)
    return total

def eng_mcas_to_nacl_thermal_state(mcas: Dict[str, float], *, flow_m3_hr: float=1.0, temperature_c: float=25.0, pressure_bar: float=1.013) -> Dict[str, Any]:
    """Collapse engineering MCAS into a NaCl-equivalent thermal-brine feed."""
    nacl_mg_l = _mcas_tds_mg_l(mcas)
    return {'flow_vol_m3_hr': flow_m3_hr, 'temperature_C': temperature_c, 'pressure_bar': pressure_bar, 'components': {'NaCl': nacl_mg_l}, 'concentration_units': 'mg/L', 'concentration_basis': 'mass', 'estimated_fields': ['NaCl'], 'available_fields': ['NaCl'], 'source_breakdown_mg_l': {key: float(value) for key, value in mcas.items() if key.endswith('_mg_l') and float(value) > 0}}

def eng_mcas_to_md_thermal_state(mcas: Dict[str, float], *, flow_m3_hr: float=1.0, temperature_c: float=25.0, pressure_bar: float=1.013) -> Dict[str, Any]:
    """Collapse engineering MCAS into a seawater/TDS membrane-distillation feed."""
    tds_mg_l = _mcas_tds_mg_l(mcas)
    return {'flow_vol_m3_hr': flow_m3_hr, 'temperature_C': temperature_c, 'pressure_bar': pressure_bar, 'components': {'TDS': tds_mg_l}, 'concentration_units': 'mg/L', 'concentration_basis': 'mass', 'estimated_fields': ['TDS'], 'available_fields': ['TDS'], 'source_breakdown_mg_l': {key: float(value) for key, value in mcas.items() if key.endswith('_mg_l') and float(value) > 0}}
