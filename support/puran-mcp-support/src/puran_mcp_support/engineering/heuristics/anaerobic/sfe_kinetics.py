from __future__ import annotations
from functools import lru_cache
from typing import Any
from puran_mcp_support.engineering.process_models.parameter_register import MADM1_SFE_KINETICS_SUBJECT
DIGESTER_MODEL_KEY = 'process_models.anaerobic_digester/v1'

def _code_register_store():
    from puran_mcp_support.engineering.process_models._core.custody import MemoryParameterStore
    from puran_mcp_support.engineering.process_models.parameter_register import madm1_sfe_kinetics_house_default
    store = MemoryParameterStore()
    for ps in madm1_sfe_kinetics_house_default():
        store.put(ps, audit_user_id='code-register', audit_reason='the house library, in memory')
    return store

def madm1_sfe_design_parameter_set(store: Any=None):
    """The ``StoredParameterSet`` an S/Fe digester DESIGN reads (house_ruled preferred, else house_default; any other
    tier refused), from ``store`` or by default the code-resident register in memory. Raises
    ``custody.NoParameterSet`` when no design set exists, never a fallback."""
    from puran_mcp_support.engineering.process_models._core.custody import design_parameter_set
    return design_parameter_set(store if store is not None else _code_register_store(), model_key=DIGESTER_MODEL_KEY, subject=MADM1_SFE_KINETICS_SUBJECT)

@lru_cache(maxsize=1)
def _default_values() -> tuple[tuple[str, float], ...]:
    stored = madm1_sfe_design_parameter_set()
    return tuple(sorted(((k, float(v.value)) for k, v in stored.parameter_set.values.items())))

def madm1_sfe_design_parameters(store: Any=None) -> dict[str, Any]:
    """{name: value} of the S/Fe design set (the S/Fe kernel's own parameter names), plus ``_identity``."""
    stored = madm1_sfe_design_parameter_set(store)
    ps = stored.parameter_set
    values = {k: float(v.value) for k, v in ps.values.items()} if store is not None else dict(_default_values())
    values['_identity'] = {'model_key': ps.model_key, 'subject': ps.subject, 'tier': ps.tier.value, 'set_hash': ps.set_hash(), 'kernel_id': ps.settings.get('kernel_id'), 'kernel_sha256': ps.settings.get('kernel_sha256')}
    return values
