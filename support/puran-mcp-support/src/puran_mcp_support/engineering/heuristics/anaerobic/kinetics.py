from __future__ import annotations
from functools import lru_cache
from typing import Any, Mapping, Sequence
DIGESTER_MODEL_KEY = 'process_models.anaerobic_digester/v1'
MADM1_KINETICS_SUBJECT = 'mADM1 core kinetics (shared biokinetic kernel)'

def _code_register_store():
    from puran_mcp_support.engineering.process_models._core.custody import MemoryParameterStore
    from puran_mcp_support.engineering.process_models.parameter_register import madm1_core_kinetics_house_default
    store = MemoryParameterStore()
    for ps in madm1_core_kinetics_house_default():
        store.put(ps, audit_user_id='code-register', audit_reason='the house library, in memory')
    return store

def madm1_design_parameter_set(store: Any=None):
    """The ``StoredParameterSet`` a digester DESIGN reads: ``custody.design_parameter_set`` (house_ruled preferred,
    else house_default; any other tier refused) for the mADM1 core kinetics, from ``store`` -- the engineering_mcp
    custody store, or by default the code-resident register in memory (``parameters verify`` proves the persisted row
    equal to it). Raises ``custody.NoParameterSet`` when no design set exists, never a fallback."""
    from puran_mcp_support.engineering.process_models._core.custody import design_parameter_set
    return design_parameter_set(store if store is not None else _code_register_store(), model_key=DIGESTER_MODEL_KEY, subject=MADM1_KINETICS_SUBJECT)

@lru_cache(maxsize=1)
def _default_values() -> tuple[tuple[str, float], ...]:
    stored = madm1_design_parameter_set()
    return tuple(sorted(((k, float(v.value)) for k, v in stored.parameter_set.values.items())))

def madm1_design_parameters(store: Any=None) -> dict[str, Any]:
    """{name: value} of the design set (the kernel's names), plus ``_identity``: the set's model key, subject, tier,
    set hash and kernel id, for the record every consumer stamps beside its answer."""
    stored = madm1_design_parameter_set(store)
    ps = stored.parameter_set
    values = {k: float(v.value) for k, v in ps.values.items()} if store is not None else dict(_default_values())
    values['_identity'] = {'model_key': ps.model_key, 'subject': ps.subject, 'tier': ps.tier.value, 'set_hash': ps.set_hash(), 'kernel_id': ps.settings.get('kernel_id'), 'kernel_sha256': ps.settings.get('kernel_sha256')}
    return values
