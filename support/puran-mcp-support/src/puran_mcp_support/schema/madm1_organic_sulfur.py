from __future__ import annotations
import hashlib
import json
from dataclasses import dataclass
from .expression import multiply
from .madm1_biokinetics import BiokineticProcessDefinition, ComponentContent, _s, _st
from .madm1_sfe_kinetics import MADM1_SFE_KERNEL
X_SORG = 'X_Sorg'
MADM1_ORGANIC_SULFUR_CONTENTS = (ComponentContent(X_SORG, 0.0, 0.0, 0.0, 0.0, 0.0, 'S'),)
ORGANIC_SULFUR_RATE_CONSTANT = 'q_pr_hyd'
ORGANIC_SULFUR_SOURCE = "PUBLIC-RULE-86fa3c55009b 'Library values + organic S release (Recommended)'; rate = protein hydrolysis constant q_pr_hyd x X_Sorg, [ASSUMED] proportionality (coordinator 2026-10-09, plan step 1)"

class OrganicSulfurProcessDefinition(BiokineticProcessDefinition):
    """A process of the opt-in organic-sulfur extension: ``index`` is its TWIN index (after the gas transfers)."""
MADM1_ORGANIC_SULFUR_HYDROLYSIS = OrganicSulfurProcessDefinition('hydrolysis_organic_S', 75, multiply(_s(ORGANIC_SULFUR_RATE_CONSTANT), _s(X_SORG)), X_SORG, {X_SORG: _st(x=-1.0)['x'], 'S_IS': _st(x=1.0)['x']})
MADM1_ORGANIC_SULFUR_PROCESSES = (MADM1_ORGANIC_SULFUR_HYDROLYSIS,)

@dataclass(frozen=True)
class OrganicSulfurExtension:
    """The opt-in extension as one object: its component, its processes, the base kernel it extends, and an identity
    hash that an engine stamps beside the base kernel's when (and only when) the extension is compiled."""
    extension_id: str
    base_kernel_id: str
    base_kernel_sha256: str
    contents: tuple[ComponentContent, ...]
    processes: tuple[OrganicSulfurProcessDefinition, ...]
    source: str

    @property
    def component_ids(self) -> tuple[str, ...]:
        return tuple((c.component_id for c in self.contents))

    @property
    def sha256(self) -> str:
        body = {'extension_id': self.extension_id, 'base': [self.base_kernel_id, self.base_kernel_sha256], 'contents': [sorted(vars(c).items()) for c in self.contents], 'processes': [[p.process_id, p.index, repr(p.rate_base), p.reference_component, sorted(((k, repr(v)) for k, v in p.stoichiometry.items())), repr(p.inhibition)] for p in self.processes], 'source': self.source}
        return hashlib.sha256(json.dumps(body, sort_keys=True, default=str).encode()).hexdigest()

def validate_organic_sulfur_extension(ext: 'OrganicSulfurExtension') -> None:
    """Refuse an extension that would not compile against its base: its component must be new to the base, every
    stoichiometric component and rate symbol must be the base's or its own, the indices must run on from 75, and the
    sulfur row must close (every coefficient is on a kg-S component)."""
    base = MADM1_SFE_KERNEL
    if ext.base_kernel_sha256 != base.sha256:
        raise ValueError('the organic-sulfur extension is defined against a different S/Fe kernel')
    clash = set(ext.component_ids) & set(base.component_ids)
    if clash:
        raise ValueError(f'organic-sulfur components already in the base kernel: {sorted(clash)}')
    known = set(base.component_ids) | set(ext.component_ids)
    params = {p.name for p in base.parameters}
    for offset, proc in enumerate(ext.processes):
        if proc.index != 75 + offset:
            raise ValueError(f'{proc.process_id}: index {proc.index}, expected {75 + offset} (after the gas transfers)')
        unknown = set(proc.stoichiometry) - known
        if unknown:
            raise ValueError(f'{proc.process_id}: stoichiometry names unknown components {sorted(unknown)}')
        if ORGANIC_SULFUR_RATE_CONSTANT not in params:
            raise ValueError(f'{ORGANIC_SULFUR_RATE_CONSTANT} is not a base-kernel parameter')
        measured = {c.component_id: c.measured_as for c in base.components + ext.contents}
        if any((measured[c] != 'S' for c in proc.stoichiometry)):
            raise ValueError(f'{proc.process_id}: every coefficient must be on a kg-S component')
        values = [getattr(v, 'value', None) for v in proc.stoichiometry.values()]
        if any((not isinstance(v, (int, float)) for v in values)):
            raise ValueError(f'{proc.process_id}: every coefficient must be a stated number (the S row is checked)')
        if abs(sum(values)) > 1e-12:
            raise ValueError(f'{proc.process_id}: the sulfur row does not close ({sum(values)!r})')
MADM1_ORGANIC_SULFUR_EXTENSION = OrganicSulfurExtension(extension_id='madm1_organic_sulfur_v1', base_kernel_id=MADM1_SFE_KERNEL.kernel_id, base_kernel_sha256=MADM1_SFE_KERNEL.sha256, contents=MADM1_ORGANIC_SULFUR_CONTENTS, processes=MADM1_ORGANIC_SULFUR_PROCESSES, source=ORGANIC_SULFUR_SOURCE)
