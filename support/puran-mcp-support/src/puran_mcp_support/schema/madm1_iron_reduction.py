from __future__ import annotations
import hashlib
import json
from dataclasses import dataclass
from .expression import monod, multiply, power
from .madm1_biokinetics import BiokineticProcessDefinition, KineticParameterDefinition, _p, _s, _st, _symbols
from .madm1_microaeration import MICROAERATION_PROCESSES, _value
from .madm1_sfe_kinetics import MADM1_SFE_KERNEL
RULING_ID = 'PUBLIC-RULE-ad9bdd59a2d4'
DONOR = 'S_ac'
COD_PER_KMOL_ELECTRON = 8.0
IRON_REDUCTION_RATE_SOURCE = "Wang, Amin, Hauduc, Graan, Guerin, Mikola & Sperandio 2026, Water Res. X 30:100482, doi:10.1016/j.wroa.2025.100482 (full text, HAL hal-05686846v1, read 2026-10-24): s2.2 'the rate of HFO reduction with organics to 0.80 d-1 according to the batch test results'; s3.2 'A first-order kinetic equation was considered for modeling reduction rate of Fe3+ into Fe2+ ... The iron reduction rates (measured at 20 C) ... were similar in both RAS and MS (0.80 +/- 0.1 d-1)' (anaerobic batches, FeCl3 added, pH held near 7; the digested-sludge batch showed no measurable rise); s4.1 'the lower end of the range reported in the literature (0.72 - 1.44 d-1)' (Prot et al. 2022, doi:10.1016/j.jece.2022.107473, and Wang R. et al. 2019, doi:10.1016/j.seppur.2019.03.024: abstracts only, the range as Wang et al. 2026 report it)"
IRON_REDUCTION_FORM_SOURCE = "Sumo2S model description (Dynamita, wiki.dynamita.com/en/biokinetic_process_models, 'Reduction of Fe3+, organic matter as electron donor', read 2026-10-24): 'this new iron reducing biomass and the associated growth is not introduced in the model. The soluble biodegradable substrate (SB) and volatile acids (SVFA) are considered as electron donor and a first order kinetic rate expression with respect to the hydrous ferric oxide (XHFO) concentration is used', X_HFO being the sum of the seven HFO states (Hauduc et al. 2019, Water Sci. Technol. 79:26, doi:10.2166/wst.2018.482, abstract); the process Wang et al. 2026 calibrated"
IRON_REDUCTION_DONOR_SOURCE = "Lovley & Phillips 1988, Appl. Environ. Microbiol. 54:1472, doi:10.1128/aem.54.6.1472-1480.1988 (abstract, cited by Sumo2S): GS-15 'oxidized acetate to carbon dioxide with the concomitant reduction of amorphic Fe(III) oxide'; 'Ethanol, butyrate, and propionate could also serve as electron donors ... hydrogen could not'. Acetate carries the reduction here (the kernel's VFA end product; the H2 route stays the kernel's abiotic one): 1 e- per Fe at 8 kg COD per kmol e- (COD_per_kmol_H2 / 2)"
IRON_REDUCTION_TEMPERATURE_SOURCE = "no temperature law: Lovley & Phillips 1988 (abstract) 'Fe(III) was reduced optimally ... at 30 to 35 degrees C', so the 20 C rate understates a mesophilic digester's; Schilling et al. 2018, Biogeochemistry, doi:10.1007/s10533-018-0520-4 (abstract) found Q10 1.5-8.9 for microbial Fe(III) reduction across soils, too wide to transfer. Carried at the measured 20 C value; bracketed to 1.44 /d (the literature top, and Q10 1.5 over 20-35 C) in lane c106's evidence"
DONOR_SWITCH_SOURCE = "a numerical donor switch on acetate so the reduction stops where no donor is left; at any digester acetate (> 0.01 kg COD/m3) it is 1 within 1 %. [ASSUMED] band 1e-5-1e-3, the micro-aeration extension's K_S0red_donor convention: Wang et al. 2026 measured the rate first order with the sludge's own donor in excess, no donor half-saturation is sourced"
IRON_REDUCTION_PARAMETERS: tuple[KineticParameterDefinition, ...] = (_p('k_Fe3t2_org', 0.8, '1/d', 'rate', IRON_REDUCTION_RATE_SOURCE + '; ' + IRON_REDUCTION_TEMPERATURE_SOURCE), _p('K_Fe3t2_org_donor', 0.0001, 'kg COD/m3', 'half_saturation', DONOR_SWITCH_SOURCE))
DONOR_SWITCH_BAND = (1e-05, 0.001)
HFO_HOSTS = (('X_HFO_H', 'HFO_H'), ('X_HFO_L', 'HFO_L'), ('X_HFO_old', 'HFO_old'), ('X_HFO_HP', 'HFO_HP'), ('X_HFO_LP', 'HFO_LP'), ('X_HFO_HP_old', 'HFO_HP_old'), ('X_HFO_LP_old', 'HFO_LP_old'))
FIRST_INDEX = MICROAERATION_PROCESSES[-1].index + 1

def _reduction(index: int, host: str, suffix: str) -> BiokineticProcessDefinition:
    """One host's reduction: first order in its Fe(III), the donor switch on acetate; kg Fe / m3 / d."""
    acetate_per_fe = multiply(-0.5, _s('COD_per_kmol_H2'), power(_s('MW_Fe'), -1.0))
    return BiokineticProcessDefinition(f'reduction_{suffix}_ac', index, multiply(_s('k_Fe3t2_org'), _s(host), monod(_s(DONOR), _s('K_Fe3t2_org_donor'))), host, {host: _st(x=-1.0)['x'], 'S_Fe2': _st(x=1.0)['x'], DONOR: acetate_per_fe})
MADM1_IRON_REDUCTION_PROCESSES: tuple[BiokineticProcessDefinition, ...] = tuple((_reduction(FIRST_INDEX + i, host, suffix) for i, (host, suffix) in enumerate(HFO_HOSTS)))

@dataclass(frozen=True)
class IronReductionExtension:
    """The extension as one object, with an identity hash an engine stamps beside the base kernel's."""
    extension_id: str
    base_kernel_id: str
    base_kernel_sha256: str
    parameters: tuple[KineticParameterDefinition, ...]
    processes: tuple[BiokineticProcessDefinition, ...]
    source: str

    def parameter_values(self) -> dict[str, float]:
        return {p.name: float(p.value) for p in self.parameters}

    @property
    def sha256(self) -> str:
        body = {'extension_id': self.extension_id, 'base': [self.base_kernel_id, self.base_kernel_sha256], 'parameters': [[p.name, p.value, p.unit, p.role, p.source] for p in self.parameters], 'processes': [[p.process_id, p.index, repr(p.rate_base), p.reference_component, sorted(((k, repr(v)) for k, v in p.stoichiometry.items())), repr(p.inhibition)] for p in self.processes], 'source': self.source}
        return hashlib.sha256(json.dumps(body, sort_keys=True, default=str).encode()).hexdigest()

def validate_iron_reduction_extension(ext: IronReductionExtension) -> None:
    """Refuse an extension that would not compile against its base: parameters disjoint from the base's, each sourced
    and the donor switch inside its band; one process per HFO host, in order from `FIRST_INDEX`; every stoichiometric
    component the base's and every symbol a parameter or a component; Fe moving host -> S_Fe2 one for one; the
    acetate giving exactly the electrons the Fe(III) takes (`COD_PER_KMOL_ELECTRON`)."""
    base = MADM1_SFE_KERNEL
    if ext.base_kernel_sha256 != base.sha256:
        raise ValueError('the iron-reduction extension is defined against a different S/Fe kernel')
    base_values = base.parameter_values()
    values = ext.parameter_values()
    overlap = set(values) & set(base_values)
    if overlap:
        raise ValueError(f'iron-reduction parameters already in the base kernel: {sorted(overlap)}')
    unsourced = [p.name for p in ext.parameters if not p.source.strip()]
    if unsourced:
        raise ValueError(f'iron-reduction parameters without a source: {unsourced}')
    low, high = DONOR_SWITCH_BAND
    if not low <= values['K_Fe3t2_org_donor'] <= high:
        raise ValueError(f'K_Fe3t2_org_donor {values['K_Fe3t2_org_donor']} is outside its [ASSUMED] band {DONOR_SWITCH_BAND}')
    hosts = [h for h, _ in HFO_HOSTS]
    if [p.reference_component for p in ext.processes] != hosts:
        raise ValueError(f'the extension reduces {[p.reference_component for p in ext.processes]}, not the hosts {hosts}')
    known = set(base.component_ids)
    params = set(base_values) | set(values)
    for offset, proc in enumerate(ext.processes):
        if proc.index != FIRST_INDEX + offset:
            raise ValueError(f'{proc.process_id}: index {proc.index}, expected {FIRST_INDEX + offset}')
        unknown = set(proc.stoichiometry) - known
        if unknown:
            raise ValueError(f'{proc.process_id}: stoichiometry names unknown components {sorted(unknown)}')
        undefined = _symbols(proc.rate) - params - known
        for coefficient in proc.stoichiometry.values():
            undefined |= _symbols(coefficient) - params
        if undefined:
            raise ValueError(f'{proc.process_id}: reads undefined symbols {sorted(undefined)}')
        host, ferrous = (proc.stoichiometry[proc.reference_component], proc.stoichiometry['S_Fe2'])
        if getattr(host, 'value', None) != -1.0 or getattr(ferrous, 'value', None) != 1.0:
            raise ValueError(f'{proc.process_id}: Fe must move from its host to S_Fe2 one for one')
        numbers = {**base_values, **values}
        taken = _value(ferrous, numbers) / numbers['MW_Fe']
        given = -_value(proc.stoichiometry[DONOR], numbers) / COD_PER_KMOL_ELECTRON
        if abs(taken - given) > 1e-12 * taken:
            raise ValueError(f'{proc.process_id}: the electron (COD) row does not close: the Fe(III) takes {taken!r} kmol e- per kg Fe, the acetate gives {given!r}')
MADM1_IRON_REDUCTION_EXTENSION = IronReductionExtension(extension_id='madm1_iron_reduction_v1', base_kernel_id=MADM1_SFE_KERNEL.kernel_id, base_kernel_sha256=MADM1_SFE_KERNEL.sha256, parameters=IRON_REDUCTION_PARAMETERS, processes=MADM1_IRON_REDUCTION_PROCESSES, source=f'{RULING_ID}; {IRON_REDUCTION_FORM_SOURCE}; {IRON_REDUCTION_RATE_SOURCE}; {IRON_REDUCTION_DONOR_SOURCE}; {IRON_REDUCTION_TEMPERATURE_SOURCE}')
