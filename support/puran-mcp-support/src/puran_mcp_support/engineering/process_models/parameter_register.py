from __future__ import annotations
from ._core.parameters import ParameterProvenance, ParameterSet, ParameterTier, ParameterValue, ValidityWindow
DIGESTER_MODEL_KEY = 'process_models.anaerobic_digester/v1'

def _v(value: float, unit: str='1', note: str='') -> ParameterValue:
    return ParameterValue(value=float(value), unit=unit, note=note)
MADM1_KINETICS_SUBJECT = 'mADM1 core kinetics (shared biokinetic kernel)'

def madm1_core_kinetics_house_default() -> list[ParameterSet]:
    """The ADM1 core's kinetic and stoichiometric parameters -> ``house_default``, house library."""
    from puran_mcp_support.engineering._sibling_libs import ensure_sibling_importable
    ensure_sibling_importable('puran_mcp_support.schema')
    from puran_mcp_support.schema.madm1_biokinetics import MADM1_BIOKINETICS_KERNEL as kernel
    values = {p.name: _v(p.value, p.unit, f'{p.role}; {p.source}') for p in kernel.parameters}
    return [ParameterSet(model_key=DIGESTER_MODEL_KEY, subject=MADM1_KINETICS_SUBJECT, tier=ParameterTier.HOUSE_DEFAULT, values=values, settings={'kernel_id': kernel.kernel_id, 'kernel_sha256': kernel.sha256, 'processes': str(len(kernel.processes)), 'reference': kernel.reference}, provenance=ParameterProvenance(method='Literature-based modified ADM1 kernel (QSDsan lineage); kinetic parameters are not site fits.', produced_by='puran_mcp_support.schema.madm1_biokinetics', evidence_refs=('Public source package; internal evidence not distributed.', 'QSDsan and WaterTAP upstream source references are retained in the kernel definitions.'), source_object='process_model_schema.madm1_biokinetics:MADM1_BIOKINETICS_KERNEL', source_sha256=kernel.sha256))]
MADM1_SFE_KINETICS_SUBJECT = 'mADM1 sulfur and iron kinetics (shared S/Fe kernel)'

def madm1_sfe_kinetics_house_default() -> list[ParameterSet]:
    """The S/Fe extension's parameters (SRB guilds, HFO reduction/ageing/P exchange, FeS) -> ``house_default``.

    Built FROM ``process_model_schema.madm1_sfe_kinetics.MADM1_SFE_PARAMETERS`` (phase 3a, owner ruling 2026-10-07
    "After SRB in WaterTAP" (HOUSE-MADM1-SRB-IN-WATERTAP-FIRST-20261007)): the twin's literature defaults for its processes 25-55 and its FeS law, the set the
    WaterTAP S/Fe digester reads through custody beside the core set. NOT a fit; the core set keeps its own subject.
    """
    from puran_mcp_support.engineering._sibling_libs import ensure_sibling_importable
    ensure_sibling_importable('puran_mcp_support.schema')
    from puran_mcp_support.schema.madm1_sfe_kinetics import MADM1_SFE_KERNEL as kernel, MADM1_SFE_PARAMETERS
    values = {p.name: _v(p.value, p.unit, f'{p.role}; {p.source}') for p in MADM1_SFE_PARAMETERS}
    return [ParameterSet(model_key=DIGESTER_MODEL_KEY, subject=MADM1_SFE_KINETICS_SUBJECT, tier=ParameterTier.HOUSE_DEFAULT, values=values, settings={'kernel_id': kernel.kernel_id, 'kernel_sha256': kernel.sha256, 'processes': str(len(kernel.processes)), 'reference': kernel.reference}, provenance=ParameterProvenance(method='Literature-based modified ADM1 kernel (QSDsan lineage); kinetic parameters are not site fits.', produced_by='puran_mcp_support.schema.madm1_sfe_kinetics', evidence_refs=('Public source package; internal evidence not distributed.', 'QSDsan and WaterTAP upstream source references are retained in the kernel definitions.'), source_object='process_model_schema.madm1_sfe_kinetics:MADM1_SFE_PARAMETERS', source_sha256=kernel.sha256))]
