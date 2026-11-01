"""Shared mADM1 mineral-extension registry: the Ca/Mg/K/Al phase set.

Companion to the Fe-S-P kernel (``vertical_slices``): the eleven non-iron
phases of the mADM1 mineral subsystem, defined ONCE so the qsdsan registry
and the water-chemistry PHREEQC benchmark consume identical constants. It is
deliberately a SEPARATE kernel — extending ``MADM1_FECL3_FE_S_P_KERNEL``
would move its sha256, which stamps every watertap chemistry payload and the
minerals-active twin fingerprint, for phases that kernel's consumers never
touch.

Convention policy (2026-08 registry closure):

- The Ca/Mg phosphate set adopts minteq.v4.dat's PROTON-inclusive dissolution
  conventions verbatim (Hydroxylapatite, Brushite/DCPD, OCP, Newberyite,
  beta-TCP), with minteq's own delta_h — including its EXPLICIT ZERO entries,
  which are a database statement ("no temperature correction"), not a missing
  value. This retires the mixed OH-/HPO4-convention table whose constants
  could not be benchmarked against any database without conversion algebra.
- Calcite/aragonite carry the Plummer & Busenberg (1982) analytical fits
  (0-90 C primary source) rather than a single van't Hoff enthalpy.
- Struvite (Ohlinger 1998 Ksp; Bhuiyan 2007 dH), K-struvite (Taylor 1963,
  held) and amorphous AlPO4 (Stumm & Morgan, held) have NO minteq.v4
  counterpart: their PHREEQC benchmark validates OUR IMPLEMENTATION of the
  chosen constants, not the constants themselves.
- ACP is represented by minteq's Ca3(PO4)2(beta): amorphous calcium phosphate
  is MORE soluble than beta-TCP, so as a P-sink screen this surrogate is
  conservative toward over-predicting capture; flagged in the source string.
"""
from __future__ import annotations
from .definitions import ChemistryKernelDefinition, PhaseDefinition
from .ontology import ParticleClass, PermeabilityClass, Phase
from .vertical_slices import _species
_MINTEQ_SOURCE = 'minteq.v4.dat bundled snapshot; verify database hash at runtime'

def _solid(component_id: str, formula: str, elements: dict[str, float]):
    return _species(component_id, formula, elements, phase=Phase.SOLID, particle=ParticleClass.PARTICULATE, permeability=PermeabilityClass.RETAINED)
MINERAL_EXTENSION_COMPONENTS = (_species('H2O', 'H2O', {'H': 2, 'O': 1}), _species('H_plus', 'H+', {'H': 1}, charge=1), _species('Ca_ii', 'Ca+2', {'Ca': 1}, charge=2), _species('Mg_ii', 'Mg+2', {'Mg': 1}, charge=2), _species('K_i', 'K+', {'K': 1}, charge=1), _species('NH4', 'NH4+', {'N': 1, 'H': 4}, charge=1), _species('Al_iii', 'Al+3', {'Al': 1}, charge=3), _species('PO4', 'PO4-3', {'P': 1, 'O': 4}, charge=-3), _species('CO3', 'CO3-2', {'C': 1, 'O': 3}, charge=-2), _solid('Calcite', 'CaCO3', {'Ca': 1, 'C': 1, 'O': 3}), _solid('Aragonite', 'CaCO3', {'Ca': 1, 'C': 1, 'O': 3}), _solid('TCP_beta', 'Ca3(PO4)2', {'Ca': 3, 'P': 2, 'O': 8}), _solid('Hydroxylapatite', 'Ca5(PO4)3OH', {'Ca': 5, 'P': 3, 'O': 13, 'H': 1}), _solid('Brushite', 'CaHPO4:2H2O', {'Ca': 1, 'P': 1, 'O': 6, 'H': 5}), _solid('OCP', 'Ca4H(PO4)3:3H2O', {'Ca': 4, 'P': 3, 'O': 15, 'H': 7}), _solid('Struvite', 'MgNH4PO4:6H2O', {'Mg': 1, 'N': 1, 'P': 1, 'O': 10, 'H': 16}), _solid('Newberyite', 'MgHPO4:3H2O', {'Mg': 1, 'P': 1, 'O': 7, 'H': 7}), _solid('Magnesite', 'MgCO3', {'Mg': 1, 'C': 1, 'O': 3}), _solid('K_struvite', 'MgKPO4:6H2O', {'Mg': 1, 'K': 1, 'P': 1, 'O': 10, 'H': 12}), _solid('AlPO4_am', 'AlPO4', {'Al': 1, 'P': 1, 'O': 4}))
_PB82_CALCITE = (-171.9065, -0.077993, 2839.319, 71.595)
_PB82_ARAGONITE = (-171.9773, -0.077993, 2903.293, 71.595)
MINERAL_EXTENSION_PHASES = (PhaseDefinition('Calcite', 'CaCO3', {'Calcite': -1, 'Ca_ii': 1, 'CO3': 1}, -8.48, 'Plummer & Busenberg 1982 analytic fit (0-90 C); minteq.v4 Calcite (-8.48, dh -8 kJ) is the coarser database counterpart', 'stable_crystalline', temperature_correction_id='analytic_plummer_busenberg_1982_calcite', analytical_expression=_PB82_CALCITE), PhaseDefinition('Aragonite', 'CaCO3', {'Aragonite': -1, 'Ca_ii': 1, 'CO3': 1}, -8.336, 'Plummer & Busenberg 1982 analytic fit (0-90 C); minteq.v4 Aragonite (-8.3, dh -12 kJ) is the coarser database counterpart', 'metastable_crystalline', temperature_correction_id='analytic_plummer_busenberg_1982_aragonite', analytical_expression=_PB82_ARAGONITE), PhaseDefinition('TCP_beta', 'Ca3(PO4)2', {'TCP_beta': -1, 'Ca_ii': 3, 'PO4': 2}, -28.92, _MINTEQ_SOURCE + '; Ca3(PO4)2(beta) used as the ACP surrogate — amorphous Ca-P is MORE soluble, so this over-predicts capture (conservative as a P-sink screen)', 'metastable_crystalline', temperature_correction_id='minteq_v4_plus_54_kj_mol', delta_h_j_mol=54000.0), PhaseDefinition('Hydroxylapatite', 'Ca5(PO4)3OH', {'Hydroxylapatite': -1, 'H_plus': -1, 'Ca_ii': 5, 'PO4': 3, 'H2O': 1}, -44.333, _MINTEQ_SOURCE, 'stable_crystalline_nucleation_inhibited', temperature_correction_id='minteq_v4_explicit_zero_delta_h', delta_h_j_mol=0.0), PhaseDefinition('Brushite', 'CaHPO4:2H2O', {'Brushite': -1, 'Ca_ii': 1, 'H_plus': 1, 'PO4': 1, 'H2O': 2}, -18.995, _MINTEQ_SOURCE, 'metastable_crystalline', temperature_correction_id='minteq_v4_plus_23_kj_mol', delta_h_j_mol=23000.0), PhaseDefinition('OCP', 'Ca4H(PO4)3:3H2O', {'OCP': -1, 'Ca_ii': 4, 'H_plus': 1, 'PO4': 3, 'H2O': 3}, -47.08, _MINTEQ_SOURCE, 'metastable_crystalline', temperature_correction_id='minteq_v4_explicit_zero_delta_h', delta_h_j_mol=0.0), PhaseDefinition('Struvite', 'MgNH4PO4:6H2O', {'Struvite': -1, 'Mg_ii': 1, 'NH4': 1, 'PO4': 1, 'H2O': 6}, -13.26, 'Ohlinger et al. 1998 (Ksp); Bhuiyan et al. 2007 (dH, endothermic dissolution). NO minteq.v4 counterpart — the PHREEQC benchmark validates our implementation of these constants, not the constants', 'candidate_p_recovery_phase', temperature_correction_id='vanthoff_bhuiyan_2007_plus_22_6_kj_mol', delta_h_j_mol=22600.0), PhaseDefinition('Newberyite', 'MgHPO4:3H2O', {'Newberyite': -1, 'Mg_ii': 1, 'H_plus': 1, 'PO4': 1, 'H2O': 3}, -18.175, _MINTEQ_SOURCE, 'metastable_crystalline', temperature_correction_id='minteq_v4_explicit_zero_delta_h', delta_h_j_mol=0.0), PhaseDefinition('Magnesite', 'MgCO3', {'Magnesite': -1, 'Mg_ii': 1, 'CO3': 1}, -7.46, _MINTEQ_SOURCE, 'stable_crystalline_nucleation_inhibited', temperature_correction_id='minteq_v4_plus_20_kj_mol', delta_h_j_mol=20000.0), PhaseDefinition('K_struvite', 'MgKPO4:6H2O', {'K_struvite': -1, 'Mg_ii': 1, 'K_i': 1, 'PO4': 1, 'H2O': 6}, -11.0, 'Taylor et al. 1963 (approximate). HELD at 25 C — no defensible dH. NO minteq.v4 counterpart — the PHREEQC benchmark validates our implementation of the chosen constants, not the constants', 'candidate_p_recovery_phase'), PhaseDefinition('AlPO4_am', 'AlPO4', {'AlPO4_am': -1, 'Al_iii': 1, 'PO4': 1}, -20.5, 'Stumm & Morgan 1996 (approximate). HELD at 25 C — no defensible dH. NO minteq.v4 counterpart — the PHREEQC benchmark validates our implementation of the chosen constants, not the constants', 'metastable_amorphous'))
MINERAL_EXTENSION_DB_COUNTERPARTS = {'TCP_beta': 'Ca3(PO4)2(beta)', 'Hydroxylapatite': 'Hydroxylapatite', 'Brushite': 'CaHPO4:2H2O', 'OCP': 'Ca4H(PO4)3:3H2O', 'Newberyite': 'MgHPO4:3H2O', 'Magnesite': 'Magnesite'}
MINERAL_EXTENSION_25C_ANCHOR_COUNTERPARTS = {'Calcite': 'Calcite (minteq log_k -8.48 matches the P&B82 25 C anchor; the minteq dH -8 kJ temperature law is NOT ours)', 'Aragonite': 'Aragonite (minteq log_k -8.3 vs our P&B82 -8.336 — near, not equal; temperature law self-rendered)'}
MADM1_MINERAL_EXTENSION_KERNEL = ChemistryKernelDefinition(schema_version='1', kernel_id='madm1_mineral_extension_v1', components=MINERAL_EXTENSION_COMPONENTS, phases=MINERAL_EXTENSION_PHASES, reactions=(), thermodynamic_dataset_id='minteq.v4+primary_literature_2026_08', phase_definition_set_id='madm1_mineral_extension_phases_v1', activity_model_id='davies_reduced_v1')
__all__ = ['MADM1_MINERAL_EXTENSION_KERNEL', 'MINERAL_EXTENSION_COMPONENTS', 'MINERAL_EXTENSION_DB_COUNTERPARTS', 'MINERAL_EXTENSION_25C_ANCHOR_COUNTERPARTS', 'MINERAL_EXTENSION_PHASES']
