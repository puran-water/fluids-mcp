from __future__ import annotations
import hashlib
import json
from dataclasses import dataclass
from typing import Mapping
from .expression import Expression, add, multiply, power
from .madm1_biokinetics import ComponentContent, KineticParameterDefinition, _p, _s, _st, _symbols
from .madm1_sfe_kinetics import ELEMENT_SOURCE, MADM1_SFE_DERIVED_SYMBOLS, MADM1_SFE_KERNEL, MW_Fe, PhysicoChemicalProcessDefinition, crystallisation_rate
from .mineral_extension import MADM1_MINERAL_EXTENSION_KERNEL
from .vertical_slices import MADM1_FECL3_FE_S_P_KERNEL
RULING_ID = 'PUBLIC-RULE-6d7f2942e54b'
MACKINAWITE_RULING_ID = 'PUBLIC-RULE-566ce0eaa22d'
X_VIVIANITE = 'X_Fe3PO42'
X_SIDERITE = 'X_Siderite'
X_MACKINAWITE = 'X_Mackinawite'
STRUVITE_RULING_ID = 'PUBLIC-RULE-ad9bdd59a2d4'
X_STRUVITE = 'X_struv'
MW_P = float(MADM1_SFE_KERNEL.parameter_values()['MW_P'])
MW_O = 15.9994
MW_C = 12.0107
MW_Mg = 24.305
MW_H = 1.00794
MW_N = float(MADM1_SFE_KERNEL.parameter_values()['MW_N'])
VIVIANITE_RATE_SOURCE = "Flores-Alsina, Solon, Kazadi Mbamba, Tait, Gernaey, Jeppsson & Batstone 2016, Water Res. 95:370-382, doi:10.1016/j.watres.2016.03.012, Supplementary Information Table 5.2 'Fe3(PO4)2 precipitation': k_cryst 0.001 /h on ((a_Fe2)^3 (a_PO4)^2 / Ksp)^(1/5) - 1, squared ('kinetic parameters were selected with reference to Kazadi Mbamba et al. (2015a,b)'; read in the SI as appended to Solon 2017, PhD thesis, Lund University, lup.lub.lu.se/search/files/24482709/Thesis_Kimberly_Solon.pdf, 2026-10-22); the twin's default (QSDsan ModifiedADM1 k_cryst, citing the same paper). A model default at 35 C, not a fitted growth rate"
SIDERITE_RATE_SOURCE = "Jiang & Tosca 2020, Geochim. Cosmochim. Acta 274:97-117, doi:10.1016/j.gca.2020.01.047 (accepted manuscript, University of Cambridge repository, read 2026-10-22), Eq. 4: surface-reaction-controlled (spiral) growth of siderite seeds, r = 10^(-13.42 +/- 0.14) (Omega - 1)^(1.868 +/- 0.307) mol/m2/s, 298.15 K, I 0.175 molal, 7 < Omega < 18, log Ksp -10.24 (the shared kernel's own). Held at 298.15 K: no activation energy is sourced for this law (the rate at digester temperature is understated); dissolution reuses the growth constant (the twin's symmetric law; the paper measured growth)"
SIDERITE_AREA_SOURCE = "Jiang & Tosca 2020 s3.2.3 / Table 3: the unreacted synthetic seed's BET area, 109.75 m2/g, to which the Eq. 4 rates are normalised (harvested solids 81.1-82.5 m2/g). [ASSUMED] transfer: digester siderite carries the seed's specific area (owner question OQ-C77-2)"
SEED_SOURCE = "the twin's nucleation-surrogate seed inventory (mineral_seed_kg_m3), as mineral_seed_FeS"
MACKINAWITE_RATE_SOURCE = "the twin's mackinawite (QSDsan ModifiedADM1 process 67: k_cryst 100 /h, n 2, nu 2, FeS's own), a surrogate for a fast precipitation: Fe(II) and sulfide in real wastewater at pH 7 reach 60 % conversion within 30 s (Kiilerich, van de Ven, Nielsen & Vollertsen 2017, Water 9:670, doi:10.3390/w9090670, abstract), and FeS precipitated from Fe(II) and H2S/HS- at pH 7.1 ages to nanocrystalline mackinawite within 47 d (Le Pape et al. 2026, Eur. J. Mineral. 38:135, doi:10.5194/ejm-38-135-2026, abstract); the kinetic study itself (Rickard 1995, Geochim. Cosmochim. Acta 59:4367, doi:10.1016/0016-7037(95)00251-t) is closed access, not read. No crystallisation-law constant is published: k HRT >> 1 at every digester residence (2,400 /d x 20-40 d), so the steady answer is the phase's saturation, not k (bracketed 10-1000 /h, lane c90 evidence)"
MACKINAWITE_SOLUBILITY_SOURCE = "minteq.v4.dat as bundled (water-chemistry-mcp databases/official): Mackinawite log_k -3.6, delta_h 0, beside FeS(ppt) -2.95 (the S/Fe kernel's amorphous FeS); phreeqc.dat and wateq4f.dat carry Mackinawite -4.648 and FeS(ppt) -3.915, so -3.6 is the soluble side of the databases' range (the solubility studies, Davison 1991 Aquat. Sci. 53:309, doi:10.1007/BF00877139, and Rickard 2006 Geochim. Cosmochim. Acta 70:5779, doi:10.1016/j.gca.2006.02.029, are closed access, not read)"
STRUVITE_RATE_SOURCE = "Flores-Alsina, Solon, Kazadi Mbamba, Tait, Gernaey, Jeppsson & Batstone 2016, Water Res. 95:370-382, doi:10.1016/j.watres.2016.03.012, Supplementary Information Table 5.2 row 7 'Struvite precipitation': k_cryst 5 /h on ((a_NH4 a_Mg a_PO4) / Ksp)^(1/3) - 1, cubed (read on the rendered table, SI p. 184 as appended to Solon 2017, PhD thesis, Lund University, 2026-10-24; the table's pKsp 13.6 is not used: the solubility is the shared registry's); the twin's default (QSDsan ModifiedADM1 k_cryst 5.0, n_cryst 3, citing the same paper). A model default, not a growth rate fitted to a digester liquor"
_PHASES = {p.phase_id: p for p in MADM1_FECL3_FE_S_P_KERNEL.phases}
_MINERAL_PHASES = {p.phase_id: p for p in MADM1_MINERAL_EXTENSION_KERNEL.phases}

def _phase_source(phase_id: str) -> str:
    if phase_id in _MINERAL_PHASES:
        return f"{_MINERAL_PHASES[phase_id].source} (the shared mineral registry's {phase_id}, process_model_schema.mineral_extension {MADM1_MINERAL_EXTENSION_KERNEL.kernel_id})"
    return f"{_PHASES[phase_id].source} (the shared Fe-S-P kernel's {phase_id}, process_model_schema.vertical_slices {MADM1_FECL3_FE_S_P_KERNEL.kernel_id})"

def _mw(fe: float, p: float, c: float, o: float) -> float:
    return fe * MW_Fe + p * MW_P + c * MW_C + o * MW_O
MW_VIVIANITE = _mw(3, 2, 0, 8)
MW_SIDERITE = _mw(1, 0, 1, 3)
MW_S = float(MADM1_SFE_KERNEL.parameter_values()['MW_S'])
MW_MACKINAWITE = MW_Fe + MW_S
MW_STRUVITE = MW_Mg + MW_N + 4.0 * MW_H + MW_P + 4.0 * MW_O
FORMULA_MASS = {X_VIVIANITE: MW_VIVIANITE, X_SIDERITE: MW_SIDERITE, X_MACKINAWITE: MW_MACKINAWITE, X_STRUVITE: MW_STRUVITE}
MADM1_IRON_MINERAL_CONTENTS = (ComponentContent(X_VIVIANITE, 0.0, 0.0, 2.0 * MW_P / MW_VIVIANITE, 0.0, 0.0, 'Fe3(PO4)2'), ComponentContent(X_SIDERITE, MW_C / MW_SIDERITE, 0.0, 0.0, 0.0, 0.0, 'FeCO3'), ComponentContent(X_MACKINAWITE, 0.0, 0.0, 0.0, 0.0, 0.0, 'FeS'), ComponentContent(X_STRUVITE, 0.0, MW_N / MW_STRUVITE, MW_P / MW_STRUVITE, 0.0, MW_Mg / MW_STRUVITE, 'MgNH4PO4'))
MADM1_IRON_MINERAL_PARAMETERS: tuple[KineticParameterDefinition, ...] = (_p('MW_O', MW_O, 'kg/kmol', 'constant', ELEMENT_SOURCE), _p('MW_C', MW_C, 'kg/kmol', 'constant', ELEMENT_SOURCE), _p('k_cryst_vivianite', 0.001, '1/h', 'rate', VIVIANITE_RATE_SOURCE), _p('n_cryst_vivianite', 2.0, '1', 'rate', VIVIANITE_RATE_SOURCE), _p('nu_vivianite', 5.0, '1', 'constant', 'five product ions (3 Fe2+, 2 PO4-3): ' + VIVIANITE_RATE_SOURCE), _p('mineral_seed_vivianite', 1e-06, 'kg/m3', 'constant', SEED_SOURCE), _p('logK_vivianite', _PHASES['Vivianite'].log_k_25c, '1', 'constant', _phase_source('Vivianite')), _p('dH_vivianite', _PHASES['Vivianite'].delta_h_j_mol, 'J/mol', 'constant', _phase_source('Vivianite')), _p('k_growth_siderite', 10.0 ** (-13.42), 'mol/m2/s', 'rate', SIDERITE_RATE_SOURCE), _p('n_cryst_siderite', 1.868, '1', 'rate', SIDERITE_RATE_SOURCE), _p('nu_siderite', 1.0, '1', 'constant', 'Eq. 4 is in (Omega - 1): root 1 renders it exactly; ' + SIDERITE_RATE_SOURCE), _p('specific_area_siderite', 109750.0, 'm2/kg', 'constant', SIDERITE_AREA_SOURCE), _p('mineral_seed_siderite', 1e-06, 'kg/m3', 'constant', SEED_SOURCE), _p('logK_siderite', _PHASES['Siderite'].log_k_25c, '1', 'constant', _phase_source('Siderite')), _p('dH_siderite', _PHASES['Siderite'].delta_h_j_mol, 'J/mol', 'constant', _phase_source('Siderite')), _p('k_cryst_mackinawite', 100.0, '1/h', 'rate', MACKINAWITE_RATE_SOURCE), _p('n_cryst_mackinawite', 2.0, '1', 'rate', MACKINAWITE_RATE_SOURCE), _p('nu_mackinawite', 2.0, '1', 'constant', "two product ions (Fe2+, HS-), FeS's convention: " + MACKINAWITE_RATE_SOURCE), _p('mineral_seed_mackinawite', 1e-06, 'kg/m3', 'constant', SEED_SOURCE), _p('logK_mackinawite', _PHASES['Mackinawite'].log_k_25c, '1', 'constant', _phase_source('Mackinawite') + '; ' + MACKINAWITE_SOLUBILITY_SOURCE), _p('dH_mackinawite', _PHASES['Mackinawite'].delta_h_j_mol, 'J/mol', 'constant', _phase_source('Mackinawite')), _p('MW_Mg', MW_Mg, 'kg/kmol', 'constant', ELEMENT_SOURCE), _p('MW_H', MW_H, 'kg/kmol', 'constant', ELEMENT_SOURCE), _p('k_cryst_struvite', 5.0, '1/h', 'rate', STRUVITE_RATE_SOURCE), _p('n_cryst_struvite', 3.0, '1', 'rate', STRUVITE_RATE_SOURCE), _p('nu_struvite', 3.0, '1', 'constant', 'three product ions (Mg+2, NH4+, PO4-3): ' + STRUVITE_RATE_SOURCE), _p('mineral_seed_struvite', 1e-06, 'kg/m3', 'constant', SEED_SOURCE), _p('logK_struvite', _MINERAL_PHASES['Struvite'].log_k_25c, '1', 'constant', _phase_source('Struvite')), _p('dH_struvite', _MINERAL_PHASES['Struvite'].delta_h_j_mol, 'J/mol', 'constant', _phase_source('Struvite')))
MADM1_IRON_MINERAL_DERIVED_SYMBOLS = {'a_PO4': "PO4-3 activity (mol/kg water) in the engine's own speciation: its free HPO4-2 by the shared third dissociation (process_model_schema.speciation 'phosphoric_3')", 'a_CO3': "CO3-2 activity (mol/kg water) in the engine's own speciation (its free CO3-2)", 'Ksp_vivianite': "vivianite solubility constant at T (Fe3(PO4)2:8H2O = 3 Fe+2 + 2 PO4-3 + 8 H2O), van't Hoff", 'Ksp_siderite': "siderite solubility constant at T (FeCO3 = Fe+2 + CO3-2), van't Hoff", 'Ksp_mackinawite': "mackinawite solubility constant at T (FeS + H+ = Fe+2 + HS-), van't Hoff", 'a_Mg': "Mg+2 activity (mol/kg water) in the engine's own speciation (its free Mg+2, after the ion pairs)", 'a_NH4': "NH4+ activity (mol/kg water) in the engine's own speciation (its free NH4+)", 'Ksp_struvite': "struvite solubility constant at T (MgNH4PO4:6H2O = Mg+2 + NH4+ + PO4-3 + 6 H2O), van't Hoff"}

def _formula_mass(fe: float, p: float, c: float, o: float) -> Expression:
    return add(multiply(fe, _s('MW_Fe')), multiply(p, _s('MW_P')), multiply(c, _s('MW_C')), multiply(o, _s('MW_O')))
MADM1_VIVIANITE_PRECIPITATION = PhysicoChemicalProcessDefinition('precipitation_Fe3PO42', 68, crystallisation_rate(multiply(power(_s('a_Fe2'), 3.0), power(_s('a_PO4'), 2.0), power(_s('Ksp_vivianite'), -1.0)), _s('k_cryst_vivianite'), X_VIVIANITE, 'mineral_seed_vivianite', 'nu_vivianite', 'n_cryst_vivianite', magnitude_floor=True), X_VIVIANITE, {'S_Fe2': multiply(-3.0, _s('MW_Fe'), power(_formula_mass(3, 2, 0, 8), -1.0)), X_VIVIANITE: _st(x=1.0)['x']})
MADM1_SIDERITE_PRECIPITATION = PhysicoChemicalProcessDefinition('precipitation_Siderite', 69, crystallisation_rate(multiply(_s('a_Fe2'), _s('a_CO3'), power(_s('Ksp_siderite'), -1.0)), multiply(3600.0, _s('k_growth_siderite'), _s('specific_area_siderite'), multiply(0.001, _formula_mass(1, 0, 1, 3))), X_SIDERITE, 'mineral_seed_siderite', 'nu_siderite', 'n_cryst_siderite', magnitude_floor=True), X_SIDERITE, {'S_Fe2': multiply(-1.0, _s('MW_Fe'), power(_formula_mass(1, 0, 1, 3), -1.0)), X_SIDERITE: _st(x=1.0)['x']})
MADM1_MACKINAWITE_PRECIPITATION = PhysicoChemicalProcessDefinition('precipitation_Mackinawite', 67, crystallisation_rate(multiply(_s('a_Fe2'), _s('a_HS'), power(multiply(_s('H_ion'), _s('Ksp_mackinawite')), -1.0)), _s('k_cryst_mackinawite'), X_MACKINAWITE, 'mineral_seed_mackinawite', 'nu_mackinawite', 'n_cryst_mackinawite'), X_MACKINAWITE, {'S_IS': multiply(-1.0, _s('MW_S'), power(add(_s('MW_S'), _s('MW_Fe')), -1.0)), 'S_Fe2': multiply(-1.0, _s('MW_Fe'), power(add(_s('MW_S'), _s('MW_Fe')), -1.0)), X_MACKINAWITE: _st(x=1.0)['x']})
MADM1_STRUVITE_PRECIPITATION = PhysicoChemicalProcessDefinition('precipitation_struv', 62, crystallisation_rate(multiply(_s('a_Mg'), _s('a_NH4'), _s('a_PO4'), power(_s('Ksp_struvite'), -1.0)), _s('k_cryst_struvite'), X_STRUVITE, 'mineral_seed_struvite', 'nu_struvite', 'n_cryst_struvite'), X_STRUVITE, {'S_Mg': multiply(-1.0, _s('MW_Mg'), power(add(_s('MW_Mg'), _s('MW_N'), multiply(4.0, _s('MW_H')), _s('MW_P'), multiply(4.0, _s('MW_O'))), -1.0)), X_STRUVITE: _st(x=1.0)['x']})
MADM1_IRON_MINERAL_PROCESSES = (MADM1_VIVIANITE_PRECIPITATION, MADM1_SIDERITE_PRECIPITATION, MADM1_MACKINAWITE_PRECIPITATION, MADM1_STRUVITE_PRECIPITATION)
IRON_PER_FORMULA = {X_VIVIANITE: 3.0, X_SIDERITE: 1.0, X_MACKINAWITE: 1.0, X_STRUVITE: 0.0}
SULFUR_PER_FORMULA = {X_MACKINAWITE: 1.0}
MAGNESIUM_PER_FORMULA = {X_STRUVITE: 1.0}

@dataclass(frozen=True)
class IronMineralExtension:
    """The extension as one object, with an identity hash an engine stamps beside the base kernel's."""
    extension_id: str
    base_kernel_id: str
    base_kernel_sha256: str
    contents: tuple[ComponentContent, ...]
    parameters: tuple[KineticParameterDefinition, ...]
    processes: tuple[PhysicoChemicalProcessDefinition, ...]
    derived_symbols: Mapping[str, str]
    source: str

    @property
    def component_ids(self) -> tuple[str, ...]:
        return tuple((c.component_id for c in self.contents))

    def content(self, component_id: str) -> ComponentContent:
        return next((c for c in self.contents if c.component_id == component_id))

    def parameter_values(self) -> dict[str, float]:
        return {p.name: float(p.value) for p in self.parameters}

    @property
    def sha256(self) -> str:
        body = {'extension_id': self.extension_id, 'base': [self.base_kernel_id, self.base_kernel_sha256], 'contents': [sorted(vars(c).items()) for c in self.contents], 'parameters': [[p.name, p.value, p.unit, p.role, p.source] for p in self.parameters], 'processes': [[p.process_id, p.index, repr(p.rate_base), p.reference_component, sorted(((k, repr(v)) for k, v in p.stoichiometry.items()))] for p in self.processes], 'derived_symbols': sorted(self.derived_symbols), 'source': self.source}
        return hashlib.sha256(json.dumps(body, sort_keys=True, default=str).encode()).hexdigest()
MADM1_IRON_MINERAL_EXTENSION = IronMineralExtension(extension_id='madm1_iron_minerals_v3', base_kernel_id=MADM1_SFE_KERNEL.kernel_id, base_kernel_sha256=MADM1_SFE_KERNEL.sha256, contents=MADM1_IRON_MINERAL_CONTENTS, parameters=MADM1_IRON_MINERAL_PARAMETERS, processes=MADM1_IRON_MINERAL_PROCESSES, derived_symbols=MADM1_IRON_MINERAL_DERIVED_SYMBOLS, source=f'{RULING_ID} OQ-C64-1; {MACKINAWITE_RULING_ID} OQ-C77-1b; {STRUVITE_RULING_ID} OQ-C90-4; {VIVIANITE_RATE_SOURCE}; {SIDERITE_RATE_SOURCE}; {MACKINAWITE_RATE_SOURCE}; {STRUVITE_RATE_SOURCE}')
