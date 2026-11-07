"""Shared mADM1 biokinetic kernel: ONE definition of the ADM1 core's kinetics, stoichiometry and parameter set.

* **QSDsan** (``qsdsan-engine-mcp/models/madm1.py``): ``ModifiedADM1`` takes its default arguments from
  :data:`MADM1_BIOKINETICS_KERNEL` (:func:`kernel_parameter_values`). Its ``rhos_madm1`` stays the twin's
  vectorised numpy implementation; a parity test proves it equals this kernel's rate expressions, evaluated by
  ``process_chemistry_reference.madm1_biokinetics``, at common states.
* **WaterTAP** (``watertap-engine-mcp/chemistry/modified_adm1_kernel.py``): the house kernel reaction block compiles
  these rate expressions into Pyomo (``chemistry.pyomo_compiler``), takes its rate-reaction stoichiometry from
  :func:`process_chemistry_reference.madm1_biokinetics.stoichiometric_matrix` and its parameter values from the
  digester parameter set read through the process-model custody API.

Scope. The kernel is the ADM1 CORE that both engines carry state for: hydrolysis, the eight uptakes, the seven
decays, PHA storage and the three lyses (the twin's process indices 0-24, the WaterTAP reactions R1-R25). The twin's
sulfate-reducing, iron/HFO and mineral subsystems are NOT here: their constants already live in the shared Fe-S-P and
mineral-extension kernels (:mod:`.vertical_slices`, :mod:`.mineral_extension`), and the steady WaterTAP digester has
no state for them. A rate law here that reads free H2S therefore reads it as a SYMBOL: the twin supplies its
sulfide, a model without sulfide supplies 0 and the factor is exactly 1.

Conventions:

* state symbols are the mADM1 component ids, in kg (measured basis) per m3: kg COD for organics and biomass, kg C /
  N / P for S_IC / S_IN / S_IP, kg P for X_PP, kg K / Mg for S_K / S_Mg;
* rates are in kg (of the process's reference component) per m3 per DAY;
* the speciation the rates read is supplied as derived symbols (:data:`DERIVED_SYMBOLS`) by each engine's own
  acid-base solve: the hydrogen-ion activity, free NH3, free molecular H2S and sodium;
* the conserved inorganic carriers (S_IC, S_IN, S_IP) close each process's C, N and P balance from the component
  contents (:data:`MADM1_CORE_CONTENTS`), exactly as QSDsan's ``conserved_for=('C', 'N', 'P')`` does.

Trace-element quarantine (fix/madm1-trace-quarantine, 2026-07-21): the former KS_Co / KS_Ni / KS_Se "defaults" were
literature dose magnitudes shipped as Monod half-saturations. They have no place in this parameter set, and
:func:`validate_madm1_biokinetic_kernel` refuses any parameter whose name a quarantine entry names
(:data:`QUARANTINED_PARAMETER_NAMES`)."""
from __future__ import annotations
import hashlib
import json
import math
from dataclasses import dataclass, field, replace
from typing import Mapping, Optional
from .expression import Expression, Symbol, add, monod, multiply, named_factor, noncompetitive_inhibition, power
SCHEMA_VERSION = '1'
TWIN_SOURCE = 'Literature-based modified ADM1 kernel (QSDsan lineage); kinetic parameters are not site fits.'
ADM1_SOURCE = 'Batstone et al. 2002, IWA STR No. 13 (ADM1) ' + '; ' + TWIN_SOURCE
BSM2_SOURCE = 'Rosen & Jeppsson 2006, BSM2 ADM1 implementation (COD-balanced splits, Hill pH form); ' + TWIN_SOURCE
PEXT_SOURCE = 'Flores-Alsina et al. 2016 (mADM1 P extension; Water Res 95:370) ' + '; ' + TWIN_SOURCE
H2S_SOURCE = 'Flores-Alsina et al. 2016 mADM1 S extension (free molecular H2S, kg S/m3 basis since the 2026-08 units fix); ' + TWIN_SOURCE
NA_SOURCE = 'house literature-held surrogate KI_Na 0.3 mol/L, Chen, Cheng & Creamer 2008 doi:10.1016/j.biortech.2007.01.057; ' + TWIN_SOURCE
GAS_SOURCE = 'ADM1 gas transfer and Henry constants (Batstone et al. 2002; Rosen & Jeppsson 2006); ' + TWIN_SOURCE
PHREEQC_CARBONATE_SOURCE = 'USGS phreeqc.dat: CO3-2 + H+ = HCO3-, log_k 10.329, delta_h -3.561 kcal/mol (association; dissociation +14899.224 J/mol), as WaterTAP house phase 1 (chemistry/modified_adm1_phase1.py)'
PHOSPHATE_SOURCE = "ADM1-P phosphate second acid H2PO4-/HPO4-2, pKa2 7.20 at 25 C, dH 3600 J/mol (the twin's _temperature_corrected_acid_constants; WaterTAP house phase 1)"
DAVIES_SOURCE = 'Davies 1962, Ion Association: log g = -A z^2 (sqrt I/(1 + sqrt I) - 0.3 I), I molal; A(T) from Kell 1975 water density and Malmberg & Maryott 1956 dielectric constant (WaterTAP house Davies block)'

@dataclass(frozen=True)
class SourcedYield:
    """A kernel yield held with its band and the measurements behind it (never a typed literal)."""
    name: str
    value: float
    band: tuple[float, float]
    unit: str
    processes: tuple[str, ...]
    ruling_id: str
    sources: tuple[str, ...]
    basis: str

    def __post_init__(self) -> None:
        low, high = self.band
        if not (math.isfinite(self.value) and low <= self.value <= high):
            raise ValueError(f'{self.name} = {self.value} lies outside its sourced band {self.band}')
        if not self.sources:
            raise ValueError(f'{self.name} carries no source')

    @property
    def source(self) -> str:
        """The parameter's source text: the ruling, the basis, the band and every measurement."""
        return f"owner ruling {self.ruling_id}'; literature screening selection; '{self.basis}; band {self.band[0]:g}-{self.band[1]:g} {self.unit}; " + '; '.join(self.sources)
SOLUBLE_INERT_UPTAKE = SourcedYield(name='Y_sI_upt', value=0.02, band=(0.01, 0.04), unit='kg COD S_I/kg COD taken up', processes=('uptake_sugars', 'uptake_amino_acids'), ruling_id='PUBLIC-RULE-9290fc54cf13', basis="measured inert soluble microbial products of anaerobic digestion, about 2 % of the COD degraded, formed in acidogenesis; decay products go to particulates (f_sI_xb 0), so this is the kernel's only S_I source", sources=('Kasapgil Ince, Ince, Sallis & Anderson 2000, Water Res. 34(16):3943-3948, doi:10.1016/S0043-1354(00)00170-6, W1964472479: inert soluble microbial products 2.2 % of influent COD, not removed with extended operation (membrane anaerobic reactor, brewery; abstract)', 'Ince, Germirli Babuna, Kasapgil & Anderson 1998, Environ. Technol. 19(4):437-442, doi:10.1080/09593331908616700, W2071886588: 0.04 kg inert soluble COD per kg soluble COD utilised (abstract)', 'Aquino & Stuckey 2003, J. Environ. Eng. 129(11):1007, doi:10.1061/(asce)0733-9372(2003)129:11(1007), W1981679070: SMP/S0 averaged 3 % at steady state (abstract)', 'Aquino & Stuckey 2004, Water Res. 38:255-266, doi:10.1016/j.watres.2003.09.031, W1970793027: SMP/S0 2 % unspiked (abstract)', 'Barker & Stuckey 2001, Water Environ. Res. 73, doi:10.2175/106143001x138831, W166321314: SMP/S0 1.4-3 % in glucose CSTRs (abstract)', 'Kuo, Sneve & Parkin 1996, Water Environ. Res. 68(3), doi:10.2175/106143096x127712, W2256771364: SMP/S0 0.6-2.5 % glucose, 0.2-1.0 % acetate (abstract)', 'Noguera, Araki & Rittmann 1994, Biotechnol. Bioeng. 44:1040, doi:10.1002/bit.260440904, W2087497752: SMP formed in glucose acidogenesis, not in acetate methanogenesis (abstract)'))
QUARANTINED_PARAMETER_NAMES = frozenset({'KS_Co', 'KS_Ni', 'KS_Se', 'trace_active', 'f_bioavail_floor', 'f_bioavail_k_sulf'})
DERIVED_SYMBOLS = {'H_ion': 'hydrogen-ion activity 10^-pH (mol/L)', 'NH3_M': 'free ammonia (mol/L)', 'H2S_free_kg_S_m3': 'free molecular H2S (kg S/m3); 0 where the model carries no sulfide', 'Na_M': 'sodium (mol/L)', 'competition_floor_kg_m3': 'floor on the butyrate/valerate and PHA substrate-share denominators (kg COD/m3)', 'pao_ratio_floor_kg_m3': 'floor on the X_PP/(X_PP + K_PP X_PAO) denominator (kg P/m3)'}
MADM1_CORE_COMPONENTS = ('S_su', 'S_aa', 'S_fa', 'S_va', 'S_bu', 'S_pro', 'S_ac', 'S_h2', 'S_ch4', 'S_IC', 'S_IN', 'S_IP', 'S_I', 'X_ch', 'X_pr', 'X_li', 'X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2', 'X_I', 'X_PHA', 'X_PP', 'X_PAO', 'S_K', 'S_Mg')
CONSERVED_CARRIERS = {'C': 'S_IC', 'N': 'S_IN', 'P': 'S_IP'}

@dataclass(frozen=True)
class ComponentContent:
    """Element contents per unit of the component's measured basis (kg element per kg COD, or 1 for its own)."""
    component_id: str
    i_C: float
    i_N: float
    i_P: float
    i_K: float
    i_Mg: float
    measured_as: str

    def content(self, element: str) -> float:
        return float(getattr(self, f'i_{element}'))

def _contents(table: Mapping[str, tuple]) -> tuple[ComponentContent, ...]:
    return tuple((ComponentContent(cid, *values) for cid, values in table.items()))
MADM1_TWIN_CONTENTS = _contents({'S_su': (0.3753484505668963, 0.0, 0.0, 0.0, 0.0, 'COD'), 'S_aa': (0.3689, 0.11065, 0.0, 0.0, 0.0, 'COD'), 'S_fa': (0.2570879798403399, 0.0, 0.0, 0.0, 0.0, 'COD'), 'S_va': (0.288729577359151, 0.0, 0.0, 0.0, 0.0, 'COD'), 'S_bu': (0.300278760453517, 0.0, 0.0, 0.0, 0.0, 'COD'), 'S_pro': (0.32172724334305397, 0.0, 0.0, 0.0, 0.0, 'COD'), 'S_ac': (0.3753484505668962, 0.0, 0.0, 0.0, 0.0, 'COD'), 'S_h2': (0.0, 0.0, 0.0, 0.0, 0.0, 'COD'), 'S_ch4': (0.18767422528344815, 0.0, 0.0, 0.0, 0.0, 'COD'), 'S_IC': (1.0, 0.0, 0.0, 0.0, 0.0, 'C'), 'S_IN': (0.0, 1.0, 0.0, 0.0, 0.0, 'N'), 'S_IP': (0.0, 0.0, 1.0, 0.0, 0.0, 'P'), 'S_I': (0.36178, 0.06003, 0.00649, 0.0, 0.0, 'COD'), 'X_ch': (0.3753484505668963, 0.0, 0.0, 0.0, 0.0, 'COD'), 'X_pr': (0.3689, 0.11065, 0.0, 0.0, 0.0, 'COD'), 'X_li': (0.26311, 0.0, 0.01067, 0.0, 0.0, 'COD'), 'X_su': (0.36612, 0.08615, 0.02154, 0.004768, 0.002964, 'COD'), 'X_aa': (0.36612, 0.08615, 0.02154, 0.004768, 0.002964, 'COD'), 'X_fa': (0.36612, 0.08615, 0.02154, 0.004768, 0.002964, 'COD'), 'X_c4': (0.36612, 0.08615, 0.02154, 0.004768, 0.002964, 'COD'), 'X_pro': (0.36612, 0.08615, 0.02154, 0.004768, 0.002964, 'COD'), 'X_ac': (0.36612, 0.08615, 0.02154, 0.004768, 0.002964, 'COD'), 'X_h2': (0.36612, 0.08615, 0.02154, 0.004768, 0.002964, 'COD'), 'X_I': (0.36178, 0.06003, 0.00649, 0.0, 0.0, 'COD'), 'X_PHA': (0.300278760453517, 0.0, 0.0, 0.0, 0.0, 'COD'), 'X_PP': (0.0, 0.0, 1.0, 0.416560294373034, 0.25894982530536087, 'P'), 'X_PAO': (0.36612, 0.08615, 0.02154, 0.004768, 0.002964, 'COD'), 'S_K': (0.0, 0.0, 0.0, 1.0, 0.0, 'K'), 'S_Mg': (0.0, 0.0, 0.0, 0.0, 1.0, 'Mg')})
DESIGN_SOLUBLE_INERT_N_CONTENT = 0.01
DESIGN_CONTENT_DEPARTURES: Mapping[tuple[str, str], str] = {('S_I', 'N'): "owner ruling PUBLIC-RULE-f116039548ff (2026-10-15): ASM2d i_N,SI 0.01 kg N/kg COD (Henze et al. 1999), the walk's content; the twin keeps 0.06003 until its ruled migration (PUBLIC-RULE-e6bb9e579ff3)"}

def _design_contents(twin: tuple[ComponentContent, ...]) -> tuple[ComponentContent, ...]:
    values = {('S_I', 'N'): DESIGN_SOLUBLE_INERT_N_CONTENT}
    if set(values) != set(DESIGN_CONTENT_DEPARTURES):
        raise ValueError('every design content departure is declared with its ruling (DESIGN_CONTENT_DEPARTURES)')
    out = []
    for content in twin:
        moved = {f'i_{element}': value for (cid, element), value in values.items() if cid == content.component_id}
        out.append(replace(content, **moved) if moved else content)
    return tuple(out)
MADM1_CORE_CONTENTS = _design_contents(MADM1_TWIN_CONTENTS)

@dataclass(frozen=True)
class KineticParameterDefinition:
    name: str
    value: float
    unit: str
    role: str
    source: str

    def __post_init__(self) -> None:
        if not math.isfinite(float(self.value)):
            raise ValueError(f'parameter {self.name!r} must be finite')

def _p(name: str, value: float, unit: str, role: str, source: str) -> KineticParameterDefinition:
    return KineticParameterDefinition(name, value, unit, role, source)
MADM1_CORE_PARAMETERS: tuple[KineticParameterDefinition, ...] = (_p('f_ch_xb', 0.275, '1', 'fraction', PEXT_SOURCE), _p('f_pr_xb', 0.275, '1', 'fraction', PEXT_SOURCE), _p('f_li_xb', 0.35, '1', 'fraction', PEXT_SOURCE), _p('f_xI_xb', 0.1, '1', 'fraction', PEXT_SOURCE), _p('f_fa_li', 0.95, '1', 'fraction', ADM1_SOURCE), _p('f_bu_su', 0.1328, '1', 'fraction', BSM2_SOURCE), _p('f_pro_su', 0.2691, '1', 'fraction', BSM2_SOURCE), _p('f_ac_su', 0.4076, '1', 'fraction', BSM2_SOURCE), _p('f_va_aa', 0.23, '1', 'fraction', ADM1_SOURCE), _p('f_bu_aa', 0.26, '1', 'fraction', ADM1_SOURCE), _p('f_pro_aa', 0.05, '1', 'fraction', ADM1_SOURCE), _p('f_ac_aa', 0.4, '1', 'fraction', ADM1_SOURCE), _p('f_ac_fa', 0.7, '1', 'fraction', ADM1_SOURCE), _p('f_pro_va', 0.54, '1', 'fraction', ADM1_SOURCE), _p('f_ac_va', 0.31, '1', 'fraction', ADM1_SOURCE), _p('f_ac_bu', 0.8, '1', 'fraction', ADM1_SOURCE), _p('f_ac_pro', 0.57, '1', 'fraction', ADM1_SOURCE), _p('Y_su', 0.1, 'kg COD/kg COD', 'yield', ADM1_SOURCE), _p('Y_aa', 0.08, 'kg COD/kg COD', 'yield', ADM1_SOURCE), _p('Y_fa', 0.06, 'kg COD/kg COD', 'yield', ADM1_SOURCE), _p('Y_c4', 0.06, 'kg COD/kg COD', 'yield', ADM1_SOURCE), _p('Y_pro', 0.04, 'kg COD/kg COD', 'yield', ADM1_SOURCE), _p('Y_ac', 0.05, 'kg COD/kg COD', 'yield', ADM1_SOURCE), _p('Y_h2', 0.06, 'kg COD/kg COD', 'yield', ADM1_SOURCE), _p(SOLUBLE_INERT_UPTAKE.name, SOLUBLE_INERT_UPTAKE.value, 'kg COD/kg COD', 'yield', SOLUBLE_INERT_UPTAKE.source), _p('f_va_pha', 0.1, '1', 'fraction', PEXT_SOURCE), _p('f_bu_pha', 0.1, '1', 'fraction', PEXT_SOURCE), _p('f_pro_pha', 0.4, '1', 'fraction', PEXT_SOURCE), _p('Y_PO4', 0.4, 'kg P/kg COD', 'yield', PEXT_SOURCE), _p('q_ch_hyd', 10.0, '1/d', 'rate_constant', ADM1_SOURCE), _p('q_pr_hyd', 10.0, '1/d', 'rate_constant', ADM1_SOURCE), _p('q_li_hyd', 10.0, '1/d', 'rate_constant', ADM1_SOURCE), _p('k_su', 30.0, '1/d', 'rate_constant', ADM1_SOURCE), _p('k_aa', 50.0, '1/d', 'rate_constant', ADM1_SOURCE), _p('k_fa', 6.0, '1/d', 'rate_constant', ADM1_SOURCE), _p('k_c4', 20.0, '1/d', 'rate_constant', ADM1_SOURCE), _p('k_pro', 13.0, '1/d', 'rate_constant', ADM1_SOURCE), _p('k_ac', 8.0, '1/d', 'rate_constant', ADM1_SOURCE), _p('k_h2', 35.0, '1/d', 'rate_constant', ADM1_SOURCE), _p('K_su', 0.5, 'kg COD/m3', 'half_saturation', ADM1_SOURCE), _p('K_aa', 0.3, 'kg COD/m3', 'half_saturation', ADM1_SOURCE), _p('K_fa', 0.4, 'kg COD/m3', 'half_saturation', ADM1_SOURCE), _p('K_c4', 0.2, 'kg COD/m3', 'half_saturation', ADM1_SOURCE), _p('K_pro', 0.1, 'kg COD/m3', 'half_saturation', ADM1_SOURCE), _p('K_ac', 0.15, 'kg COD/m3', 'half_saturation', ADM1_SOURCE), _p('K_h2', 7e-06, 'kg COD/m3', 'half_saturation', ADM1_SOURCE), _p('b_su', 0.02, '1/d', 'rate_constant', ADM1_SOURCE), _p('b_aa', 0.02, '1/d', 'rate_constant', ADM1_SOURCE), _p('b_fa', 0.02, '1/d', 'rate_constant', ADM1_SOURCE), _p('b_c4', 0.02, '1/d', 'rate_constant', ADM1_SOURCE), _p('b_pro', 0.02, '1/d', 'rate_constant', ADM1_SOURCE), _p('b_ac', 0.02, '1/d', 'rate_constant', ADM1_SOURCE), _p('b_h2', 0.02, '1/d', 'rate_constant', ADM1_SOURCE), _p('q_pha', 3.0, '1/d', 'rate_constant', PEXT_SOURCE), _p('b_pao', 0.2, '1/d', 'rate_constant', PEXT_SOURCE), _p('b_pp', 0.2, '1/d', 'rate_constant', PEXT_SOURCE), _p('b_pha', 0.2, '1/d', 'rate_constant', PEXT_SOURCE), _p('K_A', 0.004, 'kg COD/m3', 'half_saturation', PEXT_SOURCE), _p('K_PP', 0.01, 'kg P (X_PP)/kg COD (X_PAO)', 'half_saturation', PEXT_SOURCE), _p('KI_h2_fa', 5e-06, 'kg COD/m3', 'inhibition', ADM1_SOURCE), _p('KI_h2_c4', 1e-05, 'kg COD/m3', 'inhibition', ADM1_SOURCE), _p('KI_h2_pro', 3.5e-06, 'kg COD/m3', 'inhibition', ADM1_SOURCE), _p('KI_nh3', 0.0018, 'mol/L', 'inhibition', ADM1_SOURCE + '; literature sensitivity remains unqualified for a specific plant'), _p('n_nh3', 1.0, '1', 'inhibition', '1 = the ADM1 non-competitive form; ' + TWIN_SOURCE), _p('KS_IN', 0.0001, 'mol N/L', 'half_saturation', ADM1_SOURCE), _p('KS_IP', 2e-05, 'mol P/L', 'half_saturation', PEXT_SOURCE), _p('KI_h2s_c4', 0.481, 'kg S/m3', 'inhibition', H2S_SOURCE), _p('KI_h2s_pro', 0.481, 'kg S/m3', 'inhibition', H2S_SOURCE), _p('KI_h2s_ac', 0.46, 'kg S/m3', 'inhibition', H2S_SOURCE), _p('KI_h2s_h2', 0.4, 'kg S/m3', 'inhibition', H2S_SOURCE), _p('KI_Na', 0.3, 'mol/L', 'inhibition', NA_SOURCE), _p('pH_LL_aa', 4.0, 'pH', 'ph_limit', BSM2_SOURCE), _p('pH_UL_aa', 5.5, 'pH', 'ph_limit', BSM2_SOURCE), _p('pH_LL_ac', 6.0, 'pH', 'ph_limit', BSM2_SOURCE), _p('pH_UL_ac', 7.0, 'pH', 'ph_limit', BSM2_SOURCE), _p('pH_LL_h2', 5.0, 'pH', 'ph_limit', BSM2_SOURCE), _p('pH_UL_h2', 6.0, 'pH', 'ph_limit', BSM2_SOURCE), _p('MW_N', 14.0067, 'kg/kmol', 'constant', "chemicals.elements.molecular_weight({'N': 1}) as the twin imports it"), _p('MW_P', 30.973761, 'kg/kmol', 'constant', "chemicals.elements.molecular_weight({'P': 1}) as the twin imports it"), _p('kLa', 200.0, '1/d', 'gas', GAS_SOURCE), _p('KH_h2', 0.00078, 'mol/L/bar', 'gas', GAS_SOURCE), _p('KH_ch4', 0.0014, 'mol/L/bar', 'gas', GAS_SOURCE), _p('KH_co2', 0.035, 'mol/L/bar', 'gas', GAS_SOURCE), _p('KH_h2s', 0.105, 'mol/L/bar', 'gas', H2S_SOURCE), _p('dH_KH_h2', -4180.0, 'J/mol', 'gas', GAS_SOURCE), _p('dH_KH_ch4', -14240.0, 'J/mol', 'gas', GAS_SOURCE), _p('dH_KH_co2', -19410.0, 'J/mol', 'gas', GAS_SOURCE), _p('dH_KH_h2s', -19180.0, 'J/mol', 'gas', H2S_SOURCE), _p('pKa_w', 14.0, '1', 'acid_base', ADM1_SOURCE), _p('pKa_nh4', 9.25, '1', 'acid_base', ADM1_SOURCE), _p('pKa_co2', 6.35, '1', 'acid_base', ADM1_SOURCE), _p('pKa_ac', 4.76, '1', 'acid_base', ADM1_SOURCE), _p('pKa_pro', 4.88, '1', 'acid_base', ADM1_SOURCE), _p('pKa_bu', 4.82, '1', 'acid_base', ADM1_SOURCE), _p('pKa_va', 4.86, '1', 'acid_base', ADM1_SOURCE), _p('dH_Ka_w', 55900.0, 'J/mol', 'acid_base', ADM1_SOURCE), _p('dH_Ka_nh4', 51965.0, 'J/mol', 'acid_base', ADM1_SOURCE), _p('dH_Ka_co2', 7646.0, 'J/mol', 'acid_base', ADM1_SOURCE), _p('dH_Ka_ac', 0.0, 'J/mol', 'acid_base', ADM1_SOURCE), _p('dH_Ka_pro', 0.0, 'J/mol', 'acid_base', ADM1_SOURCE), _p('dH_Ka_bu', 0.0, 'J/mol', 'acid_base', ADM1_SOURCE), _p('dH_Ka_va', 0.0, 'J/mol', 'acid_base', ADM1_SOURCE), _p('pKa_co3', 10.329, '1', 'acid_base', PHREEQC_CARBONATE_SOURCE), _p('dH_Ka_co3', 14899.224, 'J/mol', 'acid_base', PHREEQC_CARBONATE_SOURCE), _p('pKa_h2po4', 7.2, '1', 'acid_base', PHOSPHATE_SOURCE), _p('dH_Ka_h2po4', 3600.0, 'J/mol', 'acid_base', PHOSPHATE_SOURCE), _p('R_gas', 8.314462618, 'J/mol/K', 'constant', 'CODATA 2018 molar gas constant (IDAES Constants.gas_constant)'), _p('davies_b', 0.3, '1', 'acid_base', DAVIES_SOURCE), _p('T_base', 298.15, 'K', 'constant', ADM1_SOURCE))

def _s(name: str) -> Symbol:
    return Symbol(name)

def _hill_ph(window: str) -> Expression:
    """BSM2 Hill pH inhibition 1/(1 + (H/K)^n), n = 3/(UL-LL), K = 10^-(UL+LL)/2, as the twin's Hill_inhibit."""
    upper, lower = (_s(f'pH_UL_{window}'), _s(f'pH_LL_{window}'))
    half_point = power(10.0, multiply(0.5, add(upper, lower)))
    exponent = multiply(3.0, power(add(upper, multiply(-1.0, lower)), -1.0))
    return named_factor(f'I_pH_{window}', noncompetitive_inhibition(power(multiply(_s('H_ion'), half_point), exponent), 1.0))

def _nutrients() -> Expression:
    return multiply(named_factor('I_IN_lim', monod(_s('S_IN'), multiply(_s('KS_IN'), _s('MW_N')))), named_factor('I_IP_lim', monod(_s('S_IP'), multiply(_s('KS_IP'), _s('MW_P')))))

def _h2(name: str) -> Expression:
    return named_factor(f'I_{name.removeprefix('KI_')}', noncompetitive_inhibition(_s('S_h2'), _s(name)))

def _h2s(name: str) -> Expression:
    return named_factor(f'I_{name.removeprefix('KI_')}', noncompetitive_inhibition(_s('H2S_free_kg_S_m3'), _s(name)))

def _nh3() -> Expression:
    ratio = multiply(_s('NH3_M'), power(_s('KI_nh3'), -1.0))
    return named_factor('I_nh3', noncompetitive_inhibition(power(ratio, _s('n_nh3')), 1.0))

def _na() -> Expression:
    return named_factor('I_Na', noncompetitive_inhibition(_s('Na_M'), _s('KI_Na')))

def _uptake(rate: str, biomass: str, substrate: str, half: str, *competition: Expression) -> Expression:
    """The uninhibited uptake: k X Monod(S) [x substrate-competition share]."""
    return multiply(_s(rate), _s(biomass), monod(_s(substrate), _s(half)), *competition)

def _inhibition(*factors: Expression) -> Expression:
    """The product of a process's inhibition and limitation factors (each in [0, 1])."""
    return multiply(*factors)
_VFA = ('S_va', 'S_bu', 'S_pro', 'S_ac')

def _pha_storage(vfa: str) -> Expression:
    others = add(*(_s(v) for v in _VFA if v != vfa), _s('competition_floor_kg_m3'))
    return multiply(_s('q_pha'), _s('X_PAO'), monod(_s(vfa), _s('K_A')), monod(_s(vfa), others), monod(_s('X_PP'), add(multiply(_s('K_PP'), _s('X_PAO')), _s('pao_ratio_floor_kg_m3'))))

def _c4_share(substrate: str, competitor: str) -> Expression:
    """The butyrate/valerate share of the c4 degraders: S / (S + S_other [+ floor])."""
    return monod(_s(substrate), add(_s(competitor), _s('competition_floor_kg_m3')))

@dataclass(frozen=True)
class BiokineticProcessDefinition:
    """One process: its id (the twin's TSV row), its twin rhos index, rate law and explicit stoichiometry.

    ``stoichiometry`` holds every coefficient the TSV states (a ``Symbol`` of a parameter, or a derived-parameter
    name: see :data:`DERIVED_PARAMETERS`); the conserved carriers (S_IC, S_IN, S_IP) are closed from the component
    contents, never stated.
    """
    process_id: str
    index: int
    rate_base: Expression
    reference_component: str
    stoichiometry: Mapping[str, Expression]
    inhibition: Optional[Expression] = None

    @property
    def rate(self) -> Expression:
        """rate_base x inhibition: the process rate, kg reference component / m3 / d."""
        return self.rate_base if self.inhibition is None else multiply(self.rate_base, self.inhibition)

def _st(**coefficients) -> dict[str, Expression]:
    from .expression import as_expression
    return {k: _s(v) if isinstance(v, str) else as_expression(v) for k, v in coefficients.items()}

def _yield_split(yield_name: str, fraction: str, inert_yield: Optional[str]=None) -> Expression:
    """(1 - Y [- Y_sI]) x f: the share of the COD taken up that goes to one product."""
    taken = [multiply(-1.0, _s(yield_name))]
    if inert_yield is not None:
        taken.append(multiply(-1.0, _s(inert_yield)))
    return multiply(add(1.0, *taken), _s(fraction))
_DECAY_PRODUCTS = dict(S_I='f_sI_xb', X_ch='f_ch_xb', X_pr='f_pr_xb', X_li='f_li_xb', X_I='f_xI_xb')

def _decay(process_id: str, index: int, rate: str, biomass: str) -> BiokineticProcessDefinition:
    return BiokineticProcessDefinition(process_id, index, multiply(_s(rate), _s(biomass)), biomass, {**_st(**_DECAY_PRODUCTS), biomass: _st(x=-1.0)['x']})

def _uptake_row(process_id: str, index: int, rate: Expression, inhibition: Expression, substrate: str, biomass: str, yield_name: str, products: Mapping[str, str]) -> BiokineticProcessDefinition:
    stoich = {substrate: _st(x=-1.0)['x'], biomass: _s(yield_name)}
    inert = SOLUBLE_INERT_UPTAKE.name if process_id in SOLUBLE_INERT_UPTAKE.processes else None
    if inert is not None:
        stoich['S_I'] = _s(inert)
    for product, fraction in products.items():
        stoich[product] = _yield_split(yield_name, fraction, inert)
    return BiokineticProcessDefinition(process_id, index, rate, substrate, stoich, inhibition)

def _pha_row(process_id: str, index: int, vfa: str) -> BiokineticProcessDefinition:
    return BiokineticProcessDefinition(process_id, index, _pha_storage(vfa), vfa, {vfa: _st(x=-1.0)['x'], 'X_PHA': _st(x=1.0)['x'], 'X_PP': multiply(-1.0, _s('Y_PO4')), 'S_K': multiply(_s('Y_PO4'), _s('K_XPP')), 'S_Mg': multiply(_s('Y_PO4'), _s('Mg_XPP'))})
_NUT = _nutrients()
MADM1_CORE_PROCESSES: tuple[BiokineticProcessDefinition, ...] = (BiokineticProcessDefinition('hydrolysis_carbs', 0, multiply(_s('q_ch_hyd'), _s('X_ch')), 'X_ch', _st(S_su=1.0, X_ch=-1.0)), BiokineticProcessDefinition('hydrolysis_proteins', 1, multiply(_s('q_pr_hyd'), _s('X_pr')), 'X_pr', _st(S_aa=1.0, X_pr=-1.0)), BiokineticProcessDefinition('hydrolysis_lipids', 2, multiply(_s('q_li_hyd'), _s('X_li')), 'X_li', {'S_su': add(1.0, multiply(-1.0, _s('f_fa_li'))), 'S_fa': _s('f_fa_li'), 'X_li': _st(x=-1.0)['x']}), _uptake_row('uptake_sugars', 3, _uptake('k_su', 'X_su', 'S_su', 'K_su'), _inhibition(_NUT, _hill_ph('aa')), 'S_su', 'X_su', 'Y_su', {'S_bu': 'f_bu_su', 'S_pro': 'f_pro_su', 'S_ac': 'f_ac_su', 'S_h2': 'f_h2_su'}), _uptake_row('uptake_amino_acids', 4, _uptake('k_aa', 'X_aa', 'S_aa', 'K_aa'), _inhibition(_NUT, _hill_ph('aa')), 'S_aa', 'X_aa', 'Y_aa', {'S_va': 'f_va_aa', 'S_bu': 'f_bu_aa', 'S_pro': 'f_pro_aa', 'S_ac': 'f_ac_aa', 'S_h2': 'f_h2_aa'}), _uptake_row('uptake_LCFA', 5, _uptake('k_fa', 'X_fa', 'S_fa', 'K_fa'), _inhibition(_NUT, _hill_ph('aa'), _h2('KI_h2_fa')), 'S_fa', 'X_fa', 'Y_fa', {'S_ac': 'f_ac_fa', 'S_h2': 'f_h2_fa'}), _uptake_row('uptake_valerate', 6, _uptake('k_c4', 'X_c4', 'S_va', 'K_c4', _c4_share('S_va', 'S_bu')), _inhibition(_NUT, _hill_ph('aa'), _h2('KI_h2_c4'), _h2s('KI_h2s_c4')), 'S_va', 'X_c4', 'Y_c4', {'S_pro': 'f_pro_va', 'S_ac': 'f_ac_va', 'S_h2': 'f_h2_va'}), _uptake_row('uptake_butyrate', 7, _uptake('k_c4', 'X_c4', 'S_bu', 'K_c4', _c4_share('S_bu', 'S_va')), _inhibition(_NUT, _hill_ph('aa'), _h2('KI_h2_c4'), _h2s('KI_h2s_c4')), 'S_bu', 'X_c4', 'Y_c4', {'S_ac': 'f_ac_bu', 'S_h2': 'f_h2_bu'}), _uptake_row('uptake_propionate', 8, _uptake('k_pro', 'X_pro', 'S_pro', 'K_pro'), _inhibition(_NUT, _hill_ph('aa'), _h2('KI_h2_pro'), _h2s('KI_h2s_pro')), 'S_pro', 'X_pro', 'Y_pro', {'S_ac': 'f_ac_pro', 'S_h2': 'f_h2_pro'}), BiokineticProcessDefinition('uptake_acetate', 9, _uptake('k_ac', 'X_ac', 'S_ac', 'K_ac'), 'S_ac', {'S_ac': _st(x=-1.0)['x'], 'S_ch4': add(1.0, multiply(-1.0, _s('Y_ac'))), 'X_ac': _s('Y_ac')}, _inhibition(_NUT, _hill_ph('ac'), _nh3(), _h2s('KI_h2s_ac'), _na())), BiokineticProcessDefinition('uptake_h2', 10, _uptake('k_h2', 'X_h2', 'S_h2', 'K_h2'), 'S_h2', {'S_h2': _st(x=-1.0)['x'], 'S_ch4': add(1.0, multiply(-1.0, _s('Y_h2'))), 'X_h2': _s('Y_h2')}, _inhibition(_NUT, _hill_ph('h2'), _h2s('KI_h2s_h2'), _na())), _decay('decay_Xsu', 11, 'b_su', 'X_su'), _decay('decay_Xaa', 12, 'b_aa', 'X_aa'), _decay('decay_Xfa', 13, 'b_fa', 'X_fa'), _decay('decay_Xc4', 14, 'b_c4', 'X_c4'), _decay('decay_Xpro', 15, 'b_pro', 'X_pro'), _decay('decay_Xac', 16, 'b_ac', 'X_ac'), _decay('decay_Xh2', 17, 'b_h2', 'X_h2'), _pha_row('storage_PHA_valerate', 18, 'S_va'), _pha_row('storage_PHA_butyrate', 19, 'S_bu'), _pha_row('storage_PHA_propionate', 20, 'S_pro'), _pha_row('storage_PHA_actate', 21, 'S_ac'), _decay('lysis_XPAO', 22, 'b_pao', 'X_PAO'), BiokineticProcessDefinition('lysis_XPP', 23, multiply(_s('b_pp'), _s('X_PP')), 'X_PP', {'X_PP': _st(x=-1.0)['x'], 'S_K': _s('K_XPP'), 'S_Mg': _s('Mg_XPP')}), BiokineticProcessDefinition('lysis_XPHA', 24, multiply(_s('b_pha'), _s('X_PHA')), 'X_PHA', {'S_va': _s('f_va_pha'), 'S_bu': _s('f_bu_pha'), 'S_pro': _s('f_pro_pha'), 'S_ac': _s('f_ac_pha'), 'X_PHA': _st(x=-1.0)['x']}))
DERIVED_PARAMETERS = ('f_h2_su', 'f_h2_aa', 'f_h2_fa', 'f_h2_va', 'f_h2_bu', 'f_h2_pro', 'f_sI_xb', 'f_ac_pha', 'K_XPP', 'Mg_XPP')

def derived_parameter_values(values: Mapping[str, float]) -> dict[str, float]:
    """The derived stoichiometric parameters, by exactly the twin's arithmetic."""
    v = values
    xpp = next((c for c in MADM1_CORE_CONTENTS if c.component_id == 'X_PP'))
    return {'f_h2_su': 1 - v['f_bu_su'] - v['f_pro_su'] - v['f_ac_su'], 'f_h2_aa': 1 - v['f_va_aa'] - v['f_bu_aa'] - v['f_pro_aa'] - v['f_ac_aa'], 'f_h2_fa': 1 - v['f_ac_fa'], 'f_h2_va': 1 - v['f_pro_va'] - v['f_ac_va'], 'f_h2_bu': 1 - v['f_ac_bu'], 'f_h2_pro': 1 - v['f_ac_pro'], 'f_sI_xb': round(1.0 - v['f_ch_xb'] - v['f_pr_xb'] - v['f_li_xb'] - v['f_xI_xb'], 4), 'f_ac_pha': 1 - v['f_va_pha'] - v['f_bu_pha'] - v['f_pro_pha'], 'K_XPP': xpp.i_K, 'Mg_XPP': xpp.i_Mg}

@dataclass(frozen=True)
class BiokineticKernelDefinition:
    schema_version: str
    kernel_id: str
    components: tuple[ComponentContent, ...]
    parameters: tuple[KineticParameterDefinition, ...]
    processes: tuple[BiokineticProcessDefinition, ...]
    conserved_carriers: Mapping[str, str] = field(default_factory=lambda: dict(CONSERVED_CARRIERS))
    derived_symbols: Mapping[str, str] = field(default_factory=lambda: dict(DERIVED_SYMBOLS))
    rate_unit: str = 'kg reference component / m3 / d'
    reference: str = 'Literature-based modified ADM1 kernel (QSDsan lineage); kinetic parameters are not site fits.'

    @property
    def component_ids(self) -> tuple[str, ...]:
        return tuple((c.component_id for c in self.components))

    @property
    def process_ids(self) -> tuple[str, ...]:
        return tuple((p.process_id for p in self.processes))

    def parameter_values(self) -> dict[str, float]:
        return {p.name: float(p.value) for p in self.parameters}

    def parameter(self, name: str) -> KineticParameterDefinition:
        for p in self.parameters:
            if p.name == name:
                return p
        raise KeyError(name)

    def content(self, component_id: str) -> ComponentContent:
        for c in self.components:
            if c.component_id == component_id:
                return c
        raise KeyError(component_id)

    @property
    def sha256(self) -> str:
        payload = {'schema_version': self.schema_version, 'kernel_id': self.kernel_id, 'components': [vars(c) for c in self.components], 'parameters': [vars(p) for p in self.parameters], 'processes': [{'process_id': p.process_id, 'index': p.index, 'reference_component': p.reference_component, 'rate_base': p.rate_base.to_dict(), 'inhibition': None if p.inhibition is None else p.inhibition.to_dict(), 'stoichiometry': {k: p.stoichiometry[k].to_dict() for k in sorted(p.stoichiometry)}} for p in self.processes], 'conserved_carriers': dict(self.conserved_carriers), 'derived_symbols': dict(self.derived_symbols), 'derived_parameters': list(DERIVED_PARAMETERS), 'rate_unit': self.rate_unit}
        canonical = json.dumps(payload, sort_keys=True, separators=(',', ':'), allow_nan=False)
        return hashlib.sha256(canonical.encode('utf-8')).hexdigest()

class BiokineticKernelError(ValueError):
    """A biokinetic kernel definition is inconsistent (unknown symbol, duplicate id, quarantined parameter)."""

def _symbols(expression: Expression) -> set[str]:
    data = expression.to_dict()
    found: set[str] = set()

    def walk(node):
        if isinstance(node, dict):
            if node.get('op') == 'symbol':
                found.add(node['name'])
            for value in node.values():
                walk(value)
        elif isinstance(node, list):
            for value in node:
                walk(value)
    walk(data)
    return found

def validate_madm1_biokinetic_kernel(kernel: BiokineticKernelDefinition) -> None:
    """Every symbol resolves, ids are unique, indices run 0..n-1, and no quarantined parameter is carried."""
    names = [p.name for p in kernel.parameters]
    if len(set(names)) != len(names):
        raise BiokineticKernelError('duplicate parameter names')
    quarantined = sorted(QUARANTINED_PARAMETER_NAMES.intersection(names))
    if quarantined:
        raise BiokineticKernelError(f'quarantined trace-element parameters {quarantined} may not enter the shared parameter set (madm1 trace quarantine: literature doses are not Monod half-saturations)')
    components = set(kernel.component_ids)
    if len(components) != len(kernel.components):
        raise BiokineticKernelError('duplicate component ids')
    ids = [p.process_id for p in kernel.processes]
    if len(set(ids)) != len(ids):
        raise BiokineticKernelError('duplicate process ids')
    if [p.index for p in kernel.processes] != list(range(len(kernel.processes))):
        raise BiokineticKernelError('process indices must run 0..n-1 in order')
    allowed = set(names) | set(DERIVED_PARAMETERS) | set(kernel.derived_symbols) | components
    for process in kernel.processes:
        unknown = _symbols(process.rate) - allowed
        if unknown:
            raise BiokineticKernelError(f'{process.process_id}: rate reads unknown symbols {sorted(unknown)}')
        for component, coefficient in process.stoichiometry.items():
            if component not in components:
                raise BiokineticKernelError(f'{process.process_id}: stoichiometry names unknown {component!r}')
            if component in kernel.conserved_carriers.values():
                raise BiokineticKernelError(f'{process.process_id}: {component} is a conserved carrier; it is closed, never stated')
            unknown = _symbols(coefficient) - set(names) - set(DERIVED_PARAMETERS)
            if unknown:
                raise BiokineticKernelError(f'{process.process_id}.{component}: coefficient reads unknown symbols {sorted(unknown)}')
        if process.reference_component not in components:
            raise BiokineticKernelError(f'{process.process_id}: unknown reference component')
MADM1_BIOKINETICS_KERNEL = BiokineticKernelDefinition(schema_version=SCHEMA_VERSION, kernel_id='madm1_biokinetics_core_v1', components=MADM1_CORE_CONTENTS, parameters=MADM1_CORE_PARAMETERS, processes=MADM1_CORE_PROCESSES)
MADM1_TWIN_BIOKINETICS_KERNEL = BiokineticKernelDefinition(schema_version=SCHEMA_VERSION, kernel_id='madm1_biokinetics_core_v1_twin_lineage', components=MADM1_TWIN_CONTENTS, parameters=MADM1_CORE_PARAMETERS, processes=MADM1_CORE_PROCESSES, reference='Literature-based modified ADM1 kernel (QSDsan lineage); kinetic parameters are not site fits.')

def kernel_parameter_values(kernel: Optional[BiokineticKernelDefinition]=None) -> dict[str, float]:
    """The kernel's parameter values by name: what both engines' defaults read."""
    return (kernel or MADM1_BIOKINETICS_KERNEL).parameter_values()
__all__ = ['BiokineticKernelDefinition', 'BiokineticKernelError', 'BiokineticProcessDefinition', 'ComponentContent', 'CONSERVED_CARRIERS', 'DERIVED_PARAMETERS', 'DERIVED_SYMBOLS', 'DESIGN_CONTENT_DEPARTURES', 'DESIGN_SOLUBLE_INERT_N_CONTENT', 'KineticParameterDefinition', 'MADM1_BIOKINETICS_KERNEL', 'MADM1_CORE_COMPONENTS', 'MADM1_CORE_CONTENTS', 'MADM1_CORE_PARAMETERS', 'MADM1_CORE_PROCESSES', 'MADM1_TWIN_BIOKINETICS_KERNEL', 'MADM1_TWIN_CONTENTS', 'QUARANTINED_PARAMETER_NAMES', 'SOLUBLE_INERT_UPTAKE', 'SourcedYield', 'derived_parameter_values', 'kernel_parameter_values', 'validate_madm1_biokinetic_kernel']
