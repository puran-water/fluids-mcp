from __future__ import annotations
import hashlib
import json
import re
from dataclasses import dataclass, field
from typing import Mapping, Optional
from .expression import Expression, add, monod, multiply, noncompetitive_inhibition, power, smooth_positive
from .madm1_biokinetics import DERIVED_PARAMETERS, BiokineticProcessDefinition, ComponentContent, KineticParameterDefinition, _h2, _hill_ph, _nutrients, _s, _st
from .madm1_organic_sulfur import MADM1_ORGANIC_SULFUR_EXTENSION
from .madm1_sfe_kinetics import MADM1_SFE_DERIVED_SYMBOLS, MADM1_SFE_KERNEL, MW_S
S_O2 = 'S_O2'
DEPOSIT_COMPONENT = 'S0_deposit'
PLAN = 'Public source package; internal evidence not distributed.'
MICROAERATION_RULING = 'PUBLIC-RULE-6dce80fc160d'
S0_REREDUCTION_RULING = 'HOUSE-MICROAERATION-S0-REREDUCTION-ZERO-BASE-20261010'
PHI_LIQ_RULING = 'PUBLIC-RULE-f700f6751d15'
BOTHEJU_2009 = 'Botheju, Lie & Bakke 2009, Model. Identif. Control 30(4):191-201, doi:10.4173/mic.2009.4.1'
BOTHEJU_2010 = 'Botheju, Lie & Bakke 2010, Model. Identif. Control 31(2):55-65, doi:10.4173/mic.2010.2.2'
BOTHEJU_2011 = 'Botheju & Bakke 2011, Open Waste Manag. J. 4:1-19, doi:10.2174/1876400201104010001'
SANDER_O2 = "Sander 2023 compilation of Henry's law constants v5, Atmos. Chem. Phys. 23:10901-12440, doi:10.5194/acp-23-10901-2023 (henrys-law.org, O2 CAS 7782-44-7): Hcp 1.3e-5 mol m-3 Pa-1, d ln H / d(1/T) 1500 K (reviews Burkholder 2019, Sander 2011; measured Rettich 2000, Wilhelm 1977)"
MULBRY_2017 = "Mulbry, Selmer & Lansing 2017, PLoS ONE 12:e0185738, doi:10.1371/journal.pone.0185738, Table 1 p. 8: maximum headspace oxidation 43, 83, 104, 101 g S m-2 d-1 of liquid/wetted surface at 3500 ppmv biogas H2S, 5 % air, 20-25 C (the 'g S' may be g H2S, about 6 % high as S: lane F4m reader A)"
VOLLERTSEN_2008 = 'Vollertsen et al. 2008, Sci. Total Environ. 394:162-170, doi:10.1016/j.scitotenv.2008.01.028 (abstract): gas-phase H2S oxidation on biologically active surfaces, nth order, n 0.45-0.75; up to 1 mg S m-2 s-1 (86 g S m-2 d-1) at 1000 ppm'
CHEMICALS = 'chemicals 1.5.1 (Bell et al., CalebBell/chemicals, MIT)'
HERLINA_2005 = 'Herlina 2005, Gas transfer at the air-water interface in a turbulent flow environment, dissertation, Universitaet Karlsruhe, doi:10.5445/KSP/1000003544, Table 5.2 p. 96'

@dataclass(frozen=True)
class MicroaerationParameter(KineticParameterDefinition):
    """A kernel parameter with its tag and band: ``low``/``high`` the bracket edges (None where the value is a sourced
    constant with no band, or a ruled rule supplies the edge), ``tag`` one of [V] [U] [DERIVED] [RULED] [ASSUMED]."""
    tag: str = '[V]'
    low: Optional[float] = None
    high: Optional[float] = None
BANDED_TAGS = ('[ASSUMED]',)
_TAGS = ('[V]', '[U]', '[DERIVED]', '[RULED]', '[ASSUMED]')

def _m(name, value, unit, role, source, tag='[V]', low=None, high=None) -> MicroaerationParameter:
    return MicroaerationParameter(name, value, unit, role, source, tag, low, high)
R_GAS = 8.314462618
MW_O2 = 31.9988
MW_N2 = 28.0134
P_REF_H2S_BAR = 0.00101325
N_DEFAULT = 0.6
_MULBRY_MAXIMA = (0.043, 0.083, 0.104, 0.101)
_MULBRY_PPMV = 3500.0
_VOLLERTSEN_RATE = 0.086
_MULBRY_MEDIAN = sorted(_MULBRY_MAXIMA)[1] / 2 + sorted(_MULBRY_MAXIMA)[2] / 2
R_A_REF_DEFAULT = _MULBRY_MEDIAN * (1000.0 / _MULBRY_PPMV) ** N_DEFAULT
R_A_REF_LOW = min(_MULBRY_MAXIMA) * (1000.0 / _MULBRY_PPMV) ** 0.75
R_A_REF_HIGH = _VOLLERTSEN_RATE
_BOTHEJU_EQUIVALENT_KLA_PER_D = (0.14, 0.27, 0.54)
_BOTHEJU_ASSUMED_DEPTH_M = (0.015, 0.05)
_HERLINA_KL_CM_S = (0.000416, 0.000558, 0.000682, 0.000793, 0.00095)
KL_O2_LOW = min(_BOTHEJU_EQUIVALENT_KLA_PER_D) * min(_BOTHEJU_ASSUMED_DEPTH_M)
KL_O2_HIGH = max(_HERLINA_KL_CM_S) * 0.01 * 86400.0
KL_O2_DEFAULT = (KL_O2_LOW * KL_O2_HIGH) ** 0.5
MICROAERATION_PARAMETERS: tuple[MicroaerationParameter, ...] = (_m('MW_O2', MW_O2, 'kg/kmol', 'constant', "chemicals.elements.molecular_weight({'O': 2}), " + CHEMICALS), _m('MW_N2', MW_N2, 'kg/kmol', 'constant', "chemicals.elements.molecular_weight({'N': 2}), " + CHEMICALS + "; argon is lumped with N2 at N2's molar mass (stated: inert in the gas, so only the partial-pressure sum moves, by Ar's share of the oxidant)"), _m('KH_o2', 0.0013, 'mol/L/bar', 'gas', SANDER_O2 + '; 1.3e-5 mol m-3 Pa-1 = 1.3e-3 mol L-1 bar-1'), _m('dH_KH_o2', -1500.0 * R_GAS, 'J/mol', 'gas', SANDER_O2 + "; the kernel's van't Hoff form log KH = log KH0 + dH/R (1/T_base - 1/T), dH = -1500 K x R"), _m('kL_O2', KL_O2_DEFAULT, 'm/d', 'gas', "liquid-side O2 transfer velocity through the wetted headspace surface: the rate is kL_O2 a_wet (S_O2* - S_O2), a_wet the digester's own wetted area per m3 of liquor. [ASSUMED] band, the default its geometric middle (sqrt(low x high) = 0.0415 m/d). LOW edge 2.1e-3 m/d, a LOWER BOUND by construction: " + BOTHEJU_2009 + " Table 3 p. 196 lists an 'equivalent kLa' of 0.14, 0.27 and 0.54 1/d for headspace air over 50 mL of liquor, which is NOT a fit: Eq. 13 (p. 196) divides the whole headspace O2 by V C* t with t the 23.5 d length of the experiment and C* at air, so it rises with the O2 load and can only understate the coefficient; the reactor's geometry is NOT reported (a 100 mL syringe, no diameter, orientation or interface area), so the edge is 0.14 1/d x an [ASSUMED] liquor depth of 0.015 m (0.015-0.05 m: 50 mL in a ~3.6 cm barrel, lying or upright). HIGH edge 0.82 m/d: " + HERLINA_2005 + ": KL,t 4.16e-4 to 9.50e-4 cm/s at 20 C for O2 into clean water through the flat surface of a grid-stirred tank with no wind and no bubbles, the top row taken ([V]). Not quantified from anything read: a digester liquor's surface film and viscosity (both lower kL), 35-38 C (raises it), and a surface broken by bursting gas bubbles (no source). An earlier draft of this constant scaled the digester's kLa and called Botheju's values 'fitted'; both were wrong", '[ASSUMED]', KL_O2_LOW, KL_O2_HIGH), _m('p_ref_h2s', P_REF_H2S_BAR, 'bar', 'constant', 'the reference H2S partial pressure: 1000 ppmv at 1 atm (definition)'), _m('n_h2s', N_DEFAULT, '1', 'rate', VOLLERTSEN_2008 + "; the default is the band's midpoint [ASSUMED]", '[ASSUMED]', 0.45, 0.75), _m('r_A_ref', R_A_REF_DEFAULT, 'kg S/m2/d', 'rate', "areal headspace-surface oxidation rate at p_ref with O2 not limiting, at 20-25 C (no temperature law is sourced: the band carries it). [DERIVED] default: Mulbry 2017's median maximum (0.092 kg S/m2/d at 3500 ppmv) referred to 1000 ppmv at n = 0.6. Band: low edge Mulbry's lowest maximum (0.043) referred at the steepest n = 0.75; high edge Vollertsen's 0.086 at 1000 ppm. " + MULBRY_2017 + '; ' + VOLLERTSEN_2008, '[ASSUMED]', R_A_REF_LOW, R_A_REF_HIGH), _m('K_pO2', 0.003, 'bar', 'half_saturation', 'no sourced O2 dependence for a headspace biofilm. [ASSUMED] band 0.1-1 % O2 (at 1 atm) with the default at 0.3 %: Mulbry 2017 p. 5 (full-scale bag digesters): biogas H2S < 500 ppmv only at residual O2 >= 0.3 % and < 100 ppmv only at >= 0.5 % -- the validation point, not a constant; the band puts the O2-saturation factor at 0.3-0.5 % O2 between 0.23 and 0.83', '[ASSUMED]', 0.001, 0.01), _m('f_SO4', 0.0, '1', 'fraction', 'share of the H2S oxidised that goes to sulfate (2 mol O2/mol S, to the liquor as acid); the design default 0 is the house product S0 alone (house_standards.microaeration SULFUR_PRODUCT; Jenicek 2014: thiosulfate and sulfate negligible in a micro-aerated digester). [ASSUMED] band 0-0.27: S0 > 98 % at limiting O2 (Pokorna-Krayzelova 2017 thesis, via lane F4m); a maximum S0 yield of 73 +/- 10 % at O2/S 0.6-1.0 (Janssen et al. 1995, doi:10.1002/bit.260470307) leaves at most ~0.27 to thiosulfate plus sulfate, a ceiling on sulfate', '[ASSUMED]', 0.0, 0.27), _m('phi_liq', 0.45, '1', 'fraction', f"share of the S0 formed that falls to the liquid S_S0; the rest is the headspace deposit. [RULED] band 0.1-0.8 ({PHI_LIQ_RULING}: 'S0 to the liquor banded 0.1-0.8'); the default is the band's midpoint [ASSUMED]. Krayzelova et al. 2014 (doi:10.1016/j.biortech.2014.09.056, liquid-dosed UASB: 33 % of influent S left as effluent solids, mostly S0, against 10 % as S0 in the gas space) against Ramos et al. 2014 (doi:10.1016/j.biortech.2014.02.001, headspace dosing: the S0 recovered from the headspace)", '[RULED]', 0.1, 0.8), _m('dH_G1_S0', -265225.0, 'J/mol', 'constant', 'H2S(g) + 1/2 O2(g) -> S(s) + H2O(l): Hf H2O(l) -285825 - Hf H2S(g) -20600 J/mol, ' + CHEMICALS), _m('dH_G2_SO4', -793400.0, 'J/mol', 'constant', 'H2S(g) + 2 O2(g) -> H2SO4(l): Hf H2SO4(l) -814000 - Hf H2S(g) -20600 J/mol, ' + CHEMICALS + " (the dilution of the acid in the liquor, about -96 kJ/mol more, is not counted: the branch's default is 0)"), _m('Y_su_aer', 0.5, 'kg COD/kg COD', 'yield', BOTHEJU_2009 + ' p. 193 (from ASM2 / Henze et al. 2002)'), _m('Y_aa_aer', 0.4, 'kg COD/kg COD', 'yield', BOTHEJU_2009 + ' p. 193'), _m('Y_fa_aer', 0.3, 'kg COD/kg COD', 'yield', BOTHEJU_2009 + ' p. 193'), _m('K_su_aer', 0.1, 'kg COD/m3', 'half_saturation', BOTHEJU_2009 + ' p. 193: one fifth of the ADM1 K_su'), _m('K_aa_aer', 0.06, 'kg COD/m3', 'half_saturation', BOTHEJU_2009 + ' p. 193: one fifth of the ADM1 K_aa'), _m('K_fa_aer', 0.08, 'kg COD/m3', 'half_saturation', BOTHEJU_2009 + ' p. 193: one fifth of the ADM1 K_fa'), _m('K_va_aer', 0.04, 'kg COD/m3', 'half_saturation', BOTHEJU_2010 + ' p. 56: 20 % of the ADM1 fermentative K; [DERIVED] 0.2 x K_c4 0.2', '[DERIVED]'), _m('K_bu_aer', 0.04, 'kg COD/m3', 'half_saturation', BOTHEJU_2010 + ' p. 56: 20 % of the ADM1 fermentative K; [DERIVED] 0.2 x K_c4 0.2', '[DERIVED]'), _m('K_pro_aer', 0.02, 'kg COD/m3', 'half_saturation', BOTHEJU_2010 + ' p. 56: 20 % of the ADM1 fermentative K; [DERIVED] 0.2 x K_pro 0.1', '[DERIVED]'), _m('K_ac_aer', 0.03, 'kg COD/m3', 'half_saturation', BOTHEJU_2010 + ' p. 56: 20 % of the ADM1 fermentative K; [DERIVED] 0.2 x K_ac 0.15', '[DERIVED]'), _m('Y_va_aer', 0.3, 'kg COD/kg COD', 'yield', BOTHEJU_2010 + ' p. 56: 5 x the fermentative yield, a third to each acidogen; [DERIVED] 5 x Y_c4 0.06', '[DERIVED]'), _m('Y_bu_aer', 0.3, 'kg COD/kg COD', 'yield', BOTHEJU_2010 + ' p. 56: 5 x the fermentative yield; [DERIVED] 5 x Y_c4 0.06', '[DERIVED]'), _m('Y_pro_aer', 0.2, 'kg COD/kg COD', 'yield', BOTHEJU_2010 + ' p. 56: 5 x the fermentative yield; [DERIVED] 5 x Y_pro 0.04', '[DERIVED]'), _m('Y_ac_aer', 0.25, 'kg COD/kg COD', 'yield', BOTHEJU_2010 + ' p. 56: 5 x the fermentative yield; [DERIVED] 5 x Y_ac 0.05', '[DERIVED]'), _m('K_O2_fac', 0.0002, 'kg O2/m3', 'half_saturation', BOTHEJU_2009 + " p. 193 takes K_O2 as ASM2's aerobic half-saturation (Henze et al. 1995) without printing it; the ASM2d heterotroph K_O2 is 0.2 g O2/m3 (Henze et al. 1999, doi:10.2166/wst.1999.0036). [ASSUMED] band 0.05-0.5 g O2/m3: no facultative-acidogen value is measured", '[ASSUMED]', 5e-05, 0.0005), _m('KI_O2', 0.0002, 'kg O2/m3', 'inhibition', BOTHEJU_2009 + " p. 193 (Eq. 4: the inhibition K equals the uptake K_O2). [ASSUMED] band: low edge 1e-7 kg/m3 (3 nM) because methane production 'stopped at below 30 nM of oxygen' in suspended rumen and digester cultures (" + BOTHEJU_2011 + ' p. 3), so a K well under 30 nM; high edge 8.1 mg/L, the influent DO at which Shen & Guiot 1996 (doi:10.1002/(SICI)1097-0290(19960320)49:6<611::AID-BIT2>3.0.CO;2-R, abstract) saw only a slight decline in granules, a tolerance threshold rather than a constant. The liquor DO is a CHECK beside it', '[ASSUMED]', 1e-07, 0.0081), _m('k_S0_red', 0.0, '1/d', 'rate', f"{S0_REREDUCTION_RULING}: zero in the base case; a bracket on liquid-phase S0 only, from zero to the [ASSUMED] bound k = 1/SRT (the digester's own SRT: the engine sets the edge, house_standards.microaeration.s0_rereduction_bracket_upper_k_per_d)", '[RULED]', 0.0, None), _m('K_S0red_donor', 0.0001, 'kg COD/m3', 'half_saturation', "a numerical donor switch on acetate so re-reduction stops where no donor is left; at any digester acetate (> 0.01 kg COD/m3) it is 1 within 1 %. [ASSUMED] band 1e-5-1e-3: the electron donor of S0 reduction in a digester is unmeasured (Xu et al. 2013 'rebound', abstract; WATS k_s0_acid, aquakin)", '[ASSUMED]', 1e-05, 0.001), _m('COD_per_kmol_S0_reduced', MW_O2 / 2.0, 'kg COD/kmol S', 'constant', "2 e- per S0 reduced to sulfide = 0.5 kmol O2-equivalent at the O2 molar mass (chemicals), the basis this extension's O2 gas balance is kept on; the S/Fe kernel's SRB basis rounds O to 16 (8 kg COD per kmol e-), 2.3e-5 apart"))
MICROAERATION_DERIVED_SYMBOLS = {'p_o2_bar': 'O2 partial pressure in the headspace (bar)', 'KH_o2_T': "Henry's constant for O2 at T on the engine's liquid basis (kg O2/m3/bar): KH_o2 x MW_O2 by van't Hoff", 'a_wet_m2_per_m3': "wetted headspace surface per m3 of liquid (m2/m3): A_wet / V_liq, A_wet from the digester's geometry (a sizing input, never invented by a rating); read by G1 and by the O2 surface transfer G4", 'gas_smoothing': "engine-bound smoothing (bar) of the H2S partial pressure G1's power law reads, 0.5 (p + sqrt(p^2 + eps^2)): the twin binds 1e-30 (its exact law), a Newton engine a value it can differentiate through at p -> 0"}
MICROAERATION_CONTENTS = (ComponentContent(S_O2, 0.0, 0.0, 0.0, 0.0, 0.0, 'O2'),)

@dataclass(frozen=True)
class GasComponent:
    component_id: str
    molar_mass: float
    measured_as: str
    liquid_partner: Optional[str]
    doc: str
BASE_GAS_COMPONENTS = (GasComponent('S_h2', 16.0, 'COD', 'S_h2', 'hydrogen (kg COD)'), GasComponent('S_ch4', 64.0, 'COD', 'S_ch4', 'methane (kg COD)'), GasComponent('S_co2', 12.0, 'C', 'S_IC', 'carbon dioxide (kg C)'), GasComponent('S_h2s', MW_S, 'S', 'S_IS', 'hydrogen sulfide (kg S)'))
MICROAERATION_GAS_COMPONENTS = (GasComponent('S_o2', MW_O2, 'O2', S_O2, 'oxygen (kg O2)'), GasComponent('S_n2', MW_N2, 'N2', None, 'nitrogen with argon lumped (kg N2; dissolved N2 omitted)'))
SULFUR_SPECIES = ('S_SO4', 'S_IS', 'S_S0', 'X_Sorg', 'S_h2s', DEPOSIT_COMPONENT)
OXYGEN_DEMAND: dict[str, float] = {**{c.component_id: 1.0 for c in MADM1_SFE_KERNEL.components if c.measured_as == 'COD'}, S_O2: -1.0, 'S_o2': -1.0, 'S_h2': 1.0, 'S_ch4': 1.0, 'S_IS': 2.0 * MW_O2 / MW_S, 'X_Sorg': 2.0 * MW_O2 / MW_S, 'S_h2s': 2.0 * MW_O2 / MW_S, 'S_S0': 1.5 * MW_O2 / MW_S, DEPOSIT_COMPONENT: 1.5 * MW_O2 / MW_S, 'S_SO4': 0.0}

@dataclass(frozen=True)
class MicroaerationProcessDefinition(BiokineticProcessDefinition):
    """A process of the micro-aeration extension. ``index`` is its TWIN index; ``stoichiometry`` is over LIQUID
    components, ``gas_stoichiometry`` over headspace components and ``deposit_stoichiometry`` the headspace S0 deposit,
    all per unit of the reference component; ``phase`` says where the rate is evaluated ("liquid", "gas_surface",
    "gas_liquid_transfer"). Rates are per m3 of liquid per day."""
    gas_stoichiometry: Mapping[str, Expression] = field(default_factory=dict)
    deposit_stoichiometry: Mapping[str, Expression] = field(default_factory=dict)
    phase: str = 'liquid'
    equations: str = ''

def _one_minus(name: str) -> Expression:
    return add(1.0, multiply(-1.0, _s(name)))

def _i1() -> Expression:
    """Botheju's I1 for the acidogens: this kernel's own nutrient x Hill pH factor (its sugar / amino-acid uptakes)."""
    return multiply(_nutrients(), _hill_ph('aa'))

def _i2() -> Expression:
    """Botheju's I2: I1 x the H2 inhibition of LCFA uptake (this kernel's KI_h2_fa)."""
    return multiply(_i1(), _h2('KI_h2_fa'))

def _o2_monod() -> Expression:
    return monod(_s(S_O2), _s('K_O2_fac'))

def _aerobic_uptake(process_id: str, index: int, substrate: str, biomass: str, rate: str, half: str, yield_name: str, inhibition: Expression, equation: str) -> MicroaerationProcessDefinition:
    """Botheju 2009 Eqs 1-3: k_m X Monod(S, K_aer) Monod(S_O2, K_O2) I; the substrate's COD to biomass and O2."""
    return MicroaerationProcessDefinition(process_id, index, multiply(_s(rate), _s(biomass), monod(_s(substrate), _s(half)), _o2_monod()), substrate, {substrate: _st(x=-1.0)['x'], biomass: _s(yield_name), S_O2: multiply(-1.0, _one_minus(yield_name))}, inhibition, equations=equation)

def _aerobic_vfa(process_id: str, index: int, substrate: str, half: str, yield_name: str, equation: str) -> MicroaerationProcessDefinition:
    """Botheju 2010 Eq. 1: (k_su X_su I1 + k_aa X_aa I1 + k_fa X_fa I2) Monod(S_vfa) Monod(S_O2); yield split in thirds."""
    capacity = add(multiply(_s('k_su'), _s('X_su'), _i1()), multiply(_s('k_aa'), _s('X_aa'), _i1()), multiply(_s('k_fa'), _s('X_fa'), _i2()))
    third = multiply(_s(yield_name), 1.0 / 3.0)
    return MicroaerationProcessDefinition(process_id, index, multiply(capacity, monod(_s(substrate), _s(half)), _o2_monod()), substrate, {substrate: _st(x=-1.0)['x'], 'X_su': third, 'X_aa': third, 'X_fa': third, S_O2: multiply(-1.0, _one_minus(yield_name))}, None, equations=equation)

def _g1_rate() -> Expression:
    """Only p_H2S is smoothed: the power law p^n needs it (its derivative is infinite at p = 0). The O2 Monod term is
    smooth at p_O2 = 0 as it stands, and a smoothed p_O2 would leave a floor rate at ZERO dose (found by the WaterTAP
    zero-dose delete test, 2026-10-10: ~1e-11 kg S/s), so zero oxygen must give exactly zero oxidation."""
    eps = _s('gas_smoothing')
    pressure = multiply(smooth_positive(_s('p_h2s_bar'), eps), power(_s('p_ref_h2s'), -1.0))
    return multiply(_s('a_wet_m2_per_m3'), _s('r_A_ref'), power(pressure, _s('n_h2s')), monod(_s('p_o2_bar'), _s('K_pO2')))
_O2_PER_S = multiply(_s('MW_O2'), power(_s('MW_S'), -1.0), add(0.5, multiply(1.5, _s('f_SO4'))))
HEADSPACE_H2S_OXIDATION = MicroaerationProcessDefinition('headspace_H2S_oxidation', 76, _g1_rate(), 'S_h2s', {'S_S0': multiply(_one_minus('f_SO4'), _s('phi_liq')), 'S_SO4': _s('f_SO4')}, None, gas_stoichiometry={'S_h2s': _st(x=-1.0)['x'], 'S_o2': multiply(-1.0, _O2_PER_S)}, deposit_stoichiometry={DEPOSIT_COMPONENT: multiply(_one_minus('f_SO4'), _one_minus('phi_liq'))}, phase='gas_surface', equations='r = a_wet r_A_ref (p_H2S/p_ref)^n p_O2/(K_pO2 + p_O2) (F4m dossier G1; Vollertsen 2008 nth order; Mulbry 2017 areal maxima); G2 sulfate share f_SO4; G3 fall-out phi_liq')
O2_TRANSFER = MicroaerationProcessDefinition('O2_transfer', 77, multiply(_s('kL_O2'), _s('a_wet_m2_per_m3'), add(_s(S_O2), multiply(-1.0, _s('KH_o2_T'), _s('p_o2_bar')))), S_O2, {S_O2: _st(x=-1.0)['x']}, None, gas_stoichiometry={'S_o2': _st(x=1.0)['x']}, phase='gas_liquid_transfer', equations="kL_O2 a_wet (S_O2 - KH_O2(T) p_O2): surface transfer through the wetted headspace area (the sign and Henry form of process 74; the coefficient is NOT the digester's kLa)")
_L1 = (_aerobic_uptake('aerobic_uptake_sugars', 78, 'S_su', 'X_su', 'k_su', 'K_su_aer', 'Y_su_aer', _i1(), 'Botheju 2009 Eq. 1'), _aerobic_uptake('aerobic_uptake_amino_acids', 79, 'S_aa', 'X_aa', 'k_aa', 'K_aa_aer', 'Y_aa_aer', _i1(), 'Botheju 2009 Eq. 2'), _aerobic_uptake('aerobic_uptake_LCFA', 80, 'S_fa', 'X_fa', 'k_fa', 'K_fa_aer', 'Y_fa_aer', _i2(), 'Botheju 2009 Eq. 3'), _aerobic_vfa('aerobic_uptake_valerate', 81, 'S_va', 'K_va_aer', 'Y_va_aer', 'Botheju 2010 Eq. 1'), _aerobic_vfa('aerobic_uptake_butyrate', 82, 'S_bu', 'K_bu_aer', 'Y_bu_aer', 'Botheju 2010 Eq. 1'), _aerobic_vfa('aerobic_uptake_propionate', 83, 'S_pro', 'K_pro_aer', 'Y_pro_aer', 'Botheju 2010 Eq. 1'), _aerobic_vfa('aerobic_uptake_acetate', 84, 'S_ac', 'K_ac_aer', 'Y_ac_aer', 'Botheju 2010 Eq. 1'))
S0_REREDUCTION = MicroaerationProcessDefinition('S0_rereduction', 85, multiply(_s('k_S0_red'), _s('S_S0'), monod(_s('S_ac'), _s('K_S0red_donor'))), 'S_S0', {'S_S0': _st(x=-1.0)['x'], 'S_IS': _st(x=1.0)['x'], 'S_ac': multiply(-1.0, _s('COD_per_kmol_S0_reduced'), power(_s('MW_S'), -1.0))}, None, equations=f'k_S0_red S_S0 Monod(S_ac); k_S0_red = 0 in the base case ({S0_REREDUCTION_RULING})')
MICROAERATION_PROCESSES: tuple[MicroaerationProcessDefinition, ...] = (HEADSPACE_H2S_OXIDATION, O2_TRANSFER, *_L1, S0_REREDUCTION)
O2_INHIBITION = noncompetitive_inhibition(_s(S_O2), _s('KI_O2'))
O2_INHIBITED_PROCESS_INDICES = (3, 4, 5, 6, 7, 8, 9, 10, 25, 27, 29, 31, 32)

@dataclass(frozen=True)
class MicroaerationExtension:
    """The opt-in extension as one object, with an identity hash an engine stamps beside the kernels it extends."""
    extension_id: str
    base_kernel_id: str
    base_kernel_sha256: str
    base_extension_id: str
    base_extension_sha256: str
    contents: tuple[ComponentContent, ...]
    gas_components: tuple[GasComponent, ...]
    parameters: tuple[MicroaerationParameter, ...]
    processes: tuple[MicroaerationProcessDefinition, ...]
    inhibition: Expression
    inhibited_process_indices: tuple[int, ...]
    derived_symbols: Mapping[str, str]
    source: str

    @property
    def component_ids(self) -> tuple[str, ...]:
        return tuple((c.component_id for c in self.contents))

    @property
    def gas_component_ids(self) -> tuple[str, ...]:
        return tuple((g.component_id for g in self.gas_components))

    def parameter(self, name: str) -> MicroaerationParameter:
        for p in self.parameters:
            if p.name == name:
                return p
        raise KeyError(name)

    def parameter_values(self) -> dict[str, float]:
        return {p.name: float(p.value) for p in self.parameters}

    @property
    def sha256(self) -> str:

        def stoich(mapping):
            return sorted(((k, v.to_dict()) for k, v in mapping.items()))
        body = {'extension_id': self.extension_id, 'base': [self.base_kernel_id, self.base_kernel_sha256, self.base_extension_id, self.base_extension_sha256], 'contents': [sorted(vars(c).items()) for c in self.contents], 'gas': [sorted(vars(g).items()) for g in self.gas_components], 'parameters': [sorted(vars(p).items()) for p in self.parameters], 'processes': [[p.process_id, p.index, p.phase, p.reference_component, p.rate.to_dict(), stoich(p.stoichiometry), stoich(p.gas_stoichiometry), stoich(p.deposit_stoichiometry)] for p in self.processes], 'inhibition': [self.inhibition.to_dict(), list(self.inhibited_process_indices)], 'derived_symbols': dict(self.derived_symbols), 'source': self.source}
        return hashlib.sha256(json.dumps(body, sort_keys=True, default=str).encode()).hexdigest()

def microaeration_parameter_values(overrides: Optional[Mapping[str, float]]=None, *, extension: Optional['MicroaerationExtension']=None) -> dict[str, float]:
    """The S/Fe kernel's parameters, its derived stoichiometric parameters and this extension's, with ``overrides``
    applied (an unknown name is refused, so a misspelt key never leaves a default in force)."""
    from .madm1_biokinetics import derived_parameter_values
    ext = extension or MADM1_MICROAERATION_EXTENSION
    values = {**MADM1_SFE_KERNEL.parameter_values(), **ext.parameter_values()}
    for name, value in (overrides or {}).items():
        if name not in values:
            raise KeyError(f'{name!r} is not a parameter of the S/Fe kernel or {ext.extension_id}')
        values[name] = float(value)
    values.update(derived_parameter_values(values))
    return values

def _symbols(expression: Expression) -> set[str]:
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
    walk(expression.to_dict())
    return found

def _value(expression: Expression, values: Mapping[str, float]) -> float:
    """A stoichiometric coefficient's number (the schema carries no evaluator; coefficients are + x ^ of numbers)."""
    data = expression.to_dict()

    def ev(node):
        op = node['op']
        if op == 'constant':
            return float(node['value'])
        if op == 'symbol':
            return float(values[node['name']])
        if op == 'add':
            return sum((ev(v) for v in node['values']))
        if op == 'multiply':
            out = 1.0
            for v in node['values']:
                out *= ev(v)
            return out
        if op == 'power':
            return ev(node['left']) ** ev(node['right'])
        raise ValueError(f'a stoichiometric coefficient may not use {op!r}')
    return ev(data)

def process_rows(process: MicroaerationProcessDefinition, values: Mapping[str, float]) -> dict[str, float]:
    """Every coefficient of a process (liquid, gas and deposit), as numbers at ``values``."""
    rows = {}
    for mapping in (process.stoichiometry, process.gas_stoichiometry, process.deposit_stoichiometry):
        for k, v in mapping.items():
            rows[k] = rows.get(k, 0.0) + _value(v, values)
    return rows
_DOI = re.compile('doi:10\\.\\d{4,9}/\\S+')

def validate_microaeration_extension(ext: 'MicroaerationExtension') -> None:
    """Refuse an extension that would not compile against its base or does not close (module docstring)."""
    base, organic = (MADM1_SFE_KERNEL, MADM1_ORGANIC_SULFUR_EXTENSION)
    if (ext.base_kernel_sha256, ext.base_extension_sha256) != (base.sha256, organic.sha256):
        raise ValueError('the micro-aeration extension is defined against a different S/Fe kernel or organic-S extension')
    known_liquid = set(base.component_ids) | set(organic.component_ids)
    clash = set(ext.component_ids) & known_liquid
    if clash:
        raise ValueError(f'micro-aeration components already in the base: {sorted(clash)}')
    base_gas = {g.component_id for g in BASE_GAS_COMPONENTS}
    if set(ext.gas_component_ids) & base_gas:
        raise ValueError('micro-aeration gas components already in the base headspace')
    liquid = known_liquid | set(ext.component_ids)
    gas = base_gas | set(ext.gas_component_ids)
    first = organic.processes[-1].index + 1
    params = {p.name for p in base.parameters} | {p.name for p in ext.parameters}
    if len(params) != len(base.parameters) + len(ext.parameters):
        raise ValueError('a micro-aeration parameter shadows a base-kernel parameter')
    for p in ext.parameters:
        if p.tag not in _TAGS:
            raise ValueError(f'{p.name}: tag {p.tag!r} is not one of {_TAGS}')
        if p.tag in BANDED_TAGS and (p.low is None or p.high is None):
            raise ValueError(f'{p.name}: an {p.tag} value carries its band (low, high), never a bare number')
        if p.low is not None and (not p.low <= p.value):
            raise ValueError(f"{p.name}: default {p.value} below its band's low edge {p.low}")
        if p.high is not None and (not p.value <= p.high):
            raise ValueError(f"{p.name}: default {p.value} above its band's high edge {p.high}")
        sourced = bool(_DOI.search(p.source)) or 'chemicals' in p.source or p.role == 'constant'
        if not sourced and p.tag not in ('[ASSUMED]', '[RULED]'):
            raise ValueError(f'{p.name}: neither a DOI nor an [ASSUMED]/[RULED] reason in its source')
    symbols = params | set(DERIVED_PARAMETERS) | set(base.derived_symbols) | set(MADM1_SFE_DERIVED_SYMBOLS) | set(ext.derived_symbols) | liquid
    values = microaeration_parameter_values(extension=ext)
    for offset, proc in enumerate(ext.processes):
        if proc.index != first + offset:
            raise ValueError(f"{proc.process_id}: index {proc.index}, expected {first + offset} (after the organic-S extension's {first - 1})")
        unknown = _symbols(proc.rate) - symbols
        if unknown:
            raise ValueError(f'{proc.process_id}: rate reads unknown symbols {sorted(unknown)}')
        if set(proc.stoichiometry) - liquid:
            raise ValueError(f'{proc.process_id}: liquid stoichiometry names unknown {sorted(set(proc.stoichiometry) - liquid)}')
        if set(proc.gas_stoichiometry) - gas:
            raise ValueError(f'{proc.process_id}: gas stoichiometry names unknown {sorted(set(proc.gas_stoichiometry) - gas)}')
        if set(proc.deposit_stoichiometry) - {DEPOSIT_COMPONENT}:
            raise ValueError(f'{proc.process_id}: the only deposit is {DEPOSIT_COMPONENT}')
        if proc.stoichiometry.keys() & set(base.conserved_carriers.values()):
            raise ValueError(f'{proc.process_id}: a conserved carrier is closed, never stated')
        for f_so4 in (values['f_SO4'], 0.0, 0.27):
            for phi in (values['phi_liq'], 0.1, 0.8):
                at = {**values, 'f_SO4': f_so4, 'phi_liq': phi}
                rows = process_rows(proc, at)
                s_row = sum((c for k, c in rows.items() if k in SULFUR_SPECIES))
                if abs(s_row) > 1e-12:
                    raise ValueError(f'{proc.process_id}: the sulfur row does not close ({s_row!r})')
                missing = sorted(set(rows) - set(OXYGEN_DEMAND) - {'S_n2'})
                if missing:
                    raise ValueError(f'{proc.process_id}: no oxygen demand stated for {missing}')
                od = sum((c * OXYGEN_DEMAND.get(k, 0.0) for k, c in rows.items()))
                if abs(od) > 1e-09:
                    raise ValueError(f'{proc.process_id}: the oxygen-demand (COD) row does not close ({od!r})')
    targets = {p.index: p for p in base.processes}
    for index in ext.inhibited_process_indices:
        if index not in targets or targets[index].process_id.startswith(('decay', 'lysis', 'hydrolysis')):
            raise ValueError(f'O2 inhibition targets {index}, which is not an uptake of the S/Fe kernel')
    unknown = _symbols(ext.inhibition) - symbols
    if unknown:
        raise ValueError(f'the O2 inhibition reads unknown symbols {sorted(unknown)}')
G1_REACTION_ENTHALPY_J_MOL = {'S0': -265225.0, 'SO4': -793400.0}
MADM1_MICROAERATION_EXTENSION = MicroaerationExtension(extension_id='madm1_microaeration_v1', base_kernel_id=MADM1_SFE_KERNEL.kernel_id, base_kernel_sha256=MADM1_SFE_KERNEL.sha256, base_extension_id=MADM1_ORGANIC_SULFUR_EXTENSION.extension_id, base_extension_sha256=MADM1_ORGANIC_SULFUR_EXTENSION.sha256, contents=MICROAERATION_CONTENTS, gas_components=MICROAERATION_GAS_COMPONENTS, parameters=MICROAERATION_PARAMETERS, processes=MICROAERATION_PROCESSES, inhibition=O2_INHIBITION, inhibited_process_indices=O2_INHIBITED_PROCESS_INDICES, derived_symbols=MICROAERATION_DERIVED_SYMBOLS, source=f"{MICROAERATION_RULING} 'In the digester, first scope (Recommended)'; {S0_REREDUCTION_RULING} 'Zero base, bracket only (Recommended)'; {PHI_LIQ_RULING} (phi_liq band); {PLAN}; rate forms {BOTHEJU_2009} Eqs 1-4, {BOTHEJU_2010} Eq. 1, Vollertsen 2008 nth order, Mulbry 2017 areal maxima")
