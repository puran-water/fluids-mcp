from __future__ import annotations
from typing import Optional
from .expression import Expression, add, multiply, power, smooth_positive
from .madm1_biokinetics import MADM1_BIOKINETICS_KERNEL, MADM1_CORE_CONTENTS, MADM1_CORE_PARAMETERS, MADM1_CORE_PROCESSES, SCHEMA_VERSION, TWIN_SOURCE, BiokineticKernelDefinition, BiokineticProcessDefinition, ComponentContent, KineticParameterDefinition, _c4_share, _decay, _h2s, _hill_ph, _inhibition, _nutrients, _p, _s, _st, validate_madm1_biokinetic_kernel
SRB_SOURCE = "Flores-Alsina et al. 2016 (mADM1 S extension; Barrera et al. 2015 SRB kinetics), the twin's defaults; " + TWIN_SOURCE
FE_SOURCE = "Flores-Alsina et al. 2016 (mADM1 Fe extension; Solon et al. 2017 HFO module), the twin's defaults; " + TWIN_SOURCE
HFO_REDUCTION_SOURCE = "Poulton, Krom & Raiswell 2004 (GCA 68:3703, 10.1016/j.gca.2004.03.012) ferrihydrite sulfidation half-lives and Kraal et al. 2022 (10.1016/j.gca.2021.12.032): screening k = 4.7e5 / MW_Fe m3/kmol/d (the twin's K_FE3_REDUCTION_SCREENING), literature-screening, NOT plant-calibrated"
FES_SOURCE = "FeS (amorphous) + H+ = Fe2+ + HS-, log K 2.95 (Davison 1991 proton-inclusive convention, the shared Fe-S-P kernel FeS_amorphous); the twin's crystallisation law and rate constant (k_cryst 100 /h)"
ELEMENT_SOURCE = 'chemicals.elements molecular weights as the twin imports them (S_mw, Fe_mw)'
MW_S = 32.065
MW_Fe = 55.845
MADM1_SFE_COMPONENTS = ('S_SO4', 'S_IS', 'X_hSRB', 'X_aSRB', 'X_pSRB', 'X_c4SRB', 'S_S0', 'S_Fe3', 'S_Fe2', 'X_HFO_H', 'X_HFO_L', 'X_HFO_old', 'X_HFO_HP', 'X_HFO_LP', 'X_HFO_HP_old', 'X_HFO_LP_old', 'X_FeS')
IP_HFO_HP = 0.6655656406124093
IP_HFO_LP = 0.17193779049153907
_BIOMASS = (0.36612, 0.08615, 0.02154, 0.004768, 0.002964, 'COD')
MADM1_SFE_CONTENTS = tuple((ComponentContent(cid, *values) for cid, values in {'S_SO4': (0.0, 0.0, 0.0, 0.0, 0.0, 'S'), 'S_IS': (0.0, 0.0, 0.0, 0.0, 0.0, 'S'), 'X_hSRB': _BIOMASS, 'X_aSRB': _BIOMASS, 'X_pSRB': _BIOMASS, 'X_c4SRB': _BIOMASS, 'S_S0': (0.0, 0.0, 0.0, 0.0, 0.0, 'S'), 'S_Fe3': (0.0, 0.0, 0.0, 0.0, 0.0, 'Fe'), 'S_Fe2': (0.0, 0.0, 0.0, 0.0, 0.0, 'Fe'), 'X_HFO_H': (0.0, 0.0, 0.0, 0.0, 0.0, 'Fe'), 'X_HFO_L': (0.0, 0.0, 0.0, 0.0, 0.0, 'Fe'), 'X_HFO_old': (0.0, 0.0, 0.0, 0.0, 0.0, 'Fe'), 'X_HFO_HP': (0.0, 0.0, IP_HFO_HP, 0.0, 0.0, 'Fe'), 'X_HFO_LP': (0.0, 0.0, IP_HFO_LP, 0.0, 0.0, 'Fe'), 'X_HFO_HP_old': (0.0, 0.0, IP_HFO_HP, 0.0, 0.0, 'Fe'), 'X_HFO_LP_old': (0.0, 0.0, IP_HFO_LP, 0.0, 0.0, 'Fe'), 'X_FeS': (0.0, 0.0, 0.0, 0.0, 0.0, 'FeS')}.items()))
MADM1_SFE_PARAMETERS: tuple[KineticParameterDefinition, ...] = (_p('Y_hSRB', 0.05, 'kg COD/kg COD', 'yield', SRB_SOURCE), _p('Y_aSRB', 0.05, 'kg COD/kg COD', 'yield', SRB_SOURCE), _p('Y_pSRB', 0.04, 'kg COD/kg COD', 'yield', SRB_SOURCE), _p('Y_c4SRB', 0.06, 'kg COD/kg COD', 'yield', SRB_SOURCE), _p('k_hSRB', 41.125, '1/d', 'rate', SRB_SOURCE), _p('k_aSRB', 10.0, '1/d', 'rate', SRB_SOURCE), _p('k_pSRB', 16.25, '1/d', 'rate', SRB_SOURCE), _p('k_c4SRB', 23.0, '1/d', 'rate', SRB_SOURCE), _p('b_hSRB', 0.02, '1/d', 'decay', SRB_SOURCE), _p('b_aSRB', 0.02, '1/d', 'decay', SRB_SOURCE), _p('b_pSRB', 0.02, '1/d', 'decay', SRB_SOURCE), _p('b_c4SRB', 0.02, '1/d', 'decay', SRB_SOURCE), _p('K_hSRB', 5.96e-06, 'kg COD/m3', 'half_saturation', SRB_SOURCE), _p('K_aSRB', 0.176, 'kg COD/m3', 'half_saturation', SRB_SOURCE), _p('K_pSRB', 0.088, 'kg COD/m3', 'half_saturation', SRB_SOURCE), _p('K_c4SRB', 0.1739, 'kg COD/m3', 'half_saturation', SRB_SOURCE), _p('K_so4_hSRB', 0.000104 * MW_S, 'kg S/m3', 'half_saturation', SRB_SOURCE), _p('K_so4_aSRB', 0.0002 * MW_S, 'kg S/m3', 'half_saturation', SRB_SOURCE), _p('K_so4_pSRB', 0.0002 * MW_S, 'kg S/m3', 'half_saturation', SRB_SOURCE), _p('K_so4_c4SRB', 0.0002 * MW_S, 'kg S/m3', 'half_saturation', SRB_SOURCE), _p('KI_h2s_hSRB', 0.499, 'kg S/m3', 'inhibition', SRB_SOURCE), _p('KI_h2s_aSRB', 0.499, 'kg S/m3', 'inhibition', SRB_SOURCE), _p('KI_h2s_pSRB', 0.52, 'kg S/m3', 'inhibition', SRB_SOURCE), _p('KI_h2s_c4SRB', 0.52, 'kg S/m3', 'inhibition', SRB_SOURCE), _p('pH_LL_h2_SRB', 5.0, '1', 'inhibition', SRB_SOURCE), _p('pH_UL_h2_SRB', 6.0, '1', 'inhibition', SRB_SOURCE), _p('pH_LL_ac_SRB', 6.0, '1', 'inhibition', SRB_SOURCE), _p('pH_UL_ac_SRB', 7.0, '1', 'inhibition', SRB_SOURCE), _p('pH_LL_aa_SRB', 6.0, '1', 'inhibition', SRB_SOURCE), _p('pH_UL_aa_SRB', 7.0, '1', 'inhibition', SRB_SOURCE), _p('k_Fe3t2_h2', 470000.0 / MW_Fe, 'm3/kmol/d', 'rate', HFO_REDUCTION_SOURCE), _p('k_Fe3t2_is', 470000.0 / MW_Fe, 'm3/kmol/d', 'rate', HFO_REDUCTION_SOURCE), _p('hfo_bound_reduction_factor', 1.0, '1', 'rate', FE_SOURCE), _p('hfo_aged_reduction_factor', 1.0, '1', 'rate', FE_SOURCE), _p('hfo_aged_bound_reduction_factor', 1.0, '1', 'rate', FE_SOURCE), _p('q_aging_H', 450.0, '1/d', 'rate', FE_SOURCE), _p('q_aging_L', 0.1, '1/d', 'rate', FE_SOURCE), _p('q_Pcoprec', 360.0, '1/d', 'rate', FE_SOURCE), _p('q_Pbinding', 0.3, '1/d', 'rate', FE_SOURCE), _p('q_diss_H', 36.0, '1/d', 'rate', FE_SOURCE), _p('q_diss_L', 36.0, '1/d', 'rate', FE_SOURCE), _p('K_Pbind', 37.2, 'kg P/m3', 'half_saturation', FE_SOURCE), _p('K_Pdiss', 0.93, 'kg P/m3', 'inhibition', FE_SOURCE), _p('MW_S', MW_S, 'kg/kmol', 'constant', ELEMENT_SOURCE), _p('MW_Fe', MW_Fe, 'kg/kmol', 'constant', ELEMENT_SOURCE), _p('COD_per_kmol_H2', 16.0, 'kg COD/kmol', 'constant', "2 e- per H2 x 8 kg COD per kmol e- (the twin's H2_COD_KG_PER_KMOL)"), _p('COD_per_kmol_S_reduced', 64.0, 'kg COD/kmol S', 'constant', "8 e- per sulfate-to-sulfide x 8 kg COD per kmol e- (the twin's _correct_srb_sulfur_mass_basis)"), _p('iP_HFO_HP', IP_HFO_HP, 'kg P/kg Fe', 'constant', TWIN_SOURCE), _p('iP_HFO_LP', IP_HFO_LP, 'kg P/kg Fe', 'constant', TWIN_SOURCE), _p('k_cryst_FeS', 100.0, '1/h', 'rate', FES_SOURCE), _p('n_cryst_FeS', 2.0, '1', 'rate', FES_SOURCE), _p('nu_FeS', 2.0, '1', 'constant', "FeS: two product ions (Fe2+, HS-); the twin's mineral nu convention"), _p('mineral_seed_FeS', 1e-06, 'kg/m3', 'constant', "the twin's nucleation-surrogate seed inventory"), _p('logK_FeS', -2.95, '1', 'constant', FES_SOURCE), _p('dH_FeS', -11000.0, 'J/mol', 'constant', FES_SOURCE))

def _one_minus(name: str) -> Expression:
    return add(1.0, multiply(-1.0, _s(name)))

def _sulfur(yield_name: str, cod_share: Expression) -> Expression:
    """kg S reduced per kg reference substrate: (1 - Y) x the oxidised share x MW_S / 64 (electron equivalents)."""
    return multiply(_one_minus(yield_name), cod_share, _s('MW_S'), power(_s('COD_per_kmol_S_reduced'), -1.0))

def _monod(s: str, k: str) -> Expression:
    from .expression import monod
    return monod(_s(s), _s(k))

def _srb_row(process_id: str, index: int, *, substrate: str, biomass: str, rate: str, half: str, so4_half: str, yield_name: str, window: str, h2s_ki: str, products: dict, oxidised: Expression, share: Optional[Expression]=None) -> BiokineticProcessDefinition:
    factors = [_s(rate), _s(biomass), _monod(substrate, half), _monod('S_SO4', so4_half)]
    if share is not None:
        factors.append(share)
    sulfur = _sulfur(yield_name, oxidised)
    stoich = {substrate: _st(x=-1.0)['x'], biomass: _s(yield_name), 'S_IS': sulfur, 'S_SO4': multiply(-1.0, sulfur)}
    for product, fraction in products.items():
        stoich[product] = multiply(_one_minus(yield_name), _s(fraction))
    return BiokineticProcessDefinition(process_id, index, multiply(*factors), substrate, stoich, _inhibition(_nutrients(), _hill_ph(window), _h2s(h2s_ki)))

def _hfo_by_h2(process_id: str, index: int, host: str, rate: Expression) -> BiokineticProcessDefinition:
    stoich = {'S_h2': multiply(-1.0, _s('COD_per_kmol_H2'), power(multiply(2.0, _s('MW_Fe')), -1.0)), 'S_Fe2': _st(x=1.0)['x'], host: _st(x=-1.0)['x']}
    return BiokineticProcessDefinition(process_id, index, multiply(rate, _s(host), _s('S_h2'), power(_s('COD_per_kmol_H2'), -1.0)), host, stoich)

def _hfo_by_sulfide(process_id: str, index: int, host: str, rate: Expression) -> BiokineticProcessDefinition:
    s_per_fe = multiply(_s('MW_S'), power(multiply(2.0, _s('MW_Fe')), -1.0))
    stoich = {'S_IS': multiply(-1.0, s_per_fe), 'S_S0': s_per_fe, 'S_Fe2': _st(x=1.0)['x'], host: _st(x=-1.0)['x']}
    return BiokineticProcessDefinition(process_id, index, multiply(rate, _s(host), _s('S_IS'), power(_s('MW_S'), -1.0)), host, stoich)

def _transfer(process_id: str, index: int, source: str, sink: str, rate: str) -> BiokineticProcessDefinition:
    return BiokineticProcessDefinition(process_id, index, multiply(_s(rate), _s(source)), source, {source: _st(x=-1.0)['x'], sink: _st(x=1.0)['x']})

def _p_exchange(process_id: str, index: int, host: str, bound: str, ip: str, rate: Expression, capture: bool) -> BiokineticProcessDefinition:
    """P co-precipitation/sorption (capture) or the bound host's dissolution, referenced to S_IP (carrier-closed)."""
    per_p = power(_s(ip), -1.0)
    sign = -1.0 if capture else 1.0
    return BiokineticProcessDefinition(process_id, index, rate, 'S_IP', {host: multiply(sign, per_p), bound: multiply(-sign, per_p)})

def _factor(history: str) -> str:
    return {'fresh_bound': 'hfo_bound_reduction_factor', 'aged_unbound': 'hfo_aged_reduction_factor', 'aged_bound': 'hfo_aged_bound_reduction_factor'}[history]
_HOST_REDUCTIONS = (('reduction_HFO_HP_h2', 'X_HFO_HP', 'h2', 'fresh_bound'), ('reduction_HFO_LP_h2', 'X_HFO_LP', 'h2', 'fresh_bound'), ('reduction_HFO_old_h2', 'X_HFO_old', 'h2', 'aged_unbound'), ('reduction_HFO_HP_old_h2', 'X_HFO_HP_old', 'h2', 'aged_bound'), ('reduction_HFO_LP_old_h2', 'X_HFO_LP_old', 'h2', 'aged_bound'), ('reduction_HFO_HP_IS', 'X_HFO_HP', 'sulfide', 'fresh_bound'), ('reduction_HFO_LP_IS', 'X_HFO_LP', 'sulfide', 'fresh_bound'), ('reduction_HFO_old_IS', 'X_HFO_old', 'sulfide', 'aged_unbound'), ('reduction_HFO_HP_old_IS', 'X_HFO_HP_old', 'sulfide', 'aged_bound'), ('reduction_HFO_LP_old_IS', 'X_HFO_LP_old', 'sulfide', 'aged_bound'))

def _host_reduction(index: int, spec: tuple) -> BiokineticProcessDefinition:
    process_id, host, donor, history = spec
    constant = 'k_Fe3t2_h2' if donor == 'h2' else 'k_Fe3t2_is'
    rate = multiply(_s(constant), _s(_factor(history)))
    build = _hfo_by_h2 if donor == 'h2' else _hfo_by_sulfide
    return build(process_id, index, host, rate)

def _noncompetitive(s: str, k: str) -> Expression:
    from .expression import noncompetitive_inhibition
    return noncompetitive_inhibition(_s(s), _s(k))
MADM1_SFE_PROCESSES: tuple[BiokineticProcessDefinition, ...] = (_srb_row('growth_SRB_h2', 25, substrate='S_h2', biomass='X_hSRB', rate='k_hSRB', half='K_hSRB', so4_half='K_so4_hSRB', yield_name='Y_hSRB', window='h2_SRB', h2s_ki='KI_h2s_hSRB', products={}, oxidised=1.0), _decay('decay_XhSRB', 26, 'b_hSRB', 'X_hSRB'), _srb_row('growth_SRB_acetate', 27, substrate='S_ac', biomass='X_aSRB', rate='k_aSRB', half='K_aSRB', so4_half='K_so4_aSRB', yield_name='Y_aSRB', window='ac_SRB', h2s_ki='KI_h2s_aSRB', products={}, oxidised=1.0), _decay('decay_XaSRB', 28, 'b_aSRB', 'X_aSRB'), _srb_row('growth_SRB_propionate', 29, substrate='S_pro', biomass='X_pSRB', rate='k_pSRB', half='K_pSRB', so4_half='K_so4_pSRB', yield_name='Y_pSRB', window='aa_SRB', h2s_ki='KI_h2s_pSRB', products={'S_ac': 'f_ac_pro'}, oxidised=_one_minus('f_ac_pro')), _decay('decay_XpSRB', 30, 'b_pSRB', 'X_pSRB'), _srb_row('growth_SRB_butyrate', 31, substrate='S_bu', biomass='X_c4SRB', rate='k_c4SRB', half='K_c4SRB', so4_half='K_so4_c4SRB', yield_name='Y_c4SRB', window='aa_SRB', h2s_ki='KI_h2s_c4SRB', products={'S_ac': 'f_ac_bu'}, oxidised=_one_minus('f_ac_bu'), share=_c4_share('S_bu', 'S_va')), _srb_row('growth_SRB_valerate', 32, substrate='S_va', biomass='X_c4SRB', rate='k_c4SRB', half='K_c4SRB', so4_half='K_so4_c4SRB', yield_name='Y_c4SRB', window='aa_SRB', h2s_ki='KI_h2s_c4SRB', products={'S_pro': 'f_pro_va', 'S_ac': 'f_ac_va'}, oxidised=add(1.0, multiply(-1.0, _s('f_pro_va')), multiply(-1.0, _s('f_ac_va'))), share=_c4_share('S_va', 'S_bu')), _decay('decay_Xc4SRB', 33, 'b_c4SRB', 'X_c4SRB'), _hfo_by_h2('reduction_HFO_H_h2', 34, 'X_HFO_H', _s('k_Fe3t2_h2')), _hfo_by_h2('reduction_HFO_L_h2', 35, 'X_HFO_L', _s('k_Fe3t2_h2')), _hfo_by_sulfide('reduction_HFO_H_IS', 36, 'X_HFO_H', _s('k_Fe3t2_is')), _hfo_by_sulfide('reduction_HFO_L_IS', 37, 'X_HFO_L', _s('k_Fe3t2_is')), _transfer('aging_HFO_H', 38, 'X_HFO_H', 'X_HFO_L', 'q_aging_H'), _transfer('aging_HFO_L', 39, 'X_HFO_L', 'X_HFO_old', 'q_aging_L'), _p_exchange('fast_P_binding', 40, 'X_HFO_H', 'X_HFO_HP', 'iP_HFO_HP', multiply(_s('q_Pcoprec'), _s('X_HFO_H'), _monod('S_IP', 'K_Pbind')), capture=True), _p_exchange('slow_P_sorption', 41, 'X_HFO_L', 'X_HFO_LP', 'iP_HFO_LP', multiply(_s('q_Pbinding'), _s('X_HFO_L'), _monod('S_IP', 'K_Pbind')), capture=True), _transfer('aging_HFO_HP', 42, 'X_HFO_HP', 'X_HFO_HP_old', 'q_aging_H'), _transfer('aging_HFO_LP', 43, 'X_HFO_LP', 'X_HFO_LP_old', 'q_aging_L'), _p_exchange('dissolution_HFO_HP', 44, 'X_HFO_H', 'X_HFO_HP', 'iP_HFO_HP', multiply(_s('q_diss_H'), _s('X_HFO_HP'), _noncompetitive('S_IP', 'K_Pdiss')), capture=False), _p_exchange('dissolution_HFO_LP', 45, 'X_HFO_L', 'X_HFO_LP', 'iP_HFO_LP', multiply(_s('q_diss_L'), _s('X_HFO_LP'), _noncompetitive('S_IP', 'K_Pdiss')), capture=False), *(_host_reduction(46 + i, spec) for i, spec in enumerate(_HOST_REDUCTIONS)))
MADM1_SFE_DERIVED_SYMBOLS = {'a_Fe2': "Fe2+ activity (mol/L) in the engine's own speciation", 'a_HS': "HS- activity (mol/L) in the engine's own speciation", 'Ksp_FeS': "FeS (amorphous) solubility constant at T (FeS + H+ = Fe2+ + HS-), the engine's own temperature law", 'mineral_smoothing': "engine-bound smoothing of the crystallisation law's max and sign (> 0; the twin binds 1e-30, numerically its exact piecewise law)", 'KH_h2s_T': "Henry's constant for H2S at T on the engine's liquid basis (kg S/m3/bar)", 'p_h2s_bar': 'H2S partial pressure in the headspace (bar)'}

class PhysicoChemicalProcessDefinition(BiokineticProcessDefinition):
    """A precipitation or gas-transfer process defined beside the biokinetic table: ``index`` is its TWIN index (it
    is not part of the kernel's contiguous process list)."""

def crystallisation_rate(saturation: Expression, k_per_hour: Expression, inventory: str, seed: str, nu: str, order: str, *, magnitude_floor: bool=False) -> Expression:
    """The twin's crystallisation law (kg mineral / m3 / d): k x X_eff x sign(D)|D|^n, D = max(SI, 0)^(1/nu) - 1, with
    X_eff = max(X, seed) when supersaturated and X when dissolving; ``k_per_hour`` is the twin's per-hour constant.

    ``magnitude_floor`` adds eps^2 under each |D|^n: in float 0.5 (x + sqrt(x^2 + eps^2)) is exactly 0 once x is a few
    eps below zero, where a non-integer order has no derivative (IPOPT's AMPL library refuses pow'(0, 1.868)). It moves
    the rate by k X eps^(2n) (about 1e-30 X /d at WaterTAP's eps); FeS (integer order) is built without it."""
    eps = _s('mineral_smoothing')
    si = smooth_positive(saturation, eps)
    drive = add(power(si, power(_s(nu), -1.0)), -1.0)
    k = multiply(24.0, k_per_hour)
    x_seeded = add(_s(inventory), smooth_positive(add(_s(seed), multiply(-1.0, _s(inventory))), eps))

    def magnitude(x: Expression) -> Expression:
        positive = smooth_positive(x, eps)
        return power(add(positive, multiply(eps, eps)) if magnitude_floor else positive, _s(order))
    growth = multiply(k, x_seeded, magnitude(drive))
    dissolution = multiply(-1.0, k, _s(inventory), magnitude(multiply(-1.0, drive)))
    return add(growth, dissolution)

def _fes_rate() -> Expression:
    return crystallisation_rate(multiply(_s('a_Fe2'), _s('a_HS'), power(multiply(_s('H_ion'), _s('Ksp_FeS')), -1.0)), _s('k_cryst_FeS'), 'X_FeS', 'mineral_seed_FeS', 'nu_FeS', 'n_cryst_FeS')
MADM1_FES_PRECIPITATION = PhysicoChemicalProcessDefinition('precipitation_FeS', 66, _fes_rate(), 'X_FeS', {'S_IS': multiply(-1.0, _s('MW_S'), power(add(_s('MW_S'), _s('MW_Fe')), -1.0)), 'S_Fe2': multiply(-1.0, _s('MW_Fe'), power(add(_s('MW_S'), _s('MW_Fe')), -1.0)), 'X_FeS': _st(x=1.0)['x']})
MADM1_H2S_TRANSFER = PhysicoChemicalProcessDefinition('IS_transfer', 74, multiply(_s('kLa'), add(_s('H2S_free_kg_S_m3'), multiply(-1.0, _s('KH_h2s_T'), _s('p_h2s_bar')))), 'S_IS', {'S_IS': _st(x=-1.0)['x']})
MADM1_SFE_KERNEL = BiokineticKernelDefinition(schema_version=SCHEMA_VERSION, kernel_id='madm1_biokinetics_sfe_v1', components=MADM1_CORE_CONTENTS + MADM1_SFE_CONTENTS, parameters=MADM1_CORE_PARAMETERS + MADM1_SFE_PARAMETERS, processes=MADM1_CORE_PROCESSES + MADM1_SFE_PROCESSES, derived_symbols={**MADM1_BIOKINETICS_KERNEL.derived_symbols, **MADM1_SFE_DERIVED_SYMBOLS}, reference='Literature-based modified ADM1 kernel (QSDsan lineage); kinetic parameters are not site fits.')

def sfe_kernel_parameter_values() -> dict[str, float]:
    """The extended kernel's parameter values by name (core + S/Fe)."""
    return MADM1_SFE_KERNEL.parameter_values()
