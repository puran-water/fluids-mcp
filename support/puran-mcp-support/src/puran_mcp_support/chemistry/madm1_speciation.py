'Davies speciation of an mADM1 liquor: the ONE acid-base / activity model both engines solve.\n\nOwner ruling HOUSE-MADM1-DAVIES-BOTH-ENGINES-20261007 (Hersh, 2026-10-07, AskUserQuestion), verbatim answer: "Davies\nfor both (Recommended)" -- the option read: the twin moves to Davies; that changes its outputs, so the change replays\nthe last 14 stored nights and stops before landing if an alarm would move.\n\nThe equations are the WaterTAP house Davies block\'s (``watertap-engine-mcp/chemistry/modified_adm1_davies.py`` over\n``modified_adm1_phase1.py``), written once as plain arithmetic so the QSDsan twin evaluates the same model:\n\n* pH is ``-log10(a_H)`` on the MOLAL basis of the engine\'s own water inventory (``water_kg_m3``): a molar\n  concentration c (mol/L) is c x 1000 / water_kg_m3 mol/kg;\n* Davies: ``log10 g_z = -A z^2 (sqrt(I)/(1 + sqrt(I)) - b I)``, I molal, A(T) from Kell (1975) water density and\n  Malmberg & Maryott (1956) dielectric constant; neutral species have g = 1;\n* weak acids (each Ka at T by van\'t Hoff from the kernel\'s pKa / dH and R): water, NH4+/NH3, CO2/HCO3-/CO3-2,\n  H2PO4-/HPO4-2, H2S/HS-, and the four VFAs; every ratio is written with activities, the balances with\n  concentrations, exactly as the WaterTAP block writes them;\n* the fully dissociated (strong) ions enter as two sums the caller forms on its own state basis: their net charge\n  ``sum z c`` and their ionic-strength weight ``sum z^2 c``. WaterTAP passes its charge-equivalent lumps (with its\n  explicit Ca/SO4 divalent corrections); the twin passes its explicit Na, K, Mg, Ca, Fe2, Fe3, Al, Cl and SO4.\n\nThe solve is deterministic: an outer fixed point on I (from the strong ions alone, no cache) around a bracketed\nBrent solve of electroneutrality in pH; it converges to ``ionic_strength_tol`` and refuses otherwise.'
from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Mapping
__all__ = ['DAVIES_MAX_IONIC_STRENGTH_MOLAL', 'ION_PAIRS', 'IonPair', 'phosphate_third_log10_ka', 'MadmLiquor', 'ion_pair_log_k', 'pair_alkalinity_M', 'acid_constants', 'davies_A', 'davies_log_gamma', 'solve_madm1_speciation', 'species_at']
LN10 = math.log(10.0)
DAVIES_MAX_IONIC_STRENGTH_MOLAL = 0.5
DAVIES_TEMPERATURE_C = (0.0, 50.0)
_ACIDS = ('w', 'nh4', 'co2', 'co3', 'h2po4', 'ac', 'pro', 'bu', 'va')
VFA_COD_PER_KMOL = {'ac': 64.0, 'pro': 112.0, 'bu': 160.0, 'va': 208.0}

@dataclass(frozen=True)
class MadmLiquor:
    """Molar totals (mol/L) of one liquor, each on the caller's own state basis."""
    strong_charge_M: float
    strong_ionic_M: float
    IN: float = 0.0
    IP: float = 0.0
    IC: float = 0.0
    IS: float = 0.0
    ac: float = 0.0
    pro: float = 0.0
    bu: float = 0.0
    va: float = 0.0
    Ca: float = 0.0
    Mg: float = 0.0
    Fe2: float = 0.0
    SO4: float = 0.0

@dataclass(frozen=True)
class IonPair:
    """One aqueous complex: ``reactants`` (species key, count) -> the complex, log K at T (``log_k``)."""
    name: str
    charge: int
    reactants: tuple
    log_k_25: float
    delta_h_kj_mol: float
    analytic: tuple | None
    source: str
    minteq_v4: str

    def log_k(self, temperature_k: float) -> float:
        t = float(temperature_k)
        if self.analytic is not None:
            a1, a2, a3, a4, a5 = self.analytic
            return a1 + a2 * t + a3 / t + a4 * math.log10(t) + a5 / (t * t)
        return self.log_k_25 - self.delta_h_kj_mol * 1000.0 / (LN10 * PHREEQC_R_J_MOL_K) * (1.0 / t - 1.0 / 298.15)
PHREEQC_R_J_MOL_K = 8.31446261815324
_KCAL = 4.184
_PQ = "phreeqc.dat (phreeqpython's bundled copy, walk_components.PHREEQC_DATABASE)"
_MQ = 'minteq.v4.dat (USGS phreeqc3 distribution)'
_NO_DH = "; phreeqc.dat states no delta_h: minteq.v4's (0 kJ) is used, which is PHREEQC's own default"
ION_PAIRS: tuple = (IonPair('CaHCO3+', 1, (('Ca', 1), ('CO3', 1), ('H', 1)), 11.435, -0.871 * _KCAL, (1317.0071, 0.34546894, -39916.84, -517.70761, 563713.9), f'{_PQ} line 336: Ca+2 + CO3-2 + H+ = CaHCO3+, -log_k 11.435, -delta_h -0.871 kcal, -analytic', f'{_MQ} line 4824: log_k 11.599, delta_h 5.4 kJ'), IonPair('CaHPO4', 0, (('Ca', 1), ('HPO4', 1)), 2.739, 3.3 * _KCAL, None, f'{_PQ} line 354: Ca+2 + HPO4-2 = CaHPO4, -log_k 2.739, -delta_h 3.3 kcal', f'{_MQ} line 4432: Ca+2 + H+ + PO4-3 = CaHPO4, log_k 15.035, delta_h -3 kJ'), IonPair('MgHCO3+', 1, (('Mg', 1), ('CO3', 1), ('H', 1)), 11.399, -2.771 * _KCAL, (48.6721, 0.03252849, -2614.335, -18.00263, 563713.9), f'{_PQ} line 376: Mg+2 + H+ + CO3-2 = MgHCO3+, -log_k 11.399, -delta_h -2.771 kcal, -analytic', f'{_MQ} line 4816: log_k 11.339, delta_h -10.6 kJ'), IonPair('MgSO4', 0, (('Mg', 1), ('SO4', 1)), 2.37, 4.55 * _KCAL, None, f'{_PQ} line 383: Mg+2 + SO4-2 = MgSO4, -log_k 2.37, -delta_h 4.550 kcal', f'{_MQ} line 3232: log_k 2.26, delta_h 5.8 kJ'), IonPair('MgHPO4', 0, (('Mg', 1), ('HPO4', 1)), 2.87, 3.3 * _KCAL, None, f'{_PQ} line 392: Mg+2 + HPO4-2 = MgHPO4, -log_k 2.87, -delta_h 3.3 kcal', f'{_MQ} line 4424: Mg+2 + H+ + PO4-3 = MgHPO4, log_k 15.175, delta_h -3 kJ'), IonPair('MgH2PO4+', 1, (('Mg', 1), ('H2PO4', 1)), 1.513, 3.4 * _KCAL, None, f'{_PQ} line 395: Mg+2 + H2PO4- = MgH2PO4+, -log_k 1.513, -delta_h 3.4 kcal', f'{_MQ} line 4416: Mg+2 + 2 H+ + PO4-3 = MgH2PO4+, log_k 21.2561, delta_h -4.6861 kJ'), IonPair('NH4SO4-', -1, (('NH4', 1), ('SO4', 1)), 1.11, 0.0, None, f'{_PQ} line 285: NH4+ + SO4-2 = NH4SO4-, -log_k 1.11{_NO_DH}', f'{_MQ} line 2968: log_k 1.03, delta_h 0 kJ'), IonPair('FeHCO3+', 1, (('Fe2', 1), ('HCO3', 1)), 2.0, 0.0, None, f'{_PQ} line 452: Fe+2 + HCO3- = FeHCO3+, -log_k 2.0{_NO_DH}', f'{_MQ} line 4760: Fe+2 + H+ + CO3-2 = FeHCO3+, log_k 11.429, delta_h 0 kJ'), IonPair('FeCO3', 0, (('Fe2', 1), ('CO3', 1)), 4.38, 0.0, None, f'{_PQ} line 450: Fe+2 + CO3-2 = FeCO3, -log_k 4.38{_NO_DH} (minteq.v4 carries no aqueous FeCO3)', f'{_MQ}: none'), IonPair('Fe(HS)2', 0, (('Fe2', 1), ('HS', 2)), 8.95, 0.0, None, f'{_PQ} line 460: Fe+2 + 2HS- = Fe(HS)2, -log_k 8.95{_NO_DH}', f'{_MQ} line 2846: log_k 8.95, delta_h 0 kJ'), IonPair('Fe(HS)3-', -1, (('Fe2', 1), ('HS', 3)), 10.987, 0.0, None, f'{_PQ} line 462: Fe+2 + 3HS- = Fe(HS)3-, -log_k 10.987{_NO_DH}', f'{_MQ} line 2854: log_k 10.987, delta_h 0 kJ'), IonPair('FeHPO4', 0, (('Fe2', 1), ('HPO4', 1)), 3.6, 0.0, None, f'{_PQ} line 464: Fe+2 + HPO4-2 = FeHPO4, -log_k 3.6{_NO_DH}', f'{_MQ} line 4296: Fe+2 + H+ + PO4-3 = FeHPO4, log_k 15.975, delta_h 0 kJ'), IonPair('FeH2PO4+', 1, (('Fe2', 1), ('H2PO4', 1)), 2.7, 0.0, None, f'{_PQ} line 466: Fe+2 + H2PO4- = FeH2PO4+, -log_k 2.7{_NO_DH}', f'{_MQ} line 4288: Fe+2 + 2 H+ + PO4-3 = FeH2PO4+, log_k 22.273, delta_h 0 kJ'))
_KEY_CHARGE = {'Ca': 2, 'Mg': 2, 'Fe2': 2, 'NH4': 1, 'SO4': -2, 'HCO3': -1, 'CO3': -2, 'HPO4': -2, 'H2PO4': -1, 'HS': -1, 'H': 1}
_TOTAL_KEYS = ('Ca', 'Mg', 'Fe2', 'SO4')
_LIGAND_TOTAL = {'HCO3': 'IC', 'CO3': 'IC', 'HPO4': 'IP', 'H2PO4': 'IP', 'HS': 'IS', 'NH4': 'IN'}

def phosphate_third_log10_ka(temperature_k: float) -> float:
    """log10 Ka3 (HPO4-2 = PO4-3 + H+) at T: the ONE shared definition (process_model_schema.speciation
    "phosphoric_3": pKa3 12.35 at 25 C with the NEA-evaluated +14.6 kJ/mol). PO4-3 is reported for mineral saturation
    only (below 1e-5 of the phosphate at digester pH, outside the balances, as both engines' phosphate always was)."""
    from puran_mcp_support.schema.speciation import ACID_DISSOCIATION_CONSTANTS
    return ACID_DISSOCIATION_CONSTANTS['phosphoric_3'].log10_ka(temperature_k)

def ion_pair_log_k(temperature_k: float) -> dict[str, float]:
    return {p.name: p.log_k(temperature_k) for p in ION_PAIRS}
_ALKALINITY_REFERENCE = {'HCO3': 1, 'CO3': 2, 'HPO4': 1, 'H2PO4': 0, 'HS': 1, 'NH4': 0, 'SO4': 0, 'Ca': 0, 'Mg': 0, 'Fe2': 0}

def pair_alkalinity_M(species: Mapping[str, float]) -> float:
    """The acid-neutralising capacity (eq/L) the ion pairs hold: each complex counts the ligand forms it is made of on
    the proton reference, less any proton its reaction consumes (CaHCO3+ = Ca + CO3 + H counts 1, as HCO3- does)."""
    total = 0.0
    for p in ION_PAIRS:
        c = species.get(p.name, 0.0)
        if c:
            total += c * sum((n * (-1 if key == 'H' else _ALKALINITY_REFERENCE[key]) for key, n in p.reactants))
    return total

def acid_constants(temperature_k: float, values: Mapping[str, float]) -> dict[str, float]:
    """pKa at ``temperature_k`` for every acid: ln K(T) = ln K(T_base) + dH/R (1/T_base - 1/T)."""
    correction = 1.0 / float(values['T_base']) - 1.0 / float(temperature_k)
    r_gas = float(values['R_gas'])
    from puran_mcp_support.schema.speciation import ACID_DISSOCIATION_CONSTANTS
    out = {'h2s': ACID_DISSOCIATION_CONSTANTS['sulfide_1'].pka(temperature_k)}
    for acid in _ACIDS:
        ln_k = -LN10 * float(values[f'pKa_{acid}']) + float(values[f'dH_Ka_{acid}']) / r_gas * correction
        out[acid] = -ln_k / LN10
    return out

def davies_A(temperature_k: float) -> float:
    """Debye-Hueckel A (kg^0.5 mol^-0.5) from water density (Kell 1975) and dielectric constant (Malmberg 1956)."""
    tc = float(temperature_k) - 273.15
    density_kg_l = (999.83952 + 16.945176 * tc - 0.0079870401 * tc ** 2 - 4.6170461e-05 * tc ** 3 + 1.0556302e-07 * tc ** 4 - 2.8054253e-10 * tc ** 5) / (1 + 0.01687985 * tc) / 1000
    dielectric = 87.74 - 0.40008 * tc + 0.0009398 * tc ** 2 - 1.41e-06 * tc ** 3
    return 1824830.0 * math.sqrt(density_kg_l) / (dielectric * float(temperature_k)) ** 1.5

def davies_log_gamma(z: int, ionic_strength_molal: float, a: float, b: float) -> float:
    root = math.sqrt(max(ionic_strength_molal, 0.0))
    return -a * z * z * (root / (1.0 + root) - b * ionic_strength_molal)

def _ligands(ph: float, free: Mapping[str, float], *, pk: Mapping[str, float], lg1: float, lg2: float) -> dict:
    """The acid-base species (mol/L) of the FREE ligand totals at activity pH ``ph``."""
    r_nh3 = 10.0 ** (ph - pk['nh4'] + lg1)
    nh4 = free['IN'] / (1.0 + r_nh3)
    r_hco3 = 10.0 ** (ph - pk['co2'] - lg1)
    r_co3 = 10.0 ** (ph - pk['co3'] + lg1 - lg2)
    co2 = free['IC'] / (1.0 + r_hco3 * (1.0 + r_co3))
    hco3 = co2 * r_hco3
    r_hpo4 = 10.0 ** (ph - pk['h2po4'] + lg1 - lg2)
    h2po4 = free['IP'] / (1.0 + r_hpo4)
    r_hs = 10.0 ** (ph - pk['h2s'] - lg1)
    h2s = free['IS'] / (1.0 + r_hs)
    return {'NH4': nh4, 'NH3': free['IN'] - nh4, 'CO2': co2, 'HCO3': hco3, 'CO3': hco3 * r_co3, 'H2PO4': h2po4, 'HPO4': h2po4 * r_hpo4, 'H2S': h2s, 'HS': h2s * r_hs}
_PAIR_TOL = 1e-15
_PAIR_MAX_ITER = 500
_PAIR_PLAIN_ITER = 60
_PAIR_TOTAL_KEYS = ('Ca', 'Mg', 'Fe2', 'SO4', 'IN', 'IC', 'IP', 'IS')

def _solve_linear(a: list[list[float]], b: list[float]) -> list[float]:
    """Gaussian elimination with partial pivoting (the pair Newton step; at most 8 unknowns)."""
    n = len(b)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(m[r][col]))
        if m[piv][col] == 0.0:
            raise ValueError('the ion-pair Newton system is singular')
        m[col], m[piv] = (m[piv], m[col])
        for r in range(col + 1, n):
            f = m[r][col] / m[col][col]
            if f:
                for c in range(col, n + 1):
                    m[r][c] -= f * m[col][c]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (m[r][n] - sum((m[r][c] * x[c] for c in range(r + 1, n)))) / m[r][r]
    return x

def _pair_structure() -> tuple:
    """Per pair, once: (name, charge, protons, the reactant keys with their counts and totals, the exponent on each
    total). A complex is ``coef * prod_t free_t ** N_t`` at fixed pH and gammas (each ligand form is a fixed fraction
    of its free total there)."""
    out = []
    for p in ION_PAIRS:
        protons, parts, expo = (0, [], {})
        for key, n in p.reactants:
            if key == 'H':
                protons += n
                continue
            total = key if key in _TOTAL_KEYS else _LIGAND_TOTAL[key]
            parts.append((key, n, total))
            expo[total] = expo.get(total, 0) + n
        out.append((p.name, p.charge, protons, tuple(parts), tuple(expo.items())))
    return tuple(out)
_PAIR_STRUCTURE = _pair_structure()
_UNIT_LIGANDS = {'IN': 1.0, 'IC': 1.0, 'IP': 1.0, 'IS': 1.0}

def _pairs(ph: float, liquor: MadmLiquor, *, pk: Mapping[str, float], log_gamma: Mapping[int, float], molar_to_molal: float, warm: dict | None=None) -> tuple[dict, dict, dict]:
    """``(ligand species, free totals, complexes)`` (mol/L) at fixed pH and gammas: the metal and ligand mass balances
    with the ION_PAIRS equilibria, by substitution on each free total ``x = T / (1 + sum n c_p / x)`` (a complex's
    concentration is proportional to its own free metal and a power of its ligand), then Newton (``_pairs_newton``)
    if it has not converged in ``_PAIR_PLAIN_ITER`` steps. No pairing total -> no iteration. The per-pair constants
    are formed once per call (``_PAIR_STRUCTURE``): the twin calls this inside every rate evaluation. ``warm``: a
    caller's dict of free totals to start from (the previous trial pH's answer); it receives this call's answer. The
    fixed point is the same from any start; only the steps to it change."""
    lg1, lg2 = (log_gamma[1], log_gamma[2])
    totals = {'IN': liquor.IN, 'IC': liquor.IC, 'IP': liquor.IP, 'IS': liquor.IS, 'Ca': liquor.Ca, 'Mg': liquor.Mg, 'Fe2': liquor.Fe2, 'SO4': liquor.SO4}
    free = dict(totals)
    if not (liquor.Ca or liquor.Mg or liquor.Fe2 or liquor.SO4):
        return (_ligands(ph, free, pk=pk, lg1=lg1, lg2=lg2), free, {p.name: 0.0 for p in ION_PAIRS})
    log_k = pk['ion_pairs']
    gamma = {z: 10.0 ** log_gamma[abs(z)] if z else 1.0 for z in (-3, -2, -1, 0, 1, 2, 3)}
    a_h = 10.0 ** (-ph)
    unit = _ligands(ph, _UNIT_LIGANDS, pk=pk, lg1=lg1, lg2=lg2)
    pairs = []
    for name, charge, protons, parts, expo in _PAIR_STRUCTURE:
        coef = 10.0 ** (log_k[name] - math.log10(gamma[charge] * molar_to_molal)) * a_h ** protons
        for key, n, total in parts:
            if totals[total] <= 0.0:
                coef = 0.0
                break
            fraction = 1.0 if key in _TOTAL_KEYS else unit[key]
            coef *= (gamma[_KEY_CHARGE[key]] * fraction * molar_to_molal) ** n
        pairs.append((name, coef, expo))
    keys = [k for k in _PAIR_TOTAL_KEYS if totals[k] > 0.0]
    if warm:
        for k in keys:
            if 0.0 < warm.get(k, 0.0) <= totals[k]:
                free[k] = warm[k]
    for step in range(_PAIR_MAX_ITER):
        complexes, held = ({}, dict.fromkeys(keys, 0.0))
        for name, coef, expo in pairs:
            c = coef
            if c:
                for total, n in expo:
                    c *= free[total] ** n
                for total, n in expo:
                    held[total] += n * c
            complexes[name] = c
        change = 0.0
        for total in keys:
            new = totals[total] / (1.0 + held[total] / free[total])
            change = max(change, abs(new - free[total]) / totals[total])
            free[total] = new
        if change <= _PAIR_TOL:
            if warm is not None:
                warm.update(free)
            return (_ligands(ph, free, pk=pk, lg1=lg1, lg2=lg2), free, complexes)
        if step + 1 >= _PAIR_PLAIN_ITER:
            out = _pairs_newton(ph, totals, free, pairs, keys, pk=pk, lg1=lg1, lg2=lg2)
            if warm is not None:
                warm.update(out[1])
            return out
    raise ValueError(f'the ion-pair mass balances did not converge at pH {ph:.6g} (last relative change {change:.3g})')

def _pairs_newton(ph, totals, free, pairs, keys, *, pk, lg1, lg2):
    """Newton on ln(free total) for the same mass balances (the substitution's fixed point): R_t = free_t +
    sum_p N_pt c_p - T_t and dR_t/d ln free_s = delta_ts free_t + sum_p N_pt N_ps c_p. Steps are limited to a factor
    e^2 and free_t <= T_t."""
    index = {t: i for i, t in enumerate(keys)}
    change = math.inf
    for _ in range(_PAIR_MAX_ITER):
        n_keys = len(keys)
        jac = [[0.0] * n_keys for _ in range(n_keys)]
        residual = [free[t] - totals[t] for t in keys]
        complexes = {}
        for name, coef, expo in pairs:
            c = coef
            if c:
                for total, n in expo:
                    c *= free[total] ** n
                for total, n in expo:
                    i = index[total]
                    residual[i] += n * c
                    for total2, n2 in expo:
                        jac[i][index[total2]] += n * n2 * c
            complexes[name] = c
        for t in keys:
            jac[index[t]][index[t]] += free[t]
        step = _solve_linear(jac, [-r for r in residual])
        change = 0.0
        for t, du in zip(keys, step):
            du = max(-2.0, min(2.0, du))
            new = min(free[t] * math.exp(du), totals[t])
            change = max(change, abs(new - free[t]) / totals[t])
            free[t] = new
        if change <= _PAIR_TOL or max((abs(r) / totals[t] for r, t in zip(residual, keys))) <= _PAIR_TOL:
            complexes = {}
            for name, coef, expo in pairs:
                c = coef
                if c:
                    for total, n in expo:
                        c *= free[total] ** n
                complexes[name] = c
            return (_ligands(ph, free, pk=pk, lg1=lg1, lg2=lg2), free, complexes)
    raise ValueError(f'the ion-pair mass balances did not converge at pH {ph:.6g} (Newton, last relative change {change:.3g})')

def species_at(ph: float, liquor: MadmLiquor, *, pk: Mapping[str, float], log_gamma: Mapping[int, float], molar_to_molal: float, warm: dict | None=None) -> dict[str, float]:
    """Every species (mol/L) at activity pH ``ph`` and fixed activity coefficients: free ligand species, the free
    pairing metals and sulfate (``Ca``, ``Mg``, ``Fe2``, ``SO4``), and each ION_PAIRS complex by its name."""
    lg1 = log_gamma[1]
    h = 10.0 ** (-ph - lg1) / molar_to_molal
    oh = 10.0 ** (ph - pk['w'] - lg1) / molar_to_molal
    lig, free, complexes = _pairs(ph, liquor, pk=pk, log_gamma=log_gamma, molar_to_molal=molar_to_molal, warm=warm)
    out = {'H': h, 'OH': oh, **lig, 'Ca': free['Ca'], 'Mg': free['Mg'], 'Fe2': free['Fe2'], 'SO4': free['SO4'], **complexes}
    if 'po4_log_k' in pk:
        out['PO4'] = 10.0 ** (pk['po4_log_k'] + ph + log_gamma[2] - log_gamma[3]) * lig['HPO4']
    for acid in ('ac', 'pro', 'bu', 'va'):
        total = getattr(liquor, acid)
        out[acid + '_ion'] = total / (1.0 + 10.0 ** (-ph + lg1 + pk[acid]))
    return out

def _pair_corrections(s: Mapping[str, float]) -> tuple[float, float]:
    """(charge, sum z^2 c) the complexes add, net of what they take from a TOTAL in the strong sums (Ca, Mg, Fe2, SO4:
    the free ligand species already exclude their share)."""
    charge = ionic = 0.0
    for p in ION_PAIRS:
        c = s.get(p.name, 0.0)
        if not c:
            continue
        from_totals = [(key, n) for key, n in p.reactants if key in _TOTAL_KEYS]
        charge += (p.charge - sum((n * _KEY_CHARGE[key] for key, n in from_totals))) * c
        ionic += (p.charge ** 2 - sum((n * _KEY_CHARGE[key] ** 2 for key, n in from_totals))) * c
    return (charge, ionic)

def _charge(s: Mapping[str, float], liquor: MadmLiquor) -> float:
    return liquor.strong_charge_M + s['H'] + s['NH4'] - s['OH'] - s['HCO3'] - 2.0 * s['CO3'] - s['H2PO4'] - 2.0 * s['HPO4'] - s['HS'] - s['ac_ion'] - s['pro_ion'] - s['bu_ion'] - s['va_ion'] + _pair_corrections(s)[0]

def _ionic_strength(s: Mapping[str, float], liquor: MadmLiquor, molar_to_molal: float) -> float:
    return 0.5 * molar_to_molal * (liquor.strong_ionic_M + s['H'] + s['OH'] + s['NH4'] + s['HCO3'] + 4.0 * s['CO3'] + s['H2PO4'] + 4.0 * s['HPO4'] + s['HS'] + s['ac_ion'] + s['pro_ion'] + s['bu_ion'] + s['va_ion'] + _pair_corrections(s)[1])

def _brent(f, a: float, b: float, *, xtol: float=1e-14, maxiter: int=200) -> float:
    """Brent's method on a bracketing interval (f(a) f(b) <= 0)."""
    fa, fb = (f(a), f(b))
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if fa * fb > 0.0:
        raise ValueError(f'electroneutrality is not bracketed on pH [{a}, {b}]: {fa:.3e}, {fb:.3e}')
    c, fc, d, e = (a, fa, b - a, b - a)
    for _ in range(maxiter):
        if fb * fc > 0.0:
            c, fc, d, e = (a, fa, b - a, b - a)
        if abs(fc) < abs(fb):
            a, b, c = (b, c, b)
            fa, fb, fc = (fb, fc, fb)
        tol = 2.0 * 2.2e-16 * abs(b) + 0.5 * xtol
        m = 0.5 * (c - b)
        if abs(m) <= tol or fb == 0.0:
            return b
        if abs(e) >= tol and abs(fa) > abs(fb):
            s = fb / fa
            if a == c:
                p, q = (2.0 * m * s, 1.0 - s)
            else:
                q, r = (fa / fc, fb / fc)
                p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)
            if p > 0.0:
                q = -q
            p = abs(p)
            if 2.0 * p < min(3.0 * m * q - abs(tol * q), abs(e * q)):
                e, d = (d, p / q)
            else:
                d, e = (m, m)
        else:
            d, e = (m, m)
        a, fa = (b, fb)
        b += d if abs(d) > tol else math.copysign(tol, m)
        fb = f(b)
    raise ValueError("Brent's method did not converge on electroneutrality")

def solve_madm1_speciation(liquor: MadmLiquor, *, temperature_k: float, water_kg_m3: float, values: Mapping[str, float], ionic_strength_tol: float=1e-13, max_outer: int=60) -> dict:
    """Solve electroneutrality with Davies activities. Returns pH (activity, molal), every species (mol/L), I
    (mol/kg water), log10 gamma for z = 1, 2, 3, the acid pKa used and the domain flag (never clipped)."""
    if not water_kg_m3 > 0:
        raise ValueError(f'the water inventory must be positive, not {water_kg_m3!r} kg/m3')
    pk = {**acid_constants(temperature_k, values), 'ion_pairs': ion_pair_log_k(temperature_k), 'po4_log_k': phosphate_third_log10_ka(temperature_k)}
    a = davies_A(temperature_k)
    b = float(values['davies_b'])
    molar_to_molal = 1000.0 / float(water_kg_m3)
    ionic = 0.5 * molar_to_molal * liquor.strong_ionic_M
    warm: dict = {}
    ph = None

    def root(log_gamma, near):

        def charge(x):
            return _charge(species_at(x, liquor, pk=pk, log_gamma=log_gamma, molar_to_molal=molar_to_molal, warm=warm), liquor)
        if near is not None:
            try:
                return _brent(charge, near - 0.05, near + 0.05)
            except ValueError:
                pass
        return _brent(charge, -2.0, 16.0)
    previous = None
    for _ in range(max_outer):
        log_gamma = {z: davies_log_gamma(z, ionic, a, b) for z in (1, 2, 3)}
        ph = root(log_gamma, ph)
        species = species_at(ph, liquor, pk=pk, log_gamma=log_gamma, molar_to_molal=molar_to_molal, warm=warm)
        updated = _ionic_strength(species, liquor, molar_to_molal)
        if abs(updated - ionic) <= ionic_strength_tol * max(1.0, updated):
            ionic = updated
            log_gamma = {z: davies_log_gamma(z, ionic, a, b) for z in (1, 2, 3)}
            ph = root(log_gamma, ph)
            species = species_at(ph, liquor, pk=pk, log_gamma=log_gamma, molar_to_molal=molar_to_molal, warm=warm)
            tc = float(temperature_k) - 273.15
            return {'pH': ph, 'species_M': species, 'ionic_strength_molal': ionic, 'log_gamma': log_gamma, 'pKa': pk, 'A': a, 'molar_to_molal': molar_to_molal, 'in_domain': ionic <= DAVIES_MAX_IONIC_STRENGTH_MOLAL and DAVIES_TEMPERATURE_C[0] <= tc <= DAVIES_TEMPERATURE_C[1], 'model': 'madm1_davies_molal_ion_pairs_v1' if liquor.Ca or liquor.Mg or liquor.Fe2 or liquor.SO4 else 'madm1_davies_molal_v1'}
        nxt = updated
        if previous is not None:
            g0, g1 = (previous[1] - previous[0], updated - ionic)
            if g1 != g0:
                cand = ionic - g1 * (ionic - previous[0]) / (g1 - g0)
                if cand > 0.0 and math.isfinite(cand):
                    nxt = cand
        previous = (ionic, updated)
        ionic = nxt
    raise ValueError(f'Davies ionic strength did not converge in {max_outer} iterations (last I {ionic:.6g})')
