from __future__ import annotations
import math
import re
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple
FINITE_DONOR_RULING_ID = 'PUBLIC-RULE-8226e460e1d8'
SOURCES = {'half_reactions': 'Rittmann B.E., McCarty P.L. (2001) Environmental Biotechnology: Principles and Applications, McGraw-Hill, Ch. 2, half reactions and the general organic donor CnHaObNc (OpenAlex W1592200047); McCarty P.L. (1975) Stoichiometry of biological reactions, Progress in Water Technology 7:157-172 (OpenAlex W3165820020)', 'acceptor_sequence': 'Froelich P.N. et al. (1979) Early oxidation of organic matter in pelagic sediments of the eastern equatorial Atlantic: suboxic diagenesis, Geochim. Cosmochim. Acta 43:1075-1090, doi 10.1016/0016-7037(79)90095-4 (OpenAlex W2043668123): nitrate, then the metal oxides, then sulfate, in order of free-energy yield', 'cells': 'C5H7O2N, 20 electron equivalents per mole (Rittmann & McCarty 2001, Table 2.4)', 'cod': '8 g O2 per electron equivalent (1/4 mol O2 accepts one electron)'}
COD_G_PER_ELECTRON_MOL = 8.0
ACCEPTOR_SEQUENCE = ('nitrate_to_n2', 'ferric_to_ferrous', 'sulfate_to_sulfide')
ACCEPTOR_ELECTRONS = {'nitrate_to_n2': 5.0, 'ferric_to_ferrous': 1.0, 'sulfate_to_sulfide': 8.0}
CELL_FORMULA = {'C': 5, 'H': 7, 'O': 2, 'N': 1}
CELL_ELECTRONS_PER_MOL = 20.0
COMPLETE_REDUCTION_RESIDUE = 1e-09

class DonorInsufficientError(ValueError):
    """The stated donor cannot reduce an acceptor that must be reduced completely."""
_FORMULA = re.compile('([A-Z][a-z]?)(\\d*)')

def parse_organic_formula(formula: str) -> Dict[str, int]:
    """CcHhOoNn as element counts; refuses any other element, a repeated one or a formula that does not parse whole."""
    text = formula.strip()
    counts: Dict[str, int] = {}
    position = 0
    for match in _FORMULA.finditer(text):
        if match.start() != position:
            break
        element, number = (match.group(1), match.group(2))
        if element not in {'C', 'H', 'O', 'N'} or element in counts:
            raise ValueError(f'donor formula {formula!r}: only C, H, O and N, each once (got {element!r})')
        counts[element] = int(number) if number else 1
        position = match.end()
    if position != len(text) or not counts:
        raise ValueError(f'donor formula {formula!r} does not parse as CcHhOoNn')
    return {element: counts.get(element, 0) for element in ('C', 'H', 'O', 'N')}

@dataclass(frozen=True)
class DonorCompound:
    """An electron donor as its neutral molecule CcHhOoNn."""
    name: str
    formula: Mapping[str, int]
    source: str

    def __post_init__(self):
        if self.electrons_per_mol <= 0:
            raise ValueError(f'donor {self.name}: {self.formula_text} releases no electrons (4C + H - 2O - 3N <= 0)')
        if self.formula.get('H', 0) < 3 * self.formula.get('N', 0):
            raise ValueError(f'donor {self.name}: {self.formula_text} has fewer than 3 H per N, so its nitrogen cannot be written as released ammonia')

    @property
    def formula_text(self) -> str:
        return ''.join((f'{e}{('' if n == 1 else n)}' for e, n in self.formula.items() if n))

    @property
    def electrons_per_mol(self) -> float:
        f = self.formula
        return float(4 * f.get('C', 0) + f.get('H', 0) - 2 * f.get('O', 0) - 3 * f.get('N', 0))

    @property
    def cod_g_per_mol(self) -> float:
        return COD_G_PER_ELECTRON_MOL * self.electrons_per_mol

    def phreeqc_formula(self, ammonia_unit: str) -> str:
        """The molecule for a PHREEQC REACTION, its nitrogen written as the released ammonia unit (``ammonia_unit``
        stands for NH3: the decoupled ammonium element of `engineering_utils.phreeqc_batch_redox`)."""
        f = self.formula
        parts = [('C', f['C']), ('H', f['H'] - 3 * f['N']), ('O', f['O']), (ammonia_unit, f['N'])]
        return ''.join((f'{e}{('' if n == 1 else n)}' for e, n in parts if n))
_RM = 'Rittmann & McCarty (2001) Table 2.4'
DONOR_COMPOUNDS: Dict[str, DonorCompound] = {'acetate': DonorCompound('acetate', parse_organic_formula('C2H4O2'), f'acetic acid, C2H4O2 ({_RM}, acetate)'), 'propionate': DonorCompound('propionate', parse_organic_formula('C3H6O2'), f'propionic acid ({_RM})'), 'butyrate': DonorCompound('butyrate', parse_organic_formula('C4H8O2'), 'butyric acid, C4H8O2'), 'valerate': DonorCompound('valerate', parse_organic_formula('C5H10O2'), 'valeric acid, C5H10O2'), 'methanol': DonorCompound('methanol', parse_organic_formula('CH4O'), f'methanol ({_RM})'), 'ethanol': DonorCompound('ethanol', parse_organic_formula('C2H6O'), f'ethanol ({_RM})'), 'glucose': DonorCompound('glucose', parse_organic_formula('C6H12O6'), f'glucose ({_RM})'), 'carbohydrate': DonorCompound('carbohydrate', parse_organic_formula('CH2O'), f'carbohydrate CH2O ({_RM})'), 'hydrogen': DonorCompound('hydrogen', parse_organic_formula('H2'), f'hydrogen ({_RM})'), 'domestic_wastewater_organics': DonorCompound('domestic_wastewater_organics', parse_organic_formula('C10H19O3N'), f'domestic wastewater organics ({_RM})'), 'protein': DonorCompound('protein', parse_organic_formula('C16H24O5N4'), f'protein ({_RM})')}

def donor_compound(compound: Optional[str], formula: Optional[str]=None) -> DonorCompound:
    """A named donor (`DONOR_COMPOUNDS`), or the caller's own formula under its own name; refuses anything else."""
    if formula:
        name = compound or formula
        if compound in DONOR_COMPOUNDS and DONOR_COMPOUNDS[compound].formula != parse_organic_formula(formula):
            raise ValueError(f'donor {compound!r} is {DONOR_COMPOUNDS[compound].formula_text}, not {formula!r}: name a different compound or give no formula')
        return DonorCompound(name, parse_organic_formula(formula), "the caller's stated formula")
    if compound not in DONOR_COMPOUNDS:
        raise ValueError(f'donor {compound!r} is not a named donor ({sorted(DONOR_COMPOUNDS)}); give its formula CcHhOoNn instead')
    return DONOR_COMPOUNDS[compound]

@dataclass(frozen=True)
class AcceptorReduction:
    acceptor: str
    demand_mol: float
    reduced_mol: float
    acceptor_electron_mol: float
    donor_electron_mol: float
    cell_mol: float

    @property
    def left_mol(self) -> float:
        return max(self.demand_mol - self.reduced_mol, 0.0)

    @property
    def complete(self) -> bool:
        return self.demand_mol <= 0 or self.reduced_mol >= self.demand_mol * (1.0 - 2 * COMPLETE_REDUCTION_RESIDUE)

@dataclass(frozen=True)
class FiniteDonorReduction:
    """One allocation of a donor to the acceptors, in moles of the system it was asked for."""
    donor: DonorCompound
    donor_electron_mol: float
    cell_synthesis_fraction: Mapping[str, float]
    acceptors: Tuple[AcceptorReduction, ...] = field(default_factory=tuple)

    def acceptor(self, name: str) -> AcceptorReduction:
        return next((a for a in self.acceptors if a.acceptor == name))

    @property
    def donor_consumed_electron_mol(self) -> float:
        return math.fsum((a.donor_electron_mol for a in self.acceptors))

    @property
    def donor_left_electron_mol(self) -> float:
        return max(self.donor_electron_mol - self.donor_consumed_electron_mol, 0.0)

    @property
    def donor_consumed_mol(self) -> float:
        return self.donor_consumed_electron_mol / self.donor.electrons_per_mol

    @property
    def cell_mol(self) -> float:
        return math.fsum((a.cell_mol for a in self.acceptors))

    def element_changes(self, names: Mapping[str, str]) -> List[Tuple[str, float]]:
        """The reaction as (formula, moles added; negative removes) for a PHREEQC REACTION on decoupled elements.
        ``names``: ``nitrate``, ``sulfate``, ``sulfide`` and ``ammonia`` (the NH3 unit) element names. N2 leaves to the
        gas, so the nitrate's element is removed alone and its oxygen stays in the water; the sulfate's sulfur moves to
        the sulfide element; the cells leave as a solid, their nitrogen from ammonium. Ferric iron is reduced by the
        batch's own equilibrium at its fixed pe, its electrons supplied by the donor added here."""
        terms: List[Tuple[str, float]] = []
        if self.donor_consumed_mol > 0:
            terms.append((self.donor.phreeqc_formula(names['ammonia']), self.donor_consumed_mol))
        nitrate = self.acceptor('nitrate_to_n2').reduced_mol
        if nitrate > 0:
            terms.append((names['nitrate'], -nitrate))
        sulfate = self.acceptor('sulfate_to_sulfide').reduced_mol
        if sulfate > 0:
            terms += [(names['sulfate'], -sulfate), (names['sulfide'], sulfate)]
        if self.cell_mol > 0:
            terms.append(('C5H4O2' + names['ammonia'], -self.cell_mol))
        return terms

    def electron_balance_mol(self) -> float:
        """Donor electrons consumed less the electrons the acceptors and the cells took: zero by construction, kept as
        the identity the tests and the ledger state."""
        accepted = math.fsum((a.acceptor_electron_mol for a in self.acceptors)) + CELL_ELECTRONS_PER_MOL * self.cell_mol
        return self.donor_consumed_electron_mol - accepted

    def ledger(self, per_litre: float=1.0, *, available_cod_mg_l: Optional[float]=None) -> Dict[str, object]:
        """The answer's record, mmol (and mg/L) per litre: moles x ``per_litre`` (litres of feed per system mole
        basis, 1 for a per-litre allocation)."""
        mmol = lambda mol: mol * per_litre * 1000.0
        cod = lambda electrons: electrons * per_litre * COD_G_PER_ELECTRON_MOL * 1000.0
        acceptors = {a.acceptor: {'demand_mmol_per_l': mmol(a.demand_mol), 'reduced_mmol_per_l': mmol(a.reduced_mol), 'left_mmol_per_l': mmol(a.left_mol), 'complete': a.complete, 'electrons_mmol_per_l': mmol(a.acceptor_electron_mol), 'donor_electrons_mmol_per_l': mmol(a.donor_electron_mol), 'cell_synthesis_fraction': self.cell_synthesis_fraction.get(a.acceptor, 0.0), 'cells_mmol_per_l': mmol(a.cell_mol)} for a in self.acceptors}
        donor_n = self.donor.formula.get('N', 0) * self.donor_consumed_mol
        alkalinity = {'nitrate_to_n2': mmol(self.acceptor('nitrate_to_n2').reduced_mol), 'sulfate_to_sulfide': 2.0 * mmol(self.acceptor('sulfate_to_sulfide').reduced_mol), 'donor_nitrogen_released': mmol(donor_n), 'cell_nitrogen_assimilated': -mmol(self.cell_mol)}
        consumed_cod = cod(self.donor_consumed_electron_mol)
        to = {a.acceptor: cod(a.acceptor_electron_mol) for a in self.acceptors}
        to_cells = cod(CELL_ELECTRONS_PER_MOL * self.cell_mol)
        supplied = cod(self.donor_electron_mol) if available_cod_mg_l is None else available_cod_mg_l
        return {'method': 'finite_donor_half_reactions', 'ruling_id': FINITE_DONOR_RULING_ID, 'sources': dict(SOURCES), 'donor': {'name': self.donor.name, 'formula': self.donor.formula_text, 'source': self.donor.source, 'electrons_per_mol': self.donor.electrons_per_mol, 'cod_g_per_mol': self.donor.cod_g_per_mol, 'form': "neutral molecule (an acid as the acid); the anion's alkalinity is the feed analysis's"}, 'acceptor_sequence': list(ACCEPTOR_SEQUENCE), 'acceptors': acceptors, 'donor_electrons_mmol_per_l': {'supplied': mmol(self.donor_electron_mol), 'consumed': mmol(self.donor_consumed_electron_mol), 'left': mmol(self.donor_left_electron_mol)}, 'products_mmol_per_l': {'n2_as_n_to_gas': mmol(self.acceptor('nitrate_to_n2').reduced_mol), 'ferrous_iron': mmol(self.acceptor('ferric_to_ferrous').reduced_mol), 'sulfide': mmol(self.acceptor('sulfate_to_sulfide').reduced_mol), 'inorganic_carbon': mmol(self.donor.formula['C'] * self.donor_consumed_mol - CELL_FORMULA['C'] * self.cell_mol), 'cells_c5h7o2n': mmol(self.cell_mol), 'ammonium_released_by_donor_as_n': mmol(donor_n), 'ammonium_assimilated_as_n': mmol(self.cell_mol)}, 'cod_mg_l': {'supplied': supplied, 'consumed': consumed_cod, 'exerted_by_nitrate': to['nitrate_to_n2'], 'exerted_by_ferric': to['ferric_to_ferrous'], 'carried_by_sulfide': to['sulfate_to_sulfide'], 'carried_by_cells': to_cells, 'left_unoxidised': supplied - consumed_cod, 'closure_mg_l': consumed_cod - math.fsum(to.values()) - to_cells}, 'alkalinity_meq_per_l': {**alkalinity, 'total': math.fsum(alkalinity.values())}, 'not_modelled': ['methanogenesis and fermentation of the donor left over: its COD is reported as left, its carbon is not added to the water', 'kinetics: each acceptor is reduced as far as the donor reaches, at equilibrium']}

def reduce_with_finite_donor(donor: DonorCompound, donor_electron_mol: float, demand_mol: Mapping[str, float], cell_synthesis_fraction: Mapping[str, float], *, complete: Sequence[str]=()) -> FiniteDonorReduction:
    """Allocate ``donor_electron_mol`` to the acceptors of ``demand_mol`` (moles of each that could be reduced) in
    `ACCEPTOR_SEQUENCE`. An acceptor named in ``complete`` must be reduced fully, else `DonorInsufficientError` names
    the shortfall in donor COD. A complete reduction takes the total less `COMPLETE_REDUCTION_RESIDUE` of it."""
    if not math.isfinite(donor_electron_mol) or donor_electron_mol < 0:
        raise ValueError(f'donor electrons must be finite and >= 0, got {donor_electron_mol!r}')
    unknown = sorted(set(demand_mol) - set(ACCEPTOR_SEQUENCE))
    if unknown:
        raise ValueError(f'unknown acceptors {unknown}; known: {list(ACCEPTOR_SEQUENCE)}')
    remaining = float(donor_electron_mol)
    reductions: List[AcceptorReduction] = []
    for acceptor in ACCEPTOR_SEQUENCE:
        demand = float(demand_mol.get(acceptor, 0.0) or 0.0)
        if not math.isfinite(demand) or demand < 0:
            raise ValueError(f'{acceptor} demand must be finite and >= 0, got {demand!r}')
        fs = float(cell_synthesis_fraction.get(acceptor, math.nan)) if demand > 0 else float(cell_synthesis_fraction.get(acceptor, 0.0))
        if not 0.0 <= fs < 1.0:
            raise ValueError(f'cell_synthesis_fraction for {acceptor} must be stated, 0 <= fs < 1 (got {fs!r})')
        per_mol = ACCEPTOR_ELECTRONS[acceptor]
        target = demand * (1.0 - COMPLETE_REDUCTION_RESIDUE)
        needed = target * per_mol / (1.0 - fs)
        if needed <= remaining:
            reduced, donor_used = (target, needed)
        else:
            donor_used = remaining
            reduced = remaining * (1.0 - fs) / per_mol
            if acceptor in complete:
                short = (needed - remaining) * COD_G_PER_ELECTRON_MOL * 1000.0
                raise DonorInsufficientError(f'the stated donor cannot reduce all the {acceptor.replace('_to_', ' to ')}: it is short by {short:.6g} mg COD in the batch ({FINITE_DONOR_RULING_ID})')
        remaining -= donor_used
        accepted = reduced * per_mol
        reductions.append(AcceptorReduction(acceptor=acceptor, demand_mol=demand, reduced_mol=reduced, acceptor_electron_mol=accepted, donor_electron_mol=donor_used, cell_mol=donor_used * fs / CELL_ELECTRONS_PER_MOL))
    return FiniteDonorReduction(donor=donor, donor_electron_mol=float(donor_electron_mol), cell_synthesis_fraction=dict(cell_synthesis_fraction), acceptors=tuple(reductions))

def donor_electron_mol(cod_mg_l: float, litres: float=1.0) -> float:
    """Electron moles a donor's COD (mg O2/L) supplies over ``litres``."""
    return float(cod_mg_l) / 1000.0 / COD_G_PER_ELECTRON_MOL * float(litres)
