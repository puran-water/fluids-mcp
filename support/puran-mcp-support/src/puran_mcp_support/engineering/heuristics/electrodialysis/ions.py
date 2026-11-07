"""The ion table the EDR model and the watertap-engine-mcp multi-ion ED profile share (one source, two consumers).

Charge and formula weight per species; tracer diffusion coefficients at 25 C are phreeqc.dat's own ``-dw`` values
(phreeqpython 1.18 bundle, the house database `phreeqc_phases.HOUSE_DATABASE`), so the transport model and PHREEQC's
specific-conductance model read one diffusivity table. Electrical mobility is the Einstein relation
``u = |z| D F / (R T)``, the relation the engine's ``MCAS_BRACKISH_5ION`` preset already uses.

Neutral solutes (SiO2, dissolved CO2, neutral organic carbon) carry z = 0: they do not migrate under the field and
cross a membrane only by diffusion (``NEUTRAL_MEMBRANE_DIFFUSIVITY_M2_S``), which is the mechanistic content of the
"organics and silica are not concentrated" claim the model tests. TOC carries the formula weight of carbon (mass as C);
its diffusivity is an [ASSUMED] small-organic value (no phreeqc.dat entry), stated in every result that carries TOC.
Anionic organics (humic/fulvic acids) are NOT represented: they migrate to and foul the AEM, flagged as a sensitivity.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Iterable, Mapping
FARADAY = 96485.33212
GAS_CONSTANT = 8.314462618
WATER_MW = 0.018015268
WATER_MOLAR_VOLUME_M3_MOL = 1.8069e-05

@dataclass(frozen=True)
class Species:
    key: str
    mcas: str
    charge: int
    mw_g_mol: float
    diffusivity_m2_s: float
    phreeqc: str
    phreeqc_as: str | None = None
    source: str = 'phreeqc.dat -dw'
SPECIES: Dict[str, Species] = {s.key: s for s in (Species('Na', 'Na_+', 1, 22.98977, 1.33e-09, 'Na'), Species('K', 'K_+', 1, 39.0983, 1.96e-09, 'K'), Species('NH4', 'NH4_+', 1, 18.03846, 1.98e-09, 'N(-3)', 'NH4'), Species('Ca', 'Ca_2+', 2, 40.078, 7.93e-10, 'Ca'), Species('Mg', 'Mg_2+', 2, 24.305, 7.05e-10, 'Mg'), Species('Sr', 'Sr_2+', 2, 87.62, 7.94e-10, 'Sr'), Species('Ba', 'Ba_2+', 2, 137.327, 8.48e-10, 'Ba'), Species('Cl', 'Cl_-', -1, 35.453, 2.03e-09, 'Cl'), Species('SO4', 'SO4_2-', -2, 96.0626, 1.07e-09, 'S(6)', 'SO4'), Species('HCO3', 'HCO3_-', -1, 61.0168, 1.18e-09, 'C(4)', 'HCO3'), Species('CO3', 'CO3_2-', -2, 60.0089, 9.55e-10, 'C(4)', 'CO3'), Species('NO3', 'NO3_-', -1, 62.0049, 1.9e-09, 'N(5)', 'NO3'), Species('F', 'F_-', -1, 18.9984, 1.46e-09, 'F'), Species('H2PO4', 'H2PO4_-', -1, 96.9872, 8.46e-10, 'P', 'H2PO4'), Species('HPO4', 'HPO4_2-', -2, 95.9793, 6.9e-10, 'P', 'HPO4'), Species('SiO2', 'SiO2', 0, 60.0843, 1.1e-09, 'Si', 'SiO2', 'phreeqc.dat -dw of H4SiO4'), Species('CO2', 'CO2', 0, 44.0095, 1.92e-09, 'C(4)', 'CO2'), Species('TOC', 'TOC', 0, 12.011, 5e-10, '', None, '[ASSUMED] small neutral organic; no phreeqc.dat entry'))}
IONS = tuple((k for k, s in SPECIES.items() if s.charge != 0))
CATIONS = tuple((k for k, s in SPECIES.items() if s.charge > 0))
ANIONS = tuple((k for k, s in SPECIES.items() if s.charge < 0))
NEUTRALS = tuple((k for k, s in SPECIES.items() if s.charge == 0))
NEUTRAL_MEMBRANE_DIFFUSIVITY_M2_S = 1e-13

def species(key: str) -> Species:
    try:
        return SPECIES[key]
    except KeyError as exc:
        raise KeyError(f'{key} is not an EDR model species; known: {', '.join(SPECIES)}') from exc

def mobility_m2_v_s(key: str, temperature_k: float=298.15) -> float:
    s = species(key)
    return abs(s.charge) * s.diffusivity_m2_s * FARADAY / (GAS_CONSTANT * temperature_k)

def mcas_package_config(keys: Iterable[str], temperature_k: float=298.15, *, phase_keys: bool=False) -> dict:
    """The MCAS ``property_package_config`` for a species set: solute list, charge (ions only), formula weight in
    kg/mol, ``diffusivity_data`` for every solute, ``elec_mobility_data`` for ions (Einstein relation). Data keys are
    the bare MCAS names (the engine's JSON-safe session form) unless ``phase_keys`` asks for ``("Liq", name)``, the form
    ``MCASParameterBlock`` takes directly."""
    keys = list(keys)
    for key in keys:
        species(key)

    def k(name):
        return ('Liq', name) if phase_keys else name
    return {'solute_list': [SPECIES[j].mcas for j in keys], 'charge': {SPECIES[j].mcas: SPECIES[j].charge for j in keys if SPECIES[j].charge}, 'mw_data': {SPECIES[j].mcas: SPECIES[j].mw_g_mol / 1000.0 for j in keys}, 'diffusivity_data': {k(SPECIES[j].mcas): SPECIES[j].diffusivity_m2_s for j in keys}, 'elec_mobility_data': {k(SPECIES[j].mcas): mobility_m2_v_s(j, temperature_k) for j in keys if SPECIES[j].charge}}

def mg_l_to_mol_m3(values: Mapping[str, float]) -> Dict[str, float]:
    return {k: float(v) / species(k).mw_g_mol for k, v in values.items()}

def mol_m3_to_mg_l(values: Mapping[str, float]) -> Dict[str, float]:
    return {k: float(v) * species(k).mw_g_mol for k, v in values.items()}

def charge_balance_meq_l(values_mg_l: Mapping[str, float]) -> float:
    """Cation minus anion equivalents (meq/L) of a mg/L composition."""
    return sum((float(v) / species(k).mw_g_mol * species(k).charge for k, v in values_mg_l.items()))

def tds_mg_l(values_mg_l: Mapping[str, float]) -> float:
    """Sum of dissolved species; dissolved CO2 (MEGA's TDS row excludes it, as a gravimetric TDS does) and the neutral
    organics (carried in COD units under TOC) are not in it."""
    return sum((float(v) for k, v in values_mg_l.items() if k not in ('CO2', 'TOC')))
MEGA_ION_NAMES = {'Na+': 'Na', 'K+': 'K', 'NH4+': 'NH4', 'Ca2+': 'Ca', 'Mg2+': 'Mg', 'Ba2+': 'Ba', 'Sr2+': 'Sr', 'F-': 'F', 'Cl-': 'Cl', 'HCO3-': 'HCO3', 'CO3 2-': 'CO3', 'NO3-': 'NO3', 'H2PO4-': 'H2PO4', 'HPO4 2-': 'HPO4', 'SO4 2-': 'SO4', 'SiO2': 'SiO2', 'CO2': 'CO2'}
MEGA_UNMODELLED = ('Br-', 'PO4 3-')

def from_mega_composition(row: Mapping[str, float]) -> Dict[str, float]:
    """A MEGA projection composition row (mg/kg) in house keys; zero rows kept only where the model carries them. A
    non-zero row the model cannot carry is refused by name."""
    out: Dict[str, float] = {}
    for name, value in row.items():
        if name in ('TDS', 'pH'):
            continue
        if name in MEGA_UNMODELLED:
            if float(value) != 0.0:
                raise ValueError(f'MEGA row {name}={value} is non-zero and the EDR model does not carry it')
            continue
        out[MEGA_ION_NAMES[name]] = float(value)
    return out
