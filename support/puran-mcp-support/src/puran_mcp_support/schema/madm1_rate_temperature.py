from __future__ import annotations
import hashlib
import json
import math
from dataclasses import dataclass
from typing import Iterable, Mapping, Optional
KERNEL_RATE_REFERENCE_K = 308.15
RATE_TEMPERATURE_GROUPS: Mapping[str, tuple[str, ...]] = {'hydrolysis': ('q_ch_hyd', 'q_pr_hyd', 'q_li_hyd'), 'acidogenesis': ('k_su', 'k_aa', 'k_fa'), 'acetogenesis': ('k_c4', 'k_pro'), 'methanogenesis': ('k_ac', 'k_h2')}
RATE_TEMPERATURE_PARAMETERS: tuple[str, ...] = tuple((p for group in RATE_TEMPERATURE_GROUPS.values() for p in group))
NOT_MOVED = {'half_saturation': 'K_S: no located mesophilic exponent (Dalby et al. 2021 Eq 13 moves K_S; exponent not stated)', 'decay': 'b: held by Bergland et al. 2015 section 2.1.1 and Dalby et al. 2021 S3', 'yield': 'Y: no located source moves it', 'KI_nh3': "no located primary source; free NH3 rises with T through the kernel's NH4+ acid constant", 'extension_guilds': 'S/Fe, iron-reduction and micro-aeration constants: no basis carries them'}

class RateTemperatureError(ValueError):
    """A rate-temperature basis is unknown, or a temperature lies outside its validity range."""

@dataclass(frozen=True)
class LinearToOptimum:
    """The BioModel steady-state temperature function (Kovalovszki et al. 2020, Eq. 2): below the optimum,
    mu(T) = mu_opt - alpha (T_opt - T); above it, mu_opt (T_max - T)/(T_max - T_opt). One source group's row."""
    group: str
    mu_opt_h: float
    alpha_per_c_h: float
    t_opt_c: float
    t_max_c: float

    def rate(self, temperature_c: float) -> float:
        if temperature_c < self.t_opt_c:
            return self.mu_opt_h - self.alpha_per_c_h * (self.t_opt_c - temperature_c)
        return self.mu_opt_h * (self.t_max_c - temperature_c) / (self.t_max_c - self.t_opt_c)

@dataclass(frozen=True)
class RateTemperatureBasis:
    """A named rate-temperature function: each moved kernel rate constant's law, the range it holds over, its sources."""
    basis_id: str
    laws: Mapping[str, LinearToOptimum]
    mapping_basis: str
    valid_range_k: tuple[float, float]
    range_basis: str
    sources: tuple[str, ...]
    reference_k: float = KERNEL_RATE_REFERENCE_K
    standing: str = 'proposed_pending_ruling'

    def __post_init__(self) -> None:
        if set(self.laws) != set(RATE_TEMPERATURE_PARAMETERS):
            raise RateTemperatureError(f'{self.basis_id}: a basis states a law for exactly {RATE_TEMPERATURE_PARAMETERS}')
        low, high = self.valid_range_k
        if not low <= self.reference_k <= high:
            raise RateTemperatureError(f"{self.basis_id}: its range {self.valid_range_k} excludes the kernel's reference")
        if not self.sources:
            raise RateTemperatureError(f'{self.basis_id} carries no source')
        for name, law in self.laws.items():
            for t_k in (low, self.reference_k, high):
                if not law.rate(t_k - 273.15) > 0:
                    raise RateTemperatureError(f'{self.basis_id}.{name}: the law is not positive at {t_k} K')

    def factors(self, temperature_k: float) -> dict[str, float]:
        """{kernel rate constant: k(T)/k(reference)}: exactly 1 at the reference; refused outside the range."""
        t = float(temperature_k)
        low, high = self.valid_range_k
        if not (math.isfinite(t) and low <= t <= high):
            raise RateTemperatureError(f"{t!r} K lies outside {self.basis_id}'s validity range {low:g}-{high:g} K ({low - 273.15:g}-{high - 273.15:g} C): {self.range_basis}")
        if t == self.reference_k:
            return {name: 1.0 for name in self.laws}
        return {name: law.rate(t - 273.15) / law.rate(self.reference_k - 273.15) for name, law in self.laws.items()}

    @property
    def sha256(self) -> str:
        payload = {'basis_id': self.basis_id, 'reference_k': self.reference_k, 'valid_range_k': list(self.valid_range_k), 'laws': {k: vars(v) for k, v in sorted(self.laws.items())}, 'mapping_basis': self.mapping_basis, 'range_basis': self.range_basis, 'sources': list(self.sources), 'standing': self.standing}
        return hashlib.sha256(json.dumps(payload, sort_keys=True, allow_nan=False).encode()).hexdigest()
_KOVALOVSZKI = 'Kovalovszki, Treu, Ellegaard, Luo & Angelidaki 2020, Applied Energy 263:114646, doi:10.1016/j.apenergy.2020.114646, OpenAlex W3008634309: Eq. 2 (the BioModel steady-state function FT_ref) and Table 2 (mu_max, alpha, T_opt, T_max per group), accepted manuscript (University of Padua repository); manure CSTRs at 37 C shifted to 25, 45 and 55 C (study S1) and at 53 C shifted to 40 and 60 C (study S2)'
_TABLE2 = {'carbohydrate degraders': LinearToOptimum('carbohydrate degraders', 0.212, 0.0039, 55.0, 65.0), 'amino acid degraders': LinearToOptimum('amino acid degraders', 0.266, 0.0049, 55.0, 65.0), 'lipid degraders': LinearToOptimum('lipid degraders', 0.022, 0.00041, 55.0, 65.0), 'long-chain fatty acid degraders': LinearToOptimum('long-chain fatty acid degraders', 0.023, 0.00042, 55.0, 65.0), 'propionic acid degraders': LinearToOptimum('propionic acid degraders', 0.02, 0.00037, 53.0, 65.0), 'butyric acid degraders': LinearToOptimum('butyric acid degraders', 0.028, 0.00052, 60.0, 70.0), 'acetic acid degraders/methanogens': LinearToOptimum('acetic acid degraders/methanogens', 0.026, 0.00047, 55.0, 65.0)}
_KOVALOVSZKI_MAPPING = {'q_ch_hyd': 'carbohydrate degraders', 'q_pr_hyd': 'amino acid degraders', 'q_li_hyd': 'lipid degraders', 'k_su': 'carbohydrate degraders', 'k_aa': 'amino acid degraders', 'k_fa': 'long-chain fatty acid degraders', 'k_c4': 'butyric acid degraders', 'k_pro': 'propionic acid degraders', 'k_ac': 'acetic acid degraders/methanogens', 'k_h2': 'acetic acid degraders/methanogens'}
KOVALOVSZKI_2020_BIOMODEL_LINEAR = RateTemperatureBasis(basis_id='kovalovszki2020_biomodel_linear', laws={name: _TABLE2[group] for name, group in _KOVALOVSZKI_MAPPING.items()}, mapping_basis="each kernel rate on the Table 2 group that carries its substrate; hydrolysis on the degraders of its polymer (BioModel models no separate hydrolysing guild); k_c4 on the butyric acid degraders (the valeric row differs by 1.6 % in alpha/mu_max); k_h2 on the acetic acid degraders/methanogens, the BioModel's only methanogen group (it carries no hydrogen state) [ASSUMED mapping]. Every row's alpha/mu_max is 0.0181-0.0186 per C, so the mapping moves a factor by under 1 % in range", valid_range_k=(298.15, 313.15), range_basis="25 C: the lowest temperature S1 ran (Kovalovszki et al. 2020 section 2.2.1). 40 C: the mesophilic upper edge; above it the community shifts, which 'cannot be reproduced by using a single model microbial group for acetic acid conversion into methane' (Kovalovszki et al. 2020 section 3.2), and transient instability appears at 42-44 C (Westerholm et al. 2018, Applied Energy 226:838, doi:10.1016/j.apenergy.2018.06.045, W2807853537, abstract); maximum biogas at 40 C on cow manure (Nandi et al. 2020, Sustainability 12:8772, doi:10.3390/su12218772, W3093903933, abstract)", sources=(_KOVALOVSZKI,))
RATE_TEMPERATURE_BASES: Mapping[str, RateTemperatureBasis] = {KOVALOVSZKI_2020_BIOMODEL_LINEAR.basis_id: KOVALOVSZKI_2020_BIOMODEL_LINEAR}
UNMOVED_STATEMENT = "the kernel's 35 C rate set, unmoved: no rate-temperature basis named (process_model_schema.madm1_rate_temperature)"

def rate_temperature_basis(basis_id: str) -> RateTemperatureBasis:
    if basis_id not in RATE_TEMPERATURE_BASES:
        raise RateTemperatureError(f'unknown rate-temperature basis {basis_id!r}; known: {sorted(RATE_TEMPERATURE_BASES)}')
    return RATE_TEMPERATURE_BASES[basis_id]

def rate_temperature_record(basis_id: Optional[str], temperature_k: Optional[float]=None, *, hold: Iterable[str]=(), hold_reason: str='') -> dict:
    """What an engine stamps beside its answer: the basis, the temperature, each rate constant's factor and sources.
    ``hold`` names rate constants kept at factor 1 (an engine rebinding a measured rate with its own temperature), with
    ``hold_reason``. ``basis_id`` None is the unmoved kernel set."""
    if basis_id is None:
        return {'basis': None, 'statement': UNMOVED_STATEMENT, 'reference_k': KERNEL_RATE_REFERENCE_K, 'factors': {name: 1.0 for name in RATE_TEMPERATURE_PARAMETERS}}
    basis = rate_temperature_basis(basis_id)
    held = tuple(hold)
    unknown = sorted(set(held) - set(RATE_TEMPERATURE_PARAMETERS))
    if unknown:
        raise RateTemperatureError(f'held rate constants {unknown} are not moved by any basis')
    if held and (not hold_reason.strip()):
        raise RateTemperatureError('a held rate constant states why it is held')
    factors = basis.factors(float(temperature_k))
    factors.update({name: 1.0 for name in held})
    return {'basis': basis.basis_id, 'basis_sha256': basis.sha256, 'standing': basis.standing, 'temperature_k': float(temperature_k), 'reference_k': basis.reference_k, 'valid_range_k': list(basis.valid_range_k), 'factors': factors, 'held': {name: hold_reason for name in held}, 'not_moved': dict(NOT_MOVED), 'mapping_basis': basis.mapping_basis, 'sources': list(basis.sources)}
