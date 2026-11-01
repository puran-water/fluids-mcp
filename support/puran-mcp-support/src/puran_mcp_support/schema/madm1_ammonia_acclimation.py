from __future__ import annotations
import hashlib
import json
import math
from dataclasses import dataclass
from typing import Mapping, Optional
from .madm1_biokinetics import MADM1_BIOKINETICS_KERNEL
_MW_N = MADM1_BIOKINETICS_KERNEL.parameter('MW_N').value

class FreeAmmoniaBasisError(ValueError):
    """An unknown free-ammonia basis, or a kernel whose ammonia form a KI basis does not fit."""

@dataclass(frozen=True)
class FreeAmmoniaBasis:
    """One located KI_nh3 (mol NH3-N/L) for the kernel's acetoclast term, with the culture and model it was found in."""
    basis_id: str
    ki_nh3_m: float
    role: str
    culture: str
    model_context: str
    sources: tuple[str, ...]
    standing: str = 'proposed_pending_ruling'

    def __post_init__(self) -> None:
        if not (math.isfinite(self.ki_nh3_m) and self.ki_nh3_m > 0):
            raise FreeAmmoniaBasisError(f'{self.basis_id}: KI_nh3 must be finite and positive')
        if not self.sources:
            raise FreeAmmoniaBasisError(f'{self.basis_id} carries no source')

    @property
    def ki_nh3_mg_n_l(self) -> float:
        return self.ki_nh3_m * _MW_N * 1000.0

    @property
    def sha256(self) -> str:
        return hashlib.sha256(json.dumps(vars(self), sort_keys=True, allow_nan=False).encode()).hexdigest()

    def record(self) -> dict:
        return {'basis': self.basis_id, 'basis_sha256': self.sha256, 'standing': self.standing, 'role': self.role, 'KI_nh3_M': self.ki_nh3_m, 'KI_nh3_mg_N_L': self.ki_nh3_mg_n_l, 'kernel_default_KI_nh3_M': kernel_default_ki_nh3_m(), 'culture': self.culture, 'model_context': self.model_context, 'sources': list(self.sources)}
_CT2021 = 'Capson-Tojo, Astals, Robles et al. 2021, Bioresource Technology 341:125802, doi:10.1016/j.biortech.2021.125802, OpenAlex W3195170166'
_CT2020 = 'Capson-Tojo, Moscoviz, Astals et al. 2020, Renewable and Sustainable Energy Reviews 117:109487, doi:10.1016/j.rser.2019.109487, OpenAlex W2980745179 (abstract read; its value as quoted by Capson-Tojo et al. 2021, the full text being behind a bot check)'
_BONK = 'Bonk, Popp, Weinrich et al. 2018, Frontiers in Microbiology 9:2921, doi:10.3389/fmicb.2018.02921, OpenAlex W2902036433'
FREE_AMMONIA_BASES: Mapping[str, FreeAmmoniaBasis] = {b.basis_id: b for b in (FreeAmmoniaBasis('capson_tojo_2021_default_adm1_fit', 0.003, 'bracket_low', 'inoculum adapted to TAN 7.27 g N/L at pH 8.10, food-waste batches at 37 C', 'the default ADM1 (acetoclasts the only acetate route) calibrated: KI_nh3 0.0030 M (up to 0.0035 M); the paper states a correct ionic-strength FAN calculation would put it considerably higher (section 3.2)', (f'{_CT2021}: Table 3 (KI,NH3 0.0030 M, control and GAC reactors), sections 3.1 and 3.3',)), FreeAmmoniaBasis('capson_tojo_2020_fan_adapted', 0.0057, 'bracket_central', 'FAN-adapted processes (1590 data points of over 50 studies recast on a modified Davies FAN)', "the 50 % inhibitory limit of acetoclastic methanogenesis for FAN-adapted processes (4.3e-4 and 0.0109 M the threshold function's start and full inhibition), applied in the authors' model with SAO", (f"{_CT2021}: section 3.3 ('estimated around 0.0057 M ... (Capson-Tojo et al., 2020)')", _CT2020)), FreeAmmoniaBasis('bonk_2018_single_population', 0.008, 'bracket_high', 'VFA-degrading community acclimating in a CSTR at 37 C, influent TAN raised 52 to 277 mM', "ADM1 with one acetoclast population (I,nh3,ac adjusted to 0.008 M; 'the value of 0.0018 leads to process failure'); the paper's two-population fit gives Methanosaeta 0.0052 M and Methanosarcina 0.3387 M", (f'{_BONK}: Table 1 and its footnote',)))}
UNMOVED_STATEMENT = "the kernel's KI_nh3, unmoved (ADM1 default, a culture unadapted to high free ammonia): no free-ammonia basis named (process_model_schema.madm1_ammonia_acclimation)"

def kernel_default_ki_nh3_m() -> float:
    return float(MADM1_BIOKINETICS_KERNEL.parameter('KI_nh3').value)

def free_ammonia_basis(basis_id: str) -> FreeAmmoniaBasis:
    """The named basis; refused unless the kernel's ammonia term is the non-competitive n = 1 form a KI basis fits."""
    if basis_id not in FREE_AMMONIA_BASES:
        raise FreeAmmoniaBasisError(f'unknown free-ammonia basis {basis_id!r}; known: {sorted(FREE_AMMONIA_BASES)}')
    if float(MADM1_BIOKINETICS_KERNEL.parameter('n_nh3').value) != 1.0:
        raise FreeAmmoniaBasisError('a KI basis is a 50 % inhibitory concentration only on the n_nh3 = 1 form')
    return FREE_AMMONIA_BASES[basis_id]

def free_ammonia_record(basis_id: Optional[str]) -> dict:
    """What an engine stamps beside its answer; ``basis_id`` None is the kernel's own KI_nh3."""
    if basis_id is None:
        return {'basis': None, 'statement': UNMOVED_STATEMENT, 'KI_nh3_M': kernel_default_ki_nh3_m()}
    return free_ammonia_basis(basis_id).record()
