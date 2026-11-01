"""Shared acid-base constants for reduced Fe-S-P speciation adapters."""
from __future__ import annotations
import math
from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping
REFERENCE_TEMPERATURE_K = 298.15
GAS_CONSTANT_J_MOL_K = 8.31446261815324

@dataclass(frozen=True)
class AcidDissociationDefinition:
    """One acid dissociation constant and its temperature law: van't Hoff from ``delta_h_j_mol`` (None: held at 25 C),
    or, where ``analytic`` is given, the published law ``log10 Ka = A1 + A2 T + A3/T + A4 log10 T + A5/T^2`` (T in K).
    ``validity_k`` is the law's stated range (refused outside it); None where the range has not been read, which
    ``validity_note`` says."""
    constant_id: str
    pka_25c: float
    delta_h_j_mol: float | None
    source: str
    analytic: tuple | None = None
    validity_k: tuple | None = None
    validity_note: str = ''

    @property
    def ka_25c(self) -> float:
        return 10.0 ** (-self.pka_25c)

    def log10_ka(self, temperature_k: float) -> float:
        temperature_k = float(temperature_k)
        if not math.isfinite(temperature_k) or temperature_k <= 0.0:
            raise ValueError('temperature_k must be finite and positive')
        if self.validity_k is not None and (not self.validity_k[0] <= temperature_k <= self.validity_k[1]):
            raise ValueError(f"{self.constant_id}: {temperature_k} K is outside the law's stated range {self.validity_k} K ({self.source}); not extrapolated")
        return self.log10_ka_law(temperature_k, math.log10)

    def log10_ka_law(self, t, log10):
        """The temperature law itself on a caller's ``t`` (K, dimensionless: a float, or a Pyomo expression with the
        caller's ``log10``), so every engine renders THIS definition rather than restating it. No range check here:
        a numeric caller goes through ``log10_ka``; an equation-oriented caller checks its temperature bounds."""
        if self.analytic is not None:
            a1, a2, a3, a4, a5 = self.analytic
            return a1 + a2 * t + a3 / t + a4 * log10(t) + a5 / (t * t)
        if self.delta_h_j_mol is None:
            return -self.pka_25c
        return -self.pka_25c + self.delta_h_j_mol / (GAS_CONSTANT_J_MOL_K * math.log(10.0)) * (1.0 / REFERENCE_TEMPERATURE_K - 1.0 / t)

    def pka(self, temperature_k: float) -> float:
        return -self.log10_ka(temperature_k)

    def ka(self, temperature_k: float) -> float:
        temperature_k = float(temperature_k)
        if not math.isfinite(temperature_k) or temperature_k <= 0.0:
            raise ValueError('temperature_k must be finite and positive')
        if self.analytic is None and self.validity_k is None:
            if self.delta_h_j_mol is None:
                return self.ka_25c
            return self.ka_25c * math.exp(self.delta_h_j_mol / GAS_CONSTANT_J_MOL_K * (1.0 / REFERENCE_TEMPERATURE_K - 1.0 / temperature_k))
        return 10.0 ** self.log10_ka(temperature_k)
_PCM_MATCH = 'Flores-Alsina/Solon mADM1 PCM convention'
_PB82 = 'Plummer and Busenberg 1982'
_STUMM_MORGAN = 'Stumm and Morgan 1996; unresolved enthalpy held at 25 C'
_ACID_RULING = 'HOUSE-MADM1-ACID-CONSTANTS-EVALUATED-LITERATURE-20261008'
_MILLERO_H2S = f"Morse, Millero, Cornwell & Rickard 1987, Earth-Science Reviews 24:1-42, doi:10.1016/0012-8252(87)90046-8, eq. 9 (Millero's evaluated curve): pKa1 = 32.55 + 1519.44/T - 15.672 log10 T + 0.02722 T, T in K [V via review]; {_ACID_RULING}"
_NEA_PO4 = f'pKa3 12.35 at 25 C with the NEA-evaluated dissociation enthalpy +14.6 +/- 3.8 kJ/mol (OECD-NEA, Chemical Thermodynamics of Uranium, phosphate auxiliary data, pp. 387-388) [V via review]; {_ACID_RULING}'
ACID_DISSOCIATION_CONSTANTS: Mapping[str, AcidDissociationDefinition] = MappingProxyType({'water': AcidDissociationDefinition('water', 14.0, 55900.0, _PCM_MATCH), 'ammonium': AcidDissociationDefinition('ammonium', 9.25, 51965.0, _PCM_MATCH), 'carbonic_1': AcidDissociationDefinition('carbonic_1', 6.35, 7646.0, _PCM_MATCH), 'carbonic_2': AcidDissociationDefinition('carbonic_2', 10.33, 14850.0, _PB82), 'sulfide_1': AcidDissociationDefinition('sulfide_1', 32.55 + 1519.44 / REFERENCE_TEMPERATURE_K - 15.672 * math.log10(REFERENCE_TEMPERATURE_K) + 0.02722 * REFERENCE_TEMPERATURE_K, None, _MILLERO_H2S, analytic=(-32.55, -0.02722, -1519.44, 15.672, 0.0), validity_k=None, validity_note="the paper's stated temperature range was not read (it was not opened here; [V via review]); applicability to thermophilic digestion (55 C) is unverified"), 'sulfide_2': AcidDissociationDefinition('sulfide_2', 14.0, None, _STUMM_MORGAN), 'phosphoric_1': AcidDissociationDefinition('phosphoric_1', 2.148, None, _STUMM_MORGAN), 'phosphoric_2': AcidDissociationDefinition('phosphoric_2', 7.2, 3600.0, _PCM_MATCH), 'phosphoric_3': AcidDissociationDefinition('phosphoric_3', 12.35, 14600.0, _NEA_PO4)})

def reduced_free_ion_fractions(*, hydrogen_activity: float, temperature_k: float) -> dict[str, float]:
    """Return free PO4, HS and CO3 fractions of their analytical totals."""
    h = float(hydrogen_activity)
    if not math.isfinite(h) or h <= 0.0:
        raise ValueError('hydrogen_activity must be finite and positive')
    constants = ACID_DISSOCIATION_CONSTANTS
    kp1 = constants['phosphoric_1'].ka(temperature_k)
    kp2 = constants['phosphoric_2'].ka(temperature_k)
    kp3 = constants['phosphoric_3'].ka(temperature_k)
    ks1 = constants['sulfide_1'].ka(temperature_k)
    ks2 = constants['sulfide_2'].ka(temperature_k)
    kc1 = constants['carbonic_1'].ka(temperature_k)
    kc2 = constants['carbonic_2'].ka(temperature_k)
    return {'PO4': kp1 * kp2 * kp3 / (h ** 3 + kp1 * h ** 2 + kp1 * kp2 * h + kp1 * kp2 * kp3), 'HS_minus': ks1 * h / (h ** 2 + ks1 * h + ks1 * ks2), 'CO3': kc1 * kc2 / (h ** 2 + kc1 * h + kc1 * kc2)}
__all__ = ['ACID_DISSOCIATION_CONSTANTS', 'GAS_CONSTANT_J_MOL_K', 'REFERENCE_TEMPERATURE_K', 'AcidDissociationDefinition', 'reduced_free_ion_fractions']
