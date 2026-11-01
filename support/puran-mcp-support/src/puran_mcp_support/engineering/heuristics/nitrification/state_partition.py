from __future__ import annotations
from dataclasses import dataclass, field
from typing import Mapping

class IonBasisError(ValueError):
    """Raised when the charge-balance basis cannot be established from measurements."""

@dataclass(frozen=True)
class IonBasis:
    """Measured (or explicitly banded) ions and inorganic carbon, in mg/L.

    ``S_IC`` is total inorganic carbon as C -- **not** alkalinity, and not derived
    from it.
    """
    S_K: float
    S_Mg: float
    S_Ca: float
    S_Na: float
    S_Cl: float
    S_IC: float
    provenance: str
    not_identifiable: tuple[str, ...] = field(default_factory=tuple)

    def as_concentrations(self) -> dict[str, float]:
        return {'S_K': self.S_K, 'S_Mg': self.S_Mg, 'S_Ca': self.S_Ca, 'S_Na': self.S_Na, 'S_Cl': self.S_Cl, 'S_IC': self.S_IC}

    @property
    def measured(self) -> bool:
        return not self.not_identifiable

def build_ion_basis(*, measured: Mapping[str, float] | None=None, assumed: Mapping[str, float] | None=None, allow_assumed: bool=False) -> IonBasis:
    """Assemble the charge-balance basis, refusing to invent anything.

    ``measured`` supplies values that were actually analysed. ``assumed`` supplies
    stated estimates, and is only accepted when ``allow_assumed=True`` -- which
    stamps every assumed species into ``not_identifiable`` so downstream output
    cannot be represented as measured.

    Raises :class:`IonBasisError` naming the missing species otherwise. That is the
    point: a twin whose pH rests on invented ions produces a confident wrong
    nitrite-lock diagnosis, and failing loudly is cheaper than that.
    """
    required = ('S_K', 'S_Mg', 'S_Ca', 'S_Na', 'S_Cl', 'S_IC')
    measured = dict(measured or {})
    assumed = dict(assumed or {})
    values: dict[str, float] = {}
    unmeasured: list[str] = []
    missing: list[str] = []

    def _check(key: str, raw: object) -> float:
        v = float(raw)
        if v != v or v in (float('inf'), float('-inf')):
            raise IonBasisError(f'{key} is not finite ({raw!r})')
        if v < 0.0:
            raise IonBasisError(f'{key} is negative ({v}). Concentrations cannot be negative, and a negative one silently biases the charge balance that sets pH.')
        return v
    for key in required:
        if key in measured and measured[key] is not None:
            values[key] = _check(key, measured[key])
        elif key in assumed and assumed[key] is not None:
            if not allow_assumed:
                missing.append(key)
            else:
                values[key] = _check(key, assumed[key])
                unmeasured.append(key)
        else:
            missing.append(key)
    if missing:
        raise IonBasisError(f'cannot build a charge-balance basis: no measurement for {sorted(missing)}. The pH solve consumes these, and pH sets free ammonia and free nitrous acid -- so substituting municipal defaults would put the nitrite-lock diagnosis on invented numbers. Measure them (major-ion panel + a closed no-headspace DIC), or pass them via `assumed=` with allow_assumed=True to proceed with an explicit not_identifiable flag.')
    provenance = 'measured' if not unmeasured else 'partly_assumed'
    return IonBasis(provenance=provenance, not_identifiable=tuple(unmeasured), **values)
