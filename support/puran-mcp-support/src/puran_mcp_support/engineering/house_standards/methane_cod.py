from __future__ import annotations
import math
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Mapping
from puran_mcp_support.engineering.thermal.gas_constants import NORMAL_MOLAR_VOLUME_M3_KMOL, NORMAL_PRESSURE_KPA, NORMAL_TEMPERATURE_K
COD_BALANCE_RULING_ID = 'PUBLIC-RULE-b60c4e378be9'
O2_KG_PER_KMOL = 2 * 16.0
COD_KG_PER_KMOL_ELECTRON = O2_KG_PER_KMOL / 4

def electron_equivalents_per_mol(c: float, h: float, o: float=0.0) -> float:
    """Electrons one mol of CcHhOo gives up when oxidised to CO2 and H2O: 4c + h - 2o (carbon to +4, hydrogen to +1,
    oxygen held at -2)."""
    return 4.0 * c + h - 2.0 * o
CH4_COD_KG_PER_KMOL = electron_equivalents_per_mol(c=1, h=4) * COD_KG_PER_KMOL_ELECTRON
METHANE_NM3_PER_KG_COD = NORMAL_MOLAR_VOLUME_M3_KMOL / CH4_COD_KG_PER_KMOL
COD_METHANE_CLOSURE_REL_TOLERANCE = 0.002

@dataclass(frozen=True)
class CodMethaneBalance:
    """One digester answer's COD balance, every term named (kg COD/d; the methane also as Nm3/d)."""
    engine: str
    cod_in_kg_d: float
    cod_out_kg_d: float
    methane_nm3_d: float
    methane_cod_kg_d: float
    other_terms_kg_d: Mapping[str, float]
    residual_kg_d: float
    relative_residual: float
    tolerance: float
    closes: bool
    engine_methane_cod_kg_d: float | None = None
    notes: tuple[str, ...] = field(default_factory=tuple)

    def as_dict(self) -> dict[str, Any]:
        return {'engine': self.engine, 'cod_in_kg_d': self.cod_in_kg_d, 'cod_out_kg_d': self.cod_out_kg_d, 'methane_nm3_d': self.methane_nm3_d, 'methane_cod_kg_d': self.methane_cod_kg_d, 'other_terms_kg_d': dict(self.other_terms_kg_d), 'residual_kg_d': self.residual_kg_d, 'relative_residual': self.relative_residual, 'tolerance': self.tolerance, 'closes': self.closes, 'engine_methane_cod_kg_d': self.engine_methane_cod_kg_d, 'notes': list(self.notes), 'methane_nm3_per_kg_cod': METHANE_NM3_PER_KG_COD, 'ruling': COD_BALANCE_RULING_ID}

class CodMethaneBalanceRefused(ValueError):
    """A digester answer whose COD does not close, or whose methane volume is not its methane COD at the equivalence."""

    def __init__(self, code: str, message: str, closure: CodMethaneBalance):
        super().__init__(f'{code}: {message}')
        self.code = code
        self.closure = closure

def _finite_nonnegative(name: str, value: float) -> float:
    v = float(value)
    if not math.isfinite(v) or v < 0:
        raise ValueError(f'{name} must be finite and non-negative, not {value!r}')
    return v

def cod_methane_balance(engine: str, *, cod_in_kg_d: float, cod_out_kg_d: float, methane_nm3_d: float, other_terms_kg_d: Mapping[str, float] | None=None, engine_methane_cod_kg_d: float | None=None, tolerance: float=COD_METHANE_CLOSURE_REL_TOLERANCE, notes: tuple[str, ...]=()) -> CodMethaneBalance:
    """The balance, computed and returned whether or not it closes (``require_cod_methane_balance`` refuses).

    The methane's COD is its normal volume over ``METHANE_NM3_PER_KG_COD``; ``engine_methane_cod_kg_d`` (the engine's
    own methane COD, where it carries one) must be the same quantity within the tolerance of the COD in, else the answer
    is refused as ``methane_volume_off_cod_equivalence``."""
    cod_in = _finite_nonnegative('cod_in_kg_d', cod_in_kg_d)
    cod_out = _finite_nonnegative('cod_out_kg_d', cod_out_kg_d)
    nm3 = _finite_nonnegative('methane_nm3_d', methane_nm3_d)
    others = {str(k): float(v) for k, v in (other_terms_kg_d or {}).items()}
    for name, value in others.items():
        if not math.isfinite(value):
            raise ValueError(f'other COD term {name!r} must be finite, not {value!r}')
    if not (math.isfinite(tolerance) and tolerance > 0):
        raise ValueError(f'the tolerance must be positive and finite, not {tolerance!r}')
    methane_cod = nm3 / METHANE_NM3_PER_KG_COD
    residual = cod_in - cod_out - methane_cod - math.fsum(others.values())
    relative = residual / cod_in if cod_in > 0 else 0.0 if residual == 0 else math.inf
    own = None if engine_methane_cod_kg_d is None else _finite_nonnegative('engine_methane_cod_kg_d', engine_methane_cod_kg_d)
    closes = abs(relative) <= tolerance and (own is None or abs(own - methane_cod) <= tolerance * max(cod_in, 1e-300))
    return CodMethaneBalance(engine=engine, cod_in_kg_d=cod_in, cod_out_kg_d=cod_out, methane_nm3_d=nm3, methane_cod_kg_d=methane_cod, other_terms_kg_d=MappingProxyType(others), residual_kg_d=residual, relative_residual=relative, tolerance=float(tolerance), closes=closes, engine_methane_cod_kg_d=own, notes=tuple(notes))

def require_cod_methane_balance(engine: str, **kwargs: Any) -> CodMethaneBalance:
    """``cod_methane_balance``, refused by name where it does not close (``CodMethaneBalanceRefused``)."""
    closure = cod_methane_balance(engine, **kwargs)
    if closure.closes:
        return closure
    if closure.engine_methane_cod_kg_d is not None and abs(closure.engine_methane_cod_kg_d - closure.methane_cod_kg_d) > closure.tolerance * max(closure.cod_in_kg_d, 1e-300):
        raise CodMethaneBalanceRefused('methane_volume_off_cod_equivalence', f'{engine}: methane {closure.methane_nm3_d:.6g} Nm3/d is {closure.methane_cod_kg_d:.6g} kg COD/d at {METHANE_NM3_PER_KG_COD:.6g} Nm3/kg COD, but the engine carries {closure.engine_methane_cod_kg_d:.6g} kg COD/d of methane: one quantity, two values ({COD_BALANCE_RULING_ID})', closure)
    others = ', '.join((f'{k} {v:.6g}' for k, v in closure.other_terms_kg_d.items())) or 'none'
    raise CodMethaneBalanceRefused('cod_methane_balance_not_closed', f'{engine}: COD does not close: in {closure.cod_in_kg_d:.6g} kg/d = out (liquid + solids) {closure.cod_out_kg_d:.6g} + methane {closure.methane_cod_kg_d:.6g} ({closure.methane_nm3_d:.6g} Nm3/d at {METHANE_NM3_PER_KG_COD:.6g} Nm3/kg COD) + other routes ({others}) leaves {closure.residual_kg_d:.6g} kg/d ({closure.relative_residual:.3%} of the COD in; tolerance {closure.tolerance:.2%}) ({COD_BALANCE_RULING_ID})', closure)
