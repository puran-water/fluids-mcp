from __future__ import annotations
import math
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Dict, List, Mapping, Optional, Sequence
MW_WATER = 0.018015
MOL_PER_NM3 = 1.0 / 0.0224138
KW_PER_TR = 3.516853
ATM_BARA = 1.01325
DEFAULT_ADP_C = 35.0
SUPERHEAT_WARN_C = 10.0
DEFAULT_DRY_COMPOSITION: Dict[str, float] = {'methane': 0.6, 'carbon dioxide': 0.4}
_CAS = {'methane': '74-82-8', 'carbon dioxide': '124-38-9', 'water': '7732-18-5', 'nitrogen': '7727-37-9', 'oxygen': '7782-44-7', 'hydrogen sulfide': '7783-06-4', 'hydrogen': '1333-74-0'}
_TRC_GAS_COEFFS: Dict[str, tuple] = {'74-82-8': (4.0, 22350000.0, 2018.0, 32.767, -31.098, 1346090000.0, 1229.0, 473.0), '124-38-9': (3.5, 1447000.0, 1029.0, 17.13, -21.542, 479500000.0, 1185.0, 57.0), '7732-18-5': (4.0, 870000.0, 1646.0, 3.111, 1.728, -54010000.0, 559.0, 304.0), '7727-37-9': (3.5, 7615000.0, 3136.0, 2.986, -2.963, 165120000.0, 1784.0, 484.0), '7782-44-7': (3.5, 312000.0, 1442.0, 3.594, -1.895, 38180000.0, 559.0, 267.0), '7783-06-4': (4.0, 1090000.0, 1553.0, 1.845, 4.965, -36440000.0, 381.0, 203.0), '1333-74-0': (4.7, 2311000.0, 1150.0, -10.157, 11.235, -1640000.0, 39.0, 41.0)}
_WATER_HVAP_TB_J_MOL = 40650.0
_WATER_TB_K = 373.15
_WATER_TC_K = 647.096
_FALLBACK_CP_J_MOL_K = 38.7
_FALLBACK_HVAP_KJ_KG = 2320.0
CP_LIQUID_WATER_J_MOL_K = 75.4

class InsufficientProcessData(ValueError):
    """Raised when the wet-gas basis is not determinable from the supplied inputs."""

@lru_cache(maxsize=1)
def _chemicals():
    """The optional ``chemicals`` package, or None -- resolved ONCE per process.

    Without the cache an absent package retried the import on every call (~90 ms each; the
    coupled train model calls this ~15x per simulate). Installing ``chemicals`` into a running
    process is not a case to support, so memoising the None is behaviour-preserving.
    """
    try:
        import chemicals
        return chemicals
    except Exception:
        return None

def water_sat_pressure_Pa(T_C: float) -> float:
    """Saturation vapour pressure of water [Pa].

    IAPWS-95 via ``chemicals.Psat_IAPWS`` when available, else Magnus (within ~0.4%
    over 0-65 C).
    """
    chem = _chemicals()
    if chem is not None:
        try:
            return float(chem.Psat_IAPWS(T_C + 273.15))
        except Exception:
            pass
    return 1000.0 * 0.61094 * math.exp(17.625 * T_C / (T_C + 243.04))

def saturated_water_mole_frac(T_C: float, P_bara: float) -> float:
    """Water mole fraction at saturation for (T, P), capped below unity."""
    return min(0.99, water_sat_pressure_Pa(T_C) / (P_bara * 100000.0))

def _dewpoint_at(partial_pressure_Pa: float) -> float:
    """Invert the saturation curve: temperature [C] at which water condenses."""
    if partial_pressure_Pa <= water_sat_pressure_Pa(-60.0):
        return float('nan')
    lo, hi = (-60.0, 370.0)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if water_sat_pressure_Pa(mid) < partial_pressure_Pa:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-09:
            break
    return 0.5 * (lo + hi)

def water_mole_frac_from_adp(adp_C: float) -> float:
    """Water mole fraction implied by an atmospheric dew point."""
    return water_sat_pressure_Pa(adp_C) / (ATM_BARA * 100000.0)

def adp_from_water_mole_frac(y_water: float) -> float:
    """Atmospheric dew point [C] for a water mole fraction."""
    return _dewpoint_at(y_water * ATM_BARA * 100000.0)

def pdp_from_water_mole_frac(y_water: float, P_bara: float) -> float:
    """Pressure dew point [C] for a water mole fraction at line pressure."""
    return _dewpoint_at(y_water * P_bara * 100000.0)

def resolve_inlet_water(*, adp_C: Optional[float]=None, pdp_C: Optional[float]=None, pdp_pressure_bara: Optional[float]=None, water_mole_frac: Optional[float]=None, allow_default: bool=True) -> float:
    """Resolve the inlet water mole fraction from whichever basis the caller supplied.

    Exactly one of ``adp_C``, (``pdp_C`` + ``pdp_pressure_bara``), or ``water_mole_frac``
    is expected. With none of them, falls back to :data:`DEFAULT_ADP_C` when
    ``allow_default``; otherwise raises so a vendor review cannot proceed on an
    unstated basis.
    """
    supplied = [adp_C is not None, pdp_C is not None, water_mole_frac is not None]
    if sum(supplied) > 1:
        raise InsufficientProcessData('Supply exactly one water basis: adp_C, pdp_C (+pdp_pressure_bara), or water_mole_frac.')
    if water_mole_frac is not None:
        return float(water_mole_frac)
    if pdp_C is not None:
        if pdp_pressure_bara is None:
            raise InsufficientProcessData('pdp_C requires pdp_pressure_bara.')
        return water_sat_pressure_Pa(pdp_C) / (pdp_pressure_bara * 100000.0)
    if adp_C is not None:
        return water_mole_frac_from_adp(adp_C)
    if not allow_default:
        raise InsufficientProcessData('No wet-gas basis supplied (adp_C / pdp_C / water_mole_frac) and defaulting is disabled.')
    return water_mole_frac_from_adp(DEFAULT_ADP_C)

def _trc_coeffs(species: str):
    chem = _chemicals()
    if chem is None:
        return None
    cas = _CAS.get(species.strip().lower())
    if cas is None:
        return None
    coeffs = _TRC_GAS_COEFFS.get(cas)
    return list(coeffs) if coeffs is not None else None

def _molar_enthalpy_drop(species: str, T_hot_C: float, T_cold_C: float) -> Optional[float]:
    """Ideal-gas enthalpy released cooling 1 mol from T_hot to T_cold [J/mol]."""
    coeffs = _trc_coeffs(species)
    if coeffs is None:
        return None
    try:
        from chemicals.heat_capacity import TRCCp_integral
        return float(TRCCp_integral(T_hot_C + 273.15, *coeffs) - TRCCp_integral(T_cold_C + 273.15, *coeffs))
    except Exception:
        return None

def dry_gas_enthalpy_drop_J_mol(composition: Mapping[str, float], T_hot_C: float, T_cold_C: float) -> tuple[float, str]:
    """Mole-weighted ideal-gas enthalpy drop for the dry gas [J/mol], and the basis used."""
    total = sum(composition.values())
    if total <= 0:
        raise InsufficientProcessData('Dry composition must sum to a positive value.')
    acc, missing = (0.0, [])
    for species, frac in composition.items():
        drop = _molar_enthalpy_drop(species, T_hot_C, T_cold_C)
        if drop is None:
            missing.append(species)
            drop = _FALLBACK_CP_J_MOL_K * (T_hot_C - T_cold_C)
        acc += frac / total * drop
    if not missing:
        return (acc, 'chemicals_TRC')
    if _chemicals() is None:
        return (acc, 'fallback_constant_cp')
    return (acc, f'fallback_constant_cp (no TRC data for {sorted(missing)})')

def water_latent_heat_J_mol(T_C: float) -> tuple[float, str]:
    """Latent heat of water at ``T_C`` [J/mol] via Watson, else a constant."""
    chem = _chemicals()
    if chem is not None:
        try:
            return (float(chem.Watson(T=T_C + 273.15, Hvap_ref=_WATER_HVAP_TB_J_MOL, T_ref=_WATER_TB_K, Tc=_WATER_TC_K)), 'chemicals_Watson')
        except Exception:
            pass
    return (_FALLBACK_HVAP_KJ_KG * 1000.0 * MW_WATER, 'fallback_constant_hvap')

@dataclass
class CoolerResult:
    """Heat-and-mass balance across one wet-gas cooler."""
    label: str
    gas_nm3h_wet: float
    pressure_bara: float
    T_in_C: float
    T_out_C: float
    sensible_kW: float
    latent_kW: float
    total_kW: float
    water_condensed_kg_h: float
    water_mole_frac_in: float
    water_mole_frac_out: float
    adp_in_C: float
    adp_out_C: float
    pdp_in_C: float
    pdp_out_C: float
    superheat_in_C: float
    condensing: bool
    property_basis: str
    warnings: List[str] = field(default_factory=list)

    @property
    def total_TR(self) -> float:
        return self.total_kW / KW_PER_TR

    def as_dict(self) -> Dict[str, object]:
        out = {k: v for k, v in self.__dict__.items()}
        out['total_TR'] = self.total_TR
        return out

def cool_wet_gas(gas_nm3h_wet: float, *, T_in_C: float, T_out_C: float, P_bara: float, water_mole_frac_in: Optional[float]=None, adp_C: Optional[float]=None, pdp_C: Optional[float]=None, pdp_pressure_bara: Optional[float]=None, dry_composition: Optional[Mapping[str, float]]=None, label: str='cooler', allow_default_water: bool=True) -> CoolerResult:
    """Sensible + latent duty across one cooler, conserving dry-gas moles.

    ``gas_nm3h_wet`` is the total (wet) normal flow at the cooler inlet. The inlet water
    content is capped at saturation for the inlet (T, P) — gas cannot carry more than
    saturation, and if the caller's basis implies it does, the excess is already liquid
    and is reported as an inlet-superheat warning.
    """
    if gas_nm3h_wet <= 0:
        raise InsufficientProcessData('gas_nm3h_wet must be positive.')
    if T_out_C > T_in_C:
        raise InsufficientProcessData('T_out_C must not exceed T_in_C for a cooler.')
    composition = dict(dry_composition or DEFAULT_DRY_COMPOSITION)
    y_basis = water_mole_frac_in if water_mole_frac_in is not None else resolve_inlet_water(adp_C=adp_C, pdp_C=pdp_C, pdp_pressure_bara=pdp_pressure_bara, allow_default=allow_default_water)
    warnings: List[str] = []
    y_sat_in = saturated_water_mole_frac(T_in_C, P_bara)
    y_in = min(y_basis, y_sat_in)
    if y_basis > y_sat_in:
        warnings.append(f'Inlet water basis ({y_basis:.4f}) exceeds saturation at {T_in_C:.1f} C / {P_bara:.2f} barA ({y_sat_in:.4f}); the stream is two-phase at the inlet and a dry-gas rating does not describe this service.')
    y_out = min(y_in, saturated_water_mole_frac(T_out_C, P_bara))
    condensing = y_out < y_in - 1e-12
    n_wet = gas_nm3h_wet * MOL_PER_NM3 / 3600.0
    n_dry = n_wet * (1.0 - y_in)
    w_in = n_dry * y_in / (1.0 - y_in)
    w_out = n_dry * y_out / (1.0 - y_out)
    cond_mol_s = max(0.0, w_in - w_out)
    dry_drop, dry_basis = dry_gas_enthalpy_drop_J_mol(composition, T_in_C, T_out_C)
    vap_drop = _molar_enthalpy_drop('water', T_in_C, T_out_C)
    vap_basis = 'chemicals_TRC'
    if vap_drop is None:
        vap_drop = 33.6 * (T_in_C - T_out_C)
        vap_basis = 'fallback_constant_cp'
    pdp_in = pdp_from_water_mole_frac(y_in, P_bara)
    T_dew = min(pdp_in, T_in_C)
    sensible_W = n_dry * dry_drop + w_out * vap_drop
    if condensing:
        vap_to_dew = _molar_enthalpy_drop('water', T_in_C, T_dew)
        if vap_to_dew is None:
            vap_to_dew = 33.6 * max(0.0, T_in_C - T_dew)
        sensible_W += cond_mol_s * vap_to_dew
        liquid_span = max(0.0, T_dew - T_out_C) * 0.5
        sensible_W += cond_mol_s * CP_LIQUID_WATER_J_MOL_K * liquid_span
    T_cond_mid = 0.5 * (T_dew + T_out_C) if condensing else T_out_C
    hvap, hvap_basis = water_latent_heat_J_mol(T_cond_mid)
    latent_W = cond_mol_s * hvap
    superheat_in = T_in_C - pdp_in
    if superheat_in < SUPERHEAT_WARN_C:
        warnings.append(f'Inlet superheat is only {superheat_in:.1f} C above the {pdp_in:.1f} C pressure dew point; a modestly wetter feed or cooler discharge puts liquid at the inlet. Confirm the inlet phase condition with the vendor.')
    basis = dry_basis
    if 'fallback' in (dry_basis, vap_basis, hvap_basis) or any(('fallback' in b for b in (dry_basis, vap_basis, hvap_basis))):
        if _chemicals() is None:
            basis = 'screening (constant-property fallback; install `chemicals` for rigorous)'
            warnings.append('Property fallback in use — results are screening-grade. Install `chemicals` for IAPWS/TRC/Watson routes.')
        else:
            named = '; '.join((b for b in (dry_basis, vap_basis, hvap_basis) if 'fallback' in b))
            basis = f'screening ({named})'
            warnings.append(f'Property fallback in use — results are screening-grade: {named}.')
    return CoolerResult(label=label, gas_nm3h_wet=gas_nm3h_wet, pressure_bara=P_bara, T_in_C=T_in_C, T_out_C=T_out_C, sensible_kW=sensible_W / 1000.0, latent_kW=latent_W / 1000.0, total_kW=(sensible_W + latent_W) / 1000.0, water_condensed_kg_h=cond_mol_s * MW_WATER * 3600.0, water_mole_frac_in=y_in, water_mole_frac_out=y_out, adp_in_C=adp_from_water_mole_frac(y_in), adp_out_C=adp_from_water_mole_frac(y_out), pdp_in_C=pdp_in, pdp_out_C=pdp_from_water_mole_frac(y_out, P_bara), superheat_in_C=superheat_in, condensing=condensing, property_basis=basis, warnings=warnings)

@dataclass
class TrainStage:
    """One cooler in a series train.

    ``gas_nm3h_wet`` is that stage's own wet normal flow, which generally differs from the
    upstream stage's because of recycle joining between them. The water *mole fraction*
    is what carries forward, not the molar flow, so the stages need not share a flow.
    """
    label: str
    gas_nm3h_wet: float
    P_bara: float
    T_in_C: float
    T_out_C: float
    knockout: bool = True

def cool_wet_gas_train(stages: Sequence[TrainStage], *, adp_C: Optional[float]=None, water_mole_frac_in: Optional[float]=None, pdp_C: Optional[float]=None, pdp_pressure_bara: Optional[float]=None, dry_composition: Optional[Mapping[str, float]]=None, allow_default_water: bool=True) -> Dict[str, object]:
    """Thread water content through consecutive coolers.

    Condensate removed at an upstream knockout lowers the water mole fraction entering
    every stage downstream, so the aftercoolers of an upgrading train must be solved
    together — rating them independently double-counts the water.

    A stage with ``knockout=False`` still condenses; the liquid simply travels on with the
    gas instead of being drained. The vapour fraction leaving is at saturation either way,
    so the downstream *duty* is unchanged, but the carried liquid is accumulated and
    reported. Re-evaporation of that liquid on recompression is **not** modelled, so a
    train with unremoved condensate is flagged rather than quietly rated.
    """
    y = water_mole_frac_in if water_mole_frac_in is not None else resolve_inlet_water(adp_C=adp_C, pdp_C=pdp_C, pdp_pressure_bara=pdp_pressure_bara, allow_default=allow_default_water)
    results: List[CoolerResult] = []
    extra_warnings: List[str] = []
    liquid_carried_kg_h = 0.0
    for stage in stages:
        res = cool_wet_gas(stage.gas_nm3h_wet, T_in_C=stage.T_in_C, T_out_C=stage.T_out_C, P_bara=stage.P_bara, water_mole_frac_in=y, dry_composition=dry_composition, label=stage.label)
        results.append(res)
        y = res.water_mole_frac_out
        if not stage.knockout and res.water_condensed_kg_h > 0:
            liquid_carried_kg_h += res.water_condensed_kg_h
            extra_warnings.append(f'{stage.label}: {res.water_condensed_kg_h:.1f} kg/h of condensate is not removed and travels on with the gas. Re-evaporation on recompression is not modelled, so downstream duties are understated and liquid carryover into the next machine is a mechanical risk in its own right.')
    total_kW = sum((r.total_kW for r in results))
    return {'stages': [r.as_dict() for r in results], 'total_kW': total_kW, 'total_TR': total_kW / KW_PER_TR, 'total_condensate_kg_h': sum((r.water_condensed_kg_h for r in results)), 'liquid_carried_forward_kg_h': liquid_carried_kg_h, 'inlet_adp_C': results[0].adp_in_C if results else None, 'outlet_adp_C': results[-1].adp_out_C if results else None, 'warnings': [w for r in results for w in r.warnings] + extra_warnings}
