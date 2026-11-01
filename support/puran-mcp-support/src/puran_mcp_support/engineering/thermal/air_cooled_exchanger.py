from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional
FACE_VELOCITY_ADVISORY_M_S = (1.5, 4.0)
FIN_DENSITY_RANGE_PER_M = (275.0, 435.0)

class InsufficientVendorData(ValueError):
    """Raised when a rating cannot be produced from the fields the vendor supplied."""

    def __init__(self, missing: List[str], context: str=''):
        self.missing = list(missing)
        detail = ', '.join(self.missing)
        msg = f'Vendor data insufficient to rate this exchanger; missing: {detail}'
        if context:
            msg = f'{msg} ({context})'
        super().__init__(msg)

@dataclass
class BundleGeometry:
    """Finned-tube bundle geometry, as it appears on a vendor datasheet."""
    tube_rows: int
    tube_passes: int
    tubes_per_row: int
    tube_length_m: float
    tube_od_m: float
    fin_thickness_m: float
    fin_density_per_m: float
    fin_height_m: float
    pitch_normal_m: Optional[float] = None
    pitch_parallel_m: Optional[float] = None
    tube_thickness_m: Optional[float] = None
    parallel_bays: int = 1
    bundles_per_bay: int = 1
    fans_per_bay: int = 1
    fan_diameter_m: Optional[float] = None
    fin_conductivity_W_mK: float = 205.0

    def required_fields(self) -> List[str]:
        missing = []
        for name in ('tube_rows', 'tube_passes', 'tubes_per_row', 'tube_length_m', 'tube_od_m', 'fin_thickness_m', 'fin_density_per_m', 'fin_height_m'):
            if getattr(self, name) in (None, 0):
                missing.append(name)
        if self.pitch_normal_m is None or self.pitch_parallel_m is None:
            missing.append('pitch_normal_m/pitch_parallel_m')
        return missing

    def build(self):
        """Construct the ``fluids`` geometry object, raising if fields are missing."""
        missing = self.required_fields()
        if missing:
            raise InsufficientVendorData(missing, 'bundle geometry')
        try:
            from fluids.geometry import AirCooledExchanger
        except Exception as exc:
            raise RuntimeError('fluids is required for ACHE rating; install the `thermal` extra.') from exc
        return AirCooledExchanger(tube_rows=self.tube_rows, tube_passes=self.tube_passes, tubes_per_row=self.tubes_per_row, tube_length=self.tube_length_m, tube_diameter=self.tube_od_m, fin_thickness=self.fin_thickness_m, fin_density=self.fin_density_per_m, fin_height=self.fin_height_m, pitch_normal=self.pitch_normal_m, pitch_parallel=self.pitch_parallel_m, tube_thickness=self.tube_thickness_m, parallel_bays=self.parallel_bays, bundles_per_bay=self.bundles_per_bay, fans_per_bay=self.fans_per_bay, fan_diameter=self.fan_diameter_m)

@dataclass
class AirState:
    """Air-side inlet condition."""
    T_in_C: float
    density_kg_m3: float
    cp_J_kgK: float
    viscosity_Pa_s: float
    conductivity_W_mK: float
    mass_flow_kg_s: Optional[float] = None
    face_velocity_m_s: Optional[float] = None
_FT_EXACT = {(1, 1), (1, 2), (1, 3), (1, 4), (2, 2), (3, 3), (4, 4), (2, 4)}

def ft_support_tier(tube_passes: int, tube_rows: int) -> str:
    """Classify how ``ht.Ft_aircooler`` will treat this (passes, rows) combination.

    ``exact`` — a dedicated coefficient set exists.
    ``substituted`` — the library substitutes a neighbouring set and calls it a reasonable
    assumption.
    ``unsupported`` — the library falls through to a default its own source describes as
    "a bad assumption"; the resulting Ft is an approximation and must be labelled one.
    """
    if (tube_passes, tube_rows) in _FT_EXACT:
        return 'exact'
    if tube_passes == 1 and tube_rows > 4:
        return 'substituted'
    if tube_passes > 4 and tube_rows > 4 and (tube_passes == tube_rows):
        return 'substituted'
    return 'unsupported'

def recover_internal_fin_efficiency(h_bare_basis: float, geom, fin_conductivity_W_mK: float) -> Optional[float]:
    """Recover the fin efficiency ``ht`` used internally, by inverting its area weighting.

    ``ht`` forms ``h_bare = (eta(h_raw)*A_fin + A_tube_showing)/A * h_raw * A_increase``
    where ``h_raw`` is the raw local coefficient and ``eta`` depends on it. Inverting that
    relation for ``h_raw`` recovers the efficiency the library actually applied — which is
    the number to report, and is not the same as calling the fin-efficiency function on the
    returned coefficient.
    """
    try:
        from ht.core import fin_efficiency_Kern_Kraus
    except Exception:
        return None

    def implied(h_raw: float) -> float:
        eta = fin_efficiency_Kern_Kraus(geom.tube_diameter, geom.fin_diameter, geom.fin_thickness, fin_conductivity_W_mK, h_raw)
        return (eta * geom.A_fin + geom.A_tube_showing) / geom.A * h_raw * geom.A_increase
    lo, hi = (0.001, max(10000.0, h_bare_basis * 10.0))
    if implied(hi) < h_bare_basis:
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if implied(mid) < h_bare_basis:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-09:
            break
    h_raw = 0.5 * (lo + hi)
    return float(fin_efficiency_Kern_Kraus(geom.tube_diameter, geom.fin_diameter, geom.fin_thickness, fin_conductivity_W_mK, h_raw))

@dataclass
class AirSideResult:
    """Air-side rating of a finned bundle."""
    h_air_effective_bare_basis: float
    h_air_effective_extended_basis: float
    h_air_briggs_young_bare_basis: float
    correlation_spread_pct: float
    fin_efficiency_internal: Optional[float]
    area_extended_m2: float
    area_bare_tube_m2: float
    area_increase_ratio: float
    area_face_m2: float
    area_min_free_flow_m2: float
    air_mass_flow_kg_s: float
    air_volumetric_flow_m3_s: float
    face_velocity_m_s: float
    bundle_dP_Pa: float
    warnings: List[str] = field(default_factory=list)

    def as_dict(self) -> Dict[str, object]:
        return dict(self.__dict__)

def rate_air_side(geometry: BundleGeometry, air: AirState, *, correlation: str='ESDU_high_fin') -> AirSideResult:
    """Air-side coefficient, area breakdown and bundle pressure drop.

    Both the ESDU and Briggs-Young correlations are evaluated every time and their spread
    reported: two independent correlations agreeing is the cheapest available confidence
    signal, and a wide spread means the geometry is outside the range where either is
    trustworthy.
    """
    try:
        from ht.air_cooler import dP_ESDU_high_fin, h_Briggs_Young, h_ESDU_high_fin
    except Exception as exc:
        raise RuntimeError('ht is required for ACHE rating; install the `thermal` extra.') from exc
    geom = geometry.build()
    warnings: List[str] = []
    if air.mass_flow_kg_s is not None:
        m_air = air.mass_flow_kg_s
    elif air.face_velocity_m_s is not None:
        m_air = air.density_kg_m3 * geom.A_face * air.face_velocity_m_s
    else:
        raise InsufficientVendorData(['air mass_flow_kg_s or face_velocity_m_s'], 'air-side flow')
    shared = dict(m=m_air, A=geom.A, A_min=geom.A_min, A_increase=geom.A_increase, A_fin=geom.A_fin, A_tube_showing=geom.A_tube_showing, tube_diameter=geom.tube_diameter, fin_diameter=geom.fin_diameter, fin_thickness=geom.fin_thickness, bare_length=geom.bare_length, rho=air.density_kg_m3, Cp=air.cp_J_kgK, mu=air.viscosity_Pa_s, k=air.conductivity_W_mK, k_fin=geometry.fin_conductivity_W_mK)
    h_esdu = float(h_ESDU_high_fin(pitch_parallel=geom.pitch_parallel, pitch_normal=geom.pitch_normal, tube_rows=geom.tube_rows, **shared))
    h_by = float(h_Briggs_Young(**shared))
    h_primary = h_esdu if correlation == 'ESDU_high_fin' else h_by
    spread = abs(h_by - h_esdu) / h_esdu * 100.0 if h_esdu else float('nan')
    if spread > 15.0:
        warnings.append(f'ESDU and Briggs-Young differ by {spread:.1f}% — this geometry is near or outside the range where either correlation is reliable.')
    dP = float(dP_ESDU_high_fin(m=m_air, A_min=geom.A_min, A_increase=geom.A_increase, flow_area_contraction_ratio=geom.flow_area_contraction_ratio, tube_diameter=geom.tube_diameter, pitch_parallel=geom.pitch_parallel, pitch_normal=geom.pitch_normal, tube_rows=geom.tube_rows, rho=air.density_kg_m3, mu=air.viscosity_Pa_s))
    q_air = m_air / air.density_kg_m3
    v_face = q_air / geom.A_face
    lo, hi = FACE_VELOCITY_ADVISORY_M_S
    if not lo <= v_face <= hi:
        warnings.append(f'Face velocity {v_face:.2f} m/s is outside the {lo}-{hi} m/s band typical of air-cooled bundles (engineering advisory, not a code requirement).')
    if not FIN_DENSITY_RANGE_PER_M[0] <= geometry.fin_density_per_m <= FIN_DENSITY_RANGE_PER_M[1]:
        warnings.append(f'Fin density {geometry.fin_density_per_m:.0f} /m is outside the {FIN_DENSITY_RANGE_PER_M[0]:.0f}-{FIN_DENSITY_RANGE_PER_M[1]:.0f} /m range the high-fin correlations were fitted over.')
    return AirSideResult(h_air_effective_bare_basis=h_primary, h_air_effective_extended_basis=h_primary / geom.A_increase, h_air_briggs_young_bare_basis=h_by, correlation_spread_pct=spread, fin_efficiency_internal=recover_internal_fin_efficiency(h_primary, geom, geometry.fin_conductivity_W_mK), area_extended_m2=geom.A, area_bare_tube_m2=geom.A_bare_tube, area_increase_ratio=geom.A_increase, area_face_m2=geom.A_face, area_min_free_flow_m2=geom.A_min, air_mass_flow_kg_s=m_air, air_volumetric_flow_m3_s=q_air, face_velocity_m_s=v_face, bundle_dP_Pa=dP, warnings=warnings)

def overall_U_extended_basis(*, h_air_extended_basis: float, h_tube_W_m2K: float, area_extended_m2: float, area_tube_inside_m2: float, fouling_air_m2K_W: float=0.0, fouling_tube_m2K_W: float=0.0, wall_resistance_m2K_W: float=0.0) -> float:
    """Overall coefficient referred to the **extended (total outside) area**.

    Tube-side and fouling resistances are scaled by the area ratio, which is the step the
    generic overall-U helper cannot do because it has no concept of extended surface.
    """
    if h_air_extended_basis <= 0 or h_tube_W_m2K <= 0:
        raise InsufficientVendorData(['h_air_extended_basis', 'h_tube_W_m2K'], 'overall U')
    ratio = area_extended_m2 / area_tube_inside_m2
    inv = 1.0 / h_air_extended_basis + fouling_air_m2K_W + wall_resistance_m2K_W + ratio * (fouling_tube_m2K_W + 1.0 / h_tube_W_m2K)
    return 1.0 / inv

@dataclass
class FanPower:
    """Fan power broken out so a vendor's motor rating can be checked against physics."""
    bundle_dP_Pa: float
    external_dP_Pa: float
    total_dP_Pa: float
    air_power_kW: float
    shaft_power_kW: float
    absorbed_electrical_kW: float
    vendor_motor_kW: Optional[float]
    reconcilable: Optional[bool]
    warnings: List[str] = field(default_factory=list)

    def as_dict(self) -> Dict[str, object]:
        return dict(self.__dict__)

def fan_power(*, air_volumetric_flow_m3_s: float, bundle_dP_Pa: float, fan_static_efficiency: Optional[float]=None, drive_efficiency: Optional[float]=None, motor_efficiency: Optional[float]=None, external_dP_Pa: float=0.0, vendor_motor_kW: Optional[float]=None) -> FanPower:
    """Air, shaft and absorbed electrical power, reported separately.

    These are three different numbers and conflating them is how a vendor's motor rating
    goes unchallenged. Efficiencies are required rather than defaulted: a plausible-looking
    shaft power built on an assumed 65% fan efficiency is exactly the kind of figure that
    cannot be defended when the vendor pushes back.
    """
    missing = [name for name, val in (('fan_static_efficiency', fan_static_efficiency), ('drive_efficiency', drive_efficiency), ('motor_efficiency', motor_efficiency)) if val is None]
    if missing:
        raise InsufficientVendorData(missing, 'fan power')
    total_dP = bundle_dP_Pa + external_dP_Pa
    air_kW = total_dP * air_volumetric_flow_m3_s / 1000.0
    shaft_kW = air_kW / fan_static_efficiency
    elec_kW = shaft_kW / drive_efficiency / motor_efficiency
    warnings: List[str] = []
    reconcilable: Optional[bool] = None
    if vendor_motor_kW is not None:
        reconcilable = elec_kW <= vendor_motor_kW * 1.05
        if vendor_motor_kW > 2.0 * elec_kW:
            warnings.append(f'Vendor motor {vendor_motor_kW:.2f} kW is more than twice the {elec_kW:.2f} kW absorbed by the stated bundle and external losses. That is not necessarily wrong — it may carry system resistance or nameplate margin — but it cannot be regenerated from the fields supplied. Ask the vendor for the fan curve, static-pressure breakdown and efficiencies.')
        elif not reconcilable:
            warnings.append(f'Absorbed electrical power {elec_kW:.2f} kW exceeds the vendor motor rating {vendor_motor_kW:.2f} kW.')
    if external_dP_Pa == 0.0:
        warnings.append('External static losses (louvres, screens, plenum, fan ring, velocity pressure) were taken as zero; shaft power is bundle-only and will understate the installed requirement.')
    return FanPower(bundle_dP_Pa=bundle_dP_Pa, external_dP_Pa=external_dP_Pa, total_dP_Pa=total_dP, air_power_kW=air_kW, shaft_power_kW=shaft_kW, absorbed_electrical_kW=elec_kW, vendor_motor_kW=vendor_motor_kW, reconcilable=reconcilable, warnings=warnings)

def lmtd_counterflow(*, T_hot_in_C: float, T_hot_out_C: float, T_cold_in_C: float, T_cold_out_C: float) -> float:
    """Counterflow log-mean temperature difference [K] -- the house one
    (``thermal.exchanger_rating.lmtd_counterflow``, Wave 2).

    Raises when either terminal difference is non-positive, which means the stated
    temperatures cross and no exchanger can deliver them.
    """
    from puran_mcp_support.engineering.thermal.exchanger_rating import TemperatureCross
    from puran_mcp_support.engineering.thermal.exchanger_rating import lmtd_counterflow as _lmtd
    try:
        return _lmtd(T_hot_in_C, T_hot_out_C, T_cold_in_C, T_cold_out_C)
    except TemperatureCross:
        d1, d2 = (T_hot_in_C - T_cold_out_C, T_hot_out_C - T_cold_in_C)
        raise InsufficientVendorData(['terminal temperature difference'], f'Temperatures cross: hot {T_hot_in_C}->{T_hot_out_C} against cold {T_cold_in_C}->{T_cold_out_C} gives approaches of {d1:.1f} and {d2:.1f} K') from None

@dataclass
class AchievableDuty:
    """What the bundle delivers, against what the service needs."""
    U_extended_basis_W_m2K: float
    area_extended_m2: float
    LMTD_K: float
    Ft: float
    MTD_K: float
    Q_achievable_kW: float
    Q_required_kW: Optional[float]
    margin_fraction: Optional[float]
    verdict: Optional[str]
    resistance_share_air_pct: float
    resistance_share_tube_pct: float
    warnings: List[str] = field(default_factory=list)

    def as_dict(self) -> Dict[str, object]:
        return dict(self.__dict__)
MARGIN_ADEQUATE = 0.2
MARGIN_MARGINAL = 0.1

def achievable_duty(*, U_extended_basis_W_m2K: float, area_extended_m2: float, T_hot_in_C: float, T_hot_out_C: float, T_air_in_C: float, T_air_out_C: float, tube_passes: int, tube_rows: int, h_air_extended_basis: float, Q_required_kW: Optional[float]=None) -> AchievableDuty:
    """Duty the bundle delivers as ``U · A · Ft · LMTD``, and the margin over the service.

    This is the question a rating exists to answer. Reporting U, area and Ft separately and
    stopping short leaves the reader to multiply them, which is exactly the step where a
    vendor review goes wrong.
    """
    lmtd = lmtd_counterflow(T_hot_in_C=T_hot_in_C, T_hot_out_C=T_hot_out_C, T_cold_in_C=T_air_in_C, T_cold_out_C=T_air_out_C)
    ft = mtd_correction(T_hot_in_C=T_hot_in_C, T_hot_out_C=T_hot_out_C, T_air_in_C=T_air_in_C, T_air_out_C=T_air_out_C, tube_passes=tube_passes, tube_rows=tube_rows)
    mtd = lmtd * ft['Ft']
    q_kW = U_extended_basis_W_m2K * area_extended_m2 * mtd / 1000.0
    warnings = list(ft['warnings'])
    margin = verdict = None
    if Q_required_kW:
        margin = q_kW / Q_required_kW - 1.0
        if margin >= MARGIN_ADEQUATE:
            verdict = 'adequate'
        elif margin >= MARGIN_MARGINAL:
            verdict = 'marginal'
        elif margin >= 0.0:
            verdict = 'tight'
        else:
            verdict = 'short'
        if margin < MARGIN_MARGINAL:
            warnings.append(f'Bundle delivers {q_kW:.1f} kW against a required {Q_required_kW:.1f} kW ({margin:+.1%}). Any fouling allowance or hot-air recirculation comes out of this margin.')
    air_share = 1.0 / h_air_extended_basis * U_extended_basis_W_m2K
    if air_share > 1.0:
        raise InsufficientVendorData(['consistent U and air-side h'], 'Air-side resistance alone exceeds the total; check the area basis of both.')
    if 1.0 - air_share > 0.7:
        warnings.append(f'Tube side carries {(1.0 - air_share) * 100:.0f}% of the total resistance, so added fin area or fan power buys very little. The design lever here is tube-side velocity and pass count.')
    return AchievableDuty(U_extended_basis_W_m2K=U_extended_basis_W_m2K, area_extended_m2=area_extended_m2, LMTD_K=lmtd, Ft=ft['Ft'], MTD_K=mtd, Q_achievable_kW=q_kW, Q_required_kW=Q_required_kW, margin_fraction=margin, verdict=verdict, resistance_share_air_pct=air_share * 100.0, resistance_share_tube_pct=(1.0 - air_share) * 100.0, warnings=warnings)

def mtd_correction(*, T_hot_in_C: float, T_hot_out_C: float, T_air_in_C: float, T_air_out_C: float, tube_passes: int, tube_rows: int) -> Dict[str, object]:
    """Crossflow MTD correction factor, with the row/pass support tier attached."""
    try:
        from ht.air_cooler import Ft_aircooler
    except Exception as exc:
        raise RuntimeError('ht is required for ACHE rating.') from exc
    tier = ft_support_tier(tube_passes, tube_rows)
    ft = float(Ft_aircooler(Thi=T_hot_in_C, Tho=T_hot_out_C, Tci=T_air_in_C, Tco=T_air_out_C, Ntp=tube_passes, rows=tube_rows))
    warnings: List[str] = []
    if tier == 'unsupported':
        warnings.append(f'No Ft coefficient set exists for {tube_passes} passes / {tube_rows} rows; the correlation falls back to a 4-row/2-pass set that its own source calls a bad assumption. Treat this Ft as an approximation, not a rated value.')
    elif tier == 'substituted':
        warnings.append(f'Ft for {tube_passes} passes / {tube_rows} rows uses a substituted coefficient set (documented as a reasonable assumption).')
    return {'Ft': ft, 'support_tier': tier, 'warnings': warnings}
