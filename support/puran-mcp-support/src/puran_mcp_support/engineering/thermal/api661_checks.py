from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional

class FindingType(str, Enum):
    STANDARD_NONCONFORMANCE = 'standard_nonconformance'
    PURCHASER_DEVIATION = 'purchaser_deviation'
    CORRELATION_RANGE_WARNING = 'correlation_range_warning'
    ENGINEERING_ADVISORY = 'engineering_advisory'

@dataclass(frozen=True)
class Rule:
    """A single check, with the provenance needed to defend it."""
    key: str
    clause: Optional[str]
    standard: Optional[str]
    summary: str
    source: str
    verified: bool = True
API661_STANDARD = 'API 661 (Air-cooled Heat Exchangers for General Refinery Service)'
API661_EDITION_NOTE = "API's current catalog lists Standard 661 at edition 8. Clause numbering below is taken from the accessible 2018 text; confirm against the contractually invoked edition before quoting a clause number to a vendor."
API661_RULES: Dict[str, Rule] = {'fan_face_area_ratio': Rule(key='fan_face_area_ratio', clause='§7.2.3.1', standard=API661_STANDARD, summary='Each fan shall cover at least 40% of the bundle face area it serves.', source='API 661 2018 text, §7.2.3.1'), 'fans_per_bay': Rule(key='fans_per_bay', clause='§7.2.3.3', standard=API661_STANDARD, summary='Two or more fans per bay unless the purchaser agrees to a single fan.', source='API 661 2018 text, §7.2.3.3'), 'face_velocity_band': Rule(key='face_velocity_band', clause=None, standard=None, summary='Bundle face velocity typically 1.5-4 m/s.', source="Engineering practice. API 661's 3.6/4 m/s figures concern approach velocity around forced/induced-draft units, NOT bundle face velocity — do not cite them for this check.", verified=False), 'fouling_allowance': Rule(key='fouling_allowance', clause=None, standard=None, summary='A non-zero air-side fouling allowance should be specified.', source='Purchaser design basis; API 661 leaves fouling to the purchaser.', verified=False)}
MIN_FAN_FACE_AREA_RATIO = 0.4
MIN_FANS_PER_BAY = 2

@dataclass
class Finding:
    """One check result, ready to render into a vendor-callout list."""
    key: str
    finding_type: FindingType
    passed: bool
    message: str
    clause: Optional[str] = None
    standard: Optional[str] = None
    source: Optional[str] = None
    observed: Optional[float] = None
    threshold: Optional[float] = None
    remedy: Optional[str] = None

    def as_dict(self) -> Dict[str, object]:
        d = dict(self.__dict__)
        d['finding_type'] = self.finding_type.value
        return d

def _finding(rule: Rule, *, passed: bool, message: str, base_type: FindingType, observed=None, threshold=None, remedy=None, api661_invoked: bool=True) -> Finding:
    ftype = base_type
    if base_type is FindingType.STANDARD_NONCONFORMANCE:
        if not api661_invoked:
            ftype = FindingType.PURCHASER_DEVIATION
        elif not rule.verified:
            ftype = FindingType.ENGINEERING_ADVISORY
    return Finding(key=rule.key, finding_type=ftype, passed=passed, message=message, clause=rule.clause, standard=rule.standard, source=rule.source, observed=observed, threshold=threshold, remedy=remedy)

def check_fan_coverage(*, fan_diameter_m: float, bundle_face_area_m2: float, fans_per_bay: int=1, parallel_bays: int=1, face_area_is_total: bool=True, api661_invoked: bool=True) -> Finding:
    """Fan area against the face area **each fan serves**, not the whole exchanger.

    The standard is written per fan, so the served area is the total face divided by the
    total number of fans. Dividing only by ``fans_per_bay`` while passing the whole
    exchanger's face area makes an N-bay unit look N times worse than it is — the bays are
    replicas, each with their own fans, so the per-fan coverage is unchanged by adding
    bays. Set ``face_area_is_total=False`` when the caller is passing a single bay's face.
    """
    import math
    if fan_diameter_m <= 0 or bundle_face_area_m2 <= 0 or fans_per_bay < 1 or (parallel_bays < 1):
        raise ValueError('fan_diameter_m, bundle_face_area_m2, fans_per_bay and parallel_bays must be positive')
    fan_area = math.pi * (fan_diameter_m / 2.0) ** 2
    total_fans = fans_per_bay * (parallel_bays if face_area_is_total else 1)
    served = bundle_face_area_m2 / total_fans
    ratio = fan_area / served
    passed = ratio >= MIN_FAN_FACE_AREA_RATIO
    needed_d = 2.0 * math.sqrt(MIN_FAN_FACE_AREA_RATIO * served / math.pi)
    rule = API661_RULES['fan_face_area_ratio']
    if passed:
        msg = f'Fan covers {ratio:.1%} of the {served:.3f} m² face area it serves, meeting the {MIN_FAN_FACE_AREA_RATIO:.0%} minimum.'
        remedy = None
    else:
        msg = f'Fan covers only {ratio:.1%} of the {served:.3f} m² face area it serves, against a {MIN_FAN_FACE_AREA_RATIO:.0%} minimum. At this coverage air short-circuits beneath the fan and much of the bundle sees little flow, so the rated air-side performance will not be realised across the face.'
        remedy = f'Increase fan diameter to at least {needed_d:.3f} m, or add fans, or reduce the face area served per fan.'
    return _finding(rule, passed=passed, message=msg, base_type=FindingType.STANDARD_NONCONFORMANCE, observed=ratio, threshold=MIN_FAN_FACE_AREA_RATIO, remedy=remedy, api661_invoked=api661_invoked)

def check_fans_per_bay(*, fans_per_bay: int, purchaser_agreed_single: bool=False, api661_invoked: bool=True) -> Finding:
    """Two fans per bay unless the purchaser has agreed otherwise."""
    rule = API661_RULES['fans_per_bay']
    passed = fans_per_bay >= MIN_FANS_PER_BAY or purchaser_agreed_single
    if passed and fans_per_bay < MIN_FANS_PER_BAY:
        msg = 'Single fan per bay, recorded as agreed by the purchaser.'
        remedy = None
    elif passed:
        msg = f'{fans_per_bay} fans per bay meets the two-fan minimum.'
        remedy = None
    else:
        msg = f'{fans_per_bay} fan per bay against a two-fan minimum. A single fan gives no airflow redundancy — a motor or belt failure takes the whole bay out rather than derating it.'
        remedy = "Quote two fans per bay, or record the purchaser's written agreement to a single fan along with the accepted loss of redundancy."
    return _finding(rule, passed=passed, message=msg, base_type=FindingType.STANDARD_NONCONFORMANCE, observed=float(fans_per_bay), threshold=float(MIN_FANS_PER_BAY), remedy=remedy, api661_invoked=api661_invoked)

def check_face_velocity(*, face_velocity_m_s: float, band=(1.5, 4.0)) -> Finding:
    """Face velocity against typical practice — advisory, deliberately not cited to API."""
    rule = API661_RULES['face_velocity_band']
    lo, hi = band
    passed = lo <= face_velocity_m_s <= hi
    msg = f'Face velocity {face_velocity_m_s:.2f} m/s is within the {lo}-{hi} m/s band typical of air-cooled bundles.' if passed else f'Face velocity {face_velocity_m_s:.2f} m/s falls outside the {lo}-{hi} m/s band typical of air-cooled bundles. This is engineering practice, not an API requirement.'
    return _finding(rule, passed=passed, message=msg, base_type=FindingType.ENGINEERING_ADVISORY, observed=face_velocity_m_s, threshold=lo)

def check_fouling_allowance(*, air_side_fouling_m2K_W: Optional[float], tube_side_fouling_m2K_W: Optional[float]=None, minimum_air_side: float=0.000176) -> Finding:
    """Fouling allowance against the purchaser's basis.

    ``minimum_air_side`` defaults to a commonly specified 0.000176 m²K/W (0.001 h·ft²·°F/Btu)
    for a clean-air site. A dusty site — a press-mud-fed CBG plant, for instance — warrants
    more, and the value belongs in the project design basis rather than here.
    """
    rule = API661_RULES['fouling_allowance']
    observed = air_side_fouling_m2K_W
    passed = observed is not None and observed > 0
    if observed is None:
        msg = 'No air-side fouling allowance stated on the datasheet.'
    elif observed <= 0:
        msg = 'Air-side fouling allowance is zero, so every margin quoted on this datasheet is against a perfectly clean bundle. Finned surfaces in dusty service lose capacity between cleans, and the quoted overdesign is consumed by it.'
    else:
        msg = f'Air-side fouling allowance {observed:.6f} m²K/W is specified.'
    remedy = None if passed else 'Specify an air-side fouling allowance from the project design basis and require the vendor to re-rate against it.'
    if tube_side_fouling_m2K_W is not None and tube_side_fouling_m2K_W <= 0:
        msg += ' Tube-side fouling is also zero.'
    return _finding(rule, passed=passed, message=msg, base_type=FindingType.PURCHASER_DEVIATION, observed=observed, threshold=minimum_air_side, remedy=remedy)
CORRELATION_RULES: Dict[str, Rule] = {'fin_density_range': Rule(key='fin_density_range', clause=None, standard=None, summary='High-fin correlations were fitted over roughly 275-435 fins/m.', source='ESDU high-fin and Briggs-Young applicability ranges.', verified=False), 'mtd_correction_support': Rule(key='mtd_correction_support', clause=None, standard=None, summary='Ft has dedicated coefficients only for certain row/pass combinations.', source='ht.air_cooler.Ft_aircooler coefficient sets.', verified=False), 'correlation_spread': Rule(key='correlation_spread', clause=None, standard=None, summary='ESDU and Briggs-Young should agree closely inside their fitted range.', source='Cross-check between two independent air-side correlations.', verified=False)}

def check_correlation_range(*, fin_density_per_m: Optional[float]=None, ft_support_tier: Optional[str]=None, correlation_spread_pct: Optional[float]=None, fin_density_range=(275.0, 435.0), max_spread_pct: float=10.0) -> List[Finding]:
    """Where our own rating is softer than it looks, as structured findings.

    These are limits on *our* numbers rather than on the vendor's design, and burying them
    in free text lets a reviewer quote a rated value that the correlation never supported.
    """
    findings: List[Finding] = []
    if fin_density_per_m is not None:
        lo, hi = fin_density_range
        inside = lo <= fin_density_per_m <= hi
        findings.append(Finding(key='fin_density_range', finding_type=FindingType.CORRELATION_RANGE_WARNING, passed=inside, message=f'Fin density {fin_density_per_m:.0f} /m is within the {lo:.0f}-{hi:.0f} /m range the high-fin correlations were fitted over.' if inside else f'Fin density {fin_density_per_m:.0f} /m lies outside the {lo:.0f}-{hi:.0f} /m range the high-fin correlations were fitted over, so the air-side coefficient is an extrapolation.', source=CORRELATION_RULES['fin_density_range'].source, observed=fin_density_per_m, threshold=hi, remedy=None if inside else 'Ask the vendor for test data at this fin density.'))
    if ft_support_tier is not None:
        passed = ft_support_tier == 'exact'
        findings.append(Finding(key='mtd_correction_support', finding_type=FindingType.CORRELATION_RANGE_WARNING, passed=passed, message={'exact': 'A dedicated Ft coefficient set exists for this row/pass count.', 'substituted': 'Ft uses a substituted coefficient set for this row/pass count, documented upstream as a reasonable assumption.'}.get(ft_support_tier, 'No Ft coefficient set exists for this row/pass count; the correlation falls back to one its own source calls a bad assumption. The MTD correction is an approximation, not a rated value.'), source=CORRELATION_RULES['mtd_correction_support'].source, remedy=None if passed else 'Ask the vendor for the MTD correction their software applied.'))
    if correlation_spread_pct is not None:
        inside = correlation_spread_pct <= max_spread_pct
        findings.append(Finding(key='correlation_spread', finding_type=FindingType.CORRELATION_RANGE_WARNING, passed=inside, message=f'ESDU and Briggs-Young agree to {correlation_spread_pct:.1f}%, which is the expected level inside their fitted range.' if inside else f'ESDU and Briggs-Young differ by {correlation_spread_pct:.1f}%, which means this geometry sits near or outside the range either was built for. Treat the air-side coefficient accordingly.', source=CORRELATION_RULES['correlation_spread'].source, observed=correlation_spread_pct, threshold=max_spread_pct))
    return findings

@dataclass
class ConformanceReport:
    findings: List[Finding] = field(default_factory=list)
    edition_note: str = API661_EDITION_NOTE

    @property
    def nonconformances(self) -> List[Finding]:
        return [f for f in self.findings if not f.passed and f.finding_type is FindingType.STANDARD_NONCONFORMANCE]

    def as_dict(self) -> Dict[str, object]:
        return {'edition_note': self.edition_note, 'findings': [f.as_dict() for f in self.findings], 'nonconformance_count': len(self.nonconformances)}

def run_checks(*, fan_diameter_m: Optional[float]=None, bundle_face_area_m2: Optional[float]=None, fans_per_bay: int=1, parallel_bays: int=1, face_velocity_m_s: Optional[float]=None, air_side_fouling_m2K_W: Optional[float]=None, tube_side_fouling_m2K_W: Optional[float]=None, purchaser_agreed_single_fan: bool=False, api661_invoked: bool=True, fin_density_per_m: Optional[float]=None, ft_support_tier: Optional[str]=None, correlation_spread_pct: Optional[float]=None) -> ConformanceReport:
    """Run every check for which the necessary vendor data is present.

    Checks whose inputs are missing are simply not run — an absent field is a data gap for
    the callout list, not a pass.
    """
    findings: List[Finding] = []
    if fan_diameter_m and bundle_face_area_m2:
        findings.append(check_fan_coverage(fan_diameter_m=fan_diameter_m, bundle_face_area_m2=bundle_face_area_m2, fans_per_bay=fans_per_bay, parallel_bays=parallel_bays, api661_invoked=api661_invoked))
    findings.append(check_fans_per_bay(fans_per_bay=fans_per_bay, purchaser_agreed_single=purchaser_agreed_single_fan, api661_invoked=api661_invoked))
    if face_velocity_m_s is not None:
        findings.append(check_face_velocity(face_velocity_m_s=face_velocity_m_s))
    findings.append(check_fouling_allowance(air_side_fouling_m2K_W=air_side_fouling_m2K_W, tube_side_fouling_m2K_W=tube_side_fouling_m2K_W))
    findings.extend(check_correlation_range(fin_density_per_m=fin_density_per_m, ft_support_tier=ft_support_tier, correlation_spread_pct=correlation_spread_pct))
    return ConformanceReport(findings=findings)
