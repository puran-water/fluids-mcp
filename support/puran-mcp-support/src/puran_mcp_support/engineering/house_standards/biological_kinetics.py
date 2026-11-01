from __future__ import annotations
import math
from dataclasses import dataclass, field
from decimal import Decimal
from typing import Literal, Mapping
from puran_mcp_support.engineering.house_standards.performance_bands import BandSource
RULINGS = {'kinetics': 'PUBLIC-RULE-462b902c614f', 'reference_15c': 'PUBLIC-RULE-5f01f35c0398', 'theta_source': 'PUBLIC-RULE-106214ea3d31', 'above_range': 'PUBLIC-RULE-a4a8099c86e2', 'y_po4': 'PUBLIC-RULE-d6ccf0230992', 'corrections_required': 'PUBLIC-RULE-36295b4c8f1e'}
NATIVE_RATE_UNIT = '1/d'
NATIVE_DIMENSIONLESS_UNIT = 'dimensionless'
_INSTALLED = 'watertap/property_models/unit_specific/activated_sludge'
_WATERTAP_VERSION = '1.8.0'

class KineticsTemperatureRangeError(ValueError):
    """A correction was requested below the fitted range, which no ruling covers."""

@dataclass(frozen=True)
class ReferenceSet:
    """One installed set of rate constants at a stated temperature, with file:line per value."""
    key: str
    source_file: str
    source_lines: Mapping[str, int]
    values: Mapping[str, float]
    temperature_c: float
    temperature_source: str

    def cite(self, name: str) -> str:
        return f'{self.values[name]:g} 1/d at {self.temperature_c:g} C, WaterTAP {_WATERTAP_VERSION} {self.source_file}:{self.source_lines[name]}'
RATE_CONSTANTS: tuple[str, ...] = ('K_H', 'mu_H', 'q_fe', 'b_H', 'q_PHA', 'q_PP', 'mu_PAO', 'b_PAO', 'b_PP', 'b_PHA', 'mu_AUT', 'b_AUT')
RULED_REFERENCE_SET = 'watertap_1.7.0_modified_asm2d_15c'
REFERENCE_SETS: dict[str, ReferenceSet] = {RULED_REFERENCE_SET: ReferenceSet(key=RULED_REFERENCE_SET, source_file=f'{_INSTALLED}/modified_asm2d_reactions.py', source_lines={'K_H': 278, 'mu_H': 351, 'q_fe': 357, 'b_H': 363, 'q_PHA': 430, 'q_PP': 436, 'mu_PAO': 442, 'b_PAO': 448, 'b_PP': 454, 'b_PHA': 460, 'mu_AUT': 497, 'b_AUT': 503}, values={'K_H': 2.46, 'mu_H': 4.23, 'q_fe': 2.11, 'b_H': 0.28, 'q_PHA': 2.46, 'q_PP': 1.23, 'mu_PAO': 0.82, 'b_PAO': 0.14, 'b_PP': 0.14, 'b_PHA': 0.14, 'mu_AUT': 0.61, 'b_AUT': 0.09}, temperature_c=15.0, temperature_source=f"{RULINGS['reference_15c']} (the WaterTAP docs' 'Value at 20°C' label, modified_ASM2D.rst:134, is overruled)")}
ANCHOR_20C = ReferenceSet(key='watertap_1.7.0_asm2d_20c', source_file=f'{_INSTALLED}/asm2d_reactions.py', source_lines={'K_H': 238, 'mu_H': 275, 'q_fe': 281, 'b_H': 287, 'q_PHA': 330, 'q_PP': 336, 'mu_PAO': 342, 'b_PAO': 348, 'b_PP': 354, 'b_PHA': 360, 'mu_AUT': 397, 'b_AUT': 403}, values={'K_H': 3.0, 'mu_H': 6.0, 'q_fe': 3.0, 'b_H': 0.4, 'q_PHA': 3.0, 'q_PP': 1.5, 'mu_PAO': 1.0, 'b_PAO': 0.2, 'b_PP': 0.2, 'b_PHA': 0.2, 'mu_AUT': 1.0, 'b_AUT': 0.15}, temperature_c=20.0, temperature_source="WaterTAP 1.7.0 docs ASM2D.rst:96 'Value at 20°C' (Henze et al. 1999 ASM2d, doi 10.2166/wst.1999.0036); the same values are QSDsan mASM2d's defaults, qsdsan process_models/_asm2d.py:919-921")
THETA_FIT_RANGE_C = (REFERENCE_SETS[RULED_REFERENCE_SET].temperature_c, ANCHOR_20C.temperature_c)
THETA_SOURCE = f"{RULINGS['theta_source']} FALLBACK theta = (k20/k15)^(1/5): the ASM2d report's 10 C column was not readable (Henze et al. 1999 doi 10.2166/wst.1999.0036 and IWA STR 9 doi 10.2166/9781780402369 closed-access)"

def _theta(k_low: float, k_high: float) -> float:
    span = THETA_FIT_RANGE_C[1] - THETA_FIT_RANGE_C[0]
    return (k_high / k_low) ** (1.0 / span)

def theta_band(name: str) -> tuple[float, float]:
    """Theta's range from k15 +/- half a unit in its last printed place (the source's own rounding)."""
    k15 = REFERENCE_SETS[RULED_REFERENCE_SET].values[name]
    half = float(Decimal(10) ** Decimal(repr(k15)).as_tuple().exponent) / 2.0
    k20 = ANCHOR_20C.values[name]
    return (_theta(k15 + half, k20), _theta(k15 - half, k20))
THETAS: dict[str, float] = {name: _theta(REFERENCE_SETS[RULED_REFERENCE_SET].values[name], ANCHOR_20C.values[name]) for name in RATE_CONSTANTS}
Guild = Literal['heterotroph', 'autotroph', 'pao']
GUILD_OF: dict[str, tuple[Guild, str]] = {'K_H': ('heterotroph', 'hydrolysis R1-R3'), 'mu_H': ('heterotroph', 'heterotrophic growth R4-R7'), 'q_fe': ('heterotroph', 'fermentation R8'), 'b_H': ('heterotroph', 'heterotrophic lysis R9'), 'q_PHA': ('pao', 'PHA storage R10'), 'q_PP': ('pao', 'polyphosphate storage R11-R12'), 'mu_PAO': ('pao', 'PAO growth R13-R14'), 'b_PAO': ('pao', 'PAO lysis R15'), 'b_PP': ('pao', 'polyphosphate lysis R16'), 'b_PHA': ('pao', 'PHA lysis R17'), 'mu_AUT': ('autotroph', "lumped nitrifier growth R18 (ASM2d's own X_AUT pair)"), 'b_AUT': ('autotroph', "lumped nitrifier lysis R19 (ASM2d's own X_AUT pair)")}
_YIELD = 'yield coefficient: held constant across temperature in the ASM family'
HENZE_1999 = BandSource(kind='peer_reviewed', doi='10.2166/wst.1999.0036', citation='Henze, Gujer, Mino, Matsuo, Wentzel, Marais & van Loosdrecht (1999) Activated Sludge Model No.2d, ASM2d, Water Sci. Technol. 39(1):165', locator=f'Y_PO4 0.40 g P/g COD, polyphosphate requirement (PO4 release) per PHA stored; read through WaterTAP {_WATERTAP_VERSION} {_INSTALLED}/asm2d_reactions.py:217-221 (the paper is closed-access)')
SMOLDERS_1994 = BandSource(kind='primary_study', doi='10.1002/bit.260430605', citation='Smolders, van der Meij, van Loosdrecht & Heijnen (1994) Model of the anaerobic metabolism of the biological phosphorus removal process: stoichiometry and pH influence, Biotechnol. Bioeng. 43', locator='abstract: observed phosphate/acetate ratio 0.25-0.75 P-mol/C-mol over pH 5.5-8.5, i.e. 0.24-0.73 g P/g COD at 30.97 g P/mol and 32 g COD per C-mol of acetate')

@dataclass(frozen=True)
class ParameterCorrection:
    """A ModifiedASM2d reaction parameter whose installed default the house corrects, at every temperature."""
    name: str
    value: float
    unit: str
    installed_value: float
    installed_line: int
    anchor_value: float
    anchor_line: int
    statement: str
    sources: tuple[BandSource, ...]

    def clause(self) -> str:
        cited = '; '.join((f'{source.citation} (doi {source.doi}): {source.locator}' for source in self.sources))
        return f'{self.name} {self.value:g} ({self.unit}): house correction of the installed {self.installed_value:g}, WaterTAP {_WATERTAP_VERSION} {REFERENCE_SETS[RULED_REFERENCE_SET].source_file}:{self.installed_line} ({RULINGS['y_po4']}); held at every temperature ({_YIELD}); {self.statement}; base ASM2d {self.anchor_value:g}, {ANCHOR_20C.source_file}:{self.anchor_line}; {cited}'
CORRECTED_PARAMETERS: dict[str, ParameterCorrection] = {'Y_PO4': ParameterCorrection(name='Y_PO4', value=0.4, unit=NATIVE_DIMENSIONLESS_UNIT, installed_value=0.0129, installed_line=240, anchor_value=0.4, anchor_line=218, statement="kg P released per kg COD of PHA stored (R10), on X_PP's kg-P basis. The installed 0.0129 is 0.40/31, the kmol value: WaterTAP's own ModifiedADM1 carries 12.903e-3 labelled kmol P/kg COD (anaerobic_digestion/modified_adm1_reactions.py:719), and the house retired the same 31-fold error on the digester side (process_model_schema/madm1_biokinetics.py:20-22, watertap-engine-mcp chemistry/modified_adm1_kernel.py:16-17)", sources=(HENZE_1999, SMOLDERS_1994))}
GRUNDITZ_2001 = BandSource(kind='primary_study', doi='10.1016/s0043-1354(00)00312-2', citation='Grunditz & Dalhammar (2001) Development of nitrification inhibition assays using pure cultures of Nitrosomonas and Nitrobacter, Water Res. 35:433', locator='abstract: optimum temperature 35 C for Nitrosomonas and 38 C for Nitrobacter')
LAPARA_2000 = BandSource(kind='primary_study', doi='10.1038/sj.jim.2900790', citation='LaPara, Konopka, Nakatsu & Alleman (2000) Effects of elevated temperature on bacterial community structure and function in bioreactors treating a synthetic wastewater, J. Ind. Microbiol. Biotechnol. 24:140-145', locator='open-access full text, Results and Fig. 3: mu_m rose 1.4x over 25-45 C then declined almost by half to 65 C; yield followed the same parabola; q_m relatively constant (0.64 ± 0.04 mg COD/mg dry cell/h) over 25-65 C; highest mu_m at 45 C')
MORGAN_SAGASTUME_2003 = BandSource(kind='primary_study', doi='10.1016/s0043-1354(03)00270-7', citation='Morgan-Sagastume & Allen (2003) Effects of temperature transient conditions on aerobic biological treatment of wastewater, Water Res. 37', locator='abstract: pulp-mill activated sludge SBRs at 35 C; shifts 35->45 C cut SCOD removal from 66-65% to 49-51% and raised effluent solids; 31.5-40 C oscillations degraded settleability')
LOPEZ_VAZQUEZ_2009 = BandSource(kind='peer_reviewed', doi='10.1016/j.watres.2008.10.032', citation='Lopez-Vazquez et al. (2009) Modeling the PAO-GAO competition: effects of carbon source, pH and temperature, Water Res. 43', locator='abstract: evaluated at 10, 20 and 30 C; PAO dominant at 10 C; at 30 C GAO tended to dominate, mixed acetate/propionate or pH 7.5 favouring PAO')
WHANG_PARK_2006 = BandSource(kind='primary_study', doi='10.2175/106143005x84459', citation='Whang & Park (2006) Competition between polyphosphate- and glycogen-accumulating organisms in EBPR systems: effect of temperature and sludge age, Water Environ. Res. 78', locator='abstract: PAO dominant at 20 C / 10 d SRT; GAO outcompeted PAO at 30 C / 10 d; unstable at 30 C / 5 d; stable EBPR at 30 C / 3 d')
PANSWAD_2003 = BandSource(kind='primary_study', doi='10.1016/s0043-1354(02)00286-5', citation='Panswad, Doungchai & Anotai (2003) Temperature effect on microbial community of enhanced biological phosphorus removal system, Water Res. 37', locator="abstract: 20 -> 30 -> 35.5 C moved the dominant group from PAO (47-70% VSS) to GAO (64-75%) to ordinary heterotrophs (90%); PAO 'lower-range mesophiles', washed out")
NIELSEN_2019 = BandSource(kind='peer_reviewed', doi='10.1016/j.copbio.2019.03.008', citation='Nielsen, McIlroy, Albertsen & Nierychlo (2019) Re-evaluating the microbiology of the enhanced biological phosphorus removal process, Curr. Opin. Biotechnol.', locator='abstract highlights: EBPR works very well under tropical conditions; GAOs do not appear to be a real problem in full-scale EBPR (counter-evidence to the lab GAO bands)')
_TOP = THETA_FIT_RANGE_C[1]
_HELD = f'above the {_TOP:g} C top of the fit the rate is HELD there as a bound ({RULINGS['above_range']})'
NITRIFIER_HIGH_T_BANDS: tuple[tuple[float, float | None, str, str], ...] = ((35.0, 38.0, 'near_optimum', 'at/near the measured optima (AOB 35 C, NOB 38 C -- Grunditz & Dalhammar 2001): rates plateau, so theta^(T-20) OVER-predicts; use the 35 C value as a ceiling, not the extrapolation'), (38.0, 42.0, 'declining', 'above both optima: measured rates DECLINE with temperature while theta^(T-20) rises -- the correction points the wrong direction. Treat rates as degraded and unquantified without site data'), (42.0, None, 'no_design_basis', 'no mesophilic design basis: beyond the measured activity range of both genera. Treat nitrification as infeasible at this temperature unless thermophilic-specific evidence is produced'))
_NITRIFIER_BANDS = ((None, _TOP, 'within_fit', 'inside the fitted theta range; theta-corrected'), (_TOP, 35.0, 'above_fit_below_optima', f'{_HELD}; the measured optima (AOB 35 C, NOB 38 C, Grunditz & Dalhammar 2001) lie above, so the held rate understates nitrifier activity (the risk the ruling accepts)'), *NITRIFIER_HIGH_T_BANDS)
HETEROTROPH_HIGH_T_BANDS: tuple[tuple[float | None, float | None, str, str], ...] = ((None, _TOP, 'within_fit', 'inside the fitted theta range; theta-corrected'), (_TOP, 45.0, 'above_fit_rising_slowly', f'{_HELD}. Acclimated mixed consortia raised mu_max only 1.4x from 25 to 45 C while q_max stayed ~constant (LaPara et al. 2000), so the held rate understates growth by at most about that factor, where theta extrapolation would have overstated it ~3x. Temperature swings in this band degrade settling (Morgan-Sagastume & Allen 2003)'), (45.0, None, 'declining_no_mesophilic_basis', 'past the measured mesophilic peak (45 C): mu_max fell by almost half from 45 to 65 C (LaPara et al. 2000), and a 35->45 C shift cut soluble COD removal from ~66% to ~49% in pulp-mill activated sludge (Morgan-Sagastume & Allen 2003). No mesophilic design basis'))
PAO_TEMPERATURE_BANDS: tuple[tuple[float | None, float | None, str, str], ...] = ((None, 20.0, 'pao_favoured', 'PAO dominant at 10 C regardless of carbon source or pH (Lopez-Vazquez et al. 2009) and at 20 C / 10 d SRT (Whang & Park 2006)'), (20.0, 30.0, 'gao_competition', f'{_HELD}. EBPR stability depends on carbon source, pH and SRT: at 30 C GAO outcompeted PAO at 10 d SRT and EBPR was stable only at 3 d (Whang & Park 2006; Lopez-Vazquez et al. 2009)'), (30.0, 35.0, 'gao_dominance', f'{_HELD}. GAO were 64-75% of VSS at 30 C (Panswad et al. 2003) and tended to dominate in the 30 C model (Lopez-Vazquez et al. 2009, not evaluated above 30 C); full-scale tropical EBPR is reported to work (Nielsen et al. 2019), so EBPR here is a risk to demonstrate, not a design basis'), (35.0, None, 'pao_washout', f'{_HELD}. PAO washed out: at 35.5 C ordinary heterotrophs were 90% of VSS (Panswad et al. 2003). No EBPR design basis'))
PAO_NO_EBPR_DESIGN_BASIS_BANDS = frozenset({'gao_dominance', 'pao_washout'})

def _band(bands, temperature_c: float, guild: str, sources) -> dict:
    t = float(temperature_c)
    for lo, hi, band, statement in bands:
        if (lo is None or t > lo) and (hi is None or t <= hi):
            return {'guild': guild, 'band': band, 'temperature_c': t, 'statement': statement, 'policy': f'theta-corrected {THETA_FIT_RANGE_C[0]:g}-{_TOP:g} C; held at the {_TOP:g} C value as a bound above ({RULINGS['above_range']}); never extrapolated', 'sources': [source.model_dump() for source in sources]}
    raise AssertionError('unreachable: bands are contiguous')

def nitrifier_high_temperature_screen(temperature_c: float) -> dict:
    """Literature-bounded statement for nitrifiers (the ported Grunditz & Dalhammar bands)."""
    return _band(_NITRIFIER_BANDS, temperature_c, 'autotroph', (GRUNDITZ_2001,))

def heterotroph_high_temperature_screen(temperature_c: float) -> dict:
    """Literature-bounded statement for heterotroph growth/decay/hydrolysis above the fit."""
    return _band(HETEROTROPH_HIGH_T_BANDS, temperature_c, 'heterotroph', (LAPARA_2000, MORGAN_SAGASTUME_2003))

def pao_temperature_screen(temperature_c: float) -> dict:
    """PAO vs GAO competition by temperature; ``ebpr_design_basis`` false in `PAO_NO_EBPR_DESIGN_BASIS_BANDS`."""
    screen = _band(PAO_TEMPERATURE_BANDS, temperature_c, 'pao', (LOPEZ_VAZQUEZ_2009, WHANG_PARK_2006, PANSWAD_2003, NIELSEN_2019))
    return {**screen, 'ebpr_design_basis': screen['band'] not in PAO_NO_EBPR_DESIGN_BASIS_BANDS}
_SCREENS = {'heterotroph': heterotroph_high_temperature_screen, 'autotroph': nitrifier_high_temperature_screen, 'pao': pao_temperature_screen}

@dataclass(frozen=True)
class ParameterRow:
    name: str
    guild: str
    process: str
    reference_value: float
    reference_ref: str
    anchor_value: float
    anchor_ref: str
    theta: float
    theta_band: tuple[float, float]
    fitted_range_c: tuple[float, float]
    value: float
    standing: Literal['corrected', 'bound']

@dataclass(frozen=True)
class NativeKineticsSet:
    temperature_k: float
    temperature_c: float
    reference_set: str
    reaction_parameter_fixes: dict[str, float]
    parameter_basis: dict
    basis_content_hash: str
    rows: tuple[ParameterRow, ...]
    bounded: tuple[str, ...]
    screens: dict[str, dict] = field(default_factory=dict)

    def report(self) -> str:
        """The set and its screens as printable text (the ruling prints the screens beside the rates)."""
        lines = [f'ModifiedASM2d rate constants at {self.temperature_c:g} C ({self.reference_set}; {THETA_SOURCE})', f'{'parameter':<8} {'k15':>6} {'k20':>6} {'theta':>7} {'band':>17} {'value':>9}  standing']
        for row in self.rows:
            lines.append(f'{row.name:<8} {row.reference_value:>6g} {row.anchor_value:>6g} {row.theta:>7.4f} [{row.theta_band[0]:.4f},{row.theta_band[1]:.4f}] {row.value:>9.4g}  {row.standing}')
        for correction in CORRECTED_PARAMETERS.values():
            lines.append(f'{correction.name:<8} installed {correction.installed_value:g} -> {correction.value:g} {correction.unit}  corrected ({RULINGS['y_po4']})')
        for guild, screen in self.screens.items():
            lines.append(f'screen {guild}: {screen['band']} -- {screen['statement']}')
        return '\n'.join(lines)

def native_kinetics_at(*, temperature_k: float, temperature_source_clause_ref: str, reference_set: str=RULED_REFERENCE_SET) -> NativeKineticsSet:
    """The ModifiedASM2d rate constants at ``temperature_k`` and the house corrections, with their native evidence.

    ``temperature_k`` has NO default: it comes from the heat balance, and
    ``temperature_source_clause_ref`` says which one. Between the fitted 15 and 20 C each rate is
    ``k15 * theta^(T - 15)``; above 20 C it is held at its 20 C value, stamped as a bound
    (``PUBLIC-RULE-a4a8099c86e2``); below 15 C the call refuses. There is no extrapolation path.

    Build the criteria with ``BiologyCriteria(kinetic_temperature_k=result.temperature_k,
    reaction_parameter_fixes=result.reaction_parameter_fixes, parameter_evidence_ids=(<id>,), ...)``;
    ``validate_native_parameter_basis(criteria, result.parameter_basis)`` then accepts the payload.
    To combine with composition fixes use :func:`merge_native_parameter_bases`.
    """
    from puran_mcp_support.engineering.process_balance.adapters.biology import COMPONENT_NAMESPACE, MODEL_KEY, NativeBiologyParameterBasis
    from puran_mcp_support.engineering.process_sizing.models import content_hash
    if isinstance(temperature_k, bool) or not isinstance(temperature_k, (int, float)) or (not math.isfinite(temperature_k)) or (temperature_k <= 0):
        raise ValueError('temperature_k must be a finite absolute temperature in K')
    if not isinstance(temperature_source_clause_ref, str) or not temperature_source_clause_ref.strip():
        raise ValueError('temperature_source_clause_ref must name where the basin temperature came from')
    if reference_set not in REFERENCE_SETS:
        raise ValueError(f'reference_set must be one of {sorted(REFERENCE_SETS)} ({RULINGS['reference_15c']})')
    reference = REFERENCE_SETS[reference_set]
    temperature_k = float(temperature_k)
    temperature_c = round(temperature_k - 273.15, 9)
    low, top = THETA_FIT_RANGE_C
    if temperature_c < low:
        raise KineticsTemperatureRangeError(f"{temperature_c:g} C is below the {low:g}-{top:g} C range the theta is fitted over, and no ruling covers holding below it ({RULINGS['above_range']} governs only above). Source the ASM2d report's 10 C column to extend the fit downward.")
    screens = {guild: screen(temperature_c) for guild, screen in _SCREENS.items()}
    rows = []
    for name in RATE_CONSTANTS:
        guild, process = GUILD_OF[name]
        k15, k20 = (reference.values[name], ANCHOR_20C.values[name])
        if temperature_c >= top:
            value, standing = (k20, 'bound' if temperature_c > top else 'corrected')
        else:
            value, standing = (k15 * THETAS[name] ** (temperature_c - low), 'corrected')
        rows.append(ParameterRow(name, guild, process, k15, reference.cite(name), k20, ANCHOR_20C.cite(name), THETAS[name], theta_band(name), THETA_FIT_RANGE_C, value, standing))
    fields = {'kinetic_temperature_k': {'value': temperature_k, 'unit': 'K', 'source_clause_ref': f'basin temperature: {temperature_source_clause_ref.strip()}; ModifiedASM2d rate constants set to it by engineering_utils.house_standards.biological_kinetics under {RULINGS['kinetics']}, {RULINGS['reference_15c']}, {RULINGS['theta_source']}, {RULINGS['above_range']}'}}
    for row in rows:
        derivation = f'k(T) = k15*theta^(T-{low:g}) with k15 = {row.reference_ref} ({reference.temperature_source}); theta {row.theta:.5f} (rounding band {row.theta_band[0]:.5f}-{row.theta_band[1]:.5f}) = (k20/k15)^(1/{top - low:g}), k20 = {row.anchor_ref}; {THETA_SOURCE}; {row.process}'
        if row.standing == 'bound':
            head = f"{row.name} at {temperature_c:g} C: BOUND -- held at the {top:g} C top of the fitted {low:g}-{top:g} C range ({RULINGS['above_range']}); {row.guild} screen band '{screens[row.guild]['band']}'"
        else:
            head = f'{row.name} at {temperature_c:g} C: corrected inside the fitted {low:g}-{top:g} C range'
        fields[f'reaction_parameter_fixes.{row.name}'] = {'value': row.value, 'unit': NATIVE_RATE_UNIT, 'source_clause_ref': f'{head}; {derivation}'}
    for correction in CORRECTED_PARAMETERS.values():
        fields[f'reaction_parameter_fixes.{correction.name}'] = {'value': correction.value, 'unit': correction.unit, 'source_clause_ref': correction.clause()}
    payload = {'schema_version': 'native_biology_parameter_basis/v1', 'model_key': MODEL_KEY, 'component_namespace': COMPONENT_NAMESPACE, 'fields': fields}
    NativeBiologyParameterBasis.model_validate(payload)
    return NativeKineticsSet(temperature_k=temperature_k, temperature_c=temperature_c, reference_set=reference.key, reaction_parameter_fixes={**{row.name: row.value for row in rows}, **{name: correction.value for name, correction in CORRECTED_PARAMETERS.items()}}, parameter_basis=payload, basis_content_hash=content_hash(payload), rows=tuple(rows), bounded=tuple((row.name for row in rows if row.standing == 'bound')), screens=screens)
