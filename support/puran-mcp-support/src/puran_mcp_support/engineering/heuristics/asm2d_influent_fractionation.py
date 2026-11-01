from __future__ import annotations
import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, ClassVar, Literal
from uuid import UUID
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator
BOD5_INPUT_ALIASES: tuple[str, ...] = ('BOD5', 'bod5', 'BOD5_mg_L', 'bod5_mg_l', 'BOD_mg_L', 'bod_mg_l')
Bod5Basis = Literal['measured', 'assumed_cod_ratio', 'class_default']
CharacteristicBasis = Literal['measured', 'derived', 'class_default']

class UnmeasuredCharacteristicError(ValueError):
    """A characteristic needed for the calculation was not supplied.

    Raised instead of substituting a generic value. A default waste-stream
    profile is a description of an industry, not of this effluent, and a design
    basis built on one cannot be checked against the plant it was written for.
    Supply the value from the project documentation, or -- only for pre-tender
    screening where no analysis exists yet -- pass
    ``allow_unmeasured_characteristics=True`` and accept that every value so
    filled is stamped and may not carry a design decision.

    Raised inside a pydantic validator, so it reaches a caller wrapped in
    ``pydantic.ValidationError``. Catch that, or read ``.errors()[0]["ctx"]
    ["error"]`` for this instance; the message is carried through intact either
    way.
    """

def _to_float(value: Any) -> float | None:
    """Coerce a possibly-string numeric to float; return None for blanks/None."""
    if value is None:
        return None
    if isinstance(value, str):
        stripped = value.strip().replace(',', '')
        if stripped == '':
            return None
        return float(stripped)
    return float(value)

class InfluentCharacter(str, Enum):
    """Canonical influent character used to select fractionation defaults."""
    MUNICIPAL = 'municipal'
    INDUSTRIAL_CETP = 'industrial_cetp'
    DAIRY = 'dairy'
    TEXTILE = 'textile'
    PHARMA = 'pharma'
    REFINERY_CHEM = 'refinery_chem'
    MIXED_INDUSTRIAL = 'mixed_industrial'

class TenderInfluentBasis(BaseModel):
    """Bulk influent basis extracted from the canonical tender artifact.

    Concentration fields are mg/L except pH, temperature_C, and flow_m3_d.
    Alkalinity is interpreted as mg/L as CaCO3, matching the Stage 5 sizing
    convention and the WaterTAP ASM2D alkalinity input conversion.
    """
    model_config = ConfigDict(extra='ignore', populate_by_name=True)
    BOD5: float | None = Field(default=None, ge=0)
    COD: float = Field(gt=0)
    TSS: float | None = Field(default=None, ge=0)
    VSS: float | None = Field(default=None, ge=0)
    TKN: float | None = Field(default=None, ge=0)
    NH4_N: float | None = Field(default=None, ge=0)
    NO3_N: float | None = Field(default=None, ge=0)
    TP: float | None = Field(default=None, ge=0)
    PO4_P: float | None = Field(default=None, ge=0)
    FOG: float | None = Field(default=None, ge=0)
    alkalinity: float | None = Field(default=None, ge=0)
    pH: float | None = Field(default=7.2, ge=0, le=14)
    temperature_C: float | None = Field(default=25.0, ge=0, le=80)
    flow_m3_d: float = Field(gt=0)
    influent_character: InfluentCharacter = InfluentCharacter.MUNICIPAL
    bod5_basis: Bod5Basis | None = None
    BOD_decay_k_per_d: float | None = Field(default=None, gt=0, le=2.0)
    source_influent_artifact_id: UUID | None = None
    TDS: float | None = Field(default=None, ge=0)
    hardness_as_CaCO3: float | None = Field(default=None, ge=0)
    chloride: float | None = Field(default=None, ge=0)
    sulphate: float | None = Field(default=None, ge=0)
    silica_as_SiO2: float | None = Field(default=None, ge=0)
    calcium: float | None = Field(default=None, ge=0)
    magnesium: float | None = Field(default=None, ge=0)
    sodium: float | None = Field(default=None, ge=0)
    allow_unmeasured_characteristics: bool = False
    characteristic_basis: dict[str, CharacteristicBasis] = Field(default_factory=dict)

    @model_validator(mode='before')
    @classmethod
    def _normalize_aliases(cls, raw: Any) -> Any:
        if not isinstance(raw, dict):
            return raw
        data = dict(raw)
        for nested_key in ('components', 'concentrations_mg_L', 'concentrations_mg_l'):
            nested = data.get(nested_key)
            if isinstance(nested, dict):
                for key, value in nested.items():
                    data.setdefault(key, value)
        aliases = {'BOD5': BOD5_INPUT_ALIASES, 'COD': ('COD', 'cod', 'COD_mg_L', 'cod_mg_l'), 'TSS': ('TSS', 'tss', 'TSS_mg_L', 'tss_mg_l'), 'VSS': ('VSS', 'vss', 'VSS_mg_L', 'vss_mg_l'), 'TKN': ('TKN', 'tkn', 'TKN_mg_L', 'tkn_mg_l'), 'NH4_N': ('NH4_N', 'nh4_n', 'NH4_N_mg_L', 'nh4_n_mg_l', 'ammonia_n_mg_l'), 'NO3_N': ('NO3_N', 'NO3_N_mg_L', 'no3_n', 'NO3N', 'NO3_N_mg_l', 'no3_n_mg_l', 'nitrate_n_mg_l'), 'TP': ('TP', 'tp', 'TP_mg_L', 'tp_mg_l', 'total_p_mg_l', 'phosphorus_mg_l'), 'PO4_P': ('PO4_P', 'po4_p', 'PO4_P_mg_L', 'po4_p_mg_l', 'orthophosphate_p_mg_l'), 'FOG': ('FOG', 'FOG_mg_L', 'fog', 'oil_grease_mg_l', 'fog_mg_l', 'oil_and_grease_mg_l'), 'alkalinity': ('alkalinity', 'alkalinity_mg_L', 'alkalinity_mg_l', 'alkalinity_mg_caco3_l', 'alkalinity_mg_L_as_CaCO3', 'alkalinity_mg_l_as_caco3'), 'pH': ('pH', 'ph'), 'temperature_C': ('temperature_C', 'temperature_c', 'temp_C', 'temp_c'), 'flow_m3_d': ('flow_m3_d', 'flow_m3_per_day', 'avg_flow_m3_d', 'design_flow_m3_d', 'feed_flow_m3d', 'feed_flow_m3_d'), 'influent_character': ('influent_character', 'wastewater_character', 'industry_character', 'industry_type'), 'source_influent_artifact_id': ('source_influent_artifact_id', 'artifact_id', 'influent_artifact_id')}
        for canonical, names in aliases.items():
            if canonical in data and data[canonical] is not None:
                continue
            for name in names:
                if name in data and data[name] is not None:
                    data[canonical] = data[name]
                    break
        if not data.get('bod5_basis'):
            data['bod5_basis'] = 'measured' if data.get('BOD5') is not None else 'assumed_cod_ratio'
        if data.get('TKN') is None:
            tn_value = None
            for name in ('TN', 'tn', 'TN_mg_L', 'TN_mg_l', 'tn_mg_l', 'total_n_mg_l', 'total_nitrogen_mg_l'):
                if name in data and data[name] is not None:
                    tn_value = _to_float(data[name])
                    break
            if tn_value is not None:
                no3_value = _to_float(data.get('NO3_N')) or 0.0
                if tn_value < no3_value:
                    raise ValueError(f'Influent basis inconsistent: TN ({tn_value}) < NO3_N ({no3_value}); total N must be >= oxidized nitrate-N')
                data['TKN'] = tn_value - no3_value
        temp_k = data.get('temperature_K') or data.get('temperature_k')
        if data.get('temperature_C') is None and temp_k is not None:
            data['temperature_C'] = float(temp_k) - 273.15
        character = data.get('influent_character')
        if isinstance(character, str):
            token = character.strip().lower().replace('-', '_').replace(' ', '_')
            aliases_by_character = {'industrial': 'mixed_industrial', 'cetp': 'industrial_cetp', 'chemical': 'refinery_chem', 'refinery': 'refinery_chem', 'domestic': 'municipal', 'stp': 'municipal'}
            data['influent_character'] = aliases_by_character.get(token, token)
        return data

    @field_validator('BOD5', 'COD', 'TSS', 'VSS', 'TKN', 'NH4_N', 'NO3_N', 'TP', 'PO4_P', 'FOG', 'alkalinity', mode='before')
    @classmethod
    def _coerce_numeric(cls, value: Any) -> Any:
        if isinstance(value, str):
            stripped = value.strip().replace(',', '')
            if stripped == '':
                return None
            return float(stripped)
        return value
    _PROFILE_DERIVATIONS: ClassVar[dict[str, str]] = {'BOD5': '0.50 x COD', 'TSS': '0.70 x COD', 'TKN': 'max(0.08 x COD, 20 mg/L)', 'TP': 'max(0.012 x COD, 4 mg/L)', 'alkalinity': '350 mg/L as CaCO3 flat'}
    _SPECIATION_DERIVATIONS: ClassVar[dict[str, tuple[str, float, str]]] = {'VSS': ('TSS', 0.75, '0.75 x TSS'), 'NH4_N': ('TKN', 0.6, '0.60 x TKN'), 'PO4_P': ('TP', 0.6, '0.60 x TP')}

    @model_validator(mode='after')
    def _resolve_characteristics(self) -> 'TenderInfluentBasis':
        """Fill what can honestly be filled; refuse the rest, and name it.

        Every characteristic here used to be back-filled from a COD ratio in
        silence, and the filled value was then indistinguishable from a
        measurement -- so a design basis could be built, checked and issued on
        numbers nobody ever observed. TKN is the case that shows the cost: a
        stream documented at 150 mg/L would be simulated at max(0.08 x COD, 20)
        = 360 mg/L instead, and the nitrogen answer -- the reason the simulation
        was run at all -- would be an artifact of the default.

        Profile characteristics therefore raise. Speciation splits off a measured
        parent are still filled, because that is how lab reports actually arrive,
        but they are stamped ``derived`` and a split off a derived parent stays
        derived rather than laundering itself into a measurement.
        """
        cod = float(self.COD)
        basis: dict[str, CharacteristicBasis] = dict(self.characteristic_basis)
        for name in ('BOD5', 'TSS', 'VSS', 'TKN', 'NH4_N', 'NO3_N', 'TP', 'PO4_P', 'alkalinity', 'TDS', 'hardness_as_CaCO3', 'chloride', 'sulphate', 'silica_as_SiO2', 'calcium', 'magnesium', 'sodium'):
            if getattr(self, name, None) is not None:
                basis.setdefault(name, 'measured')
        missing = [n for n in self._PROFILE_DERIVATIONS if getattr(self, n, None) is None]
        if missing and (not self.allow_unmeasured_characteristics):
            owed = ', '.join((f'{n} ({self._PROFILE_DERIVATIONS[n]})' for n in missing))
            raise UnmeasuredCharacteristicError(f"influent profile characteristics not supplied: {owed}. Filling these from generic ratios substitutes a description of an industry for a description of this effluent, and the filled values are indistinguishable from measurements in every artifact downstream. Supply them from the project documentation. If this is pre-tender screening with no analysis in hand, pass allow_unmeasured_characteristics=True; each value is stamped 'derived' and may not carry a design decision.")
        if self.BOD5 is None:
            self.BOD5 = cod * 0.5
            self.bod5_basis = 'assumed_cod_ratio'
            basis['BOD5'] = 'derived'
        elif self.bod5_basis is None:
            self.bod5_basis = 'measured'
        if self.TSS is None:
            self.TSS = cod * 0.7
            basis['TSS'] = 'derived'
        if self.TKN is None:
            self.TKN = max(cod * 0.08, 20.0)
            basis['TKN'] = 'derived'
        if self.TP is None:
            self.TP = max(cod * 0.012, 4.0)
            basis['TP'] = 'derived'
        if self.alkalinity is None:
            self.alkalinity = 350.0
            basis['alkalinity'] = 'derived'
        for child, (parent, ratio, _rule) in self._SPECIATION_DERIVATIONS.items():
            if getattr(self, child, None) is not None:
                continue
            parent_value = getattr(self, parent, None)
            if parent_value is None:
                continue
            setattr(self, child, float(parent_value) * ratio)
            basis[child] = 'derived'
        self.characteristic_basis = basis
        return self

    @property
    def unmeasured_characteristics(self) -> tuple[str, ...]:
        """Fields on this basis that no one measured, in artifact order.

        A design basis, a simulation report or a datasheet built from this basis
        must show these as owed rather than printing them as findings.
        """
        return tuple(sorted((n for n, b in self.characteristic_basis.items() if b != 'measured')))

    @property
    def cod_kg_m3(self) -> float:
        return self.COD / 1000.0

    @property
    def tkn_kg_m3(self) -> float:
        """Reduced (Kjeldahl) nitrogen basis in kg N/m3."""
        return float(self.TKN or 0.0) / 1000.0

    @property
    def no3_kg_m3(self) -> float:
        """Oxidized (nitrate) nitrogen basis in kg N/m3."""
        return float(self.NO3_N or 0.0) / 1000.0

    @property
    def tn_kg_m3(self) -> float:
        """Total nitrogen basis = TKN + NO3-N, in kg N/m3."""
        return self.tkn_kg_m3 + self.no3_kg_m3

    @property
    def fog_kg_m3(self) -> float:
        return float(self.FOG or 0.0) / 1000.0

    @property
    def tp_kg_m3(self) -> float:
        return float(self.TP or 0.0) / 1000.0

    @property
    def temperature_K(self) -> float:
        return float(self.temperature_C or 25.0) + 273.15

    @property
    def alkalinity_mol_per_m3(self) -> float:
        return float(self.alkalinity or 0.0) / 50.0
ASM2D_F_CV_AUT = 1.42
ASM2D_I_NBM = 0.07
ASM2D_I_NXI = 0.03
ASM2D_I_NXS = 0.04
ASM2D_I_PBM = 0.02
ASM2D_I_PXI = 0.01
ASM2D_I_PXS = 0.01
ASM2D_P_PER_XMEP = 1.0 / 4.87
ASM2D_COMPONENTS: tuple[str, ...] = ('S_O2', 'S_F', 'S_A', 'S_NH4', 'S_NO3', 'S_PO4', 'S_I', 'S_N2', 'X_I', 'X_S', 'X_H', 'X_PAO', 'X_PP', 'X_PHA', 'X_AUT', 'X_MeOH', 'X_MeP', 'X_TSS')
ASM2D_COD_COMPONENTS: tuple[str, ...] = ('S_F', 'S_A', 'S_I', 'X_S', 'X_I', 'X_H')
ASM2D_BIODEGRADABLE_COD_COMPONENTS: tuple[str, ...] = ('S_F', 'S_A', 'X_S', 'X_H')
ASM2D_INERT_COD_COMPONENTS: tuple[str, ...] = ('S_I', 'X_I')
BOD5_DECAY_K_PER_D_DEFAULT = 0.23
BOD_YIELD_LOSS_DEFAULT = 0.15
MIN_PLAUSIBLE_INERT_COD_FRACTION = 0.05
_ASM2D_PROFILES: dict[InfluentCharacter, dict[str, float]] = {InfluentCharacter.MUNICIPAL: {'S_F': 0.18, 'S_A': 0.08, 'S_I': 0.1, 'X_S': 0.43, 'X_I': 0.17, 'X_H': 0.04}, InfluentCharacter.INDUSTRIAL_CETP: {'S_F': 0.22, 'S_A': 0.09, 'S_I': 0.12, 'X_S': 0.35, 'X_I': 0.2, 'X_H': 0.02}, InfluentCharacter.DAIRY: {'S_F': 0.28, 'S_A': 0.16, 'S_I': 0.06, 'X_S': 0.38, 'X_I': 0.08, 'X_H': 0.04}, InfluentCharacter.TEXTILE: {'S_F': 0.12, 'S_A': 0.04, 'S_I': 0.36, 'X_S': 0.18, 'X_I': 0.28, 'X_H': 0.02}, InfluentCharacter.PHARMA: {'S_F': 0.15, 'S_A': 0.07, 'S_I': 0.34, 'X_S': 0.22, 'X_I': 0.2, 'X_H': 0.02}, InfluentCharacter.REFINERY_CHEM: {'S_F': 0.07, 'S_A': 0.05, 'S_I': 0.56, 'X_S': 0.12, 'X_I': 0.18, 'X_H': 0.02}, InfluentCharacter.MIXED_INDUSTRIAL: {'S_F': 0.17, 'S_A': 0.08, 'S_I': 0.3, 'X_S': 0.23, 'X_I': 0.2, 'X_H': 0.02}}
ADM1_COMPONENTS: tuple[str, ...] = ('S_su', 'S_aa', 'S_fa', 'S_va', 'S_bu', 'S_pro', 'S_ac', 'S_h2', 'S_ch4', 'S_IC', 'S_IN', 'S_IP', 'S_I', 'X_ch', 'X_pr', 'X_li', 'X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2', 'X_I')
ADM1_COD_COMPONENTS: tuple[str, ...] = ('S_su', 'S_aa', 'S_fa', 'S_va', 'S_bu', 'S_pro', 'S_ac', 'S_I', 'X_ch', 'X_pr', 'X_li', 'X_I', 'X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2')
ADM1_I_N_PROTEIN = 0.11065
ADM1_I_N_BIOMASS = 0.08615
ADM1_I_N_INERT = 0.06
ADM1_I_P_LIPID = 0.01067
ADM1_I_P_BIOMASS = 0.02
ADM1_I_P_INERT = 0.01

def coerce_tender_influent_basis(value: TenderInfluentBasis | dict[str, Any], *, allow_unmeasured_characteristics: bool=False) -> TenderInfluentBasis:
    """Return a TenderInfluentBasis from a model instance or dict-like artifact.

    ``allow_unmeasured_characteristics`` is for a stream that is DERIVED by
    construction -- a digester feed computed by upstream sizing, say -- where
    demanding a measured TKN would be demanding a measurement of something that
    does not physically exist to sample. The guard belongs at the project
    boundary, on the raw influent that comes from the project documentation, and
    that is where it stays on by default.

    Do not set it to silence the boundary check. The whole point is that a value
    filled from a generic ratio is invisible downstream.
    """
    if isinstance(value, TenderInfluentBasis):
        return value
    if allow_unmeasured_characteristics and isinstance(value, dict):
        value = dict(value)
        value.setdefault('allow_unmeasured_characteristics', True)
    return TenderInfluentBasis.model_validate(value)

def asm2d_cod_kg_m3(species: dict[str, float]) -> float:
    """Return ASM2D COD represented by a species vector in kg COD/m3."""
    total = sum((float(species.get(component, 0.0)) for component in ASM2D_COD_COMPONENTS))
    total += ASM2D_F_CV_AUT * float(species.get('X_AUT', 0.0))
    return total

def asm2d_tkn_kg_m3(species: dict[str, float]) -> float:
    """Return ASM2D reduced (Kjeldahl) nitrogen represented by a species vector.

    Reduced N = soluble ammonia (S_NH4) + dissolved-N2 closure slack (S_N2,
    ~0 at the plant boundary) + organic-bound N in biomass and particulates.
    Oxidized nitrate (S_NO3) is excluded — it is the explicit NO3 term.

    CLOSURE FUNCTION — NOT AN EFFLUENT-COMPLIANCE NUMBER. The S_N2 term makes
    nitrogen close across a denitrifying flowsheet (N leaving as N2 gas has to
    land somewhere or N appears to vanish), and S_N2 is ~0 at the *influent*
    boundary where this is used. At the *effluent* of a denitrifying plant
    S_N2 accumulates to ~15-25 mg/L, so this function will over-report
    effluent TKN/TN by that much — dissolved N2 gas is not regulatory
    TKN/TN. For an effluent vector use `asm2d_effluent_tn_kg_m3`.
    """
    soluble = float(species.get('S_NH4', 0.0)) + float(species.get('S_N2', 0.0))
    biomass = ASM2D_I_NBM * (float(species.get('X_H', 0.0)) + float(species.get('X_PAO', 0.0)) + float(species.get('X_AUT', 0.0)))
    particulate = ASM2D_I_NXI * float(species.get('X_I', 0.0))
    particulate += ASM2D_I_NXS * float(species.get('X_S', 0.0))
    return soluble + biomass + particulate

def asm2d_no3_kg_m3(species: dict[str, float]) -> float:
    """Return ASM2D oxidized (nitrate) nitrogen represented by a species vector."""
    return float(species.get('S_NO3', 0.0))

def asm2d_tn_kg_m3(species: dict[str, float]) -> float:
    """Return ASM2D total nitrogen = reduced (TKN) + oxidized (NO3).

    Carries the S_N2 closure term from `asm2d_tkn_kg_m3`; see the warning
    there. Use `asm2d_effluent_tn_kg_m3` for effluent-compliance reporting.
    """
    return asm2d_tkn_kg_m3(species) + asm2d_no3_kg_m3(species)

def asm2d_effluent_tn_kg_m3(species: dict[str, float]) -> float:
    """Return regulatory total nitrogen for an ASM2D *effluent* vector.

    TN = NH4-N + NO3-N + organic-bound N (biomass + particulate). Dissolved
    N2 gas (S_N2) is EXCLUDED: it is the product of denitrification, it
    off-gasses, and no TN method measures it. This is the canonical
    effluent-facing counterpart to the `asm2d_tn_kg_m3` closure function and
    matches the nitrogen the watertap-engine-mcp stream payload reports.
    """
    return float(species.get('S_NH4', 0.0)) + float(species.get('S_NO3', 0.0)) + ASM2D_I_NBM * (float(species.get('X_H', 0.0)) + float(species.get('X_PAO', 0.0)) + float(species.get('X_AUT', 0.0))) + ASM2D_I_NXI * float(species.get('X_I', 0.0)) + ASM2D_I_NXS * float(species.get('X_S', 0.0))

def asm2d_tp_kg_m3(species: dict[str, float]) -> float:
    """Return ASM2D total phosphorus represented by a species vector.

    Includes soluble orthophosphate (S_PO4), intracellular polyphosphate
    (X_PP), biomass/particulate-bound organic P, AND the phosphorus locked in
    precipitated metal-phosphate (X_MeP). ``X_PP`` and ``S_PO4`` are already on
    a phosphorus mass basis in ASM2d. The X_MeP term is what makes total P
    close across the chemical-P precipitation unit: P leaving S_PO4 re-appears
    as X_MeP-bound P rather than vanishing.
    """
    soluble = float(species.get('S_PO4', 0.0))
    polyphosphate = float(species.get('X_PP', 0.0))
    biomass = ASM2D_I_PBM * (float(species.get('X_H', 0.0)) + float(species.get('X_PAO', 0.0)) + float(species.get('X_AUT', 0.0)))
    particulate = ASM2D_I_PXI * float(species.get('X_I', 0.0))
    particulate += ASM2D_I_PXS * float(species.get('X_S', 0.0))
    precipitated = ASM2D_P_PER_XMEP * float(species.get('X_MeP', 0.0))
    return soluble + polyphosphate + biomass + particulate + precipitated

def _relative_error(actual: float, expected: float) -> float:
    if abs(expected) <= 1e-12:
        return abs(actual - expected)
    return abs(actual - expected) / abs(expected)

def asm2d_mass_balance_errors(species: dict[str, float], influent: TenderInfluentBasis) -> dict[str, float]:
    """Return relative COD / TKN / NO3 / TP closure errors for an ASM2D vector.

    TKN (reduced N) and NO3 (oxidized N) are validated separately so that the
    total-N basis = TKN + NO3 is not double-counted.
    """
    return {'COD': _relative_error(asm2d_cod_kg_m3(species), influent.cod_kg_m3), 'TKN': _relative_error(asm2d_tkn_kg_m3(species), influent.tkn_kg_m3), 'NO3': _relative_error(asm2d_no3_kg_m3(species), influent.no3_kg_m3), 'TP': _relative_error(asm2d_tp_kg_m3(species), influent.tp_kg_m3)}

def validate_asm2d_mass_closure(species: dict[str, float], influent: TenderInfluentBasis, *, tolerance: float=0.01) -> None:
    """Raise ValueError if ASM2D COD / TKN / NO3 / TP closure exceeds tolerance."""
    errors = asm2d_mass_balance_errors(species, influent)
    failed = {key: value for key, value in errors.items() if value > tolerance}
    if failed:
        details = ', '.join((f'{key}={value:.3%}' for key, value in sorted(failed.items())))
        raise ValueError(f'ASM2D influent fractionation does not close within {tolerance:.1%}: {details}')

def _normalize_profile(profile: dict[str, float]) -> dict[str, float]:
    total = sum((max(float(value), 0.0) for value in profile.values()))
    if total <= 0:
        raise ValueError('fractionation profile must contain a positive fraction sum')
    return {key: max(float(value), 0.0) / total for key, value in profile.items()}

def _resolve_cod_fraction_overrides(overrides: dict[str, float]) -> dict[str, float]:
    """Validate and normalize a caller-supplied COD fractionation override.

    The override REPLACES the influent-character profile. It must cover every
    ASM2D COD component (S_F / S_A / S_I / X_S / X_I / X_H), carry only those
    components, and be non-negative. The result is normalized to sum to 1.
    """
    if not isinstance(overrides, dict):
        raise TypeError('cod_fraction_overrides must be a mapping of COD component -> fraction')
    unknown = sorted(set(overrides) - set(ASM2D_COD_COMPONENTS))
    if unknown:
        raise ValueError(f'cod_fraction_overrides has unsupported components {unknown}; allowed: {list(ASM2D_COD_COMPONENTS)}')
    missing = [component for component in ASM2D_COD_COMPONENTS if component not in overrides]
    if missing:
        raise ValueError(f'cod_fraction_overrides must cover all ASM2D COD components; missing {missing}')
    cleaned: dict[str, float] = {}
    for component in ASM2D_COD_COMPONENTS:
        value = float(overrides[component])
        if not math.isfinite(value):
            raise ValueError(f'cod_fraction_overrides[{component!r}] must be finite, got {value!r}')
        if value < 0:
            raise ValueError(f'cod_fraction_overrides[{component!r}] must be non-negative')
        cleaned[component] = value
    return _normalize_profile(cleaned)

def ultimate_bod_from_bod5(bod5: float, *, decay_k_per_d: float=BOD5_DECAY_K_PER_D_DEFAULT) -> float:
    """Return ultimate carbonaceous BOD from a 5-day BOD, in the input's units.

    BOD_u = BOD5 / (1 - exp(-5k)) for first-order oxygen exertion at rate
    ``decay_k_per_d`` (base e, per day).
    """
    if not math.isfinite(bod5) or bod5 < 0:
        raise ValueError(f'bod5 must be a non-negative finite value, got {bod5!r}')
    if not math.isfinite(decay_k_per_d) or decay_k_per_d <= 0:
        raise ValueError(f'decay_k_per_d must be a positive finite value, got {decay_k_per_d!r}')
    exerted_fraction = 1.0 - math.exp(-5.0 * decay_k_per_d)
    return bod5 / exerted_fraction

def biodegradable_cod_from_bod5(bod5: float, *, decay_k_per_d: float=BOD5_DECAY_K_PER_D_DEFAULT, bod_yield_loss: float=BOD_YIELD_LOSS_DEFAULT) -> float:
    """Return biodegradable COD implied by a measured BOD5, in the input's units.

    Two corrections stand between BOD5 and biodegradable COD:
      1. BOD5 -> BOD_ultimate: only 1 - exp(-5k) of the demand is exerted by
         day 5 (~68% at the default k = 0.23 d-1).
      2. BOD_ultimate -> biodegradable COD: a fraction ``bod_yield_loss``
         (f_BOD, STOWA default 0.15) of the substrate becomes cell mass and
         endogenous residue rather than exerted oxygen demand.
    """
    if not math.isfinite(bod_yield_loss) or not 0.0 <= bod_yield_loss < 1.0:
        raise ValueError(f'bod_yield_loss must be in [0, 1), got {bod_yield_loss!r}')
    bod_ultimate = ultimate_bod_from_bod5(bod5, decay_k_per_d=decay_k_per_d)
    return bod_ultimate / (1.0 - bod_yield_loss)

def _rebalance_profile_to_measured_bod5(profile: dict[str, float], basis: TenderInfluentBasis, *, allocatable_cod: float, x_aut: float, min_inert_fraction: float=MIN_PLAUSIBLE_INERT_COD_FRACTION, bod_yield_loss: float=BOD_YIELD_LOSS_DEFAULT) -> tuple[dict[str, float], dict[str, Any]]:
    """Constrain the biodegradable/inert COD split with the measured BOD5.

    The character profile still sets the shape *within* each group (how the
    biodegradable COD divides across S_F / S_A / X_S / X_H, and how the inert
    COD divides across S_I / X_I). The measured BOD5 sets the size of the two
    groups, replacing the profile's implicit COD/BOD ratio.

    Returns ``(profile, basis_info)``. The profile is returned unchanged, with
    ``basis_info["honored"] is False`` and a warning, whenever the measurement
    cannot steer the split: no measured BOD5, or a COD/BOD5 pair whose implied
    inert COD falls below ``min_inert_fraction``.
    """
    info: dict[str, Any] = {'source': 'character_profile', 'influent_character': basis.influent_character.value, 'honored': False, 'bod5_basis': basis.bod5_basis, 'bod5_mg_L': float(basis.BOD5) if basis.BOD5 is not None else None, 'cod_mg_L': float(basis.COD), 'profile_inert_cod_fraction': round(sum((profile.get(name, 0.0) for name in ASM2D_INERT_COD_COMPONENTS)), 6), 'warnings': []}
    if basis.bod5_basis != 'measured' or basis.BOD5 is None:
        info['reason'] = "no measured BOD5 on the influent basis; the character profile's implicit COD/BOD ratio sets the biodegradable/inert split"
        return (profile, info)
    decay_k = float(basis.BOD_decay_k_per_d or BOD5_DECAY_K_PER_D_DEFAULT)
    cod_bio_mg_l = biodegradable_cod_from_bod5(float(basis.BOD5), decay_k_per_d=decay_k, bod_yield_loss=bod_yield_loss)
    info.update({'bod_decay_k_per_d': decay_k, 'bod_yield_loss': bod_yield_loss, 'measured_cod_bod5_ratio': round(float(basis.COD) / float(basis.BOD5), 4) if float(basis.BOD5) > 0 else None, 'biodegradable_cod_from_bod5_mg_L': round(cod_bio_mg_l, 4)})
    bio_target = max(cod_bio_mg_l / 1000.0 - ASM2D_F_CV_AUT * x_aut, 0.0)
    if allocatable_cod <= 0:
        info['reason'] = 'no allocatable COD'
        return (profile, info)
    bio_fraction = bio_target / allocatable_cod
    implied_inert_fraction = 1.0 - bio_fraction
    info['implied_inert_cod_fraction'] = round(implied_inert_fraction, 6)
    if implied_inert_fraction < min_inert_fraction:
        message = f'measured BOD5 {basis.BOD5:.4g} mg/L with COD {basis.COD:.4g} mg/L (COD/BOD5 = {float(basis.COD) / float(basis.BOD5):.2f}) implies {implied_inert_fraction:.1%} inert COD, below the {min_inert_fraction:.0%} physical floor: the pair is internally inconsistent (COD and BOD5 are likely unpaired design maxima, not a paired sample). Falling back to the {basis.influent_character.value} character profile ({info['profile_inert_cod_fraction']:.1%} inert COD). Resolve the influent basis with a paired COD/BOD5 sample before relying on the inert-COD-driven effluent COD floor.'
        info['reason'] = 'bod5_implies_inert_cod_below_physical_floor'
        info['warnings'].append(message)
        return (profile, info)
    bio_profile_sum = sum((profile.get(name, 0.0) for name in ASM2D_BIODEGRADABLE_COD_COMPONENTS))
    inert_profile_sum = sum((profile.get(name, 0.0) for name in ASM2D_INERT_COD_COMPONENTS))
    if bio_profile_sum <= 0 or inert_profile_sum <= 0:
        info['reason'] = 'character profile has an empty biodegradable or inert group; cannot redistribute within groups'
        return (profile, info)
    rebalanced = dict(profile)
    for name in ASM2D_BIODEGRADABLE_COD_COMPONENTS:
        rebalanced[name] = profile.get(name, 0.0) / bio_profile_sum * bio_fraction
    for name in ASM2D_INERT_COD_COMPONENTS:
        rebalanced[name] = profile.get(name, 0.0) / inert_profile_sum * implied_inert_fraction
    info.update({'source': 'measured_bod5', 'honored': True, 'reason': f'measured BOD5 constrains the biodegradable/inert COD split; the {basis.influent_character.value} profile sets the shape within each group', 'profile_implied_bod5_mg_L': round(_profile_implied_bod5_mg_l(profile, allocatable_cod=allocatable_cod, x_aut=x_aut, decay_k_per_d=decay_k, bod_yield_loss=bod_yield_loss), 4)})
    return (_normalize_profile(rebalanced), info)

def _profile_implied_bod5_mg_l(profile: dict[str, float], *, allocatable_cod: float, x_aut: float, decay_k_per_d: float=BOD5_DECAY_K_PER_D_DEFAULT, bod_yield_loss: float=BOD_YIELD_LOSS_DEFAULT) -> float:
    """Return the BOD5 a COD profile implies — the number BOD5 replaces."""
    bio_cod = allocatable_cod * sum((profile.get(name, 0.0) for name in ASM2D_BIODEGRADABLE_COD_COMPONENTS))
    bio_cod += ASM2D_F_CV_AUT * x_aut
    bod_ultimate = bio_cod * 1000.0 * (1.0 - bod_yield_loss)
    return bod_ultimate * (1.0 - math.exp(-5.0 * decay_k_per_d))

def _bound_nutrients_asm2d(species: dict[str, float]) -> tuple[float, float]:
    species_without_soluble = dict(species)
    for key in ('S_NH4', 'S_NO3', 'S_N2', 'S_PO4', 'X_MeP'):
        species_without_soluble[key] = 0.0
    return (asm2d_tn_kg_m3(species_without_soluble), asm2d_tp_kg_m3(species_without_soluble))

def _fit_asm2d_bound_nutrients(species: dict[str, float], influent: TenderInfluentBasis, profile: dict[str, float]) -> None:
    """Move particulate COD to soluble COD if organic N/P would exceed basis."""
    bound_n, bound_p = _bound_nutrients_asm2d(species)
    ratios = [1.0]
    if influent.tkn_kg_m3 > 0 and bound_n > 0:
        ratios.append(min(1.0, 0.85 * influent.tkn_kg_m3 / bound_n))
    if influent.tp_kg_m3 > 0 and bound_p > 0:
        ratios.append(min(1.0, 0.85 * influent.tp_kg_m3 / bound_p))
    scale = min(ratios)
    if scale >= 0.999:
        return
    released = 0.0
    for component in ('X_H', 'X_S', 'X_I'):
        old = species.get(component, 0.0)
        new = old * scale
        species[component] = new
        released += old - new
    soluble_total = sum((profile.get(component, 0.0) for component in ('S_F', 'S_A', 'S_I')))
    if soluble_total <= 0:
        species['S_F'] = species.get('S_F', 0.0) + released
        return
    for component in ('S_F', 'S_A', 'S_I'):
        species[component] = species.get(component, 0.0) + released * profile.get(component, 0.0) / soluble_total

@dataclass(frozen=True)
class Asm2dFractionation:
    """An ASM2D influent vector plus the provenance of its COD split."""
    species: dict[str, float]
    cod_split_basis: dict[str, Any] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()

def fractionate_to_asm2d_with_provenance(influent: TenderInfluentBasis, *, cod_fraction_overrides: dict[str, float] | None=None, honor_measured_bod5: bool=True, min_inert_cod_fraction: float=MIN_PLAUSIBLE_INERT_COD_FRACTION, bod_yield_loss: float=BOD_YIELD_LOSS_DEFAULT) -> Asm2dFractionation:
    """Return the ASM2D influent vector and the provenance of its COD split.

    See `fractionate_to_asm2d` for the species semantics. This variant also
    reports which basis set the biodegradable/inert split and any warning raised
    while deciding, so a caller can surface a re-baseline instead of absorbing it.
    """
    if not 0.0 <= min_inert_cod_fraction < 1.0:
        raise ValueError(f'min_inert_cod_fraction must be in [0, 1), got {min_inert_cod_fraction!r}')
    basis = coerce_tender_influent_basis(influent)
    cod = basis.cod_kg_m3
    x_aut = min(0.001, max(cod * 0.005, 1e-06))
    allocatable_cod = max(cod - ASM2D_F_CV_AUT * x_aut, cod * 0.95)
    if cod_fraction_overrides is not None:
        profile = _resolve_cod_fraction_overrides(cod_fraction_overrides)
        split_basis: dict[str, Any] = {'source': 'explicit_cod_fraction_overrides', 'honored': False, 'reason': 'cod_fraction_overrides replaces the COD split verbatim; measured BOD5 is not applied on top of an explicit override', 'warnings': []}
    else:
        profile = _normalize_profile(_ASM2D_PROFILES[basis.influent_character])
        if honor_measured_bod5:
            profile, split_basis = _rebalance_profile_to_measured_bod5(profile, basis, allocatable_cod=allocatable_cod, x_aut=x_aut, min_inert_fraction=min_inert_cod_fraction, bod_yield_loss=bod_yield_loss)
        else:
            split_basis = {'source': 'character_profile', 'honored': False, 'reason': 'honor_measured_bod5 disabled by caller', 'warnings': []}
    species = _allocate_asm2d_species(basis, profile, allocatable_cod=allocatable_cod, x_aut=x_aut)
    split_basis['resulting_cod_fractions'] = {name: round(profile.get(name, 0.0), 6) for name in ASM2D_COD_COMPONENTS}
    return Asm2dFractionation(species=species, cod_split_basis=split_basis, warnings=tuple(split_basis.get('warnings') or ()))

def fractionate_to_asm2d(influent: TenderInfluentBasis, *, cod_fraction_overrides: dict[str, float] | None=None, honor_measured_bod5: bool=True, min_inert_cod_fraction: float=MIN_PLAUSIBLE_INERT_COD_FRACTION, bod_yield_loss: float=BOD_YIELD_LOSS_DEFAULT) -> dict[str, float]:
    """Return ASM2D influent species in kg/m3 with COD / TKN / NO3 / TP closure.

    The COD split is decided in this precedence order:
      1. ``cod_fraction_overrides`` — REPLACES the split verbatim (must cover
         S_F / S_A / S_I / X_S / X_I / X_H, be non-negative, normalized to 1).
      2. A MEASURED BOD5 on the basis — sizes the biodegradable (S_F + S_A +
         X_S + X_H) vs inert (S_I + X_I) groups, with the influent-character
         profile setting the shape within each group. Measured plant data
         outranks a generic profile's implicit COD/BOD ratio.
      3. The influent-character profile — the fallback when BOD5 is absent,
         synthesised from COD, a generic class default, or physically
         irreconcilable with the COD (see MIN_PLAUSIBLE_INERT_COD_FRACTION).

    Reduced nitrogen (S_NH4 + organic-bound N) is closed to the TKN basis;
    oxidized nitrogen S_NO3 is set explicitly from the influent NO3-N (so total
    N = TKN + NO3 is not double-counted). Phosphorus is closed via soluble S_PO4
    after biomass/particulate-bound P. Autotrophs are seeded at ~1 g/m3 when the
    COD basis can support it, avoiding the high-SRT washout behavior caused by
    near-zero municipal defaults.

    Use `fractionate_to_asm2d_with_provenance` when the caller needs to report
    which basis set the split.
    """
    return fractionate_to_asm2d_with_provenance(influent, cod_fraction_overrides=cod_fraction_overrides, honor_measured_bod5=honor_measured_bod5, min_inert_cod_fraction=min_inert_cod_fraction, bod_yield_loss=bod_yield_loss).species

def _allocate_asm2d_species(basis: TenderInfluentBasis, profile: dict[str, float], *, allocatable_cod: float, x_aut: float) -> dict[str, float]:
    """Allocate a resolved COD profile into a closed ASM2D species vector."""
    s_no3 = max(basis.no3_kg_m3, 1e-09)
    species = {component: 0.0 for component in ASM2D_COMPONENTS}
    species.update({'S_O2': 1e-09, 'S_NO3': s_no3, 'S_N2': 1e-09, 'X_PAO': 1e-09, 'X_PP': 1e-09, 'X_PHA': 1e-09, 'X_AUT': x_aut, 'X_MeOH': 1e-09, 'X_MeP': 1e-09, 'X_TSS': float(basis.TSS or 0.0) / 1000.0})
    for component, fraction in profile.items():
        species[component] = allocatable_cod * fraction
    _fit_asm2d_bound_nutrients(species, basis, profile)
    bound_n, bound_p = _bound_nutrients_asm2d(species)
    soluble_reduced_n = max(basis.tkn_kg_m3 - bound_n - species['S_N2'], 1e-09)
    species['S_NH4'] = soluble_reduced_n
    soluble_p = max(basis.tp_kg_m3 - bound_p, 1e-09)
    species['S_PO4'] = soluble_p
    validate_asm2d_mass_closure(species, basis)
    return {key: float(value) for key, value in species.items()}

def adm1_cod_kg_m3(species: dict[str, float]) -> float:
    return sum((float(species.get(component, 0.0)) for component in ADM1_COD_COMPONENTS))

def adm1_tn_kg_m3(species: dict[str, float]) -> float:
    protein_n = ADM1_I_N_PROTEIN * (float(species.get('S_aa', 0.0)) + float(species.get('X_pr', 0.0)))
    biomass_n = ADM1_I_N_BIOMASS * sum((float(species.get(component, 0.0)) for component in ('X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2')))
    inert_n = ADM1_I_N_INERT * float(species.get('X_I', 0.0))
    return float(species.get('S_IN', 0.0)) + protein_n + biomass_n + inert_n

def adm1_tp_kg_m3(species: dict[str, float]) -> float:
    lipid_p = ADM1_I_P_LIPID * float(species.get('X_li', 0.0))
    biomass_p = ADM1_I_P_BIOMASS * sum((float(species.get(component, 0.0)) for component in ('X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2')))
    inert_p = ADM1_I_P_INERT * float(species.get('X_I', 0.0))
    return float(species.get('S_IP', 0.0)) + lipid_p + biomass_p + inert_p

def adm1_mass_balance_errors(species: dict[str, float], influent: TenderInfluentBasis) -> dict[str, float]:
    return {'COD': _relative_error(adm1_cod_kg_m3(species), influent.cod_kg_m3), 'TN': _relative_error(adm1_tn_kg_m3(species), influent.tkn_kg_m3), 'TP': _relative_error(adm1_tp_kg_m3(species), influent.tp_kg_m3)}

def validate_adm1_mass_closure(species: dict[str, float], influent: TenderInfluentBasis, *, tolerance: float=0.01) -> None:
    errors = adm1_mass_balance_errors(species, influent)
    failed = {key: value for key, value in errors.items() if value > tolerance}
    if failed:
        details = ', '.join((f'{key}={value:.3%}' for key, value in sorted(failed.items())))
        raise ValueError(f'ADM1 influent fractionation does not close within {tolerance:.1%}: {details}')

def fractionate_to_adm1(influent: TenderInfluentBasis, *, honor_measured_bod5: bool=True) -> dict[str, float]:
    """Return ADM1 influent species in kg/m3 with COD/TN/TP closure.

    The ADM1 inert pools inherit the ASM2D biodegradable/inert COD split, so a
    measured BOD5 constrains this vector for the same reason it constrains the
    ASM2D one. See `fractionate_to_asm2d` for the precedence rules.
    """
    basis = coerce_tender_influent_basis(influent)
    asm = fractionate_to_asm2d(basis, honor_measured_bod5=honor_measured_bod5)
    species = {component: 0.0 for component in ADM1_COMPONENTS}
    species.update({'S_su': asm['S_F'] * 0.55, 'S_aa': asm['S_F'] * 0.25, 'S_ac': asm['S_A'], 'S_pro': asm['S_F'] * 0.1, 'S_bu': asm['S_F'] * 0.05, 'S_va': asm['S_F'] * 0.05, 'S_I': asm['S_I'], 'X_ch': asm['X_S'] * 0.45, 'X_pr': asm['X_S'] * 0.3, 'X_li': asm['X_S'] * 0.25, 'X_I': asm['X_I'] + asm['X_H'] + ASM2D_F_CV_AUT * asm['X_AUT'], 'S_IC': basis.alkalinity_mol_per_m3 * 0.012})
    bound_n = adm1_tn_kg_m3({**species, 'S_IN': 0.0})
    species['S_IN'] = max(basis.tkn_kg_m3 - bound_n, 1e-09)
    bound_p = adm1_tp_kg_m3({**species, 'S_IP': 0.0})
    species['S_IP'] = max(basis.tp_kg_m3 - bound_p, 1e-09)
    validate_adm1_mass_closure(species, basis)
    return {key: float(value) for key, value in species.items()}
