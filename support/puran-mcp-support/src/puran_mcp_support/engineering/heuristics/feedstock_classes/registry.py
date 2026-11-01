"""Feedstock classifications. Supply the complete input basis; no dataset presets ship."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any
from puran_mcp_support.engineering.heuristics.asm2d_influent_fractionation import BOD5_INPUT_ALIASES, InfluentCharacter, TenderInfluentBasis, adm1_mass_balance_errors, asm2d_mass_balance_errors, fractionate_to_adm1, fractionate_to_asm2d_with_provenance, validate_adm1_mass_closure, validate_asm2d_mass_closure
FRACTIONATION_CONTROL_KEYS: frozenset[str] = frozenset({'cod_fraction_overrides'})

def _split_fractionation_controls(overrides: dict[str, Any] | None) -> tuple[dict[str, Any], dict[str, Any]]:
    """Split caller overrides into (bulk-basis overrides, fractionation controls)."""
    if not overrides:
        return ({}, {})
    basis_overrides = {key: value for key, value in overrides.items() if key not in FRACTIONATION_CONTROL_KEYS}
    controls = {key: value for key, value in overrides.items() if key in FRACTIONATION_CONTROL_KEYS}
    return (basis_overrides, controls)

@dataclass(frozen=True)
class FeedstockFractionation:
    """Simulator-ready fractionation vectors for one influent basis."""
    basis: TenderInfluentBasis
    asm2d_state: dict[str, float]
    adm1_state: dict[str, float]
    non_biological_state: dict[str, Any]
    closure_errors: dict[str, float]
    cod_split_basis: dict[str, Any] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()

@dataclass(frozen=True)
class FeedstockClass:
    """Feedstock-class defaults and fractionation behavior."""
    name: str
    default_basis: dict[str, Any]
    influent_character: InfluentCharacter
    compatible_template_families: tuple[str, ...]
    expected_range: dict[str, tuple[float, float]]
    non_biological_defaults: dict[str, Any] = field(default_factory=dict)
    aliases: tuple[str, ...] = ()

    def basis(self, *, flow_m3_d: float, influent_basis: dict[str, Any] | TenderInfluentBasis | None=None, overrides: dict[str, Any] | None=None) -> TenderInfluentBasis:
        """Return a canonical influent basis with class defaults and overrides."""
        if isinstance(influent_basis, TenderInfluentBasis):
            data = influent_basis.model_dump(mode='python', exclude_unset=True)
        elif isinstance(influent_basis, dict):
            data = dict(influent_basis)
        else:
            data = dict(self.default_basis)
            data.setdefault('bod5_basis', 'class_default')
        data.setdefault('influent_character', self.influent_character.value)
        data['flow_m3_d'] = flow_m3_d
        basis_overrides, _controls = _split_fractionation_controls(overrides)
        if basis_overrides:
            data.update(basis_overrides)
            data['flow_m3_d'] = flow_m3_d
            if any((name in basis_overrides for name in BOD5_INPUT_ALIASES)):
                data['bod5_basis'] = 'measured'
        data.setdefault('allow_unmeasured_characteristics', True)
        return TenderInfluentBasis.model_validate(data)

    def fractionate(self, *, flow_m3_d: float, influent_basis: dict[str, Any] | TenderInfluentBasis | None=None, overrides: dict[str, Any] | None=None, cod_fraction_overrides: dict[str, float] | None=None, include_biological: bool=True) -> FeedstockFractionation:
        """Return ASM2D, ADM1, and non-biological feed states.

        ``cod_fraction_overrides`` REPLACES the influent-character COD profile
        (S_F / S_A / S_I / X_S / X_I / X_H) and is validated + normalized by
        ``fractionate_to_asm2d``. It may be passed explicitly or nested inside
        ``overrides`` (the flat dict the ``solve_flowsheet`` MCP surface hands
        through as ``feedstock_overrides``); the explicit kwarg wins.
        """
        basis_overrides, controls = _split_fractionation_controls(overrides)
        if cod_fraction_overrides is None:
            cod_fraction_overrides = controls.get('cod_fraction_overrides')
        basis = self.basis(flow_m3_d=flow_m3_d, influent_basis=influent_basis, overrides=basis_overrides)
        asm2d_result = None
        asm2d, adm1 = ({}, {})
        if include_biological:
            asm2d_result = fractionate_to_asm2d_with_provenance(basis, cod_fraction_overrides=cod_fraction_overrides)
            asm2d = asm2d_result.species
            adm1 = fractionate_to_adm1(basis)
            validate_asm2d_mass_closure(asm2d, basis)
            validate_adm1_mass_closure(adm1, basis)
        measured_non_bio = {'TDS_mg_L': basis.TDS, 'hardness_mg_L_as_CaCO3': basis.hardness_as_CaCO3, 'chloride_mg_L': basis.chloride, 'sulphate_mg_L': basis.sulphate, 'silica_mg_L_as_SiO2': basis.silica_as_SiO2, 'calcium_mg_L': basis.calcium, 'magnesium_mg_L': basis.magnesium, 'sodium_mg_L': basis.sodium}
        non_bio: dict[str, Any] = {}
        class_defaulted: list[str] = []
        for key, default_value in self.non_biological_defaults.items():
            if measured_non_bio.get(key) is None:
                non_bio[key] = default_value
                class_defaulted.append(key)
        non_bio.update({k: v for k, v in measured_non_bio.items() if v is not None})
        non_bio.update({'flow_m3_d': basis.flow_m3_d, 'COD_mg_L': basis.COD, 'BOD5_mg_L': basis.BOD5, 'TKN_mg_L': basis.TKN, 'TP_mg_L': basis.TP, 'TSS_mg_L': basis.TSS, 'temperature_K': basis.temperature_K, 'alkalinity_mol_per_m3': basis.alkalinity_mol_per_m3})
        if class_defaulted:
            non_bio['class_defaulted_fields'] = sorted(class_defaulted)
            non_bio['class_defaulted_from'] = self.name
        return FeedstockFractionation(basis=basis, asm2d_state=asm2d, adm1_state=adm1, non_biological_state=non_bio, closure_errors={'asm2d': asm2d_mass_balance_errors(asm2d, basis), 'adm1': adm1_mass_balance_errors(adm1, basis)} if include_biological else {}, cod_split_basis=asm2d_result.cod_split_basis if asm2d_result else {'source': 'not_applicable_non_biological'}, warnings=asm2d_result.warnings if asm2d_result else ('Non-biological feed handoff: ASM2D/ADM1 vectors and biological closure were not calculated.',))

    def closure_validator(self, fractionation: FeedstockFractionation) -> None:
        """Assert COD/TN/TP closure for the biological vectors."""
        validate_asm2d_mass_closure(fractionation.asm2d_state, fractionation.basis)
        validate_adm1_mass_closure(fractionation.adm1_state, fractionation.basis)

def _basis(*, COD: float, BOD5: float | None=None, TKN: float | None=None, NH4_N: float | None=None, TP: float | None=None, TSS: float | None=None, TDS: float | None=None, FOG: float | None=None, character: InfluentCharacter) -> dict[str, Any]:
    data: dict[str, Any] = {'flow_m3_d': 1.0, 'COD': COD, 'BOD5': BOD5, 'TKN': TKN, 'NH4_N': NH4_N, 'TP': TP, 'TSS': TSS, 'influent_character': character.value}
    if TDS is not None:
        data['TDS'] = TDS
    if FOG is not None:
        data['FOG'] = FOG
    return {key: value for key, value in data.items() if value is not None}

def _ranges(cod_max: float, tn_max: float, tp_max: float, mlss: tuple[float, float]) -> dict[str, tuple[float, float]]:
    return {'effluent_COD_mg_L': (0.0, cod_max), 'effluent_TN_mg_L': (0.0, tn_max), 'effluent_TP_mg_L': (0.0, tp_max), 'MLSS_basin_mg_L': mlss, 'stream_closure_pct': (99.0, 101.0)}
BIOLOGICAL = ('asm2d', 'asm1', 'adm1', 'anaerobic')
SALTS = ('ro', 'evaporator', 'ix', 'softening', 'zld')
FEEDSTOCK_CLASSES: dict[str, FeedstockClass] = {'municipal_low_strength': FeedstockClass(name='municipal_low_strength', default_basis={}, influent_character=InfluentCharacter.MUNICIPAL, compatible_template_families=BIOLOGICAL, expected_range=_ranges(60, 15, 3, (3000, 10000)), aliases=('municipal',)), 'municipal_medium_strength': FeedstockClass(name='municipal_medium_strength', default_basis={}, influent_character=InfluentCharacter.MUNICIPAL, compatible_template_families=BIOLOGICAL, expected_range=_ranges(60, 15, 3, (3000, 10000))), 'municipal_high_strength': FeedstockClass(name='municipal_high_strength', default_basis={}, influent_character=InfluentCharacter.MUNICIPAL, compatible_template_families=BIOLOGICAL, expected_range=_ranges(80, 20, 4, (4000, 12000))), 'industrial_cetp_mixed': FeedstockClass(name='industrial_cetp_mixed', default_basis={}, influent_character=InfluentCharacter.INDUSTRIAL_CETP, compatible_template_families=BIOLOGICAL + SALTS, expected_range=_ranges(250, 30, 5, (8000, 14000)), aliases=('industrial_cetp', 'industrial_mixed'), non_biological_defaults={}), 'industrial_cetp_dairy': FeedstockClass(name='industrial_cetp_dairy', default_basis={}, influent_character=InfluentCharacter.DAIRY, compatible_template_families=BIOLOGICAL, expected_range=_ranges(150, 25, 5, (8000, 12000)), aliases=('industrial_dairy', 'dairy')), 'industrial_cetp_textile_dye': FeedstockClass(name='industrial_cetp_textile_dye', default_basis={}, influent_character=InfluentCharacter.TEXTILE, compatible_template_families=BIOLOGICAL + SALTS, expected_range=_ranges(300, 35, 5, (8000, 14000)), aliases=('industrial_textile', 'textile'), non_biological_defaults={}), 'industrial_cetp_pharma_api': FeedstockClass(name='industrial_cetp_pharma_api', default_basis={}, influent_character=InfluentCharacter.PHARMA, compatible_template_families=BIOLOGICAL + SALTS, expected_range=_ranges(350, 40, 5, (8000, 15000)), aliases=('industrial_pharma', 'pharma'), non_biological_defaults={}), 'industrial_cetp_refinery_chem': FeedstockClass(name='industrial_cetp_refinery_chem', default_basis={}, influent_character=InfluentCharacter.REFINERY_CHEM, compatible_template_families=BIOLOGICAL + SALTS, expected_range=_ranges(300, 35, 5, (7000, 14000)), aliases=('industrial_refinery_chem', 'refinery_chem'), non_biological_defaults={}), 'industrial_food_processing': FeedstockClass(name='industrial_food_processing', default_basis={}, influent_character=InfluentCharacter.DAIRY, compatible_template_families=BIOLOGICAL, expected_range=_ranges(200, 30, 5, (7000, 13000))), 'industrial_brewery': FeedstockClass(name='industrial_brewery', default_basis={}, influent_character=InfluentCharacter.DAIRY, compatible_template_families=BIOLOGICAL, expected_range=_ranges(180, 25, 5, (7000, 13000))), 'industrial_pulp_paper': FeedstockClass(name='industrial_pulp_paper', default_basis={}, influent_character=InfluentCharacter.MIXED_INDUSTRIAL, compatible_template_families=BIOLOGICAL, expected_range=_ranges(250, 25, 5, (6000, 13000))), 'agrowaste_cbg_digestate': FeedstockClass(name='agrowaste_cbg_digestate', default_basis={}, influent_character=InfluentCharacter.MIXED_INDUSTRIAL, compatible_template_families=('adm1', 'anaerobic', 'asm2d'), expected_range=_ranges(500, 160, 15, (5000, 12000)), aliases=('cbg_digestate', 'digestate')), 'agrowaste_napier_grass': FeedstockClass(name='agrowaste_napier_grass', default_basis={}, influent_character=InfluentCharacter.MIXED_INDUSTRIAL, compatible_template_families=('adm1', 'anaerobic'), expected_range=_ranges(1000, 150, 20, (10000, 30000))), 'agrowaste_cattle_manure': FeedstockClass(name='agrowaste_cattle_manure', default_basis={}, influent_character=InfluentCharacter.MIXED_INDUSTRIAL, compatible_template_families=('adm1', 'anaerobic'), expected_range={}), 'agrowaste_press_mud': FeedstockClass(name='agrowaste_press_mud', default_basis={}, influent_character=InfluentCharacter.MIXED_INDUSTRIAL, compatible_template_families=('adm1', 'anaerobic'), expected_range=_ranges(1200, 180, 25, (10000, 35000))), 'agrowaste_poultry': FeedstockClass(name='agrowaste_poultry', default_basis={}, influent_character=InfluentCharacter.MIXED_INDUSTRIAL, compatible_template_families=('adm1', 'anaerobic'), expected_range=_ranges(1200, 220, 30, (10000, 35000))), 'food_waste_press_liquor': FeedstockClass(name='food_waste_press_liquor', default_basis={}, influent_character=InfluentCharacter.MIXED_INDUSTRIAL, compatible_template_families=('adm1', 'anaerobic'), expected_range=_ranges(33000, 4000, 900, (15000, 40000)), aliases=('ofmsw_press_liquor', 'food_waste_leachate', 'compost_press_liquor', 'food_waste_recycling_wastewater')), 'saline_industrial': FeedstockClass(name='saline_industrial', default_basis={}, influent_character=InfluentCharacter.MIXED_INDUSTRIAL, compatible_template_families=SALTS, expected_range=_ranges(80, 20, 3, (0, 0)), aliases=('saline',), non_biological_defaults={}), 'zld_brine': FeedstockClass(name='zld_brine', default_basis={}, influent_character=InfluentCharacter.MIXED_INDUSTRIAL, compatible_template_families=('evaporator', 'zld'), expected_range=_ranges(80, 20, 3, (0, 0)), aliases=('brine',), non_biological_defaults={})}
_ALIASES = {alias: name for name, feedstock in FEEDSTOCK_CLASSES.items() for alias in (name, *feedstock.aliases)}

def _normalize_name(name: str) -> str:
    return name.strip().lower().replace('-', '_').replace(' ', '_')

def get_feedstock_class(name: str) -> FeedstockClass:
    """Return a feedstock class by stable name or supported alias."""
    normalized = _normalize_name(name)
    try:
        return FEEDSTOCK_CLASSES[_ALIASES[normalized]]
    except KeyError as exc:
        raise KeyError(name) from exc
