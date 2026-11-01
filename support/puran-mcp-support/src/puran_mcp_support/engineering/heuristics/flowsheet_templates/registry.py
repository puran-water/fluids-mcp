from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Callable
ExpectedRangesProvider = Callable[[str | None], dict[str, tuple[float, float]]]
SOLVE_MATURITY_STRICT_OPTIMAL = 'strict_optimal'

@dataclass(frozen=True)
class AerobicTemplateBounds:
    """Industry-practice aerobic bounds consumed by WaterTAP ASM templates.

    Bounds reflect operating envelopes from Metcalf & Eddy / WEF MOP 8 for
    each technology class. RAS bounds enforce cage/basin MLSS physics on MBR
    templates (cage MLSS = basin MLSS x (1+1/RAS); RAS >= 3.0 keeps the ratio
    <= 1.33, controlling fouling). MLR bounds enforce active pre-denit
    operation (>=2.0 internal recycle; <2.0 starves the anoxic zone of
    nitrate).
    """
    ras_recycle_ratio: tuple[float, float] = (0.0, 10.0)
    mlr_recycle_ratio: tuple[float, float] = (0.0, 10.0)
    mlss_mg_l: tuple[float, float] = (2000.0, 15000.0)
    srt_d: tuple[float, float] = (5.0, 30.0)
    kla_per_hr: tuple[float, float] = (1.0, 20.0)
CONVENTIONAL_AEROBIC_BOUNDS = AerobicTemplateBounds(ras_recycle_ratio=(0.5, 1.5), mlr_recycle_ratio=(2.0, 5.0), mlss_mg_l=(2000.0, 10000.0))
IFAS_AEROBIC_BOUNDS = AerobicTemplateBounds(ras_recycle_ratio=(0.5, 2.0), mlr_recycle_ratio=(2.0, 5.0), mlss_mg_l=(1500.0, 8000.0))
MBR_AEROBIC_BOUNDS = AerobicTemplateBounds(ras_recycle_ratio=(3.0, 5.0), mlr_recycle_ratio=(2.0, 5.0), mlss_mg_l=(3000.0, 25000.0))

@dataclass(frozen=True)
class RecycleVariableSpec:
    """One declared system-level recycle variable for a template."""
    key: str
    ratio_name: str
    lower_bound: float
    upper_bound: float
    physical_stream: str

    def bounds(self) -> tuple[float, float]:
        return (self.lower_bound, self.upper_bound)

@dataclass(frozen=True)
class TopologyTemplate:
    """Template metadata consumed by MCP server solve orchestration."""
    name: str
    family: str
    implemented: bool
    mcp_template_name: str | None
    recycle_variables: tuple[RecycleVariableSpec, ...] = ()
    bounds: dict[str, tuple[float, float]] = field(default_factory=dict)
    scaling_classes: tuple[str, ...] = ()
    initializer: str | None = None
    expected_optimal_ranges: dict[str, tuple[float, float]] | ExpectedRangesProvider = field(default_factory=dict)
    solve_maturity: str | None = None
    maturity_note: str | None = None

    def recycle_keys(self) -> tuple[str, ...]:
        return tuple((spec.key for spec in self.recycle_variables))

    @property
    def design_grade_convergence(self) -> bool:
        """True only when this template is implemented AND the engine reaches a
        strict-optimal solve for it.

        Dispatch gates the "design-grade convergence" promise on this: implemented
        structural-gap templates stay dispatchable but return ``False`` here, so
        dispatch can surface them as honest screening-grade topologies that may
        legitimately return ``infeasible_structured``.
        """
        return bool(self.implemented and self.solve_maturity == SOLVE_MATURITY_STRICT_OPTIMAL)

    @property
    def engine_template_name(self) -> str | None:
        """Return the local MCP template name for server-side builders."""
        return self.mcp_template_name

    def expected_optimal_check(self, payload: dict[str, Any]) -> dict[str, Any]:
        """Return range-check diagnostics for settled solve values."""
        expected_ranges = self.get_expected_optimal_ranges(payload.get('feedstock_class'))
        failures: dict[str, dict[str, Any]] = {}
        aliases = {'MLSS_basin_mg_L': 'MLSS_mg_L', 'MLSS_cage_mg_L': 'MLSS_mg_L', 'SRT_actual_d': 'SRT_d'}
        for key, (lower, upper) in expected_ranges.items():
            value = payload.get(key)
            if value is None and key in aliases:
                value = payload.get(aliases[key])
            try:
                numeric = float(value)
            except (TypeError, ValueError):
                continue
            if numeric < lower or numeric > upper:
                failures[key] = {'value': numeric, 'expected': [lower, upper]}
        return {'ok': not failures, 'failures': failures}

    def get_expected_optimal_ranges(self, feedstock_class: str | None=None) -> dict[str, tuple[float, float]]:
        """Return settled-value sanity ranges for a feedstock class."""
        ranges = self.expected_optimal_ranges
        if callable(ranges):
            return dict(ranges(feedstock_class))
        return dict(ranges)

def _ras(bounds: tuple[float, float], stream: str) -> RecycleVariableSpec:
    return RecycleVariableSpec(key='ras_recycle_ratio', ratio_name='asm2d_ras_recycle_ratio', lower_bound=bounds[0], upper_bound=bounds[1], physical_stream=stream)

def _mlr() -> RecycleVariableSpec:
    return RecycleVariableSpec(key='mlr_recycle_ratio', ratio_name='asm2d_mlr_recycle_ratio', lower_bound=MLR[0], upper_bound=MLR[1], physical_stream='end-of-aerobic mixed liquor to anoxic basin')
ANAEROBIC_RECIRCULATION_RATIO = (1.0, 20.0)
ANAEROBIC_WASTE_FRACTION = (0.005, 0.95)

def _anaerobic_recirculation(stream: str) -> RecycleVariableSpec:
    """EGSB recirculation:feed ratio (the high-upflow liquid recirculation)."""
    return RecycleVariableSpec(key='recirculation_ratio', ratio_name='adm1_recirculation_ratio', lower_bound=ANAEROBIC_RECIRCULATION_RATIO[0], upper_bound=ANAEROBIC_RECIRCULATION_RATIO[1], physical_stream=stream)

def _anaerobic_waste_fraction(stream: str) -> RecycleVariableSpec:
    """AnMBR / anaerobic-contact WAS split off the biomass-concentrated stream.

    The single totalFlow split that sets the solids-retention time; the biomass
    recycle (membrane retentate / sludge RAS) is the complement (``1 - waste``).
    """
    return RecycleVariableSpec(key='waste_fraction', ratio_name='adm1_waste_fraction', lower_bound=ANAEROBIC_WASTE_FRACTION[0], upper_bound=ANAEROBIC_WASTE_FRACTION[1], physical_stream=stream)
CONVENTIONAL_RAS = CONVENTIONAL_AEROBIC_BOUNDS.ras_recycle_ratio
IFAS_RAS = IFAS_AEROBIC_BOUNDS.ras_recycle_ratio
MBR_RAS = MBR_AEROBIC_BOUNDS.ras_recycle_ratio
MLR = MBR_AEROBIC_BOUNDS.mlr_recycle_ratio
MLSS_MUNICIPAL_MBR = (2500.0, 12000.0)

def _asm2d_mbr_post_expected_ranges(feedstock_class: str | None) -> dict[str, tuple[float, float]]:
    """Return MBR POST MLSS sanity ranges by feedstock class.

    Ranges are warning-only operating envelopes, not solver bounds. They follow
    typical MBR MLSS practice by wastewater class from Metcalf & Eddy 5e,
    Table 9-2, and WEF MOP 8 Section 7.
    """
    default = {'MLSS_basin_mg_L': (3000.0, 8000.0), 'MLSS_cage_mg_L': (3500.0, 10000.0)}
    by_feedstock = {'municipal_low_strength': {'MLSS_basin_mg_L': (2500.0, 6000.0), 'MLSS_cage_mg_L': (3000.0, 7500.0)}, 'municipal_medium_strength': default, 'municipal_high_strength': {'MLSS_basin_mg_L': (4000.0, 10000.0), 'MLSS_cage_mg_L': (5000.0, 13000.0)}, 'industrial_cetp_dairy': {'MLSS_basin_mg_L': (5000.0, 12000.0), 'MLSS_cage_mg_L': (6500.0, 15000.0)}, 'industrial_cetp_mixed': {'MLSS_basin_mg_L': (8000.0, 15000.0), 'MLSS_cage_mg_L': (10000.0, 19000.0)}, 'industrial_cetp_textile_dye': {'MLSS_basin_mg_L': (8000.0, 18000.0), 'MLSS_cage_mg_L': (10000.0, 23000.0)}, 'industrial_cetp_pharma_api': {'MLSS_basin_mg_L': (10000.0, 18000.0), 'MLSS_cage_mg_L': (13000.0, 23000.0)}, 'industrial_cetp_refinery_chem': {'MLSS_basin_mg_L': (8000.0, 15000.0), 'MLSS_cage_mg_L': (10000.0, 19000.0)}, 'agrowaste_cbg_digestate': {'MLSS_basin_mg_L': (5000.0, 10000.0), 'MLSS_cage_mg_L': (6000.0, 13000.0)}}
    return by_feedstock.get(str(feedstock_class or ''), default)
TOPOLOGY_TEMPLATES: dict[str, TopologyTemplate] = {'asm2d_mle': TopologyTemplate(name='asm2d_mle', family='asm2d', implemented=True, mcp_template_name='asm2d_mle', recycle_variables=(_ras(CONVENTIONAL_RAS, 'clarifier underflow to anoxic basin'), _mlr()), bounds={'RAS_q_over_qin': CONVENTIONAL_RAS, 'MLR_q_over_qin': MLR}, scaling_classes=('mass_balance', 'reaction_rates', 'hydraulic', 'retention', 'recycle_ratios'), initializer='asm2d_recycle', expected_optimal_ranges={'MLSS_basin_mg_L': (3000.0, 10000.0)}), 'asm2d_mle_post': TopologyTemplate(name='asm2d_mle_post', family='asm2d', implemented=True, mcp_template_name='asm2d_mle_post', recycle_variables=(_ras(CONVENTIONAL_RAS, 'clarifier underflow to anoxic basin'), _mlr()), bounds={'RAS_q_over_qin': CONVENTIONAL_RAS, 'MLR_q_over_qin': MLR}, scaling_classes=('mass_balance', 'reaction_rates', 'hydraulic', 'retention', 'recycle_ratios'), initializer='asm2d_recycle', expected_optimal_ranges={'MLSS_basin_mg_L': (3000.0, 10000.0)}), 'asm2d_mbr': TopologyTemplate(name='asm2d_mbr', family='asm2d', implemented=True, mcp_template_name='asm2d_mbr', recycle_variables=(_ras(MBR_RAS, 'MBR cage drain to anoxic basin'), _mlr()), bounds={'RAS_q_over_qin': MBR_RAS, 'MLR_q_over_qin': MLR, 'MLSS_cage_to_basin_ratio': (1.15, 1.35)}, scaling_classes=('mass_balance', 'reaction_rates', 'hydraulic', 'retention', 'recycle_ratios'), initializer='asm2d_recycle', expected_optimal_ranges={'MLSS_basin_mg_L': MLSS_MUNICIPAL_MBR, 'MLSS_cage_mg_L': (8000.0, 12000.0)}), 'asm2d_mbr_aerobic': TopologyTemplate(name='asm2d_mbr_aerobic', family='asm2d', implemented=True, mcp_template_name='asm2d_mbr_aerobic', recycle_variables=(_ras(MBR_RAS, 'MBR cage drain to aerobic basin'),), bounds={'RAS_q_over_qin': MBR_RAS, 'MLSS_cage_to_basin_ratio': (1.15, 1.35)}, scaling_classes=('mass_balance', 'reaction_rates', 'hydraulic', 'retention', 'recycle_ratios'), initializer='asm2d_recycle', expected_optimal_ranges={'MLSS_basin_mg_L': MLSS_MUNICIPAL_MBR, 'MLSS_cage_mg_L': (8000.0, 12000.0)}, solve_maturity=SOLVE_MATURITY_STRICT_OPTIMAL, maturity_note='Direct ModifiedASM2d fixed-volume rating at a stated WAS flow (the walk holds the SRT on it), strict-optimal from the plant-state start; the generic SRT-target ASM2d pipeline refuses this topology by name (its initializer is MLE-shaped).'), 'asm2d_mbr_post': TopologyTemplate(name='asm2d_mbr_post', family='asm2d', implemented=True, mcp_template_name='asm2d_mbr_post', recycle_variables=(_ras(MBR_RAS, 'MBR cage drain to anoxic basin'), _mlr()), bounds={'RAS_q_over_qin': MBR_RAS, 'MLR_q_over_qin': MLR, 'MLSS_cage_to_basin_ratio': (1.15, 1.35)}, scaling_classes=('mass_balance', 'reaction_rates', 'hydraulic', 'retention', 'recycle_ratios'), initializer='asm2d_recycle', expected_optimal_ranges=_asm2d_mbr_post_expected_ranges), 'asm2d_ifas': TopologyTemplate(name='asm2d_ifas', family='asm2d', implemented=False, mcp_template_name=None, recycle_variables=(_ras(IFAS_RAS, 'clarifier underflow to IFAS anoxic basin'), _mlr()), bounds={'RAS_q_over_qin': IFAS_RAS, 'MLR_q_over_qin': MLR}), 'asm2d_a2o': TopologyTemplate(name='asm2d_a2o', family='asm2d', implemented=False, mcp_template_name=None, recycle_variables=(_ras(CONVENTIONAL_RAS, 'clarifier underflow to anaerobic basin'), _mlr()), bounds={'RAS_q_over_qin': CONVENTIONAL_RAS, 'MLR_q_over_qin': MLR}), 'asm2d_a2o_mbr': TopologyTemplate(name='asm2d_a2o_mbr', family='asm2d', implemented=False, mcp_template_name=None, recycle_variables=(_ras(MBR_RAS, 'MBR cage drain to anaerobic basin'), _mlr()), bounds={'RAS_q_over_qin': MBR_RAS, 'MLR_q_over_qin': MLR}), 'asm1_cas': TopologyTemplate(name='asm1_cas', family='asm1', implemented=False, mcp_template_name=None, recycle_variables=(_ras(CONVENTIONAL_RAS, 'clarifier underflow to aeration basin'),), bounds={'RAS_q_over_qin': CONVENTIONAL_RAS}), 'asm1_eos': TopologyTemplate(name='asm1_eos', family='asm1', implemented=False, mcp_template_name=None), 'asm1_mle': TopologyTemplate(name='asm1_mle', family='asm1', implemented=False, mcp_template_name=None, recycle_variables=(_ras(CONVENTIONAL_RAS, 'clarifier underflow to anoxic basin'), _mlr()), bounds={'RAS_q_over_qin': CONVENTIONAL_RAS, 'MLR_q_over_qin': MLR}), 'anaerobic_uasb': TopologyTemplate(name='anaerobic_uasb', family='anaerobic', implemented=True, mcp_template_name='adm1_uasb', scaling_classes=('mass_balance', 'reaction_rates', 'hydraulic', 'retention'), initializer='adm1_digester', solve_maturity=SOLVE_MATURITY_STRICT_OPTIMAL, maturity_note='Single methanogenic mADM1 CSTR (combined acid/methane), biogas take-off, no recycle — the strict-optimal single-reactor baseline.'), 'anaerobic_cstr': TopologyTemplate(name='anaerobic_cstr', family='anaerobic', implemented=True, mcp_template_name='adm1_cstr', scaling_classes=('mass_balance', 'reaction_rates', 'hydraulic', 'retention'), initializer='adm1_digester', solve_maturity=SOLVE_MATURITY_STRICT_OPTIMAL, maturity_note='Literature-based modified ADM1 kernel (QSDsan lineage); kinetic parameters are not site fits.'), 'anaerobic_egsb': TopologyTemplate(name='anaerobic_egsb', family='anaerobic', implemented=True, mcp_template_name='adm1_egsb', recycle_variables=(_anaerobic_recirculation('EGSB digestate recirculation to feed mixer'),), bounds={'recirculation_ratio': ANAEROBIC_RECIRCULATION_RATIO}, scaling_classes=('mass_balance', 'reaction_rates', 'hydraulic', 'retention', 'recycle_ratios'), initializer='adm1_recycle', solve_maturity=SOLVE_MATURITY_STRICT_OPTIMAL, maturity_note='Methanogenic digester + high upflow liquid recirculation closed by the ADM1 recycle tear. Strict-optimal healthy methanogenic baseline under the corrected non-EBPR canonical feed.'), 'anaerobic_anmbr': TopologyTemplate(name='anaerobic_anmbr', family='anaerobic', implemented=True, mcp_template_name='adm1_anmbr', recycle_variables=(_anaerobic_waste_fraction('AnMBR retentate bleed (WAS); recycle is the complement'),), bounds={'waste_fraction': ANAEROBIC_WASTE_FRACTION}, scaling_classes=('mass_balance', 'reaction_rates', 'hydraulic', 'retention', 'recycle_ratios'), initializer='adm1_recycle', solve_maturity=SOLVE_MATURITY_STRICT_OPTIMAL, maturity_note='Methanogenic digester + biomass-retaining membrane with retentate recycle and WAS bleed (SRT decoupling) closed by the ADM1 recycle tear. Feasible 25 d SRT default strict-solves to a healthy methanogenic baseline.'), 'anaerobic_contact': TopologyTemplate(name='anaerobic_contact', family='anaerobic', implemented=True, mcp_template_name='adm1_contact', recycle_variables=(_anaerobic_waste_fraction('anaerobic-contact thickened sludge bleed (WAS); RAS recycle is the complement'),), bounds={'waste_fraction': ANAEROBIC_WASTE_FRACTION}, scaling_classes=('mass_balance', 'reaction_rates', 'hydraulic', 'retention', 'recycle_ratios'), initializer='adm1_recycle', solve_maturity=SOLVE_MATURITY_STRICT_OPTIMAL, maturity_note='CSTR digester + gravity thickener with sludge-return (RAS) recycle and WAS bleed closed by the ADM1 recycle tear. Feasible 25 d SRT default strict-solves to a healthy methanogenic baseline.'), 'adm1_two_stage': TopologyTemplate(name='adm1_two_stage', family='adm1', implemented=True, mcp_template_name='adm1_two_stage', scaling_classes=('mass_balance', 'reaction_rates', 'hydraulic', 'retention'), initializer='adm1_digester', solve_maturity=SOLVE_MATURITY_STRICT_OPTIMAL, maturity_note='Phase-separated acidogenic -> methanogenic forward cascade (no recycle), each stage with its own biogas headspace — strict-optimal.'), 'primary_clarifier_high_rate': TopologyTemplate(name='primary_clarifier_high_rate', family='physchem', implemented=False, mcp_template_name=None), 'primary_clarifier_lamella': TopologyTemplate(name='primary_clarifier_lamella', family='physchem', implemented=False, mcp_template_name=None), 'equalization_buffered': TopologyTemplate(name='equalization_buffered', family='hydraulic', implemented=False, mcp_template_name=None), 'equalization_continuous': TopologyTemplate(name='equalization_continuous', family='hydraulic', implemented=False, mcp_template_name=None), 'coagulation_flocculation_inline': TopologyTemplate(name='coagulation_flocculation_inline', family='physchem', implemented=False, mcp_template_name=None), 'coagulation_flocculation_plus_clarifier': TopologyTemplate(name='coagulation_flocculation_plus_clarifier', family='physchem', implemented=False, mcp_template_name=None), 'daf_inline': TopologyTemplate(name='daf_inline', family='physchem', implemented=False, mcp_template_name=None), 'daf_recycle': TopologyTemplate(name='daf_recycle', family='physchem', implemented=False, mcp_template_name=None), 'mbr_post': TopologyTemplate(name='mbr_post', family='asm2d', implemented=False, mcp_template_name=None, recycle_variables=(_ras(MBR_RAS, 'standalone MBR cage drain to upstream basin'),), bounds={'RAS_q_over_qin': MBR_RAS}), 'uf_post': TopologyTemplate(name='uf_post', family='membrane', implemented=False, mcp_template_name=None), 'ro_two_stage': TopologyTemplate(name='ro_two_stage', family='ro', implemented=True, mcp_template_name='ro_train'), 'ro_three_stage': TopologyTemplate(name='ro_three_stage', family='ro', implemented=False, mcp_template_name=None), 'ro_two_pass': TopologyTemplate(name='ro_two_pass', family='ro', implemented=False, mcp_template_name=None), 'mvr_single_effect': TopologyTemplate(name='mvr_single_effect', family='evaporator', implemented=True, mcp_template_name='mvc_evaporator'), 'mee_two_effect': TopologyTemplate(name='mee_two_effect', family='evaporator', implemented=False, mcp_template_name=None), 'mee_three_effect': TopologyTemplate(name='mee_three_effect', family='evaporator', implemented=False, mcp_template_name=None), 'atfd_horizontal': TopologyTemplate(name='atfd_horizontal', family='evaporator', implemented=False, mcp_template_name=None), 'atfd_drum': TopologyTemplate(name='atfd_drum', family='evaporator', implemented=False, mcp_template_name=None), 'ix_softener_two_train': TopologyTemplate(name='ix_softener_two_train', family='ix', implemented=False, mcp_template_name=None), 'ix_de_demineralizer': TopologyTemplate(name='ix_de_demineralizer', family='ix', implemented=False, mcp_template_name=None), 'chlorination_breakpoint': TopologyTemplate(name='chlorination_breakpoint', family='disinfection', implemented=False, mcp_template_name=None), 'ozone_combined_uv': TopologyTemplate(name='ozone_combined_uv', family='disinfection', implemented=False, mcp_template_name=None)}

def _normalize_name(name: str) -> str:
    return name.strip().lower().replace('-', '_').replace(' ', '_')

def _build_mcp_template_aliases() -> dict[str, TopologyTemplate]:
    aliases: dict[str, TopologyTemplate] = {}
    for template in TOPOLOGY_TEMPLATES.values():
        if not template.mcp_template_name:
            continue
        key = _normalize_name(template.mcp_template_name)
        if key in TOPOLOGY_TEMPLATES:
            continue
        if key not in aliases:
            aliases[key] = template
            continue
        technology_suffix = key.rsplit('_', 1)[-1]
        if _normalize_name(template.name).endswith(f'_{technology_suffix}'):
            aliases[key] = template
    return aliases
_MCP_TEMPLATE_ALIASES: dict[str, TopologyTemplate] = _build_mcp_template_aliases()

def get_topology_template(name: str) -> TopologyTemplate:
    """Return a topology template by stable registry name, or by its engine-side
    ``mcp_template_name`` alias (OBS-002 compatibility)."""
    normalized = _normalize_name(name)
    template = TOPOLOGY_TEMPLATES.get(normalized) or _MCP_TEMPLATE_ALIASES.get(normalized)
    if template is None:
        raise KeyError(name)
    return template
