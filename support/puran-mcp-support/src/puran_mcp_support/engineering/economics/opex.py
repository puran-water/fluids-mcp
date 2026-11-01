from __future__ import annotations
from enum import Enum
from typing import TYPE_CHECKING, Iterable, Literal, Optional
from pydantic import model_validator
from puran_mcp_support.engineering.pydantic_base import StrictFiniteBaseModel
from .annualize import HOURS_PER_YEAR
if TYPE_CHECKING:
    from .ledger import OpexLine
DEFAULT_SLUDGE_DISPOSAL_USD_PER_WET_TONNE = 68.5
DAYS_PER_YEAR = 365.0

class OpexCategory(str, Enum):
    """Vocabulary for ``ledger.OpexLine.category``."""
    ENERGY = 'energy'
    CHEMICAL = 'chemical'
    REPLACEMENT = 'replacement'
    SLUDGE_DISPOSAL = 'sludge_disposal'
    LABOR = 'labor'
    MAINTENANCE = 'maintenance'

class OpexItem(StrictFiniteBaseModel):
    """Common surface: ``category``, ``annual_usd``, ``to_opex_line()``."""
    category: OpexCategory
    description: Optional[str] = None

    @property
    def annual_usd(self) -> float:
        raise NotImplementedError

    def to_opex_line(self, confidence: Optional[Literal['H', 'M', 'L']]=None) -> 'OpexLine':
        from .ledger import OpexLine
        return OpexLine(category=self.category, annual_usd=self.annual_usd, basis=self.model_dump(mode='json'), confidence=confidence)

class EnergyOpex(OpexItem):
    """annual = power_kw * hours/yr * tariff."""
    category: Literal[OpexCategory.ENERGY] = OpexCategory.ENERGY
    power_kw: float
    operating_hours_per_year: float = HOURS_PER_YEAR
    tariff_usd_per_kwh: float

    @property
    def annual_usd(self) -> float:
        return self.power_kw * self.operating_hours_per_year * self.tariff_usd_per_kwh

class ChemicalOpex(OpexItem):
    """Either a direct ``kg_per_year``, or ``dose_mg_per_l`` + ``flow_m3_per_day``.

    dose [mg/L] * flow [m3/d] = g/d; /1000 -> kg/d; *365 -> kg/yr.
    """
    category: Literal[OpexCategory.CHEMICAL] = OpexCategory.CHEMICAL
    dose_mg_per_l: Optional[float] = None
    flow_m3_per_day: Optional[float] = None
    kg_per_year: Optional[float] = None
    price_usd_per_kg: float

    @model_validator(mode='after')
    def _one_basis(self) -> 'ChemicalOpex':
        dose_basis = self.dose_mg_per_l is not None and self.flow_m3_per_day is not None
        mass_basis = self.kg_per_year is not None
        if dose_basis == mass_basis:
            raise ValueError('ChemicalOpex needs exactly one basis: kg_per_year, or dose_mg_per_l together with flow_m3_per_day')
        return self

    @property
    def annual_kg(self) -> float:
        if self.kg_per_year is not None:
            return self.kg_per_year
        assert self.dose_mg_per_l is not None and self.flow_m3_per_day is not None
        return self.dose_mg_per_l * self.flow_m3_per_day / 1000.0 * DAYS_PER_YEAR

    @property
    def annual_usd(self) -> float:
        return self.annual_kg * self.price_usd_per_kg

class ReplacementOpex(OpexItem):
    """Sinking reserve: capex * fraction/yr, or capex / life_years (straight line)."""
    category: Literal[OpexCategory.REPLACEMENT] = OpexCategory.REPLACEMENT
    asset_capex_usd: float
    replacement_fraction_per_year: Optional[float] = None
    life_years: Optional[float] = None

    @model_validator(mode='after')
    def _one_basis(self) -> 'ReplacementOpex':
        if (self.replacement_fraction_per_year is None) == (self.life_years is None):
            raise ValueError('ReplacementOpex needs exactly one of replacement_fraction_per_year or life_years')
        if self.life_years is not None and self.life_years <= 0:
            raise ValueError(f'life_years must be > 0, got {self.life_years!r}')
        return self

    @property
    def annual_usd(self) -> float:
        frac = self.replacement_fraction_per_year if self.replacement_fraction_per_year is not None else 1.0 / float(self.life_years)
        return self.asset_capex_usd * frac

class SludgeDisposalOpex(OpexItem):
    """Disposal billed per WET tonne hauled (the EXPOsan WERF price basis).

    Give ``wet_tonnes_per_year`` directly, or ``dry_tonnes_per_year`` +
    ``solids_fraction`` (wet = dry / solids_fraction). The default price is the
    WERF 68.5 USD/wet-ton figure (permalink at module top); override
    ``price_usd_per_tonne`` for a local haul contract.
    """
    category: Literal[OpexCategory.SLUDGE_DISPOSAL] = OpexCategory.SLUDGE_DISPOSAL
    dry_tonnes_per_year: Optional[float] = None
    wet_tonnes_per_year: Optional[float] = None
    solids_fraction: Optional[float] = None
    price_usd_per_tonne: float = DEFAULT_SLUDGE_DISPOSAL_USD_PER_WET_TONNE

    @model_validator(mode='after')
    def _one_basis(self) -> 'SludgeDisposalOpex':
        wet_basis = self.wet_tonnes_per_year is not None
        dry_basis = self.dry_tonnes_per_year is not None and self.solids_fraction is not None
        if wet_basis == dry_basis:
            raise ValueError('SludgeDisposalOpex needs exactly one basis: wet_tonnes_per_year, or dry_tonnes_per_year together with solids_fraction')
        if self.solids_fraction is not None and (not 0.0 < self.solids_fraction <= 1.0):
            raise ValueError(f'solids_fraction must be within (0, 1], got {self.solids_fraction!r}')
        return self

    @property
    def annual_wet_tonnes(self) -> float:
        if self.wet_tonnes_per_year is not None:
            return self.wet_tonnes_per_year
        assert self.dry_tonnes_per_year is not None and self.solids_fraction is not None
        return self.dry_tonnes_per_year / self.solids_fraction

    @property
    def annual_usd(self) -> float:
        return self.annual_wet_tonnes * self.price_usd_per_tonne

class LaborOpex(OpexItem):
    """fte_count * loaded_annual_rate_usd, or a direct annual_lump_sum_usd."""
    category: Literal[OpexCategory.LABOR] = OpexCategory.LABOR
    fte_count: Optional[float] = None
    loaded_annual_rate_usd: Optional[float] = None
    annual_lump_sum_usd: Optional[float] = None

    @model_validator(mode='after')
    def _one_basis(self) -> 'LaborOpex':
        fte_basis = self.fte_count is not None and self.loaded_annual_rate_usd is not None
        lump_basis = self.annual_lump_sum_usd is not None
        if fte_basis == lump_basis:
            raise ValueError('LaborOpex needs exactly one basis: annual_lump_sum_usd, or fte_count together with loaded_annual_rate_usd')
        return self

    @property
    def annual_usd(self) -> float:
        if self.annual_lump_sum_usd is not None:
            return self.annual_lump_sum_usd
        assert self.fte_count is not None and self.loaded_annual_rate_usd is not None
        return self.fte_count * self.loaded_annual_rate_usd

class MaintenanceOpex(OpexItem):
    """annual = capex_basis_usd * pct_of_capex."""
    category: Literal[OpexCategory.MAINTENANCE] = OpexCategory.MAINTENANCE
    capex_basis_usd: float
    pct_of_capex: float

    @property
    def annual_usd(self) -> float:
        return self.capex_basis_usd * self.pct_of_capex

def build_opex_ledger(items: Iterable[OpexItem]) -> 'list[OpexLine]':
    """Flatten typed OPEX items into ledger rows (order-preserving)."""
    return [item.to_opex_line() for item in items]
