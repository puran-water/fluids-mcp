from __future__ import annotations
from typing import Any, Literal, Mapping, Optional, get_args
from pydantic import Field, model_validator
from puran_mcp_support.engineering.pydantic_base import StrictFiniteBaseModel
from .opex import OpexCategory
CostBasis = Literal['purchased', 'installed', 'package', 'bare_module', 'tic', 'service']
CostSource = Literal['parametric_model', 'vendor_quote', 'historical', 'manual', 'allowance', 'blended']
Confidence = Literal['H', 'M', 'L']
_ESTIMATE_LINE_FIELDS = ('equipment_tag', 'description', 'quantity', 'unit_cost_usd', 'cost_basis', 'cost_source', 'includes', 'installation_factor', 'sizing_params', 'upstream_artifact_hashes', 'confidence')

class CapexLine(StrictFiniteBaseModel):
    """One CAPEX line, estimate_line-shaped plus TEA provenance fields.

    ``total_cost_usd`` defaults to ``quantity * unit_cost_usd`` and
    ``installed_cost_usd`` to ``total_cost_usd * installation_factor`` (the same
    computation ``add_estimate_line`` performs server-side) when not supplied.
    """
    equipment_tag: Optional[str] = None
    description: str
    quantity: float = 1
    unit_cost_usd: float
    total_cost_usd: Optional[float] = None
    cost_basis: CostBasis = 'purchased'
    cost_source: Optional[CostSource] = None
    includes: dict[str, Any] = Field(default_factory=dict)
    installation_factor: Optional[float] = None
    installed_cost_usd: Optional[float] = None
    confidence: Optional[Confidence] = None
    upstream_artifact_hashes: list[str] = Field(default_factory=list)
    sizing_params: dict[str, Any] = Field(default_factory=dict)
    base_year: Optional[int] = None
    cost_index: Literal['CEPCI'] = 'CEPCI'
    base_index_value: Optional[float] = None
    location_factor: Optional[float] = None
    battery_limit: Optional[str] = None
    area_code: Optional[str] = None

    @model_validator(mode='after')
    def _fill_computed_costs(self) -> 'CapexLine':
        if self.total_cost_usd is None:
            object.__setattr__(self, 'total_cost_usd', self.quantity * self.unit_cost_usd)
        if self.installed_cost_usd is None and self.installation_factor is not None:
            object.__setattr__(self, 'installed_cost_usd', self.total_cost_usd * self.installation_factor)
        return self

    @property
    def effective_installed_usd(self) -> float:
        """Installed cost for rollups: explicit installed, else the total as-is."""
        if self.installed_cost_usd is not None:
            return self.installed_cost_usd
        assert self.total_cost_usd is not None
        return self.total_cost_usd

    def to_estimate_line_kwargs(self) -> dict[str, Any]:
        """kwargs for procurement-mcp ``add_estimate_line`` (TEA fields excluded).

        ``total_cost_usd`` / ``installed_cost_usd`` are excluded on purpose — the
        server recomputes both from quantity/unit_cost/installation_factor.
        """
        return {name: getattr(self, name) for name in _ESTIMATE_LINE_FIELDS}

    @classmethod
    def from_estimate_line_row(cls, row: Mapping[str, Any]) -> 'CapexLine':
        """Build from an ``estimate_line`` row mapping (extra columns ignored).

        NULL-able columns fall back to the model defaults; NUMERIC values may
        arrive as Decimal/str and are coerced to float.
        """

        def _num(key: str) -> Optional[float]:
            value = row.get(key)
            return None if value is None else float(value)
        return cls(equipment_tag=row.get('equipment_tag'), description=row['description'], quantity=_num('quantity') if row.get('quantity') is not None else 1, unit_cost_usd=float(row['unit_cost_usd']), total_cost_usd=_num('total_cost_usd'), cost_basis=row.get('cost_basis') or 'purchased', cost_source=row.get('cost_source'), includes=row.get('includes') or {}, installation_factor=_num('installation_factor'), installed_cost_usd=_num('installed_cost_usd'), confidence=row.get('confidence'), upstream_artifact_hashes=list(row.get('upstream_artifact_hashes') or []), sizing_params=row.get('sizing_params') or {})

class OpexLine(StrictFiniteBaseModel):
    """One annual operating-cost line with its serialized computation basis."""
    category: OpexCategory
    annual_usd: float
    basis: dict[str, Any] = Field(default_factory=dict)
    confidence: Optional[Confidence] = None

class CostLedger(StrictFiniteBaseModel):
    """The plant-TEA ledger: CAPEX + OPEX lines plus rollup properties."""
    capex_lines: list[CapexLine] = Field(default_factory=list)
    opex_lines: list[OpexLine] = Field(default_factory=list)
    currency: str = 'USD'
    aace_class: Optional[str] = None

    @property
    def purchased_total_usd(self) -> float:
        return sum((line.total_cost_usd or 0.0 for line in self.capex_lines))

    @property
    def installed_total_usd(self) -> float:
        return sum((line.effective_installed_usd for line in self.capex_lines))

    @property
    def annual_opex_total_usd(self) -> float:
        return sum((line.annual_usd for line in self.opex_lines))

    @property
    def opex_by_category(self) -> dict[str, float]:
        out: dict[str, float] = {}
        for line in self.opex_lines:
            key = line.category.value
            out[key] = out.get(key, 0.0) + line.annual_usd
        return out
