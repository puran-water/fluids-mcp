"""Commercial-product dose conversion shared by numeric and Pyomo adapters."""
from __future__ import annotations
import math
from dataclasses import asdict, dataclass
from typing import Any, Mapping
from puran_mcp_support.schema import PRODUCT_DOSE_UNIT_FACTORS_TO_KG_PER_DAY, CommercialProductDefinition, ProductDose, ReagentDefinition, normalized_product_dose_unit, product_dose_is_volume_based
US_GALLON_TO_LITRE = 3.785411784

@dataclass(frozen=True)
class ProductDoseConversion:
    """Resolved source terms and first-order input uncertainty for one dose."""
    product_id: str
    lot_id: str | None
    reagent_id: str
    assay_basis: str
    assay_mass_fraction_mean: float
    assay_mass_fraction_standard_uncertainty: float
    density_kg_l_mean: float | None
    density_kg_l_standard_uncertainty: float | None
    dose_basis: str
    normalized_dose_unit: str
    product_mass_kg_d: float
    product_mass_standard_uncertainty_kg_d: float
    product_volume_l_d: float | None
    active_reagent_mass_kg_d: float
    active_reagent_mol_d: float
    active_reagent_mol_standard_uncertainty_mol_d: float
    source_component_mol_d: Mapping[str, float]
    free_acidity_eq_d: float | None
    free_acidity_standard_uncertainty_eq_d: float | None
    free_acidity_status: str
    free_acidity_resolved: bool
    intrinsic_hydrolysis_acidity_treatment: str

    def public_dict(self) -> dict[str, object]:
        return asdict(self)

def _uncertain_quantity_from_mapping(value: Any, *, field: str):
    from puran_mcp_support.schema import UncertainQuantity
    if value is None or isinstance(value, UncertainQuantity):
        return value
    if not isinstance(value, Mapping):
        raise ValueError(f'{field} must be an uncertain-quantity mapping or null')
    unknown = sorted(set(value) - {'mean', 'standard_uncertainty'})
    if unknown:
        raise ValueError(f'{field} has unknown fields {unknown}')
    try:
        return UncertainQuantity(mean=float(value['mean']), standard_uncertainty=float(value['standard_uncertainty']))
    except KeyError as error:
        raise ValueError(f'{field} requires mean and standard_uncertainty') from error

def coerce_commercial_product_definition(value: Any) -> CommercialProductDefinition:
    """Coerce a JSON-like product object without accepting unknown fields."""
    if isinstance(value, CommercialProductDefinition):
        return value
    if not isinstance(value, Mapping):
        raise ValueError('commercial_product must be a mapping')
    allowed = {'product_id', 'lot_id', 'assay_mass_fraction', 'density_kg_l', 'free_acidity_eq_per_kg', 'free_acidity_status', 'assay_basis', 'free_acidity_counterion_stoichiometry'}
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise ValueError(f'commercial_product has unknown fields {unknown}')
    try:
        return CommercialProductDefinition(product_id=str(value['product_id']), lot_id=None if value.get('lot_id') is None else str(value['lot_id']), assay_mass_fraction=_uncertain_quantity_from_mapping(value['assay_mass_fraction'], field='commercial_product.assay_mass_fraction'), density_kg_l=_uncertain_quantity_from_mapping(value.get('density_kg_l'), field='commercial_product.density_kg_l'), free_acidity_eq_per_kg=_uncertain_quantity_from_mapping(value.get('free_acidity_eq_per_kg'), field='commercial_product.free_acidity_eq_per_kg'), free_acidity_status=str(value['free_acidity_status']), assay_basis=str(value.get('assay_basis', 'active_reagent_mass_fraction')), free_acidity_counterion_stoichiometry={str(component_id): float(coefficient) for component_id, coefficient in (value.get('free_acidity_counterion_stoichiometry') or {}).items()})
    except KeyError as error:
        raise ValueError(f'commercial_product is missing required field {error.args[0]!r}') from error

def coerce_product_dose(value: Any) -> ProductDose:
    """Coerce a JSON-like product-dose object without accepting unknown fields."""
    if isinstance(value, ProductDose):
        return value
    if not isinstance(value, Mapping):
        raise ValueError('product_dose must be a mapping')
    unknown = sorted(set(value) - {'value', 'unit', 'basis'})
    if unknown:
        raise ValueError(f'product_dose has unknown fields {unknown}')
    try:
        return ProductDose(value=float(value['value']), unit=str(value['unit']), basis=str(value.get('basis', 'as_delivered_product')))
    except KeyError as error:
        raise ValueError(f'product_dose is missing required field {error.args[0]!r}') from error

def _rss_product(left: float, left_uncertainty: float, right: float, right_uncertainty: float) -> float:
    """First-order standard uncertainty of a product of independent inputs."""
    return math.sqrt((right * left_uncertainty) ** 2 + (left * right_uncertainty) ** 2)

def convert_commercial_product_dose(*, product: CommercialProductDefinition, dose: ProductDose, reagent: ReagentDefinition, require_resolved_free_acidity: bool=False) -> ProductDoseConversion:
    """Convert an as-delivered product rate to conserved reagent source terms.

    The assay is explicitly the active-reagent mass fraction. Volume doses use
    the declared density; mass doses do not. Measured free acidity is emitted as
    a separate H+ source. Intrinsic FeCl3 hydrolysis acidity is deliberately not
    added here because the shared reaction network generates it when Fe(III)
    hydrolyses to HFO.
    """
    if reagent.molar_mass_kg_mol is None:
        raise ValueError(f'reagent {reagent.reagent_id!r} has no molar mass')
    if product.assay_basis != 'active_reagent_mass_fraction':
        raise ValueError('unsupported commercial-product assay basis')
    if require_resolved_free_acidity and product.free_acidity_status == 'not_measured':
        raise ValueError('commercial-product free acidity is not measured; product-dose simulation would silently impose a zero acid load')
    normalized_unit = normalized_product_dose_unit(dose.unit)
    value = float(dose.value)
    mass_uncertainty = 0.0
    if product_dose_is_volume_based(normalized_unit):
        if product.density_kg_l is None:
            raise ValueError('volume-based product dose requires measured product density')
        volume_l_d = value
        if normalized_unit == 'l_product_h':
            volume_l_d *= 24.0
        elif normalized_unit == 'us_gal_product_d':
            volume_l_d *= US_GALLON_TO_LITRE
        elif normalized_unit == 'us_gal_product_h':
            volume_l_d *= 24.0 * US_GALLON_TO_LITRE
        product_mass = volume_l_d * product.density_kg_l.mean
        mass_uncertainty = volume_l_d * product.density_kg_l.standard_uncertainty
    else:
        factor = PRODUCT_DOSE_UNIT_FACTORS_TO_KG_PER_DAY[normalized_unit]
        if factor is None:
            raise RuntimeError('mass-dose unit unexpectedly lacks a conversion factor')
        product_mass = value * factor
        volume_l_d = None if product.density_kg_l is None else product_mass / product.density_kg_l.mean
    assay = product.assay_mass_fraction
    active_mass = product_mass * assay.mean
    active_mass_uncertainty = _rss_product(product_mass, mass_uncertainty, assay.mean, assay.standard_uncertainty)
    active_moles = active_mass / reagent.molar_mass_kg_mol
    active_moles_uncertainty = active_mass_uncertainty / reagent.molar_mass_kg_mol
    sources = {component_id: active_moles * float(coefficient) for component_id, coefficient in reagent.component_stoichiometry.items()}
    free_acidity = None
    free_acidity_uncertainty = None
    if product.free_acidity_status == 'measured':
        acidity = product.free_acidity_eq_per_kg
        if acidity is None:
            raise RuntimeError('measured free acidity has no value')
        free_acidity = product_mass * acidity.mean
        free_acidity_uncertainty = _rss_product(product_mass, mass_uncertainty, acidity.mean, acidity.standard_uncertainty)
        sources['H_plus'] = sources.get('H_plus', 0.0) + free_acidity
        for component_id, coefficient in product.free_acidity_counterion_stoichiometry.items():
            sources[component_id] = sources.get(component_id, 0.0) + free_acidity * float(coefficient)
    elif product.free_acidity_status == 'not_applicable':
        free_acidity = 0.0
        free_acidity_uncertainty = 0.0
    values = (product_mass, mass_uncertainty, active_mass, active_moles, active_moles_uncertainty, *sources.values())
    if any((not math.isfinite(item) or item < 0 for item in values)):
        raise ValueError('commercial-product conversion produced a nonphysical result')
    return ProductDoseConversion(product_id=product.product_id, lot_id=product.lot_id, reagent_id=reagent.reagent_id, assay_basis=product.assay_basis, assay_mass_fraction_mean=assay.mean, assay_mass_fraction_standard_uncertainty=assay.standard_uncertainty, density_kg_l_mean=None if product.density_kg_l is None else product.density_kg_l.mean, density_kg_l_standard_uncertainty=None if product.density_kg_l is None else product.density_kg_l.standard_uncertainty, dose_basis=dose.basis, normalized_dose_unit=normalized_unit, product_mass_kg_d=product_mass, product_mass_standard_uncertainty_kg_d=mass_uncertainty, product_volume_l_d=volume_l_d, active_reagent_mass_kg_d=active_mass, active_reagent_mol_d=active_moles, active_reagent_mol_standard_uncertainty_mol_d=active_moles_uncertainty, source_component_mol_d=sources, free_acidity_eq_d=free_acidity, free_acidity_standard_uncertainty_eq_d=free_acidity_uncertainty, free_acidity_status=product.free_acidity_status, free_acidity_resolved=product.free_acidity_status != 'not_measured', intrinsic_hydrolysis_acidity_treatment='generated_by_shared_ferric_hydrolysis_reaction_not_duplicated_in_product_source')
__all__ = ['ProductDoseConversion', 'US_GALLON_TO_LITRE', 'coerce_commercial_product_definition', 'coerce_product_dose', 'convert_commercial_product_dose']
