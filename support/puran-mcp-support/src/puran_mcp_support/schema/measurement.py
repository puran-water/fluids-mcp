"""Port-specific, separation-aware measurement aggregation."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Mapping

@dataclass(frozen=True)
class MeasurementMapEntry:
    component_id: str
    phosphorus_pool: str
    kg_p_per_native_unit: float

@dataclass(frozen=True)
class PortSeparationProfile:
    profile_id: str
    component_transmission: Mapping[str, float]

    def transmission(self, component_id: str) -> float:
        if component_id not in self.component_transmission:
            raise KeyError(f'separation profile {self.profile_id!r} omits {component_id!r}')
        value = float(self.component_transmission[component_id])
        if not 0 <= value <= 1:
            raise ValueError(f'transmission for {component_id!r} must be in [0, 1]')
        return value

@dataclass(frozen=True)
class PhosphorusMeasurementAggregator:
    entries: tuple[MeasurementMapEntry, ...]

    def aggregate(self, native_state: Mapping[str, float], *, separation_profile: PortSeparationProfile) -> dict[str, float]:
        pools: dict[str, float] = {}
        for entry in self.entries:
            amount = float(native_state.get(entry.component_id, 0.0)) * entry.kg_p_per_native_unit * separation_profile.transmission(entry.component_id)
            pools[entry.phosphorus_pool] = pools.get(entry.phosphorus_pool, 0.0) + amount
        pools['total_p'] = sum(pools.values())
        pools['orthophosphate'] = pools.get('orthophosphate', 0.0)
        return pools
__all__ = ['MeasurementMapEntry', 'PhosphorusMeasurementAggregator', 'PortSeparationProfile']
