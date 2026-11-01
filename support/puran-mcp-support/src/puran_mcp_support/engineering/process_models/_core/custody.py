from __future__ import annotations
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Iterable, Mapping, Protocol, Sequence
from .parameters import ParameterSet, ParameterTier
HOUSE_PROJECT_REF = 'house-process-library'
DESIGN_TIERS = frozenset({ParameterTier.HOUSE_DEFAULT, ParameterTier.HOUSE_RULED})
_DESIGN_REFUSAL_REASON = {ParameterTier.SITE_CALIBRATED: "one asset's calibration belongs to that site's twin; fleet evidence reaches design only as a family prior the owner has ruled (house_ruled)", ParameterTier.FAMILY_PRIOR_CANDIDATE: 'a candidate standard is put to the owner (priors.ruling_proposal), never read by a design: required is ours', ParameterTier.CATALOG: 'vendor catalog data is an input to a house default, not one', ParameterTier.VENDOR_PROJECTION_FIT: "a fit to a vendor's projection software is an input to a house default, not one, until the house adopts it"}

class ParameterCustodyError(ValueError):
    """A parameter set used, stored or asked for outside its custody."""

class TierRefusal(ParameterCustodyError):
    """A set of a tier the reader may not use."""

class DesignTierRefusal(TierRefusal):
    """A design offered anything but a house_default or house_ruled set."""

class NoParameterSet(LookupError):
    """No current set of the asked tier exists: the caller keeps its own (code) path and says so."""

def custody_scope(ps: ParameterSet) -> dict[str, str]:
    """The supersession key of a set, as the ``scope_refs`` of its row."""
    if not (ps.subject or '').strip():
        raise ParameterCustodyError(f'a stored {ps.model_key} set names its subject (the element, product or feed the values describe): it is part of the custody key')
    scope = {'model_key': ps.model_key, 'subject': ps.subject.strip(), 'tier': ps.tier.value}
    if ps.tier is ParameterTier.SITE_CALIBRATED:
        scope['asset_uid'] = ps.asset_uid
        scope['site'] = ps.site
    return scope

def project_ref_for(ps: ParameterSet) -> str:
    """A site set lives under its site's project; every other tier under the house project."""
    return ps.site if ps.tier is ParameterTier.SITE_CALIBRATED else HOUSE_PROJECT_REF

def persistence_boundary(content: Mapping[str, Any], *, project_ref: str, scope_refs: Mapping[str, Any] | None) -> tuple[dict[str, Any], dict[str, Any]]:
    """Validate a ``process_model_parameter_set`` row at the persistence boundary (both house writers call this).

    Returns the canonical content and the scope_refs carrying the custody key. Refuses a malformed set, a set filed
    under any project but its own, and scope_refs that contradict the content (they would scope supersession to
    somebody else's row)."""
    ps = ParameterSet.model_validate(dict(content))
    scope = custody_scope(ps)
    expected = project_ref_for(ps)
    if project_ref != expected:
        raise ParameterCustodyError(f'a {ps.tier.value} set of {ps.model_key} [{ps.subject}] lives under project {expected!r}, not {project_ref!r}')
    refs = dict(scope_refs or {})
    clash = {k: refs[k] for k in scope if k in refs and str(refs[k]) != scope[k]}
    if clash:
        raise ParameterCustodyError(f"scope_refs {clash} contradict the set's own custody key {scope}")
    refs.update(scope)
    return (ps.model_dump(mode='json'), refs)

@dataclass(frozen=True)
class StoredParameterSet:
    parameter_set: ParameterSet
    project_ref: str
    content_hash: str
    is_current: bool
    artifact_id: str | None = None
    created_at: datetime | None = None

class ParameterStore(Protocol):

    def candidates(self, *, model_key: str, tiers: Iterable[ParameterTier], subject: str | None=None, asset_uid: str | None=None, project_ref: str | None=None, include_history: bool=False) -> list[StoredParameterSet]:
        ...

    def put(self, ps: ParameterSet, *, audit_user_id: str, audit_reason: str) -> StoredParameterSet:
        ...

def _require_audit(audit_user_id: str, audit_reason: str) -> None:
    if not (audit_user_id or '').strip() or not (audit_reason or '').strip():
        raise ParameterCustodyError('a parameter-set write carries its audit context (audit_user_id, audit_reason)')

@dataclass
class MemoryParameterStore:
    """The same custody semantics in memory: tests, dry runs, and the code-resident register."""
    rows: list[StoredParameterSet] = field(default_factory=list)

    def put(self, ps: ParameterSet, *, audit_user_id: str='memory', audit_reason: str='memory store') -> StoredParameterSet:
        _require_audit(audit_user_id, audit_reason)
        project = project_ref_for(ps)
        content, scope = persistence_boundary(ps.model_dump(mode='json'), project_ref=project, scope_refs=None)
        self.rows = [StoredParameterSet(r.parameter_set, r.project_ref, r.content_hash, False, r.artifact_id, r.created_at) if r.is_current and r.project_ref == project and (custody_scope(r.parameter_set) == scope) else r for r in self.rows]
        stored = StoredParameterSet(ParameterSet.model_validate(content), project, ps.set_hash(), True, f'memory-{len(self.rows) + 1}', datetime.now(timezone.utc))
        self.rows.append(stored)
        return stored

    def candidates(self, *, model_key: str, tiers: Iterable[ParameterTier], subject: str | None=None, asset_uid: str | None=None, project_ref: str | None=None, include_history: bool=False) -> list[StoredParameterSet]:
        tiers = frozenset(tiers)
        out = [r for r in self.rows if r.parameter_set.model_key == model_key and r.parameter_set.tier in tiers and (subject is None or r.parameter_set.subject == subject) and (asset_uid is None or r.parameter_set.asset_uid == asset_uid) and (project_ref is None or r.project_ref == project_ref) and (include_history or r.is_current)]
        return sorted(out, key=lambda r: r.created_at or datetime.min.replace(tzinfo=timezone.utc), reverse=True)

def admit_for_design(ps: ParameterSet, *, model_key: str | None=None) -> ParameterSet:
    """A design reads house_default or house_ruled only."""
    if ps.tier not in DESIGN_TIERS:
        raise DesignTierRefusal(f'a design reads house_default or house_ruled parameters only; {ps.model_key} [{ps.subject}] is {ps.tier.value}: {_DESIGN_REFUSAL_REASON[ps.tier]}')
    if model_key is not None and ps.model_key != model_key:
        raise ParameterCustodyError(f'asked for {model_key}, offered {ps.model_key}')
    return ps

def _single(rows: Sequence[StoredParameterSet], what: str) -> StoredParameterSet | None:
    if len(rows) > 1:
        raise ParameterCustodyError(f'{len(rows)} current sets for {what}: a duplicate is_current to repair')
    return rows[0] if rows else None

def design_parameter_set(store: ParameterStore, *, model_key: str, subject: str) -> StoredParameterSet:
    """The set a design reads: the current house_ruled set, else the current house_default set.

    Raises :class:`NoParameterSet` when there is neither (the design keeps its code path), and
    :class:`DesignTierRefusal` when the store offers any other tier."""
    rows = store.candidates(model_key=model_key, tiers=DESIGN_TIERS, subject=subject, project_ref=HOUSE_PROJECT_REF)
    for row in rows:
        admit_for_design(row.parameter_set, model_key=model_key)
    current = [r for r in rows if r.is_current and r.parameter_set.subject == subject]
    for tier in (ParameterTier.HOUSE_RULED, ParameterTier.HOUSE_DEFAULT):
        chosen = _single([r for r in current if r.parameter_set.tier is tier], f'{model_key} [{subject}] {tier.value}')
        if chosen is not None:
            return chosen
    raise NoParameterSet(f'no house_ruled or house_default set of {model_key} [{subject}]')
