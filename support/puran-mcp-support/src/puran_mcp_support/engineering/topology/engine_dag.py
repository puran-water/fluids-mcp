from __future__ import annotations
from collections import defaultdict, deque
from dataclasses import dataclass
from typing import Any, Mapping, Sequence
from puran_mcp_support.engineering.models.simulation_seed import EngineKey, PortRole, SimulationSeed, SimulationSeedConnection

class EngineTopologyError(ValueError):
    """Raised when Stage-7b engine topology is invalid."""

@dataclass(frozen=True)
class ExternalStreamSpec:
    stream_id: str
    role: PortRole
    composition_hint: Mapping[str, Any] | None = None

@dataclass(frozen=True)
class ResolvedSeedConnection:
    source_seed_id: str
    sink_seed_id: str
    stream_id: str
    downstream_inlet_role: PortRole
    tear_stream: bool
    tear_rationale: str | None
    explicit: bool

@dataclass(frozen=True)
class EngineNode:
    engine_key: EngineKey
    seed_ids: tuple[str, ...]
    upstream_engines: frozenset[str]
    downstream_engines: frozenset[str]
    terminal_outlet_streams: tuple[str, ...]

@dataclass(frozen=True)
class TearEdge:
    source_seed_id: str
    sink_seed_id: str
    stream_id: str
    tear_rationale: str

@dataclass(frozen=True)
class EngineDAG:
    nodes: tuple[EngineNode, ...]
    layers: tuple[tuple[EngineNode, ...], ...]
    edges: tuple[tuple[str, str, str], ...]
    tear_edges: tuple[TearEdge, ...]

@dataclass(frozen=True)
class _InterEngineEdge:
    source_node_id: str
    sink_node_id: str
    source_engine: str
    sink_engine: str
    source_seed_id: str
    sink_seed_id: str
    stream_id: str
    tear_stream: bool
    tear_rationale: str | None

def _coerce_seeds(seeds: Sequence[SimulationSeed | Mapping[str, Any]]) -> list[SimulationSeed]:
    parsed = [seed if isinstance(seed, SimulationSeed) else SimulationSeed.model_validate(seed) for seed in seeds]
    seen: set[str] = set()
    for seed in parsed:
        if seed.seed_id in seen:
            raise EngineTopologyError(f'E9001 duplicate simulation seed_id {seed.seed_id!r}.')
        seen.add(seed.seed_id)
    return parsed

def _matching_inlets(*, stream_id: str, source_seed_id: str, inlets_by_stream: Mapping[str, Sequence[tuple[SimulationSeed, Any]]]) -> list[tuple[SimulationSeed, Any]]:
    return [(seed, inlet) for seed, inlet in inlets_by_stream.get(stream_id, ()) if seed.seed_id != source_seed_id]

def _explicit_target(*, source_seed: SimulationSeed, stream_id: str, connection: SimulationSeedConnection, seed_by_id: Mapping[str, SimulationSeed]) -> tuple[SimulationSeed, Any]:
    downstream = seed_by_id.get(connection.downstream_seed_id)
    if downstream is None:
        raise EngineTopologyError(f'E9002 simulation seed {source_seed.seed_id!r} outlet {stream_id!r} targets unknown downstream seed {connection.downstream_seed_id!r}.')
    stream_matches = [inlet for inlet in downstream.stream_topology.inlets if inlet.stream_id == stream_id]
    if len(stream_matches) == 1:
        return (downstream, stream_matches[0])
    role_matches = [inlet for inlet in downstream.stream_topology.inlets if inlet.role == connection.downstream_inlet_role]
    if len(role_matches) != 1:
        raise EngineTopologyError(f'E9002 simulation seed {source_seed.seed_id!r} outlet {stream_id!r} targets {downstream.seed_id!r}/{connection.downstream_inlet_role!r}, but that downstream inlet is not uniquely present.')
    return (downstream, role_matches[0])

def resolve_seed_connections(seeds: Sequence[SimulationSeed | Mapping[str, Any]]) -> tuple[ResolvedSeedConnection, ...]:
    """Resolve explicit and stream-id stitched seed-to-seed outlet connections."""
    parsed = _coerce_seeds(seeds)
    seed_by_id = {seed.seed_id: seed for seed in parsed}
    inlets_by_stream: dict[str, list[tuple[SimulationSeed, Any]]] = defaultdict(list)
    outlets_by_stream: dict[str, list[tuple[SimulationSeed, Any]]] = defaultdict(list)
    for seed in parsed:
        for inlet in seed.stream_topology.inlets:
            inlets_by_stream[inlet.stream_id].append((seed, inlet))
        for outlet in seed.stream_topology.outlets:
            outlets_by_stream[outlet.stream_id].append((seed, outlet))
    duplicate_outlets = {stream_id: [seed.seed_id for seed, _ in rows] for stream_id, rows in outlets_by_stream.items() if len(rows) > 1}
    if duplicate_outlets:
        labels = ', '.join((f'{stream_id} -> {seed_ids}' for stream_id, seed_ids in sorted(duplicate_outlets.items())))
        raise EngineTopologyError(f'E9003 duplicate outlet stream_id values in simulation seeds: {labels}.')
    resolved: list[ResolvedSeedConnection] = []
    for source_seed in parsed:
        for outlet in source_seed.stream_topology.outlets:
            connection = outlet.connects_to
            if connection is not None:
                downstream, inlet = _explicit_target(source_seed=source_seed, stream_id=outlet.stream_id, connection=connection, seed_by_id=seed_by_id)
                resolved.append(ResolvedSeedConnection(source_seed_id=source_seed.seed_id, sink_seed_id=downstream.seed_id, stream_id=outlet.stream_id, downstream_inlet_role=inlet.role, tear_stream=connection.tear_stream, tear_rationale=connection.tear_rationale, explicit=True))
                continue
            matches = _matching_inlets(stream_id=outlet.stream_id, source_seed_id=source_seed.seed_id, inlets_by_stream=inlets_by_stream)
            if not matches:
                continue
            if len(matches) > 1:
                labels = ', '.join((f'{seed.seed_id}/{inlet.role}' for seed, inlet in matches))
                raise EngineTopologyError(f'E9003 ambiguous outlet stream_id {outlet.stream_id!r} from seed {source_seed.seed_id!r}; matched multiple downstream inlets: {labels}.')
            downstream, inlet = matches[0]
            resolved.append(ResolvedSeedConnection(source_seed_id=source_seed.seed_id, sink_seed_id=downstream.seed_id, stream_id=outlet.stream_id, downstream_inlet_role=inlet.role, tear_stream=False, tear_rationale=None, explicit=False))
    return tuple(resolved)

def _external_stream_ids(external_streams: Mapping[str, ExternalStreamSpec] | None) -> set[str]:
    if external_streams is None:
        return set()
    ids: set[str] = set()
    for key, spec in external_streams.items():
        if key != spec.stream_id:
            raise EngineTopologyError(f'E9001 external_streams keys must match ExternalStreamSpec.stream_id: {key!r} != {spec.stream_id!r}.')
        ids.add(spec.stream_id)
    return ids

def _validate_unfed_inlets(*, seeds: Sequence[SimulationSeed], external_streams: Mapping[str, ExternalStreamSpec] | None) -> None:
    if external_streams is None:
        return
    produced = {outlet.stream_id for seed in seeds for outlet in seed.stream_topology.outlets}
    external_ids = _external_stream_ids(external_streams)
    missing = sorted({inlet.stream_id for seed in seeds for inlet in seed.stream_topology.inlets if inlet.stream_id not in produced and inlet.stream_id not in external_ids})
    if missing:
        raise EngineTopologyError(f'E9001 simulation seed inlets are not produced by another seed and are not declared external streams: {missing}.')

def _seed_components(seeds: Sequence[SimulationSeed], connections: Sequence[ResolvedSeedConnection]) -> list[tuple[str, ...]]:
    adjacency: dict[str, set[str]] = {seed.seed_id: set() for seed in seeds}
    for edge in connections:
        if edge.source_seed_id in adjacency and edge.sink_seed_id in adjacency:
            adjacency[edge.source_seed_id].add(edge.sink_seed_id)
            adjacency[edge.sink_seed_id].add(edge.source_seed_id)
    seen: set[str] = set()
    components: list[tuple[str, ...]] = []
    for seed_id in sorted(adjacency):
        if seed_id in seen:
            continue
        queue = deque([seed_id])
        seen.add(seed_id)
        current: list[str] = []
        while queue:
            node = queue.popleft()
            current.append(node)
            for neighbor in sorted(adjacency[node]):
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
        components.append(tuple(sorted(current)))
    return sorted(components, key=lambda item: item[0])

def _node_sort_key(node: EngineNode) -> tuple[str, tuple[str, ...]]:
    return (str(node.engine_key), node.seed_ids)

def _path_between(*, start: str, target: str, adjacency: Mapping[str, Sequence[_InterEngineEdge]]) -> list[_InterEngineEdge] | None:
    stack: list[tuple[str, list[_InterEngineEdge]]] = [(start, [])]
    seen: set[str] = set()
    while stack:
        node, path = stack.pop()
        if node == target:
            return path
        if node in seen:
            continue
        seen.add(node)
        for edge in sorted(adjacency.get(node, ()), key=lambda e: (e.sink_engine, e.stream_id, e.sink_seed_id), reverse=True):
            stack.append((edge.sink_node_id, [*path, edge]))
    return None

def _cycle_edges(edges: Sequence[_InterEngineEdge]) -> dict[_InterEngineEdge, list[_InterEngineEdge]]:
    adjacency: dict[str, list[_InterEngineEdge]] = defaultdict(list)
    for edge in edges:
        adjacency[edge.source_node_id].append(edge)
    out: dict[_InterEngineEdge, list[_InterEngineEdge]] = {}
    for edge in edges:
        path = _path_between(start=edge.sink_node_id, target=edge.source_node_id, adjacency=adjacency)
        if path is not None:
            out[edge] = [edge, *path]
    return out

def _raise_for_implicit_cycles(edges: Sequence[_InterEngineEdge]) -> None:
    participating = _cycle_edges(edges)
    implicit = [edge for edge in participating if not edge.tear_stream]
    if not implicit:
        return
    edge = sorted(implicit, key=lambda e: (e.source_engine, e.sink_engine, e.stream_id))[0]
    cycle = participating[edge]
    names = [cycle[0].source_engine]
    names.extend((item.sink_engine for item in cycle))
    raise EngineTopologyError(f'E9004 implicit inter-engine cycle requires every participating edge to declare tear_stream=True with tear_rationale: {' -> '.join(names)} via stream {edge.stream_id!r}.')

def _topological_layers(*, node_ids: Sequence[str], nodes_by_id: Mapping[str, EngineNode], edges: Sequence[_InterEngineEdge]) -> tuple[tuple[EngineNode, ...], ...]:
    adjacency: dict[str, set[str]] = {node_id: set() for node_id in node_ids}
    indegree: dict[str, int] = {node_id: 0 for node_id in node_ids}
    for edge in edges:
        if edge.sink_node_id not in adjacency[edge.source_node_id]:
            adjacency[edge.source_node_id].add(edge.sink_node_id)
            indegree[edge.sink_node_id] += 1
    ready = sorted((node_id for node_id, degree in indegree.items() if degree == 0), key=lambda nid: _node_sort_key(nodes_by_id[nid]))
    layers: list[tuple[EngineNode, ...]] = []
    emitted = 0
    while ready:
        current = ready
        emitted += len(current)
        layers.append(tuple((nodes_by_id[node_id] for node_id in current)))
        next_ready: list[str] = []
        for node_id in current:
            for downstream in sorted(adjacency[node_id], key=lambda nid: _node_sort_key(nodes_by_id[nid])):
                indegree[downstream] -= 1
                if indegree[downstream] == 0:
                    next_ready.append(downstream)
        ready = sorted(next_ready, key=lambda nid: _node_sort_key(nodes_by_id[nid]))
    if emitted != len(node_ids):
        remaining = sorted((nodes_by_id[node_id].engine_key for node_id, degree in indegree.items() if degree > 0))
        raise EngineTopologyError(f'E9005 inter-engine graph remains cyclic after tear removal: {remaining}.')
    return tuple(layers)

def build_engine_dag(seeds: list[SimulationSeed], *, external_streams: dict[str, ExternalStreamSpec] | None=None) -> EngineDAG:
    parsed = _coerce_seeds(seeds)
    _validate_unfed_inlets(seeds=parsed, external_streams=external_streams)
    engine_seeds = [seed for seed in parsed if seed.engine_key != 'none']
    if not engine_seeds:
        return EngineDAG(nodes=(), layers=(), edges=(), tear_edges=())
    engine_seed_ids = {seed.seed_id for seed in engine_seeds}
    seed_by_id = {seed.seed_id: seed for seed in engine_seeds}
    connections = resolve_seed_connections(parsed)
    engine_connections = [edge for edge in connections if edge.source_seed_id in engine_seed_ids and edge.sink_seed_id in engine_seed_ids]
    components = _seed_components(engine_seeds, engine_connections)
    seed_to_component: dict[str, int] = {}
    for idx, component in enumerate(components):
        for seed_id in component:
            seed_to_component[seed_id] = idx
    seed_to_node_id: dict[str, str] = {}
    grouped: dict[str, list[str]] = defaultdict(list)
    for seed in engine_seeds:
        node_id = f'{seed.engine_key}::{seed_to_component[seed.seed_id]}'
        seed_to_node_id[seed.seed_id] = node_id
        grouped[node_id].append(seed.seed_id)
    raw_inter_edges: list[_InterEngineEdge] = []
    for edge in engine_connections:
        source_node = seed_to_node_id[edge.source_seed_id]
        sink_node = seed_to_node_id[edge.sink_seed_id]
        if source_node == sink_node:
            continue
        source_seed = seed_by_id[edge.source_seed_id]
        sink_seed = seed_by_id[edge.sink_seed_id]
        raw_inter_edges.append(_InterEngineEdge(source_node_id=source_node, sink_node_id=sink_node, source_engine=str(source_seed.engine_key), sink_engine=str(sink_seed.engine_key), source_seed_id=edge.source_seed_id, sink_seed_id=edge.sink_seed_id, stream_id=edge.stream_id, tear_stream=edge.tear_stream, tear_rationale=edge.tear_rationale))
    deduped_edges = tuple({(edge.source_node_id, edge.sink_node_id, edge.stream_id, edge.source_seed_id, edge.sink_seed_id): edge for edge in raw_inter_edges}.values())
    _raise_for_implicit_cycles(deduped_edges)
    tear_edges = tuple((TearEdge(source_seed_id=edge.source_seed_id, sink_seed_id=edge.sink_seed_id, stream_id=edge.stream_id, tear_rationale=edge.tear_rationale or '') for edge in sorted(deduped_edges, key=lambda e: (e.source_engine, e.sink_engine, e.stream_id)) if edge.tear_stream))
    dag_edges = tuple((edge for edge in deduped_edges if not edge.tear_stream))
    downstream_by_node: dict[str, set[str]] = defaultdict(set)
    upstream_by_node: dict[str, set[str]] = defaultdict(set)
    outgoing_streams_by_node: dict[str, set[str]] = defaultdict(set)
    consumed_inter_streams = {edge.stream_id for edge in engine_connections}
    for edge in deduped_edges:
        outgoing_streams_by_node[edge.source_node_id].add(edge.stream_id)
    for edge in dag_edges:
        downstream_by_node[edge.source_node_id].add(edge.sink_engine)
        upstream_by_node[edge.sink_node_id].add(edge.source_engine)
    for seed in engine_seeds:
        node_id = seed_to_node_id[seed.seed_id]
        for outlet in seed.stream_topology.outlets:
            if outlet.stream_id not in consumed_inter_streams:
                outgoing_streams_by_node[node_id].add(outlet.stream_id)
    nodes_by_id: dict[str, EngineNode] = {}
    for node_id, seed_ids in grouped.items():
        first_seed = seed_by_id[sorted(seed_ids)[0]]
        nodes_by_id[node_id] = EngineNode(engine_key=first_seed.engine_key, seed_ids=tuple(sorted(seed_ids)), upstream_engines=frozenset(sorted(upstream_by_node[node_id])), downstream_engines=frozenset(sorted(downstream_by_node[node_id])), terminal_outlet_streams=tuple(sorted(outgoing_streams_by_node[node_id])))
    sorted_node_ids = sorted(nodes_by_id, key=lambda nid: _node_sort_key(nodes_by_id[nid]))
    layers = _topological_layers(node_ids=sorted_node_ids, nodes_by_id=nodes_by_id, edges=dag_edges)
    public_edges = tuple(sorted({(edge.source_engine, edge.sink_engine, edge.stream_id) for edge in dag_edges}))
    return EngineDAG(nodes=tuple((nodes_by_id[node_id] for node_id in sorted_node_ids)), layers=layers, edges=public_edges, tear_edges=tear_edges)
