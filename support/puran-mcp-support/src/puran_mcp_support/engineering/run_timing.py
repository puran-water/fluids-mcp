"""Timing spans for the proposal workflow's critical path: a record only (speed plan P0).

THE RULE (`TIMING_RULE`): a span is an observation. Nothing reads it to decide anything, it never changes a result, a
refusal or an exit code, and a failure to write it is swallowed (counted, never raised). It never records source
contents: names, ids, counts, hashes already computed, seconds and an exception's type with a short message only.

Where the lines go: nowhere, unless a sink is named. ``PURANOS_TIMING_LOG`` names a JSON-lines file (opened append-only,
one ``os.write`` per line so concurrent processes interleave whole lines), and ``PURANOS_TIMING_PARENT`` the span a
process's top-level spans hang from. `chain_run` sets both for every stage it launches (the file beside the attempt's
own records, OUTSIDE its ``out`` directory, so no output hash moves); `child_env` sets them for any launch. A standalone
command is timed by exporting ``PURANOS_TIMING_LOG`` yourself.

Events (``run_timing_event/v1``): ``start`` (span id, parent id, name, category, kind, pid, process id, thread, boot id,
monotonic and epoch start, attributes) and ``end`` (wall, process CPU, waited-children CPU, peak RSS of the process and
its children, outcome, attributes added while open). Launch accounting (`LAUNCH_RULE`) adds two: ``launch``, written by
`child_env` in the launching process just before the launch (its id, the spawning span, pid, instant), and ``process``,
written once per process before its first span (process id, pid, parent pid, the launch id it inherited, its parent
span, the kernel's process start instant). A span opened by a thread that carries no span of its own hangs from its
process's outermost open span (``parent_via: process_root``), never from the launching process's span; a thread
started through `carried` runs in a copy of its starter's context, so its spans hang from the starter's (`THREAD_RULE`).
A span still open when its parent ends is ended with it (``closed_by_parent``), so every refusal that unwinds through a
span leaves an end event; a span whose process died leaves a start with no end, which the report lists as open, never
as finished.

Standard library only, importable from any interpreter in the repository (the vendor-projection venv carries its own
writer of the same schema, `vendor_projection.timing`, because it cannot import this package)."""
from __future__ import annotations
import contextlib
import contextvars
import functools
import itertools
import json
import os
import threading
import time
SCHEMA = 'run_timing_event/v1'
LOG_ENV = 'PURANOS_TIMING_LOG'
PARENT_ENV = 'PURANOS_TIMING_PARENT'
LAUNCH_ENV = 'PURANOS_TIMING_LAUNCH'
LOG_NAME = 'timing.jsonl'
TIMING_RULE = 'timing spans are a record only: nothing reads them to decide anything, they never change a result, a refusal or an exit code, a write failure is swallowed, and no source contents are recorded'
LAUNCH_RULE = 'a launch is one process, counted once: `child_env` writes a `launch` event (spawning span, instant) and the launched process a `process` event (launch id, parent pid, kernel start instant) before its first span; launch latency is spawn -> process start -> first instrumented span, only where both events are present and the parent pid matches. Every further span of that process is a worker span, never another launch. A log without them has no causal launch latency: its parent offset includes everything the launching span did before the launch and is labelled non-causal'
THREAD_RULE = "a worker thread started through `carried` runs in a copy of its starter's context: its spans (a lock wait included) and the processes it launches (`child_env`) hang from the starter's current span; a thread started bare hangs from its process's outermost span (parent_via process_root)"
CATEGORIES: dict[str, str] = {'a': 'walks/solves', 'b': 'RO option search', 'c': 'chain/render/recalculation', 'd': 'agent reading/reasoning/orchestration', 'e': 'shared models/fixes/integration', 'f': 'owner/stop/relaunch', 'g': 'self-inflicted regeneration', 'h': 'serialisation (queue and lock wait)'}
MESSAGE_CHARS = 160
_current: contextvars.ContextVar = contextvars.ContextVar('puranos_timing_span', default=None)
_lock = threading.Lock()
_sinks: dict[str, int] = {}
_override: list[str] = []
_boot_id: list[str] = []
_NONCE = os.urandom(4).hex()
_counter = itertools.count(1)
_announced: set = set()
_roots: dict[int, list] = {}
write_failures = 0

def _proc() -> str:
    """This process's identity in the log: the span-id prefix (random per interpreter, plus the pid, so a reused pid
    or a fork is a different process)."""
    return f'{_NONCE}{os.getpid():x}'

def _process_start() -> tuple[float | None, float | None, float | None]:
    """The kernel's start instant of this process on the monotonic and epoch clocks (``/proc/self/stat`` field 22,
    truncated to a clock tick) and that resolution, or Nones where it cannot be read."""
    try:
        with open('/proc/self/stat') as handle:
            ticks = int(handle.read().rsplit(')', 1)[1].split()[19])
        tick = 1.0 / os.sysconf('SC_CLK_TCK')
        age = time.clock_gettime(time.CLOCK_BOOTTIME) - ticks * tick
        if age < 0:
            return (None, None, None)
        return (round(time.monotonic() - age, 6), round(time.time() - age, 6), tick)
    except (OSError, ValueError, IndexError, AttributeError):
        return (None, None, None)

def _announce() -> None:
    """Write this process's `process` event once per sink, before its first span (`LAUNCH_RULE`)."""
    key = (os.getpid(), log_path())
    if key in _announced:
        return
    _announced.add(key)
    start_mono, start_epoch, tick = _process_start()
    _write({'schema': SCHEMA, 'event': 'process', 'proc': _proc(), 'pid': os.getpid(), 'ppid': os.getppid(), 'boot': _boot(), 'launch': os.environ.get(LAUNCH_ENV) or None, 'parent': os.environ.get(PARENT_ENV) or None, 'mono': round(time.monotonic(), 6), 'epoch': round(time.time(), 6), 'process_start_mono': start_mono, 'process_start_epoch': start_epoch, 'process_start_resolution_s': tick})

def _boot() -> str | None:
    if not _boot_id:
        try:
            with open('/proc/sys/kernel/random/boot_id') as handle:
                _boot_id.append(handle.read().strip())
        except OSError:
            _boot_id.append('')
    return _boot_id[0] or None

def log_path() -> str | None:
    """The sink this process writes to now, or None (timing off)."""
    if _override:
        return _override[-1]
    return os.environ.get(LOG_ENV) or None

def enabled() -> bool:
    return bool(_override) or bool(os.environ.get(LOG_ENV))

def _write(event: dict) -> None:
    global write_failures
    path = log_path()
    if not path:
        return
    try:
        line = (json.dumps(event, separators=(',', ':'), default=str, allow_nan=False) + '\n').encode()
        with _lock:
            fd = _sinks.get(path)
            if fd is None:
                directory = os.path.dirname(path)
                if directory:
                    os.makedirs(directory, exist_ok=True)
                fd = _sinks[path] = os.open(path, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 420)
            os.write(fd, line)
    except (OSError, ValueError, TypeError):
        write_failures += 1

@contextlib.contextmanager
def use_log(path):
    """Write this process's spans to ``path`` while the block runs (the environment's sink otherwise)."""
    _override.append(str(path))
    try:
        yield str(path)
    finally:
        _override.remove(str(path))

def _usage() -> tuple[float, float, int | None]:
    times = os.times()
    try:
        import resource
        rss = max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    except (ImportError, OSError, ValueError):
        rss = None
    return (times.user + times.system, times.children_user + times.children_system, rss)

class Span:
    """One open span. `annotate` adds attributes to its end event; `end` closes it (any open child first)."""
    __slots__ = ('id', 'parent_id', 'name', 'category', 'kind', 'attrs', '_mono', '_cpu', '_child_cpu', '_children', '_parent', '_pid', 'closed')

    def __init__(self, name, category, kind, parent, attrs, root=None):
        _announce()
        self._pid = os.getpid()
        self.id = f'{_proc()}-{next(_counter):x}'
        self._parent = parent
        via = None
        if parent is not None:
            self.parent_id = parent.id
        elif root is not None:
            self.parent_id, via = (root.id, 'process_root')
        else:
            self.parent_id = os.environ.get(PARENT_ENV) or None
            via = 'environment' if self.parent_id else None
            if threading.current_thread() is threading.main_thread():
                _roots.setdefault(self._pid, []).append(self)
        self.name = name
        self.category = category if category is not None else parent.category if parent is not None else None
        self.kind = kind
        self.attrs: dict = {}
        self._children: list[Span] = []
        self.closed = False
        self._cpu, self._child_cpu, _ = _usage()
        self._mono = time.monotonic()
        thread = threading.current_thread()
        _write({'schema': SCHEMA, 'event': 'start', 'span': self.id, 'parent': self.parent_id, 'name': name, 'category': self.category, 'kind': kind, 'pid': self._pid, 'proc': _proc(), 'boot': _boot(), 'mono': round(self._mono, 6), 'epoch': round(time.time(), 6), 'attrs': attrs or {}, **({'parent_via': via} if via else {}), **({'thread': thread.name} if thread is not threading.main_thread() else {})})
        if parent is not None:
            parent._children.append(self)

    def annotate(self, **attrs) -> None:
        self.attrs.update(attrs)

    def end(self, outcome: str='ok', *, exception: BaseException | None=None, **attrs) -> None:
        if self.closed:
            return
        for child in list(self._children):
            child.end('closed_by_parent', parent_outcome=outcome)
        self.closed = True
        mono = time.monotonic()
        cpu, child_cpu, rss = _usage()
        self.attrs.update(attrs)
        event = {'schema': SCHEMA, 'event': 'end', 'span': self.id, 'mono': round(mono, 6), 'wall_s': round(mono - self._mono, 6), 'cpu_s': round(cpu - self._cpu, 6), 'child_cpu_s': round(child_cpu - self._child_cpu, 6), 'peak_rss_kb': rss, 'outcome': self.attrs.pop('outcome', outcome), 'attrs': self.attrs}
        if exception is not None:
            event['exception'] = {'type': type(exception).__name__, 'message': str(exception)[:MESSAGE_CHARS]}
        _write(event)
        if self._parent is not None and self in self._parent._children:
            self._parent._children.remove(self)
        roots = _roots.get(self._pid)
        if roots and self in roots:
            roots.remove(self)
        if _current.get() is self:
            _current.set(self._parent if self._parent is not None and (not self._parent.closed) else None)

class _NullSpan:
    """What a span is when timing is off: every call does nothing."""
    id = None
    closed = True

    def annotate(self, **attrs) -> None:
        pass

    def end(self, outcome: str='ok', **attrs) -> None:
        pass
NULL_SPAN = _NullSpan()

def _outcome_of(exc: BaseException) -> str:
    if isinstance(exc, SystemExit):
        return 'ok' if exc.code in (0, None) else 'refused'
    if isinstance(exc, (KeyboardInterrupt, GeneratorExit)):
        return 'interrupted'
    return 'raised'

def start(name: str, *, category: str | None=None, kind: str | None=None, **attrs):
    """Open a span under the current one and make it current; close it with ``.end()``. `NULL_SPAN` when timing is off."""
    if not enabled():
        return NULL_SPAN
    parent = _current.get()
    if parent is not None and parent.closed:
        parent = None
    root = None
    if parent is None and threading.current_thread() is not threading.main_thread():
        root = next((r for r in _roots.get(os.getpid(), ()) if not r.closed), None)
    sp = Span(name, category, kind, parent, attrs, root)
    _current.set(sp)
    return sp

@contextlib.contextmanager
def span(name: str, *, category: str | None=None, kind: str | None=None, **attrs):
    """A span around a block. An exception ends it with its outcome and type, then propagates unchanged."""
    if not enabled():
        yield NULL_SPAN
        return
    sp = start(name, category=category, kind=kind, **attrs)
    try:
        yield sp
    except BaseException as exc:
        sp.end(_outcome_of(exc), exception=None if isinstance(exc, SystemExit) and exc.code in (0, None) else exc)
        raise
    sp.end()

def timed(name: str, *, category: str | None=None, kind: str | None=None):
    """Decorator form of `span`: costs one environment lookup per call when timing is off."""

    def wrap(function):

        @functools.wraps(function)
        def call(*args, **kwargs):
            if not enabled():
                return function(*args, **kwargs)
            with span(name, category=category, kind=kind):
                return function(*args, **kwargs)
        return call
    return wrap

def annotate(**attrs) -> None:
    """Add attributes to the end event of the current span (nothing when timing is off or no span is open)."""
    sp = _current.get()
    if sp is not None and (not sp.closed):
        sp.annotate(**attrs)

def carried(function):
    """``function`` bound to a copy of the context in force NOW, for a worker thread (`THREAD_RULE`). Each call runs in a
    fresh copy, so one bound function may run in several threads at once."""
    context = contextvars.copy_context()

    @functools.wraps(function)
    def run(*args, **kwargs):
        return context.copy().run(function, *args, **kwargs)
    return run

def current_id() -> str | None:
    sp = _current.get()
    if sp is not None and (not sp.closed):
        return sp.id
    return os.environ.get(PARENT_ENV) or None

class Steps:
    """A run of sibling spans with no block to hang them on (a pass loop): `next` ends the previous step and opens the
    next; `close` (or the enclosing span ending) ends the last."""

    def __init__(self, name: str, *, category: str | None=None, kind: str | None=None):
        self.name, self.category, self.kind = (name, category, kind)
        self._open = NULL_SPAN
        self._parent = _current.get()

    def next(self, **attrs):
        self._open.end()
        if self._parent is not None and _current.get() is not self._parent:
            _current.set(self._parent)
        self._open = start(self.name, category=self.category, kind=self.kind, **attrs)
        return self._open

    def close(self, outcome: str='ok') -> None:
        self._open.end(outcome)
        self._open = NULL_SPAN

class Phases:
    """Consecutive named phases of one block (a native model's build, scaling, initialisation, solve, extraction):
    `enter` ends the open phase and opens the next as its sibling, named ``<prefix>.<phase>`` with kind
    ``<prefix>_<phase>`` so `timing_report` splits them; `close` ends the last. A phase still open when the block
    raises is ended with its parent (``closed_by_parent``). Nothing at all when timing is off."""

    def __init__(self, prefix: str, *, category: str | None=None):
        self.prefix, self.category = (prefix, category)
        self._open = NULL_SPAN
        self._parent = _current.get()

    def enter(self, phase: str, /, **attrs):
        global write_failures
        try:
            self._open.end()
            if self._parent is not None and _current.get() is not self._parent:
                _current.set(self._parent)
            self._open = start(f'{self.prefix}.{phase}', category=self.category, kind=f'{self.prefix}_{phase}', **attrs)
        except Exception:
            write_failures += 1
            self._open = NULL_SPAN
        return self._open

    def close(self, outcome: str='ok') -> None:
        global write_failures
        try:
            self._open.end(outcome)
        except Exception:
            write_failures += 1
        self._open = NULL_SPAN

def child_env(env: dict | None=None, *, log=None) -> dict | None:
    """The environment a launched process inherits so its spans hang from the current one and go to the same sink
    (``log``: another file for them), with a `launch` event written here first (`LAUNCH_RULE`). Call it once per
    launch, just before it.

    With timing off, returns ``env`` unchanged (None stays None, so ``subprocess`` inherits exactly as before)."""
    path = log_path()
    if not path:
        return env
    out = dict(os.environ if env is None else env)
    out[LOG_ENV] = str(log) if log is not None else path
    parent = current_id()
    if parent:
        out[PARENT_ENV] = parent
    launch = f'{_proc()}-L{next(_counter):x}'
    _announce()
    _write({'schema': SCHEMA, 'event': 'launch', 'launch': launch, 'span': parent, 'pid': os.getpid(), 'proc': _proc(), 'boot': _boot(), 'mono': round(time.monotonic(), 6), 'epoch': round(time.time(), 6)})
    out[LAUNCH_ENV] = launch
    return out
