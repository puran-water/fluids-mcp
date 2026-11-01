from __future__ import annotations
import math
import os
import re
import threading
from typing import Dict, Iterable, Optional, Tuple
from .phreeqc_phases import _KEYWORDS
RULING_ID = 'PUBLIC-RULE-e5ab57ace3c5'
NO_RUN_RULING_ID = 'PUBLIC-RULE-5164c01473a4'
OSMOTIC_MODEL_KEYWORDS = ('PITZER', 'SIT')
WATER_ACTIVITY_MOLALITY_COEFFICIENT = 0.017
REBUILD_SPECIES_FLOOR_MOL_KGW = 1e-12
REBUILD_RELATIVE_TOLERANCE = 1e-05
_REACTION_CACHE: Dict[Tuple[str, float, int], Dict[str, Tuple[str, float]]] = {}
_MODEL_CACHE: Dict[Tuple[str, float, int], str] = {}
_PRIVATE: Dict[Tuple[str, float, int, str], Tuple[object, threading.Lock]] = {}
_PRIVATE_LOCK = threading.Lock()
_TERM = re.compile('^\\s*((?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][-+]?\\d+)?)?\\s*(\\S+)\\s*$')

class SaturationIndexError(ValueError):
    """A saturation index that cannot be put on PHREEQC's own definition: refused by name, never the uncorrected one."""

def _water_in(side: str, *, skip_first: bool) -> float:
    """The water coefficient on one side of a reaction: terms joined by `` + `` or `` - `` (PHREEQC writes a negative
    coefficient as ``- H2O``, e.g. akermanite's ``... + 2 H4SiO4 - H2O``); a species' own charge (``Ca+2``) carries no
    spaces."""
    total = 0.0
    parts = re.split('\\s+([+-])\\s+', ' ' + side.strip())
    sign, terms = (1.0, [])
    first = parts[0].strip()
    if first:
        terms.append((1.0, first))
    for op, term in zip(parts[1::2], parts[2::2]):
        terms.append((1.0 if op == '+' else -1.0, term.strip()))
    for i, (sign, term) in enumerate(terms):
        if skip_first and i == 0:
            continue
        m = _TERM.match(term)
        if m and m.group(2) == 'H2O':
            total += sign * (float(m.group(1)) if m.group(1) else 1.0)
    return total

def phase_water_coefficients(database_text: str) -> Dict[str, Tuple[str, float]]:
    """``{lower-case phase name: (its spelling, nu_w)}`` for every phase a PHREEQC database text defines: ``nu_w`` the
    net water coefficient of the phase's dissolution reaction, products positive (gypsum +2, amorphous silica -2,
    halite 0). A later definition of a name replaces an earlier one, as PHREEQC reads it."""
    out: Dict[str, Tuple[str, float]] = {}
    in_phases = False
    current: Optional[str] = None
    for raw in database_text.splitlines():
        line = raw.split('#', 1)[0].rstrip()
        if not line.strip():
            continue
        token = line.split()[0]
        if not line[:1].isspace() and token.upper() in _KEYWORDS and (token.upper() == token):
            in_phases = token.upper() == 'PHASES'
            current = None
            continue
        if not in_phases:
            continue
        if not line[:1].isspace() and '=' not in line and (not token.startswith('-')):
            current = token
            continue
        if current is not None and '=' in line and (not token.startswith('-')):
            left, right = line.split('=', 1)
            out[current.lower()] = (current, _water_in(right, skip_first=False) - _water_in(left, skip_first=True))
            current = None
    return out

def database_water_coefficients(database_path: str) -> Dict[str, Tuple[str, float]]:
    """`phase_water_coefficients` of a database file, cached by its path, mtime and size."""
    key = _file_key(database_path)
    if key not in _REACTION_CACHE:
        with open(database_path, 'r', errors='ignore') as handle:
            _REACTION_CACHE[key] = phase_water_coefficients(handle.read())
    return _REACTION_CACHE[key]

def _file_key(database_path: str) -> Tuple[str, float, int]:
    stat = os.stat(database_path)
    return (os.path.realpath(database_path), stat.st_mtime, stat.st_size)

def activity_model(database_path: str) -> str:
    """``"ion_association"``, ``"pitzer"`` or ``"sit"``: how PHREEQC computes the activity of water on this database
    file, read from the file itself (a ``PITZER`` or ``SIT`` block at the start of a line), cached by path, mtime and
    size."""
    key = _file_key(database_path)
    if key not in _MODEL_CACHE:
        model = 'ion_association'
        with open(database_path, 'r', errors='ignore') as handle:
            for raw in handle:
                code = raw.split('#', 1)[0]
                if code[:1].isspace() or not code.strip():
                    continue
                token = code.split()[0].upper()
                if token in OSMOTIC_MODEL_KEYWORDS:
                    model = token.lower()
                    break
        _MODEL_CACHE[key] = model
    return _MODEL_CACHE[key]

def _dumped_water_activity(ip, solution_number: int) -> float:
    """-ah2o of a stored solution by ``DUMP``: run ONLY on a private instance (a run replaces its selected output)."""
    ip.set_dump_string_on()
    try:
        ip.run_string(f'DUMP\n -solution {int(solution_number)}\nEND\n')
        dump = ip.get_dump_string()
    finally:
        ip.set_dump_string_off()
    if isinstance(dump, bytes):
        dump = dump.decode(errors='ignore')
    for line in str(dump).splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0] == '-ah2o':
            value = float(parts[1])
            if math.isfinite(value) and value > 0:
                return value
    raise SaturationIndexError(f'PHREEQC reports no water activity for solution {solution_number} ({RULING_ID})')

def _water_kg(ip, solution_number: int) -> float:
    water = float(ip.get_mass(int(solution_number)))
    if not (math.isfinite(water) and water > 0):
        raise SaturationIndexError(f'solution {solution_number} reports {water!r} kg of water: no water activity ({NO_RUN_RULING_ID})')
    return water

def _ion_association_water_activity(ip, solution_number: int) -> float:
    moles = ip.get_species_moles(int(solution_number))
    if not moles:
        raise SaturationIndexError(f'solution {solution_number} reports no aqueous species: no water activity ({NO_RUN_RULING_ID})')
    total = math.fsum((float(v) for k, v in moles.items() if k not in ('H2O', 'e-')))
    value = 1.0 - WATER_ACTIVITY_MOLALITY_COEFFICIENT * total / _water_kg(ip, solution_number)
    if not (math.isfinite(value) and value > 0):
        raise SaturationIndexError(f"solution {solution_number}'s aqueous species sum to {total!r} mol: PHREEQC's ion-association water activity is {value!r} ({NO_RUN_RULING_ID})")
    return value
SPECIATION_DEFINITION_KEYWORDS = frozenset({'SOLUTION_MASTER_SPECIES', 'SOLUTION_SPECIES', 'PITZER', 'SIT', 'NAMED_EXPRESSIONS', 'NAMED_ANALYTICAL_EXPRESSION', 'NAMED_ANALYTICAL_EXPRESSIONS', 'NAMED_LOG_K', 'LLNL_AQUEOUS_MODEL_PARAMETERS', 'LLNL_AQUEOUS_MODEL'})
_BLOCK_KEYWORDS = frozenset({'END', 'TITLE', 'COMMENT', 'DATABASE', 'SOLUTION', 'SOLUTION_SPREAD', 'SPREAD_SOLUTION', 'SOLUTION_MASTER_SPECIES', 'SOLUTION_SPECIES', 'PHASES', 'EXCHANGE_MASTER_SPECIES', 'EXCHANGE_SPECIES', 'SURFACE_MASTER_SPECIES', 'SURFACE_SPECIES', 'RATES', 'PITZER', 'SIT', 'LLNL_AQUEOUS_MODEL_PARAMETERS', 'LLNL_AQUEOUS_MODEL', 'NAMED_EXPRESSIONS', 'NAMED_ANALYTICAL_EXPRESSION', 'NAMED_ANALYTICAL_EXPRESSIONS', 'NAMED_LOG_K', 'ISOTOPES', 'ISOTOPE_ALPHAS', 'ISOTOPE_RATIOS', 'CALCULATE_VALUES', 'EQUILIBRIUM_PHASES', 'EQUILIBRIUM_PHASE', 'EQUILIBRIA', 'EQUILIBRIUM', 'PURE_PHASES', 'PURE', 'REACTION', 'REACTIONS', 'REACTION_TEMPERATURE', 'REACTION_TEMPERATURES', 'REACTION_PRESSURE', 'REACTION_PRESSURES', 'MIX', 'USE', 'SAVE', 'COPY', 'DELETE', 'DUMP', 'RUN_CELLS', 'EXCHANGE', 'SURFACE', 'GAS_PHASE', 'KINETICS', 'SOLID_SOLUTIONS', 'SOLID_SOLUTION', 'INVERSE_MODELING', 'TRANSPORT', 'ADVECTION', 'SELECTED_OUTPUT', 'SELECT_OUTPUT', 'SELECTED_OUT', 'SELECT_OUT', 'USER_PUNCH', 'USER_PRINT', 'USER_GRAPH', 'PRINT', 'KNOBS', 'INCREMENTAL_REACTIONS', 'INCREMENTAL', 'GAS_BINARY_PARAMETERS'})

def _speciation_definitions(extra_text: str) -> str:
    """The run-time blocks of ``extra_text`` that change a solution's speciation (`SPECIATION_DEFINITION_KEYWORDS`), in
    order; every other block (a SOLUTION, a reaction, PHASES, output) is left out."""
    out, keep = ([], False)
    for raw in (extra_text or '').split('\n'):
        code = raw.split('#', 1)[0].strip()
        if code:
            word = code.split()[0].upper()
            if word in _BLOCK_KEYWORDS or word.endswith(('_RAW', '_MODIFY')):
                keep = word in SPECIATION_DEFINITION_KEYWORDS
        if keep:
            out.append(raw)
    return '\n'.join(out)

def _private_instance(database_path: str, definitions: str=''):
    """A private PhreeqPython instance of the database file (plus a run's own speciation definitions), cached by the
    file's path, mtime and size and the definitions' text."""
    key = (*_file_key(database_path), definitions)
    with _PRIVATE_LOCK:
        if key not in _PRIVATE:
            from pathlib import Path
            from phreeqpython import PhreeqPython
            pp = PhreeqPython(database=os.path.basename(database_path), database_directory=Path(os.path.dirname(os.path.realpath(database_path))))
            if definitions:
                try:
                    pp.ip.run_string(definitions + '\nEND\n')
                except Exception as exc:
                    raise SaturationIndexError(f"the run's own speciation definitions do not load on {os.path.basename(database_path)}: {str(exc)[:300]} ({NO_RUN_RULING_ID})") from exc
            _PRIVATE[key] = (pp.ip, threading.Lock())
        return _PRIVATE[key]

def _rebuild_disagreement(caller: Dict[str, float], rebuilt: Dict[str, float]) -> Optional[str]:
    for species in sorted(set(caller) | set(rebuilt)):
        if species == 'H2O':
            continue
        a, b = (float(caller.get(species, 0.0)), float(rebuilt.get(species, 0.0)))
        if max(abs(a), abs(b)) <= REBUILD_SPECIES_FLOOR_MOL_KGW:
            continue
        if a == 0.0 or abs(b / a - 1.0) > REBUILD_RELATIVE_TOLERANCE:
            return f'{species} {a!r} against {b!r}'
    return None

def _rebuilt_water_activity(ip, solution_number: int, database_path: str, model: str, extra_text: str='') -> float:
    n = int(solution_number)
    water = _water_kg(ip, n)
    totals = ip.get_elements_totals(n)
    body = ''.join((f' {name} {float(moles) / water!r}\n' for name, moles in sorted(totals.items()) if float(moles) > 0))
    deck = f'SOLUTION 1\n units mol/kgw\n temp {float(ip.get_temperature(n))!r}\n pH {float(ip.get_ph(n))!r}\n pe {float(ip.get_pe(n))!r}\n -water {water!r}\n{body}END\n'
    private, lock = _private_instance(database_path, _speciation_definitions(extra_text))
    with lock:
        try:
            private.run_string(deck)
        except Exception as exc:
            raise SaturationIndexError(f'solution {n} on {os.path.basename(database_path)} ({model}) cannot be reproduced from its totals, pH, pe, temperature and water to read its water activity: {str(exc)[:300]} ({NO_RUN_RULING_ID})') from exc
        for label, read in (('molality', 'get_species_molalities'), ('activity', 'get_species_activities')):
            off = _rebuild_disagreement(getattr(ip, read)(n), getattr(private, read)(1))
            if off is not None:
                raise SaturationIndexError(f"solution {n} on {os.path.basename(database_path)} ({model}) cannot be reproduced from its totals, pH, pe, temperature and water (a pressure, or a run-time {model.upper()} parameter or species the database file does not carry): species {label} {off}, beyond {REBUILD_RELATIVE_TOLERANCE:g} relative; its water activity is not read at the file's parameters ({NO_RUN_RULING_ID})")
        return _dumped_water_activity(private, 1)

def water_activity(ip, solution_number: int, *, database_path: str, extra_text: str='') -> float:
    'Public engineering calculation interface.'
    model = activity_model(database_path)
    for line in _speciation_definitions(extra_text).split('\n'):
        code = line.split('#', 1)[0].strip()
        if code and code.split()[0].upper() in OSMOTIC_MODEL_KEYWORDS:
            model = code.split()[0].lower()
    if model == 'ion_association':
        return _ion_association_water_activity(ip, solution_number)
    return _rebuilt_water_activity(ip, solution_number, database_path, model, extra_text)

def _coefficients(database_path: str, extra_text: str) -> Dict[str, Tuple[str, float]]:
    table = dict(database_water_coefficients(database_path))
    if extra_text:
        table.update(phase_water_coefficients(extra_text))
    return table

def corrected_si(ip, solution_number: int, phase: str, *, database_path: str, extra_text: str='', log_water_activity: Optional[float]=None) -> float:
    """PHREEQC's ``SI(phase)`` of a stored solution: VIPhreeqc's ``get_si`` plus ``nu_w log10(a_w)``. ``database_path``
    is the database the instance loaded (`phreeqc_phases.phreeqpython_database`); ``extra_text`` any PHASES the caller
    defined at run time. A phase neither defines is refused (VIPhreeqc answers -999 for it without an error). A defined
    phase whose elements the solution does not carry has no SI: VIPhreeqc's own sentinel (-999) is returned unchanged,
    as every house caller reads it (``> -900``), never corrected into a number."""
    table = _coefficients(database_path, extra_text)
    hit = table.get(phase.lower())
    if hit is None:
        raise SaturationIndexError(f'{phase} is not a phase of {os.path.basename(database_path)}{(' or the run-time PHASES' if extra_text else '')}: no saturation index ({RULING_ID})')
    raw = float(ip.get_si(int(solution_number), hit[0]))
    if raw <= -998.0:
        return raw
    nu = hit[1]
    if nu == 0.0:
        return raw
    lw = math.log10(water_activity(ip, solution_number, database_path=database_path, extra_text=extra_text)) if log_water_activity is None else float(log_water_activity)
    return raw + nu * lw

def saturation_index(solution, phase: str, *, database_path: str, extra_text: str='') -> float:
    """`corrected_si` for a phreeqpython ``Solution`` (the replacement for ``Solution.si(phase)``)."""
    return corrected_si(solution.pp.ip, solution.number, phase, database_path=database_path, extra_text=extra_text)

def saturation_indices(solution, *, database_path: str, extra_text: str='', phases: Optional[Iterable[str]]=None) -> Dict[str, float]:
    """Every phase the database defines whose elements the phreeqpython ``Solution`` carries (the replacement for
    ``Solution.phases``), or the named ``phases``, on PHREEQC's own definition. The water activity is read once."""
    ip, number = (solution.pp.ip, solution.number)
    table = _coefficients(database_path, extra_text)
    lw = math.log10(water_activity(ip, number, database_path=database_path, extra_text=extra_text))
    if phases is None:
        out = {}
        for spelling, nu in table.values():
            value = float(ip.get_si(int(number), spelling))
            if value > -998.0:
                out[spelling] = value + nu * lw
        return out
    return {p: corrected_si(ip, number, p, database_path=database_path, extra_text=extra_text, log_water_activity=lw) for p in phases}
