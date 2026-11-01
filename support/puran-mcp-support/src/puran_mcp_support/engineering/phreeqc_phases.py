'D3 of the ruling: "refuse a missing mineral unless the substitute is the same phase under another name". D5: the\nno-fallback rule extends to every PHREEQC caller in engineering_utils and the servers ("a dependency must never choose a\ndesign basis"). This module is the one declared table of same-phase names (water-chemistry-mcp\'s\n`waterchem.phase_choice` reads it from here) and the reader every caller uses to learn the phases a database file\ndefines. Standard library only; phreeqpython is imported inside the one function that needs it.\n\nA database defines a phase by name in its PHASES block, and PHREEQC matches the name without regard to case. A name a\ndatabase does not define is answered only by a member of the same `SAME_PHASE` entry that it does define; each entry\ncarries its evidence (read from the shipped databases 2026-10-16: enc-e/runs/phase_defs_alias_candidates.txt, or the\nmineral\'s nomenclature where no shipped database carries both names). `TWO_PHASES` lists the swaps earlier tables made\nthat are two phases, each proved by a shipped database defining both names with different data.\n\nAn unnamed ``PhreeqPython()`` loads phreeqpython\'s own default, vitens.dat (a 2015 derivative of phreeqc.dat whose\nCO2(g) differs from phreeqc.dat\'s by 0.005-0.010 log units: enc-e/runs/probe_vitens_vs_phreeqc.txt). Every caller names\nits database (`phreeqpython_database`).'
from __future__ import annotations
import os
from dataclasses import dataclass
from typing import Dict, Optional, Tuple
RULING_ID = 'PUBLIC-RULE-762de2e5f9c5'
PHREEQPYTHON_UNNAMED_DEFAULT = 'vitens.dat'

@dataclass(frozen=True)
class SamePhase:
    """Names under which databases define ONE phase, and the evidence that they are one."""
    names: Tuple[str, ...]
    phase: str
    evidence: str
SAME_PHASE: Tuple[SamePhase, ...] = (SamePhase(('Fe(OH)3(a)', 'Ferrihydrite'), 'amorphous ferric hydroxide (ferrihydrite)', 'minteq.dat defines Ferrihydrite with the reaction and log K that phreeqc.dat and wateq4f.dat give Fe(OH)3(a): Fe(OH)3 + 3 H+ = Fe+3 + 3 H2O, log K 4.891'), SamePhase(('Al(OH)3(a)', 'Al(OH)3(am)'), 'amorphous aluminium hydroxide', 'minteq.v4.dat defines Al(OH)3(am) with the reaction and log K that phreeqc.dat and wateq4f.dat give Al(OH)3(a): Al(OH)3 + 3 H+ = Al+3 + 3 H2O, log K 10.8'), SamePhase(('Dolomite', 'Dolomite-ord'), 'ordered dolomite, CaMg(CO3)2', "llnl.dat defines both names with one reaction and one log K (2.5135); wateq4f.dat's Dolomite is its ordered dolomite beside Dolomite(d)"), SamePhase(('Hydroxyapatite', 'Hydroxylapatite'), 'hydroxylapatite, Ca5(PO4)3OH', 'two spellings of one mineral (hydroxylapatite is the IMA-CNMNC name: Pasero et al. 2010, Eur. J. Mineral. 22, 163-179, doi 10.1127/0935-1221/2010/0022-2022); phreeqc.dat, wateq4f.dat, minteq.dat and sit.dat spell it Hydroxyapatite, llnl.dat and minteq.v4.dat Hydroxylapatite, and no shipped database defines both'), SamePhase(('Brushite', 'CaHPO4:2H2O'), 'brushite, CaHPO4.2H2O', 'the mineral name and the formula name of dicalcium phosphate dihydrate; llnl.dat and sit.dat define Brushite, minteq.v4.dat CaHPO4:2H2O, and no shipped database defines both'), SamePhase(('Brucite', 'Mg(OH)2'), 'brucite, Mg(OH)2', 'the mineral name and the formula name of crystalline magnesium hydroxide; no shipped database defines a phase named Mg(OH)2 (in phreeqc.dat it is only an aqueous species, and phreeqc.dat defines no brucite)'))
TWO_PHASES: Tuple[Tuple[str, str, str], ...] = (('Strengite', 'Vivianite', 'an Fe(III) phosphate and an Fe(II) phosphate; minteq.v4.dat defines both'), ('Gibbsite', 'Al(OH)3(a)', 'crystalline and amorphous Al(OH)3; phreeqc.dat defines both (log K 8.11 and 10.8)'), ('Gibbsite', 'Al(OH)3(am)', 'crystalline and amorphous Al(OH)3; minteq.v4.dat defines both (8.291 and 10.8)'), ('FeS(ppt)', 'Mackinawite', 'amorphous FeS and mackinawite; phreeqc.dat defines both (log K -3.915 and -4.648)'), ('SiO2(a)', 'SiO2(am)', 'minteq.dat defines both with different log K (-3.018 and -2.71)'), ('CaHPO4', 'CaHPO4:2H2O', 'monetite and brushite; minteq.v4.dat defines both (log K -19.275 and -18.995)'), ('Dolomite', 'Dolomite(d)', 'ordered and disordered dolomite; wateq4f.dat defines both (log K -17.09 and -16.54)'))
_KEYWORDS = {'SOLUTION_MASTER_SPECIES', 'SOLUTION_SPECIES', 'PHASES', 'EXCHANGE_MASTER_SPECIES', 'EXCHANGE_SPECIES', 'SURFACE_MASTER_SPECIES', 'SURFACE_SPECIES', 'RATES', 'END', 'PITZER', 'SIT', 'LLNL_AQUEOUS_MODEL_PARAMETERS', 'NAMED_EXPRESSIONS', 'ISOTOPES', 'ISOTOPE_ALPHAS', 'ISOTOPE_RATIOS', 'CALCULATE_VALUES', 'SPECIES', 'MASTER_SPECIES', 'SOLUTION_SPREAD', 'SOLUTION', 'EQUILIBRIUM_PHASES', 'KNOBS', 'DATABASE'}
_PHASES_CACHE: Dict[Tuple[str, float, int], Dict[str, str]] = {}

def phase_names(database_text: str) -> Dict[str, str]:
    """The phases a PHREEQC database text defines, keyed by lower-case name (PHREEQC's match), valued by its spelling.

    A phase is a line at column 0 inside a PHASES block whose first token is not an option (``-log_k``) and carries no
    reaction (``=``); an indented line is the phase's reaction or an option."""
    names: Dict[str, str] = {}
    in_phases = False
    for raw in database_text.splitlines():
        line = raw.split('#', 1)[0].rstrip()
        if not line.strip():
            continue
        token = line.split()[0]
        if not line[:1].isspace() and token.upper() in _KEYWORDS and (token.upper() == token):
            in_phases = token.upper() == 'PHASES'
            continue
        if in_phases and (not line[:1].isspace()) and ('=' not in line) and (not token.startswith('-')):
            names.setdefault(token.lower(), token)
    return names

def database_phase_names(database_path: str) -> Dict[str, str]:
    """`phase_names` of a database file, cached by its path, mtime and size."""
    stat = os.stat(database_path)
    key = (os.path.realpath(database_path), stat.st_mtime, stat.st_size)
    if key not in _PHASES_CACHE:
        with open(database_path, 'r', errors='ignore') as handle:
            _PHASES_CACHE[key] = phase_names(handle.read())
    return _PHASES_CACHE[key]

def same_phase_as(name: str) -> Optional[SamePhase]:
    """The `SAME_PHASE` entry a name belongs to, or None."""
    lowered = name.lower()
    return next((entry for entry in SAME_PHASE if lowered in {member.lower() for member in entry.names}), None)

def resolve_phase(name: str, defined: Dict[str, str]) -> Tuple[Optional[str], Optional[SamePhase]]:
    """``(the name the database defines for this phase, the SAME_PHASE entry used or None)``; ``(None, None)`` when the
    database defines neither the name nor a same-phase name. ``defined`` is `phase_names` of the database."""
    if name.lower() in defined:
        return (defined[name.lower()], None)
    entry = same_phase_as(name)
    for member in entry.names if entry else ():
        if member.lower() in defined:
            return (defined[member.lower()], entry)
    return (None, None)

class MissingPhaseError(ValueError):
    """A mineral the database does not define, with no same-phase name there: refused, never replaced or dropped."""

def require_phases(names, database_path: str, *, purpose: str) -> Dict[str, str]:
    """``{name asked: the name the database defines}`` for every mineral a model needs, or `MissingPhaseError` naming
    each one the database does not define under its own or a same-phase name (D3 of the ruling, applied to every
    PHREEQC caller by D5). phreeqpython's ``si()`` of an undefined phase returns -999 without an error, so a model that
    only reads SI would skip the mineral silently."""
    defined = database_phase_names(database_path)
    chosen: Dict[str, str] = {}
    missing = []
    for name in names:
        hit, _entry = resolve_phase(name, defined)
        if hit is None:
            missing.append(name)
        else:
            chosen[name] = hit
    if missing:
        raise MissingPhaseError(f'{', '.join(missing)} {('is' if len(missing) == 1 else 'are')} not a phase of {os.path.basename(database_path)} under its own name or a same-phase name, and {purpose} needs {('it' if len(missing) == 1 else 'them')}: refused rather than skipped ({RULING_ID}, D3/D5). Name a database that defines it.')
    return chosen

def phreeqpython_database(database: str) -> Tuple[Dict[str, object], str]:
    """``(PhreeqPython keyword arguments, the database file)`` for a named database: a bare name is phreeqpython's own
    file of that name (the file the walk's ``PhreeqcModel(database="phreeqc.dat")`` loads); an absolute path is that
    file; a relative path with a folder refuses (it would be read against the working directory). An
    empty name refuses: an unnamed ``PhreeqPython()`` would load phreeqpython's default, vitens.dat. A file that does
    not exist refuses by name."""
    if not database or not str(database).strip():
        raise ValueError(f"No PHREEQC database named: an unnamed PhreeqPython() loads phreeqpython's default {PHREEQPYTHON_UNNAMED_DEFAULT}, so every caller names one ({RULING_ID}, D5).")
    database = str(database)
    if os.path.dirname(database):
        if not os.path.isabs(database):
            raise ValueError(f"PHREEQC database path {database!r} is relative: it would be read against the working directory, so a path must be absolute (a bare name is phreeqpython's own file) ({RULING_ID}, D5).")
        path = database
        from pathlib import Path
        kwargs: Dict[str, object] = {'database': os.path.basename(path), 'database_directory': Path(os.path.dirname(path))}
    else:
        import phreeqpython
        path = os.path.join(os.path.dirname(phreeqpython.__file__), 'database', database)
        kwargs = {'database': database}
    if not os.path.isfile(path):
        raise FileNotFoundError(f'PHREEQC database {database!r} does not exist ({path}); nothing loads in its place ({RULING_ID}, D5).')
    return (kwargs, path)
HOUSE_IRON_SULFIDE_PHASE = 'Mackinawite'
AMORPHOUS_IRON_SULFIDE_PHASE = 'FeS(ppt)'
IRON_SULFIDE_RULING_ID = 'PUBLIC-RULE-1ddf4c4cc59b'
HOUSE_DATABASE = 'phreeqc.dat'
HOUSE_DATABASE_RULING_ID = 'PUBLIC-RULE-8226e460e1d8'

def database_file_record(path: str) -> Dict[str, str]:
    'Public engineering calculation interface.'
    import hashlib
    real = os.path.realpath(path)
    cached = _FILE_RECORDS.get(real)
    if cached is None:
        with open(real, 'rb') as handle:
            cached = _FILE_RECORDS[real] = {'database': os.path.basename(real), 'sha256': hashlib.sha256(handle.read()).hexdigest()}
    return dict(cached)
_FILE_RECORDS: Dict[str, Dict[str, str]] = {}

def house_database() -> Dict[str, object]:
    """``{"database", "path", "sha256", "kwargs", "ruling"}`` of the house phreeqc.dat (`HOUSE_DATABASE`): the file a
    named ``phreeqpython_database("phreeqc.dat")`` loads, hashed once per process. Every walk-chemistry PHREEQC consumer
    takes its database from here (never constructs a path of its own) and records the sha256 on its answers."""
    import hashlib
    kwargs, path = phreeqpython_database(HOUSE_DATABASE)
    real = os.path.realpath(path)
    cached = _HOUSE_DATABASE.get(real)
    if cached is None:
        with open(real, 'rb') as handle:
            digest = hashlib.sha256(handle.read()).hexdigest()
        cached = _HOUSE_DATABASE[real] = {'database': HOUSE_DATABASE, 'path': real, 'sha256': digest, 'kwargs': dict(kwargs), 'ruling': HOUSE_DATABASE_RULING_ID}
    return {**cached, 'kwargs': dict(cached['kwargs'])}
_HOUSE_DATABASE: Dict[str, Dict[str, object]] = {}
__all__ = ['AMORPHOUS_IRON_SULFIDE_PHASE', 'HOUSE_IRON_SULFIDE_PHASE', 'IRON_SULFIDE_RULING_ID', 'database_file_record', 'HOUSE_DATABASE', 'HOUSE_DATABASE_RULING_ID', 'house_database', 'MissingPhaseError', 'PHREEQPYTHON_UNNAMED_DEFAULT', 'require_phases', 'RULING_ID', 'SAME_PHASE', 'SamePhase', 'TWO_PHASES', 'database_phase_names', 'phase_names', 'phreeqpython_database', 'resolve_phase', 'same_phase_as']
