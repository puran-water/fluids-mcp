"A PHREEQC BATCH step (USE solution + EQUILIBRIUM_PHASES / REACTION; phreeqpython's equalize, add, change_ph,\nchange_temperature) re-equilibrates every redox couple of the solution at one pe, where an initial SOLUTION keeps each\nvalence state it was given. No biology acts in an equilibrium step, so a batch step carries ammonium, nitrate and\nsulfide on redox-uncoupled elements derived from the database's own definitions (`decoupled_redox_definitions`, never\ntyped), renders the phases that hold them on those elements (`decoupled_phase_definitions`), and refuses where any redox\nstate moved beyond what a phase took (`check_batch_redox`).\n\nThe definitions are read from any database that writes its blocks the way phreeqc.dat and minteq.v4.dat do: every\nSOLUTION_SPECIES block (minteq.v4.dat has two, and defines NH4+ and HS- in the second), valence states with or without\nthe sign (phreeqc.dat's N(+5), minteq.v4.dat's N(5)), CRLF or LF."
from __future__ import annotations
import functools
import math
import re
from puran_mcp_support.engineering.house_standards import soluble_inert_and_batch_redox as _redox_ruling
REDOX_DECOUPLING_RULING_ID = _redox_ruling.RULING_ID
DECOUPLED_AMMONIUM, DECOUPLED_AMMONIUM_ION = ('Amm', 'AmmH+')
DECOUPLED_SULFIDE = 'Sulf'
DECOUPLED_NITRATE = 'Nit'

def _ammonium_rename(formula: str) -> str:
    'Public engineering calculation interface.'

    def amm(match):
        h = int(match.group(1) or 1)
        if h < 3:
            return match.group(0)
        rest = h - 3
        return DECOUPLED_AMMONIUM + ('' if rest == 0 else 'H' if rest == 1 else f'H{rest}')
    return re.sub('NH(\\d*)(?![0-9])', amm, formula)
DECOUPLED_STATES = {'N(-3)': (DECOUPLED_AMMONIUM, 'N', _ammonium_rename, lambda rows: float(rows['N'][4]) + 3.0 * float(rows['H'][4])), 'N(+5)': (DECOUPLED_NITRATE, 'N', lambda f: re.sub('N(?![a-z])', DECOUPLED_NITRATE, f), lambda rows: float(rows['N'][4])), 'S(-2)': (DECOUPLED_SULFIDE, 'S', lambda f: re.sub('S(?![a-z])', DECOUPLED_SULFIDE, f), lambda rows: float(rows['S'][4]))}
DECOUPLED_ANALYSIS_KEYS = {'N(-3)': DECOUPLED_AMMONIUM, 'N(5)': DECOUPLED_NITRATE, 'N(+5)': DECOUPLED_NITRATE, 'S(-2)': DECOUPLED_SULFIDE}
DECOUPLED_SULFATE = 'Sulfat'
DECOUPLED_NITRITE = 'Nitr'
DECOUPLED_FERROUS = 'Ferr'
OPTIONAL_DECOUPLED_STATES = {'S(+6)': (DECOUPLED_SULFATE, None, lambda f: re.sub('S(?![a-z])', DECOUPLED_SULFATE, f), lambda rows: float(rows['S'][4])), 'N(+3)': (DECOUPLED_NITRITE, 'N', lambda f: re.sub('N(?![a-z])', DECOUPLED_NITRITE, f), lambda rows: float(rows['N'][4])), 'Fe(+2)': (DECOUPLED_FERROUS, 'Fe', lambda f: re.sub('Fe(?![a-z])', DECOUPLED_FERROUS, f), lambda rows: float(rows['Fe'][4]))}
DECOUPLABLE_STATES = {**DECOUPLED_STATES, **OPTIONAL_DECOUPLED_STATES}
BATCH_REDOX_ELEMENTS = {'N': (DECOUPLED_AMMONIUM, DECOUPLED_NITRATE, 'N(-3)', 'N(0)', 'N(3)', 'N(5)'), 'S': (DECOUPLED_SULFIDE, 'S(6)', 'S(-2)'), 'C': ('C(4)', 'C(-4)')}
BATCH_REDOX_STATES = tuple((state for states in BATCH_REDOX_ELEMENTS.values() for state in states))
_UNCHECKED_STATES = frozenset({'C(4)'})
BATCH_REDOX_RELATIVE_TOLERANCE = 0.0001
BATCH_REDOX_ABSOLUTE_TOLERANCE_MOL = 1e-12
_RATE_OPTIONS = {'log_k', 'logk', 'delta_h', 'deltah', 'analytic', 'analytical', 'analytical_expression', 'a_e'}
_TERM = re.compile('^([0-9.]+)?\\s*(\\S+)$')
_BLOCK_KEYWORDS = frozenset({'SOLUTION_MASTER_SPECIES', 'SOLUTION_SPECIES', 'PHASES', 'SURFACE_MASTER_SPECIES', 'SURFACE_SPECIES', 'EXCHANGE_MASTER_SPECIES', 'EXCHANGE_SPECIES', 'RATES', 'END', 'PITZER', 'SIT', 'LLNL_AQUEOUS_MODEL_PARAMETERS', 'NAMED_EXPRESSIONS', 'CALCULATE_VALUES', 'ISOTOPES', 'ISOTOPE_RATIOS', 'ISOTOPE_ALPHAS', 'GAS_BINARY_PARAMETERS', 'SOLUTION', 'SOLUTION_SPREAD', 'SELECTED_OUTPUT', 'USER_PUNCH', 'USER_PRINT', 'KNOBS', 'TITLE', 'DATABASE', 'EQUILIBRIUM_PHASES', 'REACTION', 'SURFACE', 'EXCHANGE', 'GAS_PHASE', 'KINETICS', 'SAVE', 'USE', 'PRINT'})

def _selection(states) -> dict[str, tuple]:
    """The held states, in `DECOUPLABLE_STATES` order whatever the caller's (so the definitions are deterministic);
    ``N(5)`` and ``N(+5)`` name one state. None is `DECOUPLED_STATES`. Refuses a state with no decoupling."""
    if states is None:
        return dict(DECOUPLED_STATES)
    asked = {_valence_key(state) for state in states}
    unknown = sorted(asked - {_valence_key(state) for state in DECOUPLABLE_STATES})
    if unknown:
        raise ValueError(f'no decoupling is defined for {unknown}; decouplable: {sorted(DECOUPLABLE_STATES)}')
    return {state: entry for state, entry in DECOUPLABLE_STATES.items() if _valence_key(state) in asked}

def _species_terms(side: str) -> list[tuple[str, str]]:
    """The (coefficient, species) terms of one side of a reaction. Terms are joined by a "+" with whitespace both sides,
    spaces or tabs (phreeqc.dat writes ``Hfo_sOH	+ H+`` in its SURFACE_SPECIES)."""
    terms = []
    for raw in re.split('\\s+\\+\\s+', side.strip()):
        match = _TERM.match(raw.strip())
        if match is None:
            raise ValueError(f'PHREEQC reaction term {raw!r} not read')
        terms.append((match.group(1) or '', match.group(2)))
    return terms

def _valence_key(name: str) -> str:
    """``N(+5)`` and ``N(5)`` name one state: phreeqc.dat writes the sign, minteq.v4.dat does not."""
    return re.sub('\\(\\+', '(', name)

def _sections(text: str, keyword: str) -> list[list[str]]:
    """Every ``keyword`` block of a PHREEQC text, each as its lines with comments and trailing space stripped (blank
    lines dropped). A database may repeat a block: minteq.v4.dat has two SOLUTION_SPECIES and two PHASES blocks."""
    sections: list[list[str]] = []
    current: list[str] | None = None
    for raw in text.splitlines():
        line = raw.split('#', 1)[0].rstrip()
        if not line.strip():
            continue
        head = line.split()[0] if not line[0].isspace() else None
        if head is not None and head in _BLOCK_KEYWORDS:
            current = [] if head == keyword else None
            if current is not None:
                sections.append(current)
            continue
        if current is not None:
            current.append(line)
    return sections

def _species_entries(text: str, keyword: str) -> list[tuple[str, list[str]]]:
    """(reaction, option lines) of every entry of every ``keyword`` block: an entry starts at its equation, the one line
    carrying ``=`` (phreeqc.dat writes it at column 0, wateq4f.dat indents it like its options). Options written on the
    equation's own line after a ``;`` (wateq4f.dat's ``... + H2O; log_K 4.28``) are its first options."""
    entries: list[tuple[str, list[str]]] = []
    for section in _sections(text, keyword):
        for line in section:
            if '=' in line:
                equation, *inline = line.split(';')
                entries.append((equation.strip(), [option.strip() for option in inline if option.strip()]))
            elif entries:
                entries[-1][1].append(line.strip())
    return entries

def _master_rows(text: str) -> dict[str, list[str]]:
    rows: dict[str, list[str]] = {}
    for section in _sections(text, 'SOLUTION_MASTER_SPECIES'):
        for line in section:
            fields = line.split()
            if len(fields) >= 4:
                rows.setdefault(_valence_key(fields[0]), fields)
    return rows

def decoupled_members(database_text: str, states=None) -> dict[str, set[str]]:
    """The species of each held state (`DECOUPLED_STATES`, or the ``states`` selected): its master species and every
    species the database builds from one, anywhere in a non-redox reaction."""
    return {state: set(built['members']) for state, built in _decoupling(database_text, states).items()}

def _decoupling(database_text: str, states=None) -> dict[str, dict]:
    """The decoupling of ``states`` read from the database (cached: water-chemistry-mcp's batch tools read it on every
    run; the entries are never mutated by a caller)."""
    return _decoupling_cached(database_text, None if states is None else tuple(states))

@functools.lru_cache(maxsize=32)
def _decoupling_cached(database_text: str, states) -> dict[str, dict]:
    if 'SOLUTION_MASTER_SPECIES' not in database_text or 'SOLUTION_SPECIES' not in database_text:
        raise ValueError('the PHREEQC database has no SOLUTION_MASTER_SPECIES / SOLUTION_SPECIES blocks')
    selected = _selection(states)
    rows = _master_rows(database_text)
    entries = _species_entries(database_text, 'SOLUTION_SPECIES')
    built: dict[str, dict] = {}
    claimed: dict[int, dict[str, set[str]]] = {}
    for state, (element, gfw_formula, rename, element_gfw) in selected.items():
        row = rows.get(_valence_key(state))
        if row is None or any((e not in rows or len(rows[e]) < 5 for e in ('N', 'H', 'S'))):
            raise ValueError(f'the PHREEQC database does not define {state} beside the N, H and S element weights')
        members, identity, others = ({row[1]}, None, [])
        for index, (reaction, options) in enumerate(entries):
            left, right = reaction.split('=', 1)
            lhs, rhs = (_species_terms(left), _species_terms(right))
            defined = rhs[0][1]
            named = {name for _c, name in lhs + rhs[1:]} & members
            if defined == row[1] and ('e-' in {name for _c, name in lhs + rhs} or reaction.replace(' ', '') == f'{row[1]}={row[1]}'):
                identity = [o for o in options if o.split()[0].lstrip('-').lower() not in _RATE_OPTIONS]
                continue
            if 'e-' in {name for _c, name in lhs + rhs} or _implicit_redox(options):
                continue
            if not named and defined not in members:
                continue
            members.add(defined)
            claimed.setdefault(index, {})[state] = {name for _c, name in lhs + rhs} & members
            others.append(index)
        if identity is None:
            raise ValueError(f"the PHREEQC database gives {state}'s master species {row[1]} no definition")
        built[state] = {'row': row, 'rows': rows, 'members': members, 'identity': identity, 'others': others, 'element': element, 'gfw_formula': row[3] if gfw_formula is None else gfw_formula, 'rename': rename, 'element_gfw': element_gfw}
    written: set[int] = set()
    for state, b in built.items():
        texts = []
        for index in b['others']:
            if index in written:
                continue
            written.add(index)
            reaction, options = entries[index]
            left, right = reaction.split('=', 1)

            def name(species, takers=claimed[index]):
                renamed = species
                for taker, members_named in takers.items():
                    if species in members_named:
                        renamed = built[taker]['rename'](renamed)
                return renamed
            side = lambda terms: ' + '.join((f'{c}{name(n)}' for c, n in terms))
            texts.append('\n'.join([f'{side(_species_terms(left))} = {side(_species_terms(right))}', *('    ' + o for o in options)]))
        b['others'] = texts
    return built
_MASS_BALANCE_OPTIONS = frozenset({'-mass_balance', '-mb', 'mass_balance', 'mb'})

def _implicit_redox(options: list[str]) -> bool:
    return any((option.split()[0].lower() in _MASS_BALANCE_OPTIONS for option in options if option.split()))

def decoupled_species_names(database_text: str, states=None) -> dict[str, str]:
    """Every species of a held state (`decoupled_members`) -> its name on the decoupled elements, as
    `decoupled_redox_definitions` writes it (an ion pair of two held states renamed by both: NH4SO4- -> AmmHSulfatO4-),
    so an answer read on the held elements can be mapped back to the database's own names."""
    built = _decoupling(database_text, states)
    names: dict[str, str] = {}
    for species in set().union(*(b['members'] for b in built.values())):
        name = species
        for b in built.values():
            if species in b['members']:
                name = b['rename'](name)
        names[species] = name
    return names

def defined_valence_states(database_text: str) -> dict[str, str]:
    """Each valence state the database's SOLUTION_MASTER_SPECIES defines, keyed without a sign (N(5)) -> as the
    database spells it (phreeqc.dat's N(+5), minteq.v4.dat's N(5))."""
    spelled: dict[str, str] = {}
    for section in _sections(database_text, 'SOLUTION_MASTER_SPECIES'):
        for line in section:
            fields = line.split()
            if len(fields) >= 4 and '(' in fields[0]:
                spelled.setdefault(_valence_key(fields[0]), fields[0])
    return spelled

def decoupled_master_species(database_text: str, states=None) -> dict[str, str]:
    """Each held state's master species on its decoupled element (NO3- -> NitO3-, SO4-2 -> SulfatO4-2)."""
    return {state: built['rename'](built['row'][1]) for state, built in _decoupling(database_text, states).items()}

def decoupled_analysis_keys(states=None) -> dict[str, str]:
    """Analysis key -> decoupled element for the held states: each state under both spellings (N(5) and N(+5)). None is
    `DECOUPLED_ANALYSIS_KEYS`."""
    if states is None:
        return dict(DECOUPLED_ANALYSIS_KEYS)
    return {key: entry[0] for state, entry in _selection(states).items() for key in (state, _valence_key(state))}

def decoupled_redox_definitions(database_text: str, states=None) -> str:
    """The walk's ammonium, nitrate and sulfide as redox-uncoupled elements (`DECOUPLED_STATES`): a PHREEQC input
    (SOLUTION_MASTER_SPECIES + SOLUTION_SPECIES) re-writing every species the database builds from each state's master
    species, every option and constant copied as the database writes it. The master species' activity, diffusion and
    volume options (on its identity, or on its redox definition from another state) become the identity's. A species
    built from a member, anywhere in a non-redox reaction, is a member. A redox reaction is never copied: it is the
    coupling being cut. Refuses, never guesses, where a state is missing or gives its master species no definition.
    ``states``: the held states (None: `DECOUPLED_STATES`, the walk's)."""
    masters, species = (['SOLUTION_MASTER_SPECIES'], ['SOLUTION_SPECIES'])
    for built in _decoupling(database_text, states).values():
        row, rename = (built['row'], built['rename'])
        decoupled_master = rename(row[1])
        masters.append(f'{built['element']}    {decoupled_master}    {row[2]}    {built['gfw_formula']}    {built['element_gfw'](built['rows']):.6f}')
        species += [f'{decoupled_master} = {decoupled_master}', *('    ' + o for o in built['identity']), *built['others']]
    return '\n'.join(masters + species)

def _phase_entries(text: str) -> dict[str, tuple[str, list[str]]]:
    """name -> (reaction, option lines) of every phase of every PHASES block of ``text``; a later definition wins, as
    it does in PHREEQC."""
    phases: dict[str, tuple[str, list[str]]] = {}
    for section in _sections(text, 'PHASES'):
        name: str | None = None
        for line in section:
            if not line[0].isspace():
                name = line.split()[0]
                phases[name] = ('', [])
            elif name is not None:
                reaction, options = phases[name]
                if not reaction and '=' in line:
                    phases[name] = (line.strip(), options)
                else:
                    options.append(line.strip())
    return phases

def decoupled_phase_definitions(database_text: str, texts: list[str], phase_names, states=None) -> tuple[str, dict[str, str]]:
    """A PHASES block re-defining, under its own name, every phase of ``phase_names`` whose dissolution reaction names a
    decoupled state's species (`decoupled_members` of ``database_text``): the phase's formula and each member species on
    the decoupled element (struvite's MgNH4PO4:6H2O = ... + NH4+ becomes MgAmmHPO4:6H2O = ... + AmmH+), every constant
    copied. Definitions are read from ``texts`` in order (the database, then any input that defines phases), the last
    one winning as in PHREEQC. Returns the block ("" when no phase needs one) and each re-defined phase's formula. A
    phase defined nowhere is left to PHREEQC, which refuses it by name. Refuses a redox reaction that names a member (it
    cannot be re-written without the coupling being cut). ``states``: the held states (None: `DECOUPLED_STATES`)."""
    built = _decoupling(database_text, states)
    known: dict[str, tuple[str, list[str]]] = {}
    for text in texts:
        known.update(_phase_entries(text))
    block, formulas = (['PHASES'], {})
    for name in dict.fromkeys(phase_names):
        if name not in known or not known[name][0]:
            continue
        reaction, options = known[name]
        left, right = reaction.split('=', 1)
        lhs, rhs = (_species_terms(left), _species_terms(right))
        names = {species for _c, species in lhs[1:] + rhs}
        hits = [state for state, b in built.items() if names & b['members']]
        if not hits:
            continue
        if 'e-' in names:
            raise ValueError(f'phase {name!r} ({reaction}) reacts a decoupled redox state in a redox reaction: it cannot be carried onto the decoupled elements ({REDOX_DECOUPLING_RULING_ID})')
        formula = lhs[0][1]
        for state in hits:
            formula = built[state]['rename'](formula)

        def side(terms):
            out = []
            for coefficient, species in terms:
                for state in hits:
                    if species in built[state]['members']:
                        species = built[state]['rename'](species)
                out.append(f'{coefficient}{species}')
            return ' + '.join(out)
        lhs_text = ' + '.join([f'{lhs[0][0]}{formula}', side(lhs[1:])]) if len(lhs) > 1 else f'{lhs[0][0]}{formula}'
        block += [name, f'    {lhs_text} = {side(rhs)}', *('    ' + o for o in options)]
        formulas[name] = formula
    return ('\n'.join(block) if formulas else '', formulas)

def phase_formulas(texts: list[str]) -> dict[str, str]:
    """name -> formula (the first term of its dissolution reaction) of every phase of every PHASES block of ``texts``, a
    later definition winning as in PHREEQC; a phase with no readable reaction is left out."""
    known: dict[str, tuple[str, list[str]]] = {}
    for text in texts:
        known.update(_phase_entries(text))
    formulas = {}
    for name, (reaction, _options) in known.items():
        if '=' not in reaction:
            continue
        try:
            formulas[name] = _species_terms(reaction.split('=', 1)[0])[0][1]
        except ValueError:
            continue
    return formulas

def phases_naming_held_states(database_text: str, texts: list[str], states=None) -> tuple[list[str], list[str]]:
    """The phases of ``texts`` (the database, then any input that defines phases) whose dissolution reaction names a
    held state's species (`decoupled_members` of ``database_text``), split into those `decoupled_phase_definitions`
    carries onto the decoupled elements and those it cannot: a redox reaction (FeS2 + 2 H+ + 2 e- = ..., S + 2 H+ +
    2 e- = H2S), or one the term reader cannot read while it names a member. A batch that lets one of the second kind
    react cannot hold its states."""
    members = set().union(*decoupled_members(database_text, states).values())
    known: dict[str, tuple[str, list[str]]] = {}
    for text in texts:
        known.update(_phase_entries(text))
    holdable, coupled = ([], [])
    for name, (reaction, _options) in known.items():
        if '=' not in reaction:
            continue
        try:
            left, right = reaction.split('=', 1)
            lhs, rhs = (_species_terms(left), _species_terms(right))
        except ValueError:
            if any((re.search(f'(?<![A-Za-z0-9]){re.escape(m)}(?![A-Za-z0-9+-])', reaction) for m in members)):
                coupled.append(name)
            continue
        names = {species for _c, species in lhs[1:] + rhs}
        if not names & members:
            continue
        (coupled if 'e-' in names else holdable).append(name)
    return (holdable, coupled)

def reactions_naming_decoupled_members(database_text: str, texts: list[str] | None=None, keywords=('SURFACE_SPECIES', 'EXCHANGE_SPECIES'), states=None) -> list[str]:
    """The surface and exchange reactions of ``texts`` (default: the database itself; pass any input that defines its
    own) that name a decoupled state's species (`decoupled_members` of ``database_text``). Left as written they bind a
    coupled species, so a caller whose batch has an exchanger refuses on any (the adsorption would vanish silently); a
    surface is carried onto the decoupled elements by `decoupled_surface_definitions`."""
    members = set().union(*decoupled_members(database_text, states).values())
    found = []
    for keyword, text in ((k, x) for k in keywords for x in (texts if texts is not None else [database_text])):
        for reaction, _options in _species_entries(text, keyword):
            if '=' in reaction and {s for _c, s in _species_terms(reaction.replace('=', ' + '))} & members:
                found.append(reaction)
    return found

def decoupled_surface_definitions(database_text: str, texts: list[str], states=None) -> str:
    """A SURFACE_SPECIES block re-defining every surface reaction of ``texts`` (the database, then any input that defines
    its own; a later definition of one surface species wins, as in PHREEQC) that names a held state's species: the
    member species AND the surface species the reaction defines renamed onto the decoupled element (PHREEQC reads a
    species' composition from its name, so Hfo_wSO4- formed from SulfatO4-2 must be Hfo_wSulfatO4-), every option
    copied. The coupled definition stays and binds nothing (its species is not in the batch). "" when no reaction names
    a member. Refuses a redox surface reaction that names one (it cannot be re-written without the coupling being
    cut)."""
    built = _decoupling(database_text, states)
    found: dict[str, str] = {}
    for text in texts:
        for reaction, options in _species_entries(text, 'SURFACE_SPECIES'):
            left, right = reaction.split('=', 1)
            lhs, rhs = (_species_terms(left), _species_terms(right))
            names = {species for _c, species in lhs + rhs}
            hits = [state for state, b in built.items() if names & b['members']]
            if not hits:
                continue
            if 'e-' in names:
                raise ValueError(f'surface reaction {reaction!r} reacts a held redox state in a redox reaction: it cannot be carried onto the decoupled elements ({REDOX_DECOUPLING_RULING_ID})')
            defined = rhs[0][1]

            def name(species, hits=hits, defined=defined):
                renamed = species
                for state in hits:
                    if species == defined or species in built[state]['members']:
                        renamed = built[state]['rename'](renamed)
                return renamed
            side = lambda terms: ' + '.join((f'{c}{name(n)}' for c, n in terms))
            found[defined] = '\n'.join([f'{side(lhs)} = {side(rhs)}', *('    ' + o for o in options)])
    return '\n'.join(['SURFACE_SPECIES', *found.values()]) if found else ''
_FORMULA_TOKEN = re.compile('([A-Z][a-z]*)(\\d*\\.?\\d*)|(\\()|(\\))(\\d*\\.?\\d*)')

def formula_element_counts(formula: str) -> dict[str, float]:
    """Atoms of each element in a phase formula (parentheses and ``:`` hydrate parts read): ``MgAmmHPO4:6H2O`` ->
    Mg 1, Amm 1, H 13, P 1, O 10."""
    counts: dict[str, float] = {}
    for index, part in enumerate(formula.split(':')):
        multiplier = 1.0
        if index:
            match = re.match('^(\\d*\\.?\\d+)', part)
            if match:
                multiplier, part = (float(match.group(1)), part[match.end():])
        stack: list[dict[str, float]] = [{}]
        position = 0
        while position < len(part):
            match = _FORMULA_TOKEN.match(part, position)
            if match is None:
                raise ValueError(f'phase formula {formula!r} not read at {part[position:]!r}')
            position = match.end()
            if match.group(1):
                stack[-1][match.group(1)] = stack[-1].get(match.group(1), 0.0) + float(match.group(2) or 1)
            elif match.group(3):
                stack.append({})
            else:
                group, factor = (stack.pop(), float(match.group(5) or 1))
                for element, n in group.items():
                    stack[-1][element] = stack[-1].get(element, 0.0) + n * factor
        if len(stack) != 1:
            raise ValueError(f'phase formula {formula!r} has unbalanced parentheses')
        for element, n in stack[0].items():
            counts[element] = counts.get(element, 0.0) + n * multiplier
    return counts
_REAGENT_GROUPS = (('SO4', 'S(+6)', '{}O4'), ('NO3', 'N(+5)', '{}O3'), ('NO2', 'N(+3)', '{}O2'), ('NH4', 'N(-3)', '{}H'), ('NH3', 'N(-3)', '{}'))
_SULFIDE_S = re.compile('S(?![a-z0-9])')

def held_reagent_formula(formula: str, states=None) -> str:
    """The formula a REACTION adds a reagent with when ``states`` are held (`DECOUPLABLE_STATES`; None:
    `DECOUPLED_STATES`): each N and S group put on its held state's element, as its own groups name it (SO4 -> SulfatO4,
    NO3 -> NitO3, NO2 -> NitrO2, NH4 -> AmmH, NH3 -> Amm, a sulfide S -> Sulf: FeSO4 -> FeSulfatO4, (NH4)2SO4 ->
    (AmmH)2SulfatO4, Na2S:9H2O -> Na2Sulf:9H2O). A sulfide S stands alone in a formula with no oxygen outside its
    hydrate water and with another element beside it (Na2S, NaHS, FeS; not S, S8 or FeS2). An element no held state
    belongs to is left as written. Refuses (ValueError naming the reagent and `REDOX_DECOUPLING_RULING_ID`) any N or S
    the groups do not place on a held state: the valence would be guessed, or the reagent's own redox would move a held
    state."""
    selected = _selection(states)
    held_elements = {state.split('(', 1)[0] for state in selected}
    if not held_elements & {'N', 'S'}:
        return formula
    out = formula
    for group, state, template in _REAGENT_GROUPS:
        if state[0] not in held_elements:
            continue
        pattern = re.compile(re.escape(group) + '(?![0-9])')
        if not pattern.search(out):
            continue
        entry = next((e for k, e in selected.items() if _valence_key(k) == _valence_key(state)), None)
        if entry is None:
            raise ValueError(f'reagent {formula} carries {group} ({state}), a state this batch has not held: it cannot be added on a held element ({REDOX_DECOUPLING_RULING_ID})')
        out = pattern.sub(template.format(entry[0]), out)
    body = out.split(':', 1)[0]
    if 'S' in held_elements and _SULFIDE_S.search(body):
        counts = formula_element_counts(body)
        others = {e for e, n in counts.items() if e != 'S' and n}
        sulfide = next((e for k, e in selected.items() if _valence_key(k) == 'S(-2)'), None)
        if sulfide is None or 'O' in others or (not others) or (len(_SULFIDE_S.findall(body)) != 1):
            raise ValueError(f'reagent {formula} carries sulfur in a valence no held state carries (only SO4 and a lone sulfide S are placed): it cannot be added on a held element ({REDOX_DECOUPLING_RULING_ID})')
        out = _SULFIDE_S.sub(sulfide[0], body, count=1) + out[len(body):]
    left = formula_element_counts(out)
    unplaced = sorted((e for e in ('N', 'S') if e in held_elements and left.get(e)))
    if unplaced:
        raise ValueError(f'reagent {formula} carries {' and '.join(unplaced)} in a valence no held state carries (the groups placed: {', '.join((g for g, *_ in _REAGENT_GROUPS))} and a lone sulfide S): it cannot be added on a held element ({REDOX_DECOUPLING_RULING_ID})')
    return out

def batch_redox_elements(states=None) -> dict[str, tuple[str, ...]]:
    """`BATCH_REDOX_ELEMENTS` for the held ``states``: by element, the decoupled element of each held state, then every
    coupled state (a re-equilibration would move a held state into one). None is `BATCH_REDOX_ELEMENTS` itself."""
    if states is None:
        return BATCH_REDOX_ELEMENTS
    selected = _selection(states)
    return {element: tuple((entry[0] for state, entry in selected.items() if state.split('(', 1)[0] == element)) + tuple((state for state in names if '(' in state)) for element, names in BATCH_REDOX_ELEMENTS.items()}

def batch_redox_states(database_text: str, states=None) -> tuple[str, ...]:
    """The states of `batch_redox_elements` a database defines (a decoupled element is defined by
    `decoupled_redox_definitions`); the rest cannot occur on it (minteq.v4.dat defines no N(0) and no C(-4)), and a
    SELECTED_OUTPUT naming one fails."""
    rows = _master_rows(database_text)
    decoupled = {built[0] for built in _selection(states).values()}
    return tuple((state for names in batch_redox_elements(states).values() for state in names if state in decoupled or _valence_key(state) in rows))

def batch_redox_moles(ip, solution: int) -> dict[str, float]:
    """Moles of every `BATCH_REDOX_STATES` state in one PHREEQC solution. phreeqpython's ``get_total`` is the
    solution's MOLES (TOTMOLE), not mol/kgw: measured 2026-10-15, Na after a CO2 batch reads 0.1304927 from get_total
    and from TOTMOLE, 0.1297535 mol/kgw on 1.0057 kg of water."""
    return {state: float(ip.get_total(solution, state)) for state in BATCH_REDOX_STATES}

def check_batch_redox(before: dict, after: dict, *, step: str, taken: dict | None=None, unchecked=(), elements=None) -> None:
    """Refuses (ValueError naming each state) where a batch step moved a redox state's moles by more than the phases
    ``taken`` from it (`REDOX_DECOUPLING_RULING_ID`): an equilibrium step re-equilibrating redox is a defect, never an
    answer. ``unchecked``: further states read only to scale their element's tolerance, for a caller whose batch moves
    one legitimately and has proved it cannot move by redox (as C(4) leaves as CO2 gas or a carbonate). ``elements``:
    the states checked, by element (`batch_redox_elements` of the held states; None: `BATCH_REDOX_ELEMENTS`)."""
    taken = taken or {}
    skip = _UNCHECKED_STATES | frozenset(unchecked)
    moved = {}
    for states in (BATCH_REDOX_ELEMENTS if elements is None else elements).values():
        scale = math.fsum((abs(before[state]) for state in states))
        for state in states:
            if state in skip:
                continue
            expected = before[state] - taken.get(state, 0.0)
            if abs(after[state] - expected) <= BATCH_REDOX_RELATIVE_TOLERANCE * scale + BATCH_REDOX_ABSOLUTE_TOLERANCE_MOL:
                continue
            moved[state] = f'{before[state]:.6g} -> {after[state]:.6g} mol (phases took {taken.get(state, 0.0):.6g})'
    if moved:
        raise ValueError(f'{step}: the PHREEQC batch step moved redox state(s) {moved} ({REDOX_DECOUPLING_RULING_ID}: an equilibrium step re-equilibrates no redox couple; carry the reduced states decoupled)')
WATER_REDOX_RULING_ID = 'PUBLIC-RULE-d6cee394f32b'
WATER_REDOX_STATES = ('O(0)', 'H(0)')
WATER_REDOX_NOISE_MOL_PER_KGW = 1e-07
WATER_REDOX_STRICT_RULING_ID = 'PUBLIC-RULE-2938127a7d60'

def water_redox_tolerance_mol(before: dict, *, held_mol: float, water_kg: float) -> float:
    """The moles of O(0) or H(0) a held step may move (`check_water_redox`): the shared relative tolerance
    (`BATCH_REDOX_RELATIVE_TOLERANCE`) of the redox inventory the step holds back from reacting, the held states' moles
    (``held_mol``: the electrons a re-equilibration would have put there are the ones that would land in H2 or O2) and the
    feed's own O(0) and H(0), plus PHREEQC's own noise on ``water_kg`` of water (`WATER_REDOX_NOISE_MOL_PER_KGW`)."""
    scale = math.fsum((abs(float(before.get(state, 0.0))) for state in WATER_REDOX_STATES)) + abs(float(held_mol))
    return BATCH_REDOX_RELATIVE_TOLERANCE * scale + max(WATER_REDOX_NOISE_MOL_PER_KGW * abs(float(water_kg)), BATCH_REDOX_ABSOLUTE_TOLERANCE_MOL)

def check_water_redox(before: dict, after: dict, *, step: str, held_mol: float, water_kg: float) -> None:
    """Refuses (ValueError naming each state) where a held batch step moved dissolved O2 (O(0)) or H2 (H(0)), in moles,
    from the feed's by more than `water_redox_tolerance_mol` (`WATER_REDOX_RULING_ID`, C15-2). The caller exempts a step
    where a gas phase, a Fix_pe or an O2(g)/H2(g) phase sets them, and names the exemption; it never calls this there."""
    tolerance = water_redox_tolerance_mol(before, held_mol=held_mol, water_kg=water_kg)
    moved = {state: f'{float(before.get(state, 0.0)):.6g} -> {float(after.get(state, 0.0)):.6g} mol' for state in WATER_REDOX_STATES if abs(float(after.get(state, 0.0)) - float(before.get(state, 0.0))) > tolerance}
    if moved:
        raise ValueError(f"{step}: the held PHREEQC batch step moved dissolved O2/H2 {moved} from the feed's, beyond {tolerance:.3g} mol ({WATER_REDOX_RULING_ID}, C15-2: with the N and S couples held, the electrons a reagent handed the water, or took from it, ended in H2 or O2; no gas phase or Fix_pe sets them here)")
__all__ = ['BATCH_REDOX_ABSOLUTE_TOLERANCE_MOL', 'BATCH_REDOX_ELEMENTS', 'BATCH_REDOX_RELATIVE_TOLERANCE', 'WATER_REDOX_NOISE_MOL_PER_KGW', 'WATER_REDOX_RULING_ID', 'WATER_REDOX_STATES', 'WATER_REDOX_STRICT_RULING_ID', 'check_water_redox', 'water_redox_tolerance_mol', 'BATCH_REDOX_STATES', 'DECOUPLABLE_STATES', 'DECOUPLED_AMMONIUM', 'DECOUPLED_AMMONIUM_ION', 'DECOUPLED_ANALYSIS_KEYS', 'DECOUPLED_FERROUS', 'DECOUPLED_NITRATE', 'DECOUPLED_NITRITE', 'DECOUPLED_STATES', 'DECOUPLED_SULFATE', 'DECOUPLED_SULFIDE', 'OPTIONAL_DECOUPLED_STATES', 'REDOX_DECOUPLING_RULING_ID', 'batch_redox_elements', 'batch_redox_moles', 'batch_redox_states', 'check_batch_redox', 'decoupled_analysis_keys', 'decoupled_master_species', 'decoupled_members', 'decoupled_phase_definitions', 'decoupled_redox_definitions', 'decoupled_species_names', 'decoupled_surface_definitions', 'defined_valence_states', 'formula_element_counts', 'held_reagent_formula', 'phase_formulas', 'phases_naming_held_states', 'reactions_naming_decoupled_members']
