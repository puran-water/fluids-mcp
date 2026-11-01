from __future__ import annotations
import ctypes
import re
CONVERGENCE_RULING_ID = 'PUBLIC-RULE-d6cee394f32b'
ADD_FORMULA_NOT_CONVERGED = 'Pure phase with add formula has not converged'
_WARNING = re.compile('(?:WARNING:\\s*)?(?P<phase>[^\\s,][^,\\n]*?),\\s*' + re.escape(ADD_FORMULA_NOT_CONVERGED) + '\\.?(?P<detail>[^\\n]*(?:\\n[ \\t]+[^\\n]*)?)')
_RESIDUAL = re.compile('Residual:\\s*(?P<value>[-+0-9.eE]+)')

class AddFormulaNotConvergedError(ValueError):
    """PHREEQC did not converge a pure phase with an add formula: the answer is at the wrong pe, pH or SI."""

    def __init__(self, message: str, phases: tuple[str, ...]):
        super().__init__(message)
        self.phases = phases

def unconverged_add_formula_phases(text: str) -> list[dict]:
    """Every pure phase with an add formula PHREEQC's warnings say did not converge, in order, once each: ``phase``,
    ``residual`` (float, where PHREEQC printed it, else None) and the ``warning`` line. Any other warning is not this."""
    found: list[dict] = []
    seen: set[str] = set()
    for match in _WARNING.finditer(text or ''):
        phase = match.group('phase').strip()
        if phase in seen:
            continue
        seen.add(phase)
        residual = _RESIDUAL.search(match.group('detail') or '')
        found.append({'phase': phase, 'residual': float(residual.group('value')) if residual else None, 'warning': f'{phase}, {ADD_FORMULA_NOT_CONVERGED}.'})
    return found

def refuse_unconverged_add_formula(text: str, *, step: str) -> None:
    """Refuses (`AddFormulaNotConvergedError`, a ValueError) where ``text`` (the run's PHREEQC warnings) says a pure
    phase with an add formula did not converge: the run answered at the wrong value of what that phase holds."""
    found = unconverged_add_formula_phases(text)
    if not found:
        return
    phases = tuple((item['phase'] for item in found))
    residuals = ', '.join((f'{item['phase']} (residual {item['residual']:.4g})' if item['residual'] is not None else item['phase'] for item in found))
    raise AddFormulaNotConvergedError(f'{step}: PHREEQC did not converge the pure phase(s) with an add formula {residuals} ("{ADD_FORMULA_NOT_CONVERGED}"), so the run answered at the wrong value of what each holds (a pe, a pH), with a warning only: a Fix_pe in a water nothing poises answers at the wrong pe. Refused, never answered ({CONVERGENCE_RULING_ID}, C15-3): give the water the couple that sets the pe, or a reagent that can reach the pH', phases)

def phreeqpython_warning_text(ip) -> str:
    """The warnings of the last run on a phreeqpython VIPhreeqc instance (``PhreeqPython().ip``): IPhreeqc's own
    ``GetWarningString``, which phreeqpython does not wrap. Each ``run_string`` replaces it. Refuses (RuntimeError)
    where the loaded library has no such symbol: the check cannot run, and a check that cannot run is not a pass."""
    function = getattr(getattr(ip, 'dll', None), 'GetWarningString', None)
    if function is None:
        raise RuntimeError(f"the loaded phreeqpython library exposes no GetWarningString: PHREEQC's warnings cannot be read, so the add-formula convergence check cannot run ({CONVERGENCE_RULING_ID}, C15-3)")
    function.restype = ctypes.c_char_p
    function.argtypes = [ctypes.c_int]
    raw = function(ip.id_)
    return raw.decode('utf-8', 'replace') if raw else ''
