"""Numeric evaluation of the shared mADM1 biokinetic kernel (``process_model_schema.madm1_biokinetics``).

The reference arithmetic both engines are held to: :func:`process_rates` evaluates the kernel's rate expressions at a
state, :func:`stoichiometric_matrix` evaluates its stoichiometric table and closes the conserved carriers (S_IC,
S_IN, S_IP) from the component contents. The QSDsan twin's ``rhos_madm1`` and the WaterTAP kernel reaction block are
each tested against these functions at common states, so the two engines agree by construction of the tests, not by
assertion.
"""
from __future__ import annotations
from typing import Mapping, Optional
from puran_mcp_support.schema.madm1_biokinetics import DERIVED_SYMBOLS, MADM1_BIOKINETICS_KERNEL, BiokineticKernelDefinition, derived_parameter_values
from .evaluator import evaluate
__all__ = ['madm1_parameter_values', 'madm1_process_rates', 'madm1_stoichiometric_matrix']

def madm1_parameter_values(overrides: Optional[Mapping[str, float]]=None, *, kernel: BiokineticKernelDefinition=MADM1_BIOKINETICS_KERNEL) -> dict[str, float]:
    """The kernel's parameter values with ``overrides`` applied, plus the derived stoichiometric parameters.

    An override must name a kernel parameter: an unknown name is refused, so a misspelt key never silently leaves
    the default in force.
    """
    values = kernel.parameter_values()
    for name, value in (overrides or {}).items():
        if name not in values:
            raise KeyError(f'{name!r} is not a parameter of {kernel.kernel_id}')
        values[name] = float(value)
    values.update(derived_parameter_values(values))
    return values

def madm1_process_rates(state: Mapping[str, float], derived: Mapping[str, float], *, parameters: Optional[Mapping[str, float]]=None, kernel: BiokineticKernelDefinition=MADM1_BIOKINETICS_KERNEL) -> dict[str, float]:
    """Each process's rate (kg reference component / m3 / d) at ``state`` (kg/m3 on the measured basis) with the
    engine's speciation ``derived`` (:data:`DERIVED_SYMBOLS`)."""
    missing = sorted(set(DERIVED_SYMBOLS) - set(derived))
    if missing:
        raise KeyError(f'derived speciation symbols missing: {missing}')
    values = madm1_parameter_values(parameters, kernel=kernel)
    symbols = {**values, **{k: float(v) for k, v in state.items()}, **{k: float(v) for k, v in derived.items()}}
    return {p.process_id: evaluate(p.rate, symbols) for p in kernel.processes}

def madm1_stoichiometric_matrix(parameters: Optional[Mapping[str, float]]=None, *, kernel: BiokineticKernelDefinition=MADM1_BIOKINETICS_KERNEL) -> dict[str, dict[str, float]]:
    """{process_id: {component_id: coefficient}} over every kernel component; carriers closed by C, N, P balance.

    The closure is the twin's ``conserved_for=('C', 'N', 'P')``: carrier coefficient = - sum over the stated
    coefficients of coefficient x element content (each carrier holds only its own element).
    """
    values = madm1_parameter_values(parameters, kernel=kernel)
    carriers = dict(kernel.conserved_carriers)
    matrix: dict[str, dict[str, float]] = {}
    for process in kernel.processes:
        row = {cid: 0.0 for cid in kernel.component_ids}
        for component, coefficient in process.stoichiometry.items():
            row[component] = evaluate(coefficient, values)
        for element, carrier in carriers.items():
            row[carrier] = -sum((row[cid] * kernel.content(cid).content(element) for cid in kernel.component_ids if cid not in carriers.values()))
        matrix[process.process_id] = row
    return matrix
