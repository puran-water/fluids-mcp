"""Mechanical tools delegate to the shared deterministic calculation library."""
from puran_mcp_support.fluid_mechanics import EductorInput, NoiseInput, ReliefValveInput

def relief_valve_size(request: ReliefValveInput) -> dict:
    """API 520 gas/steam/liquid sizing; API 521 case/load supplied. SI units.

    Missing set pressure, back pressure, case/load/factors or fire basis return
    named holds. Set pressure is Pa gauge; back pressure is Pa absolute.
    """
    from puran_mcp_support.fluid_mechanics import relief_valve_size as calculate
    return calculate(request)

def valve_noise(request: NoiseInput) -> dict:
    """IEC 60534 gas/liquid valve noise, or API 521 unweighted flare-tip noise.

    Explicit SI properties, geometry, correction factors and basis are required.
    Select the correlation; API 521 is not an IEC valve-noise calculation.
    """
    from puran_mcp_support.fluid_mechanics import valve_noise as calculate
    return calculate(request)

def eductor_size(request: EductorInput) -> dict:
    """Single-phase liquid jet eductor; exactly two unknowns, SI inputs.

    Explicit material, diffuser geometry and loss coefficients; no gas or
    cavitation-qualified ejector design. Uses fluids.jet_pump.
    """
    from puran_mcp_support.fluid_mechanics import eductor_size as calculate
    return calculate(request)
