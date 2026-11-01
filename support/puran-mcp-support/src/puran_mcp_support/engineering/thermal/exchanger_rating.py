from __future__ import annotations
import math

class TemperatureCross(ValueError):
    """A terminal temperature difference is not positive: no counterflow exchanger delivers these temperatures."""

def lmtd_counterflow(t_hot_in_C: float, t_hot_out_C: float, t_cold_in_C: float, t_cold_out_C: float) -> float:
    """Counterflow log-mean temperature difference [K].

    ``d1 = T_hot,in - T_cold,out`` (hot end), ``d2 = T_hot,out - T_cold,in`` (cold end);
    ``LMTD = (d1 - d2) / ln(d1 / d2)``, ``d1`` when the two ends are equal. Raises :class:`TemperatureCross` when either
    end is not positive."""
    d1 = t_hot_in_C - t_cold_out_C
    d2 = t_hot_out_C - t_cold_in_C
    if d1 <= 0 or d2 <= 0:
        raise TemperatureCross(f'temperatures cross: hot {t_hot_in_C}->{t_hot_out_C} against cold {t_cold_in_C}->{t_cold_out_C} gives terminal differences {d1:.1f} and {d2:.1f} K')
    if abs(d1 - d2) < 1e-09:
        return d1
    return (d1 - d2) / math.log(d1 / d2)
