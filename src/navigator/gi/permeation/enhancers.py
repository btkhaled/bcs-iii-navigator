"""Enhancers: Hill / TEER / Bliss + dynamiques tau_on / tau_rec."""

from __future__ import annotations

from typing import Dict, List


def hill_drop(c_mm: float, ec50_mm: float, dmax: float, n: float) -> float:
    """Fraction de chute TEER (0..dmax)."""
    if c_mm <= 0 or ec50_mm <= 0:
        return 0.0
    e = dmax * (c_mm ** n) / (ec50_mm ** n + c_mm ** n)
    return max(0.0, min(dmax, e))


def bliss_total(drops: List[float]) -> float:
    """Indépendance de Bliss: E_tot = 1 - prod(1 - E_i)."""
    p = 1.0
    for e in drops:
        p *= (1.0 - max(0.0, min(1.0, e)))
    return 1.0 - p


def teer_target(teer0: float, e_tot: float) -> float:
    return teer0 * (1.0 - e_tot)


def teer_step(teer: float, target: float, dt_h: float, tau_on_h: float, tau_rec_h: float) -> float:
    tau = tau_on_h if target < teer else tau_rec_h
    if tau <= 0:
        return target
    alpha = min(1.0, dt_h / tau)
    return teer + alpha * (target - teer)


def teer_multiplier(teer0: float, teer: float, teer_max_drop: float = 10.0) -> float:
    if teer <= 0:
        return teer_max_drop
    return max(1.0, min(teer_max_drop, teer0 / max(teer, 1e-9)))
