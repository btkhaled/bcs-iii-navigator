"""Perméation passive: Fick + Goldman complet + UWL + Renkin."""

from __future__ import annotations

import math

from ...chemistry.descriptors import bohley_radius_nm
from ...core.units import FARADAY, GAS_R, BODY_TEMP_K


def goldman_factor(z: float, vm_mv: float = -60.0) -> float:
    """Facteur Goldman-Hodgkin-Katz pour le flux passif.

    M = xi / (exp(xi) - 1), xi = z*F*Vm/(R*T), Vm = V_inside - V_outside.
    z=0 -> 1 (limite).
    """
    if abs(z) < 1e-9:
        return 1.0
    vm_v = vm_mv / 1000.0
    xi = z * FARADAY * vm_v / (GAS_R * BODY_TEMP_K)
    # garde-fou numérique
    if xi > 50:
        return xi * math.exp(-xi)  # ~0
    if xi < -50:
        return -xi
    denom = math.exp(xi) - 1.0
    if abs(denom) < 1e-12:
        return 1.0
    return xi / denom


def renkin_factor(radius_nm: float, pore_nm: float) -> float:
    a = radius_nm / max(pore_nm, 1e-9)
    if a >= 1.0:
        return 0.0
    return (1 - a) ** 2 * (1 - 2.104 * a + 2.09 * a ** 3 - 0.95 * a ** 5)


def two_pore_factor(radius_nm: float, pore_small_nm: float, pore_large_nm: float, frac_large: float) -> float:
    return (1 - frac_large) * renkin_factor(radius_nm, pore_small_nm) + frac_large * renkin_factor(radius_nm, pore_large_nm)


def uwl_permeability_cm_s(mw: float, h_um: float = 30.0) -> float:
    # Stokes-Einstein simplifié
    r_nm = bohley_radius_nm(mw)
    kB = 1.380649e-23
    T = 310.15
    eta = 0.69e-3
    D_m2_s = kB * T / (6 * math.pi * eta * max(r_nm * 1e-9, 1e-10))
    D_cm2_s = D_m2_s * 1e4
    h_cm = max(h_um * 1e-4, 1e-6)
    return D_cm2_s / h_cm


def series_permeability(p_mem: float, p_uwl: float) -> float:
    if p_mem <= 0:
        return 0.0
    if p_uwl <= 0:
        return p_mem
    return 1.0 / (1.0 / p_mem + 1.0 / p_uwl)
