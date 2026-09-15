"""Dissolution Noyes-Whitney + Johnson (plugin)."""

from __future__ import annotations

import math

from ..core.registry import DISSOLUTION_REGISTRY


def _stokes_diffusion_cm2_s(radius_nm: float) -> float:
    # Stokes-Einstein à 37C, eta eau 0.69e-3 Pa.s
    kB = 1.380649e-23
    T = 310.15
    eta = 0.69e-3
    r_m = max(radius_nm, 0.1) * 1e-9
    D_m2_s = kB * T / (6.0 * math.pi * eta * r_m)
    return D_m2_s * 1e4


@DISSOLUTION_REGISTRY.register("noyes_whitney")
def noyes_whitney_dMdt_mg_h(
    m_solid_mg: float,
    c_mg_l: float,
    cs_mg_l: float,
    radius_nm: float,
    h_um: float = 30.0,
    rho_mg_cm3: float = 1300.0,
    mw: float = 300.0,
) -> float:
    """dM/dt = D*A_s/h * (Cs - C). A_s = 3*M/(rho*r_particule).

    D = Stokes-Einstein sur le rayon MOLÉCULAIRE (diffusion des molécules
    dissoutes), A_s sur le rayon PARTICULAIRE (surface du solide).
    """
    from ..chemistry.descriptors import bohley_radius_nm

    if m_solid_mg <= 0:
        return 0.0
    if cs_mg_l <= c_mg_l:
        return 0.0
    D = _stokes_diffusion_cm2_s(bohley_radius_nm(mw))
    r_cm = radius_nm * 1e-7
    h_cm = h_um * 1e-4
    A_s = 3.0 * m_solid_mg / (rho_mg_cm3 * max(r_cm, 1e-7))  # cm2
    # D*A/h en cm3/s ; (Cs-C) mg/L = mg/1000cm3
    dMdt_mg_s = D * A_s / max(h_cm, 1e-6) * (cs_mg_l - c_mg_l) / 1000.0
    return max(0.0, dMdt_mg_s * 3600.0)


@DISSOLUTION_REGISTRY.register("wang_flanagan")
def wang_flanagan_dMdt_mg_h(m_solid_mg: float, c_mg_l: float, cs_mg_l: float, radius_nm: float, h_um: float = 30.0, rho_mg_cm3: float = 1300.0, mw: float = 300.0) -> float:
    # même base, facteur de forme sphérique (identique ici, hook pour extension)
    return noyes_whitney_dMdt_mg_h(m_solid_mg, c_mg_l, cs_mg_l, radius_nm, h_um, rho_mg_cm3, mw)
