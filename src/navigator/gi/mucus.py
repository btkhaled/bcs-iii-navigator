"""Mucus: obstruction Amsden + zêta + dynamique (influx k + clearance)."""

from __future__ import annotations


def amsden_factor(particle_nm: float, d_crit_nm: float = 150.0) -> float:
    if particle_nm <= 0:
        return 1.0
    return 1.0 / (1.0 + (particle_nm / max(d_crit_nm, 1e-9)) ** 2)


def zeta_factor(zeta_mv: float) -> float:
    # piégeage si très chargé
    z = abs(zeta_mv)
    if z <= 10:
        return 1.0
    if z >= 40:
        return 0.2
    return 1.0 - 0.8 * (z - 10.0) / 30.0


def f_mucus(particle_nm: float, thickness_um: float, zeta_mv: float, d_crit_nm: float = 150.0) -> float:
    import math

    f = amsden_factor(particle_nm, d_crit_nm)
    f *= math.exp(-0.01 * (thickness_um / 100.0))
    f *= zeta_factor(zeta_mv)
    return max(0.01, min(1.0, f))
