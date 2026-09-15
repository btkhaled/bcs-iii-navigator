"""Endocytose: MM sur C_mucus, pondérée taille."""

from __future__ import annotations

from ...core.units import pmol_cm2_min_to_mg_h, mg_l_to_um
from .carriers import mm_rate


def endo_flux_mg_h(c_mucus_mg_l: float, area_cm2: float, mw: float, vmax_scale: float = 1.0) -> float:
    # nulle sous ~500 Da, pleine au-delà ; Vmax basal faible
    if mw < 500:
        w = 0.0
    elif mw < 1000:
        w = (mw - 500) / 500.0
    else:
        w = 1.0
    if w <= 0:
        return 0.0
    vmax = 50.0 * vmax_scale * w  # pmol/cm2/min basal
    km = 500.0  # uM
    v = mm_rate(vmax, km, mg_l_to_um(max(c_mucus_mg_l, 0.0), mw))
    return max(0.0, pmol_cm2_min_to_mg_h(v, area_cm2, mw))
