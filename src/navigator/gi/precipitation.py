"""Précipitation: cinétique supersaturation -> solide."""

from __future__ import annotations


def precipitation_rate_mg_h(c_mg_l: float, cs_mg_l: float, volume_l: float, k_precip_h: float = 5.0) -> float:
    if c_mg_l <= cs_mg_l:
        return 0.0
    return k_precip_h * (c_mg_l - cs_mg_l) * volume_l
