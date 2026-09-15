"""Premier passage Fg x Fh (well-stirred hépatique)."""

from __future__ import annotations


def fh_liver(fu: float, clint_liver_ul_min_mg: float, qh_l_h: float) -> float:
    if clint_liver_ul_min_mg <= 0:
        return 1.0
    cl_l_h = clint_liver_ul_min_mg * 1e-6 * 60.0 * 40.0 * 1500.0 / 1000.0  # ~ foie 1500g
    # borne
    cl_l_h = min(cl_l_h, qh_l_h * 5)
    return qh_l_h / (qh_l_h + fu * cl_l_h)
