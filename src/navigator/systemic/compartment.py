"""Compartiment systémique 1-compartment."""

from __future__ import annotations


def elimination_rate_h(cl_l_h: float, vd_l: float) -> float:
    if vd_l <= 0:
        return 0.0
    return cl_l_h / vd_l
