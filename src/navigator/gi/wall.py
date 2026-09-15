"""Paroi entérocytaire: Fg + transit fini."""

from __future__ import annotations


def fg_gut(fu: float, clint_gut_ul_min_mg: float, tgut_h: float, vgut_l: float) -> float:
    """Fg = 1 / (1 + fu*CLint*Tgut/Vgut). CLint en uL/min/mg -> L/h (facteur microsomes)."""
    if clint_gut_ul_min_mg <= 0:
        return 1.0
    # conversion simplifiée: 1 uL/min/mg * 40 mg prot/g foie... pour gut on utilise facteur 10
    cl_l_h = clint_gut_ul_min_mg * 1e-6 * 60.0 * 10.0
    x = fu * cl_l_h * tgut_h / max(vgut_l, 1e-9)
    return 1.0 / (1.0 + x)
