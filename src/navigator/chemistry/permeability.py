"""QSPR P_int0 (propriété physico-chimique, pas prédiction F).

P_int0 = 10^(a + b*logP_clamp) ; a,b depuis configs/permeation.json
"""

from __future__ import annotations


def pint0_cm_s(logp: float, a: float, b: float) -> float:
    lp = max(-6.0, min(5.0, logp))
    return 10.0 ** (a + b * lp)
