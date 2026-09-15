"""PK: C(t), AUC, Cmax."""

from __future__ import annotations

import numpy as np


def pk_metrics(t_h: np.ndarray, a_central_mg: np.ndarray, vd_l: float):
    c = a_central_mg / max(vd_l, 1e-9)
    auc = float(np.trapezoid(c, t_h)) if len(t_h) > 1 else 0.0
    cmax = float(np.max(c)) if len(c) else 0.0
    return auc, cmax, c
