"""Métriques validation: RMSE, R2, AFE, AAFE."""

from __future__ import annotations

import math
from typing import Dict


def rmse(pred: Dict[str, float], obs: Dict[str, float]) -> float:
    ks = [k for k in pred if k in obs]
    if not ks:
        return float("nan")
    return math.sqrt(sum((pred[k] - obs[k]) ** 2 for k in ks) / len(ks))


def r2(pred: Dict[str, float], obs: Dict[str, float]) -> float:
    ks = [k for k in pred if k in obs]
    if len(ks) < 2:
        return float("nan")
    m = sum(obs[k] for k in ks) / len(ks)
    ss_tot = sum((obs[k] - m) ** 2 for k in ks)
    ss_res = sum((obs[k] - pred[k]) ** 2 for k in ks)
    return 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
