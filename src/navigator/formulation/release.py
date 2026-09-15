"""Release kinetics + enteric pH-dissolution."""

from __future__ import annotations

import math


def release_fraction(t_h: float, k_release_h: float) -> float:
    return 1.0 - math.exp(-max(k_release_h, 1e-9) * max(t_h, 0.0))


def enteric_open(ph: float, threshold: float = 5.5) -> float:
    # sigmoïde pH: fermé à pH bas, ouvert au-delà du seuil
    return 1.0 / (1.0 + math.exp(-8.0 * (ph - threshold)))


def effective_release(t_h: float, k_release_h: float, ph: float, enteric: float, threshold: float = 5.5) -> float:
    f = release_fraction(t_h, k_release_h)
    if enteric <= 0:
        return f
    gate = enteric_open(ph, threshold)
    # enteric=1 -> strictement gaté ; enteric<1 -> partiellement retardé
    return f * (1.0 - enteric * (1.0 - gate))
