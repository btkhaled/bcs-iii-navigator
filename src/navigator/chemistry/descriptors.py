"""Descripteurs: rayon Bohley, volume. Pur mécanistique."""

from __future__ import annotations

import math


def bohley_radius_nm(mw: float) -> float:
    """rayon = 0.091 * MW^0.33 (nm)."""
    if mw <= 0:
        raise ValueError("MW must be > 0")
    return 0.091 * (mw ** 0.33)


def molecular_volume_nm3(radius_nm: float) -> float:
    return (4.0 / 3.0) * math.pi * radius_nm ** 3
