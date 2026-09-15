"""Bilan masse exact."""

from __future__ import annotations


def mass_balance(dose_mg: float, remaining_mg: float, absorbed_mg: float, eliminated_mg: float, fecal_mg: float, biliary_mg: float = 0.0) -> float:
    return abs(dose_mg - (remaining_mg + eliminated_mg + fecal_mg + biliary_mg))
    # absorbed est inclus dans eliminated+central ; on vérifie via remaining+elim+fecal
