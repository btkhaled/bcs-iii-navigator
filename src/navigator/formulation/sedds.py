"""SEDDS: microémulsion dynamique avec sels biliaires."""

from __future__ import annotations


def sedds_solubilization(lipid_frac: float, bile_mm: float) -> float:
    # facteur multiplicatif de Cs effective
    return 1.0 + 2.0 * lipid_frac * (1.0 + 0.1 * bile_mm)
