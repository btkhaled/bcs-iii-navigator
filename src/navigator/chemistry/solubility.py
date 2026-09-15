"""Solubilité pH-dépendante Cs(pH) + solubilisation micellaire bile."""

from __future__ import annotations

from ..core.schemas import Compound


def cs_at_ph(compound: Compound, ph: float, bile_mm: float = 0.0, bile_sol_factor: float = 0.05) -> float:
    """Cs(pH) = Cs0 * (1 + sum 10^|pH-pKa|) + effet micellaire bile.

    Cs0 = solubilité intrinsèque (forme neutre). bile augmente Cs effective.
    """
    mult = 1.0
    for g in compound.pkas:
        if g.type == "acid":
            mult += 10.0 ** (ph - g.pka)
        else:
            mult += 10.0 ** (g.pka - ph)
    # borne anti-explosion (polycations)
    mult = min(mult, 1e6)
    cs = compound.cs_mg_l * mult / max(1.0, len(compound.pkas) if compound.pkas else 1.0)
    # micelles: +5% par mM de bile (paramètre formulation_effects -> ici défaut, surchargé par caller)
    cs = cs * (1.0 + bile_sol_factor * bile_mm)
    return max(cs, 1e-6)
