"""Ionisation Henderson-Hasselbalch polyprotique, charge nette z, f_unionized."""

from __future__ import annotations

from typing import List

from ..core.schemas import Compound


def ionization_at_ph(compound: Compound, ph: float):
    """Retourne (z, f_unionized, f_ionized). Sites indépendants."""
    z = 0.0
    f_union = 1.0
    for g in compound.pkas:
        if g.type == "acid":
            # HA <-> A- + H+ ; f_union = 1/(1+10^(pH-pKa))
            f_u = 1.0 / (1.0 + 10.0 ** (ph - g.pka))
            f_ion = 1.0 - f_u
            z += -f_ion
            f_union *= f_u
        else:
            # B+H+ <-> BH+ ; f_union = 1/(1+10^(pKa-pH))
            f_u = 1.0 / (1.0 + 10.0 ** (g.pka - ph))
            f_ion = 1.0 - f_u
            z += f_ion
            f_union *= f_u
    f_union = max(0.0, min(1.0, f_union))
    return z, f_union, 1.0 - f_union
