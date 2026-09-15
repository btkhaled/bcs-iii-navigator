"""Effets repas: modulation vidange / bile / pH (simplifié mécanistique)."""

from __future__ import annotations


def fed_factors(fed: bool) -> dict:
    if not fed:
        return {"bile_mult": 1.0, "transit_mult": 1.0, "volume_add_l": 0.0}
    return {"bile_mult": 3.0, "transit_mult": 1.3, "volume_add_l": 0.3}
