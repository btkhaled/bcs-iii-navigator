"""Composition: mapping formulation JSON -> effets quantitatifs."""

from __future__ import annotations

from dataclasses import dataclass

from ..core.schemas import Formulation


@dataclass
class FormulationEffectsValues:
    m_hip: float
    f_hip: float
    f_trans_extra: float
    lipid_factor: float
    snac_factor: float
    chitosan: float
    k_release_h: float
    enteric: float
    f_free: float


def hip_fraction(m_hip: float) -> float:
    # f_HIP croît avec M_HIP, saturée
    if m_hip <= 1.0:
        return 0.0
    return min(0.95, (m_hip - 1.0) / (m_hip + 4.0))


def evaluate(formulation: Formulation, snac_per_mm: float = 0.02, lipid_per_frac: float = 1.0) -> FormulationEffectsValues:
    f_hip = hip_fraction(formulation.m_hip)
    lipid_factor = 1.0 + lipid_per_frac * formulation.lipid_frac
    snac_factor = 1.0 + snac_per_mm * formulation.snac_mm
    return FormulationEffectsValues(
        m_hip=formulation.m_hip,
        f_hip=f_hip,
        f_trans_extra=f_hip * 0.5,
        lipid_factor=lipid_factor,
        snac_factor=snac_factor,
        chitosan=formulation.chitosan,
        k_release_h=formulation.k_release_h,
        enteric=formulation.enteric,
        f_free=formulation.f_free,
    )
