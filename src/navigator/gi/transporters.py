"""Transporters: façade expression segmentaire."""

from __future__ import annotations

from typing import Dict

from ..core.settings import Settings


def expr_for(settings: Settings, segment_name: str, transporter: str) -> float:
    for s in settings.physiology.segments:
        if s.name == segment_name:
            return float(s.expr.get(transporter, 1.0))
    return 1.0


def vmax_km(settings: Settings, transporter: str):
    d = settings.transporters.get(transporter, {"vmax_pmol_cm2_min": 0.0, "km_um": 1.0})
    ref = float(d.get("ref", 1.0))
    return float(d["vmax_pmol_cm2_min"]) * ref, float(d["km_um"])
