"""Report: mise en forme résultat."""

from __future__ import annotations

from ..core.schemas import SimResult


def summarize(r: SimResult) -> dict:
    return {
        "F_oral": round(r.f_oral, 4),
        "F_oral_pct": round(r.f_oral * 100, 2),
        "Fg": round(r.f_g, 4),
        "Fh": round(r.f_h, 4),
        "AUC": round(r.auc_mg_h_l, 3),
        "Cmax": round(r.cmax_mg_l, 4),
        "dominant": r.dominant,
        "breakdown": r.breakdown.model_dump(),
        "z": round(r.z, 2),
        "mass_balance_mg": r.mass_balance_mg,
        "warnings": r.warnings,
        "run_id": r.run_id,
    }
