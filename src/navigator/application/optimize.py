"""Optimize: NSGA-II sur formulation."""

from __future__ import annotations

from ..core.loader import load_compound
from ..core.schemas import Formulation
from ..core.settings import Settings
from ..engine.runner import run
from ..optimization.nsga2 import optimize


def optimize_formulation(compound_path: str, config_dir: str, space: dict, n_gen: int = 8, pop: int = 10, seed: int = 0):
    s = Settings(config_dir)
    c = load_compound(compound_path)

    def obj(g: dict) -> float:
        f = Formulation(id="opt", dose_mg=g.get("dose_mg", 500.0), m_hip=g.get("hip", 1.0),
                        c10_mm=g.get("c10_mm", 0.0), snac_mm=g.get("snac_mm", 0.0),
                        edta_mm=g.get("edta_mm", 0.0), k_release_h=g.get("k_release_h", 2.0))
        try:
            return run(c, f, s).f_oral
        except Exception:
            return 0.0

    bounds = {k: tuple(v) for k, v in space.items()}
    return optimize(bounds, obj, n_gen=n_gen, pop=pop, seed=seed)
