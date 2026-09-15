"""Formulation NSGA-II — parallèle, multi-objectif F oral tobramycine."""

from __future__ import annotations

import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

from ..core.loader import load_compound, load_formulation, load_json
from ..core.schemas import Formulation
from ..core.settings import Settings
from ..engine.runner import run
from ..optimization.moo import optimize_moo

try:
    from threadpoolctl import threadpool_limits
    threadpool_limits(limits=1)
except ImportError:
    pass

_FMTX: dict = {}


def _finit(config_dir: str, compound_path: str, bounds=None, keys=None):
    _FMTX["settings"] = Settings(config_dir)
    _FMTX["compound"] = load_compound(compound_path)
    if bounds is not None:
        _FMTX["bounds"] = bounds
    if keys is not None:
        _FMTX["keys"] = keys


def _feval(args):
    x, bounds, keys = args
    # thread limit per worker
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(limits=1)
    except ImportError:
        pass
    # clip to bounds (garde-fou SBX/mutation numérique)
    # space not in worker; bounds passed, reconstruct per key order
    # keys order == bounds order from caller, so clip via bounds
    g = {}
    for k, v, (lo, hi) in zip(keys, x, _FMTX.get("bounds", [(1,32)]*len(keys))):
        g[k] = min(hi, max(lo, v))
    f = Formulation(
        id="opt",
        dose_mg=g.get("dose_mg", 500.0),
        m_hip=g.get("hip", 1.0),
        hip_counterion="SDS" if g.get("hip", 1.0) > 2 else None,
        lipid_frac=g.get("lipid_frac", 0.0),
        particle_nm=g.get("particle_nm", 150.0),
        c10_mm=g.get("c10_mm", 0.0),
        snac_mm=g.get("snac_mm", 0.0),
        edta_mm=g.get("edta_mm", 0.0),
        chitosan=g.get("chitosan", 0.0),
        k_release_h=g.get("k_release_h", 2.0),
        enteric=g.get("enteric", 0.0),
        f_free=0.95 if g.get("hip", 1.0) > 5 else 1.0,
        labrasol_mm=g.get("labrasol_mm", 0.0),
        nacdc_mm=g.get("nacdc_mm", 0.0),
        pip_mm=g.get("pip_mm", 0.0),
        tpgs_mm=g.get("tpgs_mm", 0.0),
        denovo_p_trans=g.get("denovo_p_trans", 1.0),
        denovo_p_para=g.get("denovo_p_para", 1.0),
        denovo_mucus=g.get("denovo_mucus", 1.0),
    )
    try:
        r = run(_FMTX["compound"], f, _FMTX["settings"])
        # 3 objectifs à minimiser: -F, -Cmax/CMI-like (approx), +dose (on veut F/dose)
        # pour Pareto: on maximise F et F/dose
        f_oral = r.f_oral
        eff = f_oral * 500.0 / max(f.dose_mg, 1.0)  # F normalisé dose
        # minimize: -F, -eff  (2 objectifs)
        return ([-f_oral, -eff], 0, {"genes": g, "f": f_oral, "eff": eff, "dom": r.dominant, "form": f.model_dump()})
    except Exception as e:
        return ([0.0, 0.0], 1, {"genes": g, "error": str(e)})


def run_formulation_optim(config_dir: str = "/Users/kalo/NAVIGATOR/configs",
                          compound_path: str = "/Users/kalo/NAVIGATOR/data/compounds/tobramycin.json",
                          space_path: str = "/Users/kalo/NAVIGATOR/data/optimization_space/tobramycin_full.json",
                          pop: int = 36, gen: int = 30, seed: int = 0, workers: int = 6):
    import random, time
    space = load_json(space_path)
    keys = list(space.keys())
    bounds = [tuple(space[k]) for k in keys]
    t0 = time.time()

    with ProcessPoolExecutor(max_workers=workers, initializer=_finit, initargs=(config_dir, compound_path, bounds, keys)) as pool:
        def evaluate_many(xs):
            futs = {pool.submit(_feval, (x, bounds, keys)): i for i, x in enumerate(xs)}
            res = [None] * len(xs)
            for fut in as_completed(futs):
                i = futs[fut]
                res[i] = fut.result()
            return res

        def log_fn(g, info):
            # info has best_rmse etc for moo; adapt
            feas = info.get("n_feas", 0)
            # find best F among feas
            print(f"[form-opt] gen {g:3d} t={time.time()-t0:5.0f}s feas={feas} best_F={-info.get('best_rmse',0)*100:.2f}%", flush=True)

        # Wrap evaluate_many to match moo signature: returns [(objs, viol, aux)]
        # moo expects best_rmse = min objs[0]; we use -F as objs[0]
        out = optimize_moo(evaluate_many, bounds, pop=pop, gen=gen, seed=seed, seed_with=None, log_fn=log_fn, log_every=2)

    # Pareto front: rank 0, viol 0
    front = [p for p in out["population"] if p["rank"] == 0 and p["viol"] == 0]
    # also re-evaluate front precisely (already precise, same numerics)
    front_sorted = sorted(front, key=lambda p: -p["aux"]["f"])
    return {"front": front_sorted, "population": out["population"], "keys": keys, "bounds": bounds}
