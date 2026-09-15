"""Orchestration calibration: NSGA-II -> front -> knee -> candidats + rapport."""

from __future__ import annotations

import json
import os
import shutil
from datetime import datetime, timezone

from ..core.loader import load_json, write_json
from ..core.provenance import run_id
from ..inference.calibrate import Calibrator
from ..optimization.moo import optimize_moo

ACCEPT = {"rmse_train_pp": 12.0, "rmse_holdout_pp": 15.0}


def _margin_ok(x, bounds) -> bool:
    for v, (lo, hi) in zip(x, bounds):
        if not (lo + 0.2 * (hi - lo) <= v <= hi - 0.2 * (hi - lo)):
            return False
    return True


def calibrate(config_dir: str, compounds_dir: str, formulation_path: str, space_path: str,
              constraints_path: str, flit_path: str, pop: int = 40, gen: int = 40,
              seed: int = 0, workers: int | None = None, apply: bool = False) -> dict:
    t0 = datetime.now(timezone.utc).isoformat()
    import time as _t
    _start = _t.time()
    with Calibrator(config_dir, compounds_dir, formulation_path, space_path,
                    constraints_path, flit_path, workers) as cal:
        def log_fn(g, info):
            el = _t.time() - _start
            print(f"[nsga2] gen {g:3d} t={el:6.0f}s best_rmse={info['best_rmse']*100:6.2f}pp viol={info['viol']} feas={info['n_feas']}", flush=True)

        out = optimize_moo(cal.evaluate_many, cal.bounds, pop=pop, gen=gen, seed=seed,
                           seed_with=[cal.current_x], log_fn=log_fn, log_every=1)
        pop_state = out["population"]
        feas = [p for p in pop_state if p["viol"] == 0 and p["rank"] == 0]
        pool = feas or [p for p in pop_state if p["rank"] == 0]
        # knee: L2 normalisée vers utopie
        keys = range(3)
        cols = [[p["objs"][k] for p in pool] for k in keys]
        lo = [min(c) for c in cols]
        hi = [max(c) for c in cols]
        def dist(p):
            s = 0.0
            for k in keys:
                r = hi[k] - lo[k]
                v = 0.0 if r <= 0 else (p["objs"][k] - lo[k]) / r
                s += v * v
            return s
        pool_sorted = sorted(pool, key=lambda p: (dist(p), p["objs"][0]))
        # réévaluation précise des 5 meilleurs (corrige biais fast): 5*10 runs séquentiels local
        cands = []
        for p in pool_sorted[:5]:
            from ..inference.calibrate import apply_genes
            from ..core.settings import Settings as _S
            import math as _m
            g_prec = p["aux"]["genes"]
            s_prec = apply_genes(_S(config_dir), g_prec, fast=False)
            # panel précis local
            from ..core.loader import load_compound as _lc2, load_formulation as _lf2, load_json as _lj2
            import os as _os
            space = _lj2(os.path.join(os.path.dirname(config_dir.rstrip("/")), "data", "calibration_space", "search_space.json"))
            comps_prec = {m: _lc2(_os.path.join(os.path.dirname(config_dir.rstrip("/")), "data", "compounds", m + ".json")) for m in space["train"]}
            form_prec = _lf2(_os.path.join(os.path.dirname(config_dir.rstrip("/")), "data", "formulations", "native.json"))
            flit_prec = _lj2(_os.path.join(os.path.dirname(config_dir.rstrip("/")), "data", "validation", "f_lit.json"))
            cons_prec = _lj2(_os.path.join(os.path.dirname(config_dir.rstrip("/")), "data", "calibration_space", "mechanism_constraints.json"))
            from ..engine.runner import run as _run
            panel = {}
            for m in space["train"]:
                try:
                    r = _run(comps_prec[m], form_prec, s_prec)
                    panel[m] = (r.f_oral, r.dominant, r.mass_balance_mg, None)
                except Exception as e:
                    panel[m] = (0.0, "FAILED", 1e9, str(e))
            errs = [(panel[m][0] - flit_prec[m])**2 for m in space["train"]]
            rmse = _m.sqrt(sum(errs)/len(errs))
            maxerr = _m.sqrt(max(errs))
            obs = [flit_prec[m] for m in space["train"]]
            mean = sum(obs)/len(obs)
            r2 = 1 - sum(errs)/ (sum((o-mean)**2 for o in obs) or 1e-12)
            tol = max(s_prec.numerics.mass_tol_mg * 500, 0.5)
            viol = 0
            for m in space["train"]:
                f, dom, mb, err = panel[m]
                if err is not None or dom == "FAILED": viol += 3; continue
                if mb > tol: viol += 1
                if dom not in cons_prec["dominant"].get(m, [dom]): viol += 1
                band = cons_prec["bands"].get(m)
                if band and not (band[0] <= f <= band[1]): viol += 1
            aux_prec2 = {"genes": g_prec, "rmse_pp": rmse*100, "maxerr_pp": maxerr*100, "r2": r2,
                         "pred": {m: panel[m][0] for m in space["train"]}, "dom": {m: panel[m][1] for m in space["train"]}}
            objs_prec = [min(rmse,2.0), min(maxerr,2.0), min(1.0-r2,3.0)]
            h = cal.holdout(p["x"], fast=False)
            cands.append({"x": p["x"], "objs": objs_prec, "viol": viol, "aux": aux_prec2, "holdout": h})
        best = min(cands, key=lambda c: (c["holdout"]["viol"], c["holdout"]["rmse_pp"], c["aux"]["rmse_pp"]))
        gene_names = [d["name"] for d in cal.space["genes"]]
        genes = {n: v for n, v in zip(gene_names, best["x"])}
        decoded = best["aux"]["genes"]
        accept = (
            best["aux"]["rmse_pp"] < ACCEPT["rmse_train_pp"]
            and best["holdout"]["rmse_pp"] < ACCEPT["rmse_holdout_pp"]
            and best["viol"] == 0
            and best["holdout"]["viol"] == 0
            and _margin_ok(best["x"], cal.bounds)
        )
        rid = run_id({"genes": genes, "seed": seed, "pop": pop, "gen": gen})
        cdir = os.path.join(os.path.dirname(config_dir.rstrip("/")), "configs", "candidates", rid)
        os.makedirs(cdir, exist_ok=True)
        for fn in ("permeation.json", "transporters.json", "gi_physiology.json"):
            shutil.copy(os.path.join(config_dir, fn), os.path.join(cdir, fn + ".bak"))
        # candidats patchés
        perm = load_json(os.path.join(config_dir, "permeation.json"))
        perm.update({k: decoded[k] for k in ("p_scale", "t_scale", "vmax_scale", "frac_large_pore")})
        write_json(os.path.join(cdir, "permeation.json"), perm)
        tr = load_json(os.path.join(config_dir, "transporters.json"))
        tr["OCT3"]["ref"] = decoded["ref_OCT3"]
        tr["ENT"]["ref"] = decoded["ref_ENT"]
        write_json(os.path.join(cdir, "transporters.json"), tr)
        phys = load_json(os.path.join(config_dir, "gi_physiology.json"))
        for s in phys["segments"]:
            if s["name"] == "colon":
                s["fecal_bound"] = decoded["fecal_bound_colon"]
        write_json(os.path.join(cdir, "gi_physiology.json"), phys)
        report = {
            "run_id": rid, "t0": t0, "seed": seed, "pop": pop, "gen": gen,
            "accept": accept, "applied": bool(apply and accept),
            "train": {"rmse_pp": best["aux"]["rmse_pp"], "maxerr_pp": best["aux"].get("maxerr_pp"),
                      "r2": best["aux"]["r2"], "pred": best["aux"]["pred"], "dom": best["aux"]["dom"]},
            "holdout": best["holdout"],
            "genes_search_units": genes, "genes_decoded": decoded,
            "margin_ok": _margin_ok(best["x"], cal.bounds),
            "front_size": len(pool),
            "front": [{"x": p["x"], "objs": p["objs"], "viol": p["viol"]} for p in pool_sorted],
        }
        write_json(os.path.join(cdir, "report.json"), report)
        if apply and accept:
            for fn in ("permeation.json", "transporters.json", "gi_physiology.json"):
                shutil.copy(os.path.join(cdir, fn), os.path.join(config_dir, fn))
            calib = load_json(os.path.join(config_dir, "calibration.json"))
            calib.update({k: decoded[k] for k in
                          ("p_scale", "t_scale", "vmax_scale", "frac_large_pore",
                           "ref_OCT3", "ref_ENT", "fecal_bound_colon")})
            calib["notes"] = f"NSGA-II contrainte {rid} appliquée."
            write_json(os.path.join(config_dir, "calibration.json"), calib)
        return {"candidate_dir": cdir, "report": report}
