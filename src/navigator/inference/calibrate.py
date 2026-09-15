"""Calibration NSGA-II contrainte du moteur (train/hold-out, mécanismes éliminatoires)."""

from __future__ import annotations

import copy
import math
import os
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List, Tuple

from ..core.loader import load_compound, load_formulation, load_json
from ..core.settings import Settings
from ..engine.runner import run

import os as _os

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    _os.environ.setdefault(_k, "1")
try:
    from threadpoolctl import threadpool_limits as _tpl
    _tpl(limits=1)
except ImportError:
    pass

_CTX: dict = {}
_WCACHE: dict = {}
_FAST = {"t_end_h": 24.0, "rtol": 1e-4, "atol": 1e-7, "n_eval": 120}


def _init(config_dir: str, compounds_dir: str, formulation_path: str, space_path: str,
          constraints_path: str, flit_path: str):
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(limits=1)
    except ImportError:
        pass
    from ..core.loader import load_compound as _lc, load_formulation as _lf, load_json as _lj
    from ..core.settings import Settings as _S

    space = _lj(space_path)
    compounds = {}
    for mol in space["train"] + space["holdout"]:
        compounds[mol] = _lc(os.path.join(compounds_dir, mol + ".json"))
    _CTX.update(
        settings=_S(config_dir),
        compounds=compounds,
        formulation=_lf(formulation_path),
        space=space,
        constraints=_lj(constraints_path),
        flit=_lj(flit_path),
    )


def gene_defs() -> List[dict]:
    return _CTX["space"]["genes"]


def decode(x: List[float]) -> Dict[str, float]:
    g = {}
    for v, d in zip(x, _CTX["space"]["genes"]):
        g[d["name"]] = (10.0 ** v) if d["scale"] == "log10" else v
    return g


def encode(current: Dict[str, float]) -> List[float]:
    out = []
    for d in _CTX["space"]["genes"]:
        v = current[d["name"]]
        out.append(math.log10(v) if d["scale"] == "log10" else v)
    return out


def apply_genes(settings: Settings, g: Dict[str, float], fast: bool = False) -> Settings:
    s = copy.deepcopy(settings)
    s.permeation.p_scale = g["p_scale"]
    s.permeation.t_scale = g["t_scale"]
    s.permeation.vmax_scale = g["vmax_scale"]
    s.permeation.frac_large_pore = g["frac_large_pore"]
    s.transporters["OCT3"]["ref"] = g["ref_OCT3"]
    s.transporters["ENT"]["ref"] = g["ref_ENT"]
    for seg in s.physiology.segments:
        if seg.name == "colon":
            seg.fecal_bound = g["fecal_bound_colon"]
    if fast:
        s.numerics.t_end_h = _FAST["t_end_h"]
        s.numerics.rtol = _FAST["rtol"]
        s.numerics.atol = _FAST["atol"]
        s.numerics.n_eval = _FAST["n_eval"]
    return s


def current_genes_dict(settings: Settings) -> Dict[str, float]:
    colon_fb = next(s.fecal_bound for s in settings.physiology.segments if s.name == "colon")
    return {
        "p_scale": settings.permeation.p_scale,
        "t_scale": settings.permeation.t_scale,
        "vmax_scale": settings.permeation.vmax_scale,
        "ref_OCT3": settings.transporters["OCT3"]["ref"],
        "ref_ENT": settings.transporters["ENT"]["ref"],
        "frac_large_pore": settings.permeation.frac_large_pore,
        "fecal_bound_colon": colon_fb,
    }


def _panel(settings: Settings, mols: List[str]):
    """Retourne {id: (F, dominant, mass_bal)} avec garde-fous."""
    out = {}
    for mol in mols:
        try:
            r = run(_CTX["compounds"][mol], _CTX["formulation"], settings)
            out[mol] = (r.f_oral, r.dominant, r.mass_balance_mg, None)
        except Exception as e:  # noqa: BLE001
            out[mol] = (0.0, "FAILED", 1e9, f"{type(e).__name__}: {e}")
    return out


def _metrics(panel: dict, mols: List[str]):
    flit = _CTX["flit"]
    errs = [(panel[m][0] - flit[m]) ** 2 for m in mols]
    rmse = math.sqrt(sum(errs) / len(errs))
    maxerr = math.sqrt(max(errs))
    obs = [flit[m] for m in mols]
    mean = sum(obs) / len(obs)
    ss_tot = sum((o - mean) ** 2 for o in obs)
    ss_res = sum(errs)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return rmse, maxerr, r2


def _violations(panel: dict, mols: List[str], settings: Settings) -> int:
    cons = _CTX["constraints"]
    n = 0
    tol = max(settings.numerics.mass_tol_mg * 500.0, 0.5)
    for m in mols:
        f, dom, mb, err = panel[m]
        if err is not None or dom == "FAILED":
            n += 3
            continue
        if mb > tol:
            n += 1
        allowed = cons["dominant"].get(m)
        if allowed is not None and dom not in allowed:
            n += 1
        band = cons["bands"].get(m)
        if band is not None and not (band[0] <= f <= band[1]):
            n += 1
    return n


def _eval_key(x: List[float]) -> tuple:
    return tuple(round(v, 6) for v in x)


def eval_task(x: List[float], fast: bool = True):
    key = (_eval_key(x), fast)
    if key in _WCACHE:
        return _WCACHE[key]
    g = decode(x)
    s = apply_genes(_CTX["settings"], g, fast=fast)
    train = _CTX["space"]["train"]
    panel = _panel(s, train)
    rmse, maxerr, r2 = _metrics(panel, train)
    viol = _violations(panel, train, s)
    if not math.isfinite(r2):
        r2 = -1.0
    objs = [min(rmse, 2.0), min(maxerr, 2.0), min(1.0 - r2, 3.0)]
    aux = {"genes": g, "rmse_pp": rmse * 100, "maxerr_pp": maxerr * 100, "r2": r2,
           "pred": {m: panel[m][0] for m in train}, "dom": {m: panel[m][1] for m in train}}
    res = (objs, viol, aux)
    _WCACHE[key] = res
    return res


def reeval_precise(x: List[float]):
    return eval_task(x, fast=False)


class Calibrator:
    def __init__(self, config_dir: str, compounds_dir: str, formulation_path: str,
                 space_path: str, constraints_path: str, flit_path: str, workers: int | None = None):
        self.args = (config_dir, compounds_dir, formulation_path, space_path, constraints_path, flit_path)
        self.workers = workers or max(1, (os.cpu_count() or 4) - 1)
        self.space = load_json(space_path)
        self.bounds = [(d["lo"], d["hi"]) for d in self.space["genes"]]
        self._cache: dict = {}
        self._pool: ProcessPoolExecutor | None = None

    def _new_pool(self):
        if self._pool is not None:
            self._pool.shutdown(wait=False, cancel_futures=True)
        self._pool = ProcessPoolExecutor(max_workers=self.workers, initializer=_init, initargs=self.args)

    def __enter__(self):
        self._new_pool()
        # encode solution courante (seed) — charge settings en local
        s = Settings(self.args[0])
        cur = current_genes_dict(s)
        self.current_x = [math.log10(cur[d["name"]]) if d["scale"] == "log10" else cur[d["name"]]
                          for d in self.space["genes"]]
        return self

    def __exit__(self, *a):
        if self._pool is not None:
            self._pool.shutdown(wait=True)

    TASK_TIMEOUT_S = 300.0

    def evaluate_many(self, xs: List[List[float]]):
        from concurrent.futures import as_completed

        keys = [_eval_key(x) for x in xs]
        pending = {}
        futs = {}
        for i, (x, k) in enumerate(zip(xs, keys)):
            if k not in self._cache:
                fut = self._pool.submit(eval_task, x)
                futs[fut] = (i, k, x)
        try:
            for fut in as_completed(futs, timeout=self.TASK_TIMEOUT_S * max(1, len(futs) // max(1, self.workers) + 1)):
                i, k, x = futs.pop(fut)
                try:
                    self._cache[k] = fut.result(timeout=1)
                except Exception:  # noqa: BLE001 — timeout ou crash worker: échec pénalisé
                    self._cache[k] = ([2.0, 2.0, 3.0], 99, {"genes": {}, "error": "task-timeout"})
        except Exception:  # noqa: BLE001 — global timeout: recycle le pool, marque le reste en échec
            for fut, (i, k, x) in futs.items():
                fut.cancel()
                self._cache[k] = ([2.0, 2.0, 3.0], 99, {"genes": {}, "error": "pool-timeout"})
            self._new_pool()
        return [self._cache[k] for k in keys]

    def holdout(self, x: List[float], fast: bool = False):
        """Évalue le hold-out pour un individu (processus local, 4 runs)."""
        s = Settings(self.args[0])
        argentinos = load_json(self.args[3])
        g = {}
        for v, d in zip(x, argentinos["genes"]):
            g[d["name"]] = (10.0 ** v) if d["scale"] == "log10" else v
        s = apply_genes(s, g, fast=fast)
        compounds = {m: load_compound(os.path.join(self.args[1], m + ".json")) for m in argentinos["holdout"]}
        f = load_formulation(self.args[2])
        flit = load_json(self.args[5])
        cons = load_json(self.args[4])
        pred, dom = {}, {}
        for m in argentinos["holdout"]:
            try:
                r = run(compounds[m], f, s)
                pred[m], dom[m] = r.f_oral, r.dominant
            except Exception:  # noqa: BLE001
                pred[m], dom[m] = 0.0, "FAILED"
        errs = [(pred[m] - flit[m]) ** 2 for m in argentinos["holdout"]]
        rmse = math.sqrt(sum(errs) / len(errs))
        obs = [flit[m] for m in argentinos["holdout"]]
        mean = sum(obs) / len(obs)
        ss_tot = sum((o - mean) ** 2 for o in obs) or 1e-12
        r2 = 1 - sum(errs) / ss_tot
        viol = sum(
            0 if (dom.get(m) in cons["dominant"].get(m, [dom.get(m)])
                  and cons["bands"].get(m, [-1, 2])[0] <= pred[m] <= cons["bands"].get(m, [-1, 2])[1])
            else 1
            for m in argentinos["holdout"]
        )
        return {"rmse_pp": rmse * 100, "r2": r2, "pred": pred, "dom": dom, "viol": viol}
