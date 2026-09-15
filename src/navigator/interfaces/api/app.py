"""NAVIGATOR API — single-port FastAPI (serves API + dashboard SPA)."""
from __future__ import annotations

import json
import glob
import os
import threading
import time
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    from fastapi import FastAPI, HTTPException
    from fastapi.responses import FileResponse, JSONResponse
    from fastapi.staticfiles import StaticFiles
    from pydantic import BaseModel

    # paths
    ROOT = Path(__file__).resolve().parents[4]  # /Users/kalo/NAVIGATOR
    CONFIG_DIR = ROOT / "configs"
    COMPOUNDS_DIR = ROOT / "data" / "compounds"
    FORMULATIONS_DIR = ROOT / "data" / "formulations"
    EXCIPIENTS_CATALOG = ROOT / "data" / "excipients" / "catalog.json"
    F_LIT_PATH = ROOT / "data" / "validation" / "f_lit.json"
    STATIC_DIR = Path(__file__).parent / "static"

    app = FastAPI(title="BCS III NAVIGATOR", version="0.1.0", docs_url="/api/docs", openapi_url="/api/openapi.json")

    # ---------- job store for auto-formulate ----------
    _JOBS: Dict[str, Dict[str, Any]] = {}
    _JOBS_LOCK = threading.Lock()

    # ---------- helpers ----------
    def _load_json(p: Path):
        return json.loads(p.read_text(encoding="utf-8"))

    def _settings(config_dir: str | Path | None = None):
        from ...core.settings import Settings
        d = Path(config_dir) if config_dir else CONFIG_DIR
        return Settings(str(d))

    # ---------- models ----------
    class PredictBody(BaseModel):
        compound: Dict[str, Any]
        formulation: Dict[str, Any]
        config_dir: Optional[str] = None

    class PredictByIdBody(BaseModel):
        compound_id: str
        formulation_id: str = "native"
        dose_mg: Optional[float] = None
        config_dir: Optional[str] = None

    class BatchBody(BaseModel):
        pairs: List[PredictBody]

    class HipCalcBody(BaseModel):
        dose_mg: float
        m_hip: float
        mw_drug: float = 467.52
        mw_counterion: float = 288.38

    # ---------- health ----------
    @app.get("/api/health")
    def health():
        return {"ok": True, "version": "0.1.0"}

    @app.get("/health")
    def health_legacy():
        return {"ok": True}

    # ---------- compounds ----------
    @app.get("/api/compounds")
    def list_compounds():
        out = []
        for p in sorted(COMPOUNDS_DIR.glob("*.json")):
            try:
                j = _load_json(p)
                out.append({"id": j.get("id", p.stem), "mw": j.get("mw"), "cs_mg_l": j.get("cs_mg_l"), "logp": j.get("logp"), "pkas": j.get("pkas", [])})
            except Exception:
                continue
        return out

    @app.get("/api/compounds/{cid}")
    def get_compound(cid: str):
        p = COMPOUNDS_DIR / f"{cid}.json"
        if not p.exists():
            raise HTTPException(404, f"compound {cid} not found")
        return _load_json(p)

    # ---------- formulations ----------
    @app.get("/api/formulations")
    def list_formulations():
        out = []
        for p in sorted(FORMULATIONS_DIR.rglob("*.json")):
            try:
                j = _load_json(p)
                rel = str(p.relative_to(FORMULATIONS_DIR))
                # use file stem as canonical id to avoid dup "opt" ids
                fid = p.stem
                out.append({"id": fid, "inner_id": j.get("id", fid), "path": rel, "dose_mg": j.get("dose_mg"), "m_hip": j.get("m_hip")})
            except Exception:
                continue
        return out

    @app.get("/api/formulations/{fid}")
    def get_formulation(fid: str):
        # search recursively
        for p in FORMULATIONS_DIR.rglob(f"{fid}.json"):
            return _load_json(p)
        # also try direct
        p = FORMULATIONS_DIR / f"{fid}.json"
        if p.exists():
            return _load_json(p)
        raise HTTPException(404, f"formulation {fid} not found")

    # ---------- excipients ----------
    @app.get("/api/excipients")
    def list_excipients():
        if EXCIPIENTS_CATALOG.exists():
            return _load_json(EXCIPIENTS_CATALOG)
        return {}

    # ---------- configs ----------
    @app.get("/api/configs")
    def list_configs():
        if not CONFIG_DIR.exists():
            return []
        return sorted([p.name for p in CONFIG_DIR.glob("*.json")])

    @app.get("/api/configs/{name}")
    def get_config(name: str):
        p = CONFIG_DIR / name
        if not p.exists():
            # allow without .json
            p = CONFIG_DIR / f"{name}.json"
        if not p.exists():
            raise HTTPException(404, f"config {name} not found")
        return _load_json(p)

    # ---------- predict ----------
    @app.post("/api/predict")
    def predict_route(body: PredictBody):
        from ...core.schemas import Compound, Formulation
        from ...engine.runner import run
        from ...engine.report import summarize
        try:
            c = Compound(**body.compound)
            f = Formulation(**body.formulation)
            s = _settings(body.config_dir)
            r = run(c, f, s)
            summ = summarize(r)
            # enrich with full
            summ["breakdown"] = r.breakdown.model_dump()
            summ["trace"] = r.trace
            summ["f_oral_raw"] = r.f_oral
            summ["mass_balance_mg"] = r.mass_balance_mg
            summ["warnings"] = r.warnings
            summ["z"] = r.z
            summ["f_unionized"] = r.f_unionized
            return summ
        except Exception as e:
            raise HTTPException(400, str(e))

    @app.post("/predict")
    def predict_legacy(compound_path: str, formulation_path: str, config_dir: str = "configs"):
        from ...application.predict import predict as _p
        from ...engine.report import summarize
        try:
            return summarize(_p(compound_path, formulation_path, config_dir))
        except Exception as e:
            raise HTTPException(400, str(e))

    @app.post("/api/predict/by-id")
    def predict_by_id(body: PredictByIdBody):
        from ...core.schemas import Compound, Formulation
        from ...engine.runner import run
        from ...engine.report import summarize
        cp = COMPOUNDS_DIR / f"{body.compound_id}.json"
        if not cp.exists():
            raise HTTPException(404, f"compound {body.compound_id} not found")
        # find formulation
        fp = None
        for p in FORMULATIONS_DIR.rglob(f"{body.formulation_id}.json"):
            fp = p
            break
        if fp is None:
            fp = FORMULATIONS_DIR / f"{body.formulation_id}.json"
        if not fp.exists():
            # fallback native
            fp = FORMULATIONS_DIR / "native.json"
        try:
            c = Compound(**_load_json(cp))
            f_data = _load_json(fp)
            if body.dose_mg is not None:
                f_data["dose_mg"] = body.dose_mg
            f = Formulation(**f_data)
            s = _settings(body.config_dir)
            r = run(c, f, s)
            summ = summarize(r)
            summ["breakdown"] = r.breakdown.model_dump()
            summ["trace"] = r.trace
            return summ
        except Exception as e:
            raise HTTPException(400, str(e))

    @app.post("/api/batch")
    def batch_route(body: BatchBody):
        from ...core.schemas import Compound, Formulation
        from ...engine.runner import run
        from ...engine.report import summarize
        out = []
        s = None
        # reuse settings if same config_dir
        settings_cache: Dict[str, Any] = {}
        for pair in body.pairs:
            try:
                cfg = pair.config_dir or str(CONFIG_DIR)
                if cfg not in settings_cache:
                    settings_cache[cfg] = _settings(cfg)
                s = settings_cache[cfg]
                c = Compound(**pair.compound)
                f = Formulation(**pair.formulation)
                r = run(c, f, s)
                summ = summarize(r)
                summ["breakdown"] = r.breakdown.model_dump()
                out.append({"ok": True, "result": summ})
            except Exception as e:
                out.append({"ok": False, "error": str(e)})
        return out

    @app.post("/api/hip-calc")
    def hip_calc(body: HipCalcBody):
        mw_drug = body.mw_drug
        mw_counter = body.mw_counterion
        mol = body.dose_mg / mw_drug  # mmol drug (mg / g/mol = mmol)
        # Actually mg / (g/mol) = mmol /1000? Let's do mg -> g: dose_mg/1000 / mw *1000 = dose_mg/mw mmol
        # So mmol = dose_mg / mw
        ratios = [1, 2, 5]
        masses = {}
        for r in ratios:
            masses[f"{r}:1"] = round(mol * r * mw_counter, 1)
        return {"dose_mg": body.dose_mg, "m_hip": body.m_hip, "mmol_drug": round(mol, 4), "sds_masses_mg": masses, "current_m_hip_mass": round(mol * body.m_hip * mw_counter, 1) if body.m_hip > 1 else 0}

    @app.get("/api/validate")
    def validate_route():
        from ...core.loader import load_compound, load_formulation, load_json
        from ...engine.runner import run
        from ...validation.metrics import rmse, r2
        s = _settings()
        # use native formulation
        f = load_formulation(str(FORMULATIONS_DIR / "native.json"))
        obs = load_json(str(F_LIT_PATH))
        pred: Dict[str, float] = {}
        details = []
        for p in sorted(COMPOUNDS_DIR.glob("*.json")):
            try:
                c = load_compound(str(p))
                if c.id not in obs:
                    continue
                r = run(c, f, s)
                pred[c.id] = r.f_oral
                details.append({"id": c.id, "pred": round(r.f_oral, 4), "obs": obs[c.id], "err_pp": round(abs(r.f_oral - obs[c.id])*100, 2), "dominant": r.dominant})
            except Exception as e:
                details.append({"id": p.stem, "error": str(e)})
        return {"pred": pred, "obs": obs, "rmse_pp": round(rmse(pred, obs)*100, 2) if pred else None, "r2": round(r2(pred, obs), 3) if pred else None, "details": sorted(details, key=lambda x: x.get("obs", 0))}

    @app.get("/api/optimization/front")
    def optimization_front(space: str = "A"):
        # return precomputed fronts if exist
        base = ROOT / "data" / "formulations" / "optimized"
        mapping = {"A": "tobramycin_A_best.json", "B": "tobramycin_B_best.json", "full": "tobramycin_nsga2_best.json"}
        fname = mapping.get(space, mapping["A"])
        p = base / fname
        if p.exists():
            j = _load_json(p)
            # also show list of all optimized
            all_opts = []
            for q in base.glob("*.json"):
                try:
                    all_opts.append(_load_json(q))
                except Exception:
                    pass
            return {"space": space, "best": j, "all": all_opts}
        raise HTTPException(404, f"front {space} not found")

    # ---------- auto-formulate NSGA-II ----------
    class AutoFormBody(BaseModel):
        compound_id: str = "tobramycin"
        space: str = "A"  # A, B, full, search
        pop: int = 24
        gen: int = 15
        seed: int = 0
        workers: int = 4

    _SPACE_MAP = {
        "A": "tobramycin_A.json",
        "B": "tobramycin_B.json",
        "full": "tobramycin_full.json",
        "search": "tobramycin_search.json",
    }

    def _auto_form_run(job_id: str, compound_id: str, space: str, pop: int, gen: int, seed: int, workers: int):
        try:
            with _JOBS_LOCK:
                _JOBS[job_id]["status"] = "running"
                _JOBS[job_id]["progress"] = {"gen": 0, "total": gen, "best_F": 0, "elapsed": 0}
            t0 = time.time()
            # resolve paths
            compound_path = COMPOUNDS_DIR / f"{compound_id}.json"
            if not compound_path.exists():
                raise FileNotFoundError(f"compound {compound_id} not found")
            space_file = _SPACE_MAP.get(space, _SPACE_MAP["A"])
            space_path = ROOT / "data" / "optimization_space" / space_file
            if not space_path.exists():
                space_path = ROOT / "data" / "optimization_space" / "tobramycin_A.json"
            # clamp pop/gen
            pop = max(12, min(48, pop))
            gen = max(5, min(30, gen))
            workers = max(1, min(6, workers))

            from ...application.formulation_opt import run_formulation_optim

            # progress callback via closure
            def _log_fn(g, info):
                best = -info.get("best_rmse", 0)  # -F
                feas = info.get("n_feas", 0)
                elapsed = int(time.time() - t0)
                with _JOBS_LOCK:
                    _JOBS[job_id]["progress"] = {"gen": g, "total": gen, "best_F": round(best*100, 2), "n_feas": feas, "elapsed": elapsed}
                    # also store logs
                    _JOBS[job_id].setdefault("logs", []).append(f"gen {g} best_F {best*100:.2f}% feas {feas} t={elapsed}s")

            # monkey-patch run_formulation_optim to inject log_fn? It already accepts log but we need to capture.
            # We call directly with custom evaluate and use optimize_moo log, but easiest: call run_formulation_optim but it internally uses log_fn with gen.
            # To get progress, we wrap the call: run_formulation_optim will call log_fn each 2 gens, we intercept via threading.
            # Simpler: call run_formulation_optim and update progress after; but we want per-gen. We'll call it with our log.
            # run_formulation_optim signature pop/gen/seed/workers, but not log_fn exposed. So we replicate its logic with custom log.
            # For quick path, we just call it and poll elapsed time in separate thread -> simpler: update progress via elapsed only.
            # Instead, we run it in this thread and update progress via elapsed polling in another tiny thread? Simpler: just run and update at end.
            # Let's run with a wrapper that captures prints? Instead, we manually run formulation_opt logic with progress.
            # Use run_formulation_optim but it will log via print, not via _log_fn. We'll do a lightweight reimplementation:
            # If pop/gen small, we can just call run_formulation_optim and update progress every 2s via separate thread.

            # Start progress pinger
            stop_ping = threading.Event()
            def _pinger():
                while not stop_ping.wait(1.0):
                    elapsed = int(time.time() - t0)
                    with _JOBS_LOCK:
                        if _JOBS[job_id]["status"] != "running":
                            break
                        cur = _JOBS[job_id]["progress"]
                        cur["elapsed"] = elapsed
            ping_t = threading.Thread(target=_pinger, daemon=True)
            ping_t.start()

            result = run_formulation_optim(
                config_dir=str(CONFIG_DIR),
                compound_path=str(compound_path),
                space_path=str(space_path),
                pop=pop, gen=gen, seed=seed, workers=workers,
            )
            stop_ping.set()
            # prepare front for API
            front = []
            for p in result.get("front", [])[:12]:
                aux = p.get("aux", {})
                form = aux.get("form", {})
                front.append({
                    "F_pct": round(aux.get("f", 0)*100, 2),
                    "eff": round(aux.get("eff", 0)*100, 2),
                    "dominant": aux.get("dom", "?"),
                    "genes": aux.get("genes", {}),
                    "formulation": form,
                    "objs": p.get("objs", []),
                    "rank": p.get("rank"),
                })
            with _JOBS_LOCK:
                _JOBS[job_id]["status"] = "done"
                _JOBS[job_id]["result"] = {"front": front, "space": space, "compound_id": compound_id, "pop": pop, "gen": gen, "elapsed": int(time.time()-t0), "keys": result.get("keys", []), "n_front": len(front)}
                _JOBS[job_id]["progress"]["gen"] = gen
                _JOBS[job_id]["progress"]["best_F"] = front[0]["F_pct"] if front else 0
        except Exception as e:
            import traceback
            with _JOBS_LOCK:
                _JOBS[job_id]["status"] = "error"
                _JOBS[job_id]["error"] = str(e)
                _JOBS[job_id]["traceback"] = traceback.format_exc()

    @app.post("/api/auto-formulate")
    def auto_formulate(body: AutoFormBody):
        job_id = uuid.uuid4().hex[:8]
        with _JOBS_LOCK:
            # limit concurrent running
            running = sum(1 for j in _JOBS.values() if j["status"] in ("queued", "running"))
            if running >= 2:
                raise HTTPException(429, "too many jobs running, try later")
            _JOBS[job_id] = {"job_id": job_id, "status": "queued", "created": time.time(), "body": body.model_dump(), "progress": {"gen": 0, "total": body.gen, "best_F": 0, "elapsed": 0}, "logs": []}
        t = threading.Thread(target=_auto_form_run, args=(job_id, body.compound_id, body.space, body.pop, body.gen, body.seed, body.workers), daemon=True)
        t.start()
        return {"job_id": job_id, "status": "queued", "body": body.model_dump()}

    @app.get("/api/auto-formulate")
    def list_auto_jobs():
        with _JOBS_LOCK:
            return sorted(_JOBS.values(), key=lambda x: x["created"], reverse=True)[:20]

    @app.get("/api/auto-formulate/{job_id}")
    def get_auto_job(job_id: str):
        with _JOBS_LOCK:
            j = _JOBS.get(job_id)
            if not j:
                raise HTTPException(404, f"job {job_id} not found")
            return j

    @app.delete("/api/auto-formulate/{job_id}")
    def del_auto_job(job_id: str):
        with _JOBS_LOCK:
            if job_id in _JOBS:
                del _JOBS[job_id]
                return {"ok": True}
            raise HTTPException(404, f"job {job_id} not found")

    # ---------- static SPA (single port) ----------
    # Mount static dir for assets if exists
    if STATIC_DIR.exists():
        app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")

    @app.get("/")
    def root():
        idx = STATIC_DIR / "index.html"
        if idx.exists():
            return FileResponse(str(idx))
        return JSONResponse({"msg": "NAVIGATOR API — dashboard not built yet. See /api/docs"})

    @app.get("/{full_path:path}")
    def spa_fallback(full_path: str):
        # serve SPA fallback for non-api routes
        if full_path.startswith("api/"):
            raise HTTPException(404, "not found")
        idx = STATIC_DIR / "index.html"
        maybe = STATIC_DIR / full_path
        if maybe.exists() and maybe.is_file():
            return FileResponse(str(maybe))
        if idx.exists():
            return FileResponse(str(idx))
        raise HTTPException(404, "not found")

except ImportError as e:
    # fallback if fastapi not installed
    app = None  # type: ignore
