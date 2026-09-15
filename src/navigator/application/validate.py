"""Validate: run panel vs F_lit."""

from __future__ import annotations

from ..core.loader import load_compound, load_formulation, load_json
from ..core.settings import Settings
from ..engine.runner import run
from ..validation.metrics import r2, rmse


def validate(config_dir: str, compounds_dir: str, formulation_path: str, f_lit_path: str):
    import glob, os

    s = Settings(config_dir)
    f = load_formulation(formulation_path)
    obs = load_json(f_lit_path)
    pred = {}
    for cp in sorted(glob.glob(os.path.join(compounds_dir, "*.json"))):
        c = load_compound(cp)
        if c.id not in obs:
            continue
        pred[c.id] = run(c, f, s).f_oral
    return {"pred": pred, "obs": obs, "rmse": rmse(pred, obs), "r2": r2(pred, obs)}
