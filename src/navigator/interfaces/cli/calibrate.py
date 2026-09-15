"""CLI calibration NSGA-II contrainte."""
from __future__ import annotations

import json
import sys


def main():
    from ...application.calibrate import calibrate

    base = "/Users/kalo/NAVIGATOR"
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--pop", type=int, default=40)
    ap.add_argument("--gen", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    out = calibrate(f"{base}/configs", f"{base}/data/compounds",
                    f"{base}/data/formulations/native.json",
                    f"{base}/data/calibration_space/search_space.json",
                    f"{base}/data/calibration_space/mechanism_constraints.json",
                    f"{base}/data/validation/f_lit.json",
                    pop=a.pop, gen=a.gen, seed=a.seed, workers=a.workers, apply=a.apply)
    r = out["report"]
    print(json.dumps({
        "candidate_dir": out["candidate_dir"],
        "accept": r["accept"], "applied": r["applied"],
        "train_rmse_pp": round(r["train"]["rmse_pp"], 2),
        "train_r2": round(r["train"]["r2"], 3),
        "holdout_rmse_pp": round(r["holdout"]["rmse_pp"], 2),
        "holdout_viol": r["holdout"]["viol"],
        "genes": {k: round(v, 4) for k, v in r["genes_decoded"].items()},
    }, indent=2))


if __name__ == "__main__":
    main()
