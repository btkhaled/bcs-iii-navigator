# BCS III NAVIGATOR

**A mechanistic PBPK framework for exploring oral delivery of poorly permeable drugs.**

> Predicts `F_oral` (fraction reaching systemic circulation) from chemistry, formulation and gut physiology — parameter-driven, fully reproducible, mass-balance residual < 10⁻¹⁰.

[![Version](https://img.shields.io/badge/version-0.1.0-0B1426)](https://github.com/btkhaled/bcs-iii-navigator)
[![Python](https://img.shields.io/badge/python-3.11%2B-0E7A6E)](https://www.python.org)
[![PBPK](https://img.shields.io/badge/model-mechanistic-informational)](https://github.com/btkhaled/bcs-iii-navigator)
[![Private](https://img.shields.io/badge/visibility-private-lightgrey)](https://github.com/btkhaled/bcs-iii-navigator)

**Author:** Khaled Ben Taieb — bt.khaled@gmail.com  
**Poster:** [`poster.html`](./poster.html) — A4 scientific poster (print-ready)

---

## 1. Problem

BCS Class III: **high solubility / low intestinal permeability** (FDA). Archetypal case study: **tobramycin** (MW 467.52 g·mol⁻¹, XlogP3-AA = −6.2, pKₐ 6.7–9.1, mean charge ≈ +4.43 at pH 6.5) — high solubility alone does not establish high permeability nor BCS classification.

The framework separates **literature-derived / calculated / assumed / predicted** quantities with explicit uncertainty — *mechanism first, quantification second*.

## 2. What it does

```
Compound + Formulation + Physiology  →  SimResult
  compound.json   formulation.json     configs/      →  { F_oral, Fg, Fh, AUC, Cmax, breakdown, trace, run_id }
```

**Case study — tobramycin (500 mg):**

| Formulation | F_oral | Cmax (mg·L⁻¹) | Dominant path |
|-------------|--------|---------------|---------------|
| Native | **1.38%** | 0.13 | para |
| **Optimized A** — excipients with regulatory precedents (SDS HIP m=32, C10 47 mM, Labrasol 18 mM, NaCDC 14 mM, PIP 9.6 mM, TPGS 4.7 mM, chitosan 0.8, enteric 0.66 @ pH 5.5) | **18.62%** | **3.83** | trans 68% |

> Within the modeled formulation space, conventional excipient combinations increase predicted F to 18.6%, but systemic exposure remains limited (see dose–Cmax table in poster). SDS:tobramycin 2:1 → 617 mg SDS per 500 mg drug.

```python
from navigator.application.predict import predict  # src/navigator/application/predict.py:12

r = predict("data/compounds/tobramycin.json", "data/formulations/native.json", "configs/")
# r.f_oral == 0.0138, r.breakdown, r.mass_balance_mg < 1e-10, r.trace, r.run_id

r_opt = predict("data/compounds/tobramycin.json", "data/formulations/optimized/tobramycin_A_best.json", "configs/")
# r_opt.f_oral == 0.1862, r_opt.dominant == "trans"
```

## 3. Mechanistic model

**4 GI segments** duodenum → colon (`configs/gi_physiology.json:4`) — pH 6.0→6.8, TEER 60→120 Ω·cm², bile 0.5–8 mM, UWL 30 µm, 5N+13 ODE states.

```
INPUT (SMILES) → DISSOLUTION (Noyes-Whitney) → MUCUS (Amsden) → ENTEROCYTE (Fg·P-gp/BCRP) → BLOOD (Fh·AUC)
                   src/navigator/gi/dissolution.py   src/navigator/gi/mucus.py   src/navigator/engine/odes.py
```

**Net epithelial flux** (`poster.html:136`):

```
J_total = J_trans + J_para + J_endo + J_OCT + J_PEPT1 + J_ENT + J_SGLT + J_OATP
J_net   = max(0, J_total − J_pgp)      · F_oral = A_abs / dose
```

| Pathway | Model | File |
|---------|-------|------|
| **Transcellular** | Overton + Goldman (−60 mV) + UWL | `src/navigator/gi/permeation/passive.py`, `src/navigator/chemistry/permeability.py` |
| **Paracellular** | Renkin 2 pores (0.6 / 5 nm) + √Goldman | `src/navigator/gi/wall.py` |
| **Transporters (candidate)** | OCT, PEPT1, ENT, SGLT, OATP (Michaelis–Menten) | `src/navigator/gi/transporters.py`, `src/navigator/gi/permeation/carriers.py` |
| **Efflux** | P-gp / BCRP | `src/navigator/gi/permeation/efflux.py` |
| **Endocytosis** | Transcytosis fraction | `src/navigator/gi/permeation/endocytosis.py` |
| **Enhancers** | Hill EC₅₀/Dmax/n, TEER modulation (Bliss) | `src/navigator/gi/permeation/enhancers.py` |
| **First-pass** | `Fg = 1/(1+fu·CLint·Tgut/Vgut)`, `Fh = Qh/(Qh+fu·CLint)` | `src/navigator/systemic/first_pass.py` |

- **Solver:** SciPy LSODA (`src/navigator/engine/solver.py`) — `rtol 1e-6, atol 1e-9, t_end 96h` (`configs/numerics.json`)
- **Mass balance:** `|dose − (GI_rest + A_central + A_elim + A_fecal)| < 1e-10` (`src/navigator/engine/runner.py:42`)
- **Schemas strict** `extra=forbid` (`src/navigator/core/schemas.py:10`) — every JSON validated, `run_id` deterministic hash (`src/navigator/core/provenance.py`)

## 4. Tech stack

- **Core:** Python 3.11+, NumPy 2.x, SciPy 1.15, Pydantic 2.10 (strict schemas)
- **API:** FastAPI 0.136 + Uvicorn (single-port API + dashboard SPA) — `src/navigator/interfaces/api/app.py:28`
- **Optimization:** NSGA-II MOO (`src/navigator/optimization/nsga2.py`)
- **No external PBPK black box** — fully mechanistic, auditable `trace[]` (`src/navigator/core/trace.py`)

## 5. Project structure

```
.
├── configs/                      # Validated physiology & numerics
│   ├── gi_physiology.json        # 4 segments, pH/transit/volume/area/TEER/mucus/bile
│   ├── permeation.json           # pore sizes, QSPR, UWL, Goldman, vmax/p scales
│   ├── transporters.json         # Vmax/Km per transporter
│   ├── formulation_effects.json  # Hill EC50/Dmax/n, TEER, chitosan
│   ├── first_pass.json / systemic.json
│   ├── numerics.json             # LSODA, mass_tol, n_eval
│   ├── calibration.json          # 7 calibrated params (RMSE 18.2→13.9 pp)
│   └── candidates/               # gitignored — NSGA-II trials
├── schemas/                      # JSON Schema v1 (compound, formulation, excipient, result)
├── data/
│   ├── compounds/                # 14 BCS III drugs (tobramycin, metformin, atenolol...)
│   ├── compounds_raw/            # raw vendor sheets
│   ├── formulations/             # native.json + optimized/ (A/B/full)
│   ├── excipients/               # catalog.json + detailed.json + NAV-880
│   ├── optimization_space/       # tobramycin_A/B/full/search.json
│   ├── validation/f_lit.json     # literature F_oral benchmark
│   └── calibration_space/
├── src/navigator/
│   ├── chemistry/                # ionization, lipophilicity, solubility, descriptors, permeability
│   ├── gi/                       # dissolution, mucus, wall, precipitation, transporters, permeation/*
│   ├── formulation/              # hip, sedds, composition, release
│   ├── systemic/                 # first_pass, pk, compartment, mass_balance
│   ├── engine/                   # odes, fluxes, solver, runner, state, report
│   ├── optimization/             # nsga2, moo
│   ├── validation/               # metrics (RMSE, R²)
│   ├── application/              # predict, batch, calibrate, validate, formulation_opt
│   ├── core/                     # schemas, settings, loader, trace, provenance, units, errors
│   └── interfaces/               # api (FastAPI), cli, worker
├── tests/
│   ├── unit/test_core.py         # Bohley, Goldman, Renkin, ionization
│   ├── integration/test_engine.py
│   ├── contract/ / golden/ / property/
├── poster.html                   # A4 scientific poster
├── pyproject.toml
└── README.md
```

## 6. Quick start

### Install

```bash
git clone https://github.com/btkhaled/bcs-iii-navigator.git
cd bcs-iii-navigator
python -m venv .venv && source .venv/bin/activate
pip install -e .              # or pip install -e ".[dev]"
# pip install fastapi uvicorn scipy pydantic  # minimal
```

### CLI — predict

```bash
python -m navigator.interfaces.cli.main data/compounds/tobramycin.json data/formulations/native.json configs/
# → {"F_oral": 0.0138, "F_oral_pct": 1.38, "Fg": ..., "Fh": ..., "dominant": "para", ...}
# see src/navigator/interfaces/cli/main.py:12
```

### Python API

```python
from navigator.application.predict import predict
from navigator.engine.report import summarize
import json

r = predict("data/compounds/tobramycin.json", "data/formulations/native.json", "configs/")
print(json.dumps(summarize(r), indent=2))
# {"F_oral": 0.0138, "F_oral_pct": 1.38, "AUC": ..., "Cmax": ..., "breakdown": {...}}

# Direct runner with full trace
from navigator.core.schemas import Compound, Formulation
from navigator.core.settings import Settings
from navigator.engine.runner import run  # src/navigator/engine/runner.py:18

s = Settings("configs/")
c = Compound(**json.load(open("data/compounds/tobramycin.json")))
f = Formulation(**json.load(open("data/formulations/optimized/tobramycin_A_best.json")))
r = run(c, f, s)
print(r.f_oral, r.dominant, r.mass_balance_mg, r.trace[:2])
```

### API + Dashboard (single port)

```bash
uvicorn navigator.interfaces.api.app:app --reload --port 8000
# http://localhost:8000/            → dashboard SPA (src/navigator/interfaces/api/static/index.html)
# http://localhost:8000/api/docs    → OpenAPI (Swagger)
# http://localhost:8000/api/health  → {"ok": true, "version": "0.1.0"}
```

```bash
# Predict via JSON body
curl -X POST http://localhost:8000/api/predict \
  -H "Content-Type: application/json" \
  -d '{"compound": {"id":"tobramycin","mw":467.52,"cs_mg_l":100000,"logp":-5.8,"pkas":[{"type":"base","pka":8.6}]}, "formulation": {"id":"native","dose_mg":500}}' | jq

# By ID (uses data/compounds/*.json + data/formulations/*.json)
curl -X POST http://localhost:8000/api/predict/by-id \
  -H "Content-Type: application/json" \
  -d '{"compound_id":"tobramycin","formulation_id":"native","dose_mg":500}' | jq

# Batch
curl -X POST http://localhost:8000/api/batch -H "Content-Type: application/json" -d @- <<'JSON' | jq
{"pairs": [{"compound": {"id":"tobramycin","mw":467.52,"cs_mg_l":100000}, "formulation": {"dose_mg":500}}]}
JSON

# HIP calculator (SDS mass for given m_hip)
curl -X POST http://localhost:8000/api/hip-calc -H "Content-Type: application/json" \
  -d '{"dose_mg":500,"m_hip":32}' | jq

# Validation benchmark
curl http://localhost:8000/api/validate | jq
# {"rmse_pp": 13.9, "r2": 0.65, "details": [...]}

# Precomputed optimization fronts
curl "http://localhost:8000/api/optimization/front?space=A" | jq
```

### Auto-formulate — NSGA-II (async job)

```bash
curl -X POST http://localhost:8000/api/auto-formulate \
  -H "Content-Type: application/json" \
  -d '{"compound_id":"tobramycin","space":"A","pop":24,"gen":15,"seed":0}' | jq
# → {"job_id":"a1b2c3d4","status":"queued"}

curl http://localhost:8000/api/auto-formulate/a1b2c3d4 | jq
# → {"status":"running","progress":{"gen":8,"total":15,"best_F":18.2,"elapsed":42}}
# → {"status":"done","result":{"front":[{"F_pct":18.62,"dominant":"trans",...}]}}
```

Direct Python:

```python
from navigator.application.formulation_opt import run_formulation_optim

res = run_formulation_optim(
    config_dir="configs/",
    compound_path="data/compounds/tobramycin.json",
    space_path="data/optimization_space/tobramycin_A.json",
    pop=24, gen=15, seed=0, workers=4
)
print(res["front"][0])  # best formulation
```

## 7. Configuration

Every `configs/*.json` loaded via `Settings` (`src/navigator/core/settings.py:20`) and validated against Pydantic schemas (`src/navigator/core/schemas.py`):

| File | Key fields |
|------|------------|
| `gi_physiology.json` | `segments[].ph, transit_h, volume_l, area_cm2, radius_cm, teer0, mucus_thickness_um, bile_mm, expr{PGP,BCRP,OCT3,PEPT1...}` |
| `permeation.json` | `pore_nm, pore_large_nm, frac_large_pore, p_base_cm_s, qspr_a/b, uwl_thickness_um, membrane_potential_mv, vmax_scale, p_scale, t_scale` |
| `transporters.json` | `{transporter: {vmax_pmol_cm2_min, km_um}}` |
| `formulation_effects.json` | `enhancers{ec50_mm,dmax,n}, snac_factor, chitosan_teer_factor, teer_max_drop` |
| `calibration.json` | `p_scale, t_scale, vmax_scale, frac_large_pore, ref_OCT3/ENT, fecal_bound_colon` |
| `numerics.json` | `method: LSODA, rtol, atol, mass_tol_mg, n_eval, t_end_h` |
| `first_pass.json` | `qh_l_h, vgut_l, tgut_h` |
| `systemic.json` | `vd_l, cl_l_h` |

Compound schema (`schemas/compound.v1.schema.json`, `src/navigator/core/schemas.py:19`):
```json
{ "id": "tobramycin", "mw": 467.52, "logp": -5.8, "tpsa": 268, "cs_mg_l": 100000, "pkas": [{"type":"base","pka":8.6}], "fu_plasma": 0.9 }
```

Formulation schema (`src/navigator/core/schemas.py:45`): `dose_mg, m_hip, c10_mm, labrasol_mm, tpgs_mm, pip_mm, nacdc_mm, chitosan, enteric, lipid_frac, ...`

## 8. Validation — preliminary qualification

External benchmark: **14 BCS III drugs** (`data/validation/f_lit.json:2`) — 10 train / 4 holdout, 7 calibrated parameters.

| Metric | Before calibration | After constrained calibration |
|--------|-------------------|-------------------------------|
| **RMSE** | 18.2 pp | **13.9 pp** |
| **R²** | 0.36 | **0.65** |

Mechanistic constraints reduce the risk of compensating structural misspecification through parameter fitting.

```bash
curl http://localhost:8000/api/validate | jq '.rmse_pp, .r2'
python -m navigator.application.validate  # src/navigator/application/validate.py
```

## 9. Optimization

MOO NSGA-II (`src/navigator/optimization/nsga2.py:optimize`, `src/navigator/optimization/moo.py`):

- **Space A:** excipients with regulatory precedents (`data/optimization_space/tobramycin_A.json`)
- **Space B / full:** de-novo theoretical enhancer `denovo_p_trans / denovo_p_para / denovo_mucus` (`data/optimization_space/tobramycin_full.json`)

```bash
python -m navigator.interfaces.cli.calibrate --help  # src/navigator/interfaces/cli/calibrate.py
```

## 10. Reproducibility & Audit

- **Deterministic `run_id`:** SHA hash of compound + formulation + physiology (`src/navigator/core/provenance.py`, `src/navigator/engine/runner.py:57`)
- **Auditable `trace[]`:** every calculation step logged (`src/navigator/core/trace.py`, `src/navigator/engine/runner.py:61`)
- **Versioned inputs:** `compound.json + formulation.json + configs/ → result.json` (mass_balance, breakdown, warnings)
- **Schemas:** `schemas/result.v1.schema.json` — `F_oral, Fg, Fh, AUC, Cmax, breakdown{trans,para,endo,oct...pgp}`

```python
from navigator.engine.report import summarize  # src/navigator/engine/report.py:8
summarize(r)
# {"F_oral":0.1862, "F_oral_pct":18.62, "Fg":..., "Fh":..., "AUC":..., "Cmax":..., "dominant":"trans", ...}
```

## 11. Evidence levels

| Level | Examples |
|-------|----------|
| **Literature-derived / Measured** | MW 467.52, pKₐ series, solubility 94 mg·mL⁻¹ @25°C |
| **Calculated** | XlogP −6.2, charge +4.43 @pH6.5, Bohley radius, SDS masses (617 mg @2:1) |
| **Assumed** | Hill EC₅₀/Dmax/n, HIP m=32, hypothetical Vmax/Km |
| **Predicted** | F 1.38%→18.62%, AUC, Cmax, breakdown — depends on assumptions |

## 12. Tests

```bash
pip install -e ".[dev]"
pytest tests/unit/test_core.py -v          # Bohley, Goldman, Renkin, ionization — src: tests/unit/test_core.py:12
pytest tests/integration/test_engine.py -v
pytest
```

## 13. Roadmap

- [ ] Quantify parameter + structural uncertainty (value-of-information)
- [ ] Identify highest-value experiments
- [ ] CI + golden/property tests (`tests/golden/`, `tests/property/`)
- [ ] Docker image + GHCR
- [ ] Calibration under stricter mechanistic constraints

## 14. Citation

```bibtex
@software{navigator2026,
  title   = {BCS III NAVIGATOR: A mechanistic PBPK framework for poorly permeable drugs},
  author  = {Ben Taieb, Khaled},
  year    = {2026},
  version = {0.1.0},
  url     = {https://github.com/btkhaled/bcs-iii-navigator}
}
```

## 15. Disclaimer

Mechanistically structured, selectively calibrated under mechanistic constraints. Predictions depend on assumptions above — **for research/exploration, not clinical decision-making**. No claim of BCS classification beyond the tobramycin case study.

## 16. License

Proprietary — All rights reserved. Private repository. Contact the author for access, collaboration or licensing: **bt.khaled@gmail.com**

---

*Résumé FR — Framework PBPK mécanistique pour prédire la biodisponibilité orale (F_oral) des molécules BCS III peu perméables. Séparation explicite littérature / calculé / supposé / prédit, traçabilité complète et bilan de masse < 1e-10. Cas d’étude tobramycine : F natif 1,38% → 18,62% optimisé (excipients à précédents réglementaires).*
