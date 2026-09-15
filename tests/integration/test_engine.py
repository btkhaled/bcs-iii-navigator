"""Tests intégration: bilan masse exact + directions (formulé > natif)."""

import sys
sys.path.insert(0, "/Users/kalo/NAVIGATOR/src")

from navigator.core.loader import load_compound, load_formulation
from navigator.core.settings import Settings
from navigator.engine.runner import run

CFG = "/Users/kalo/NAVIGATOR/configs"
DATA = "/Users/kalo/NAVIGATOR/data"


def test_mass_balance():
    s = Settings(CFG)
    c = load_compound(f"{DATA}/compounds/tobramycin.json")
    f = load_formulation(f"{DATA}/formulations/native.json")
    r = run(c, f, s)
    assert r.mass_balance_mg < 1.0, r.mass_balance_mg
    assert 0.0 <= r.f_oral <= 1.0


def test_tobra_native_below_2pct():
    s = Settings(CFG)
    c = load_compound(f"{DATA}/compounds/tobramycin.json")
    f = load_formulation(f"{DATA}/formulations/native.json")
    r = run(c, f, s)
    assert r.f_oral < 0.02, r.f_oral
    assert r.dominant == "para", r.dominant


def test_formulated_beats_native():
    s = Settings(CFG)
    c = load_compound(f"{DATA}/compounds/tobramycin.json")
    nat = load_formulation(f"{DATA}/formulations/native.json")
    opt = load_formulation(f"{DATA}/formulations/tobramycin_ode_01.json")
    r_nat = run(c, nat, s)
    r_opt = run(c, opt, s)
    assert r_opt.f_oral > r_nat.f_oral, (r_opt.f_oral, r_nat.f_oral)


if __name__ == "__main__":
    test_mass_balance(); test_tobra_native_below_2pct(); test_formulated_beats_native()
    print("integration OK")
