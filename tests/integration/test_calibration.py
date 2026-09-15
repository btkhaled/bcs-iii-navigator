"""Test calibration: un bon-RMSE à mécanisme faux est REJETÉ."""

import sys
sys.path.insert(0, "/Users/kalo/NAVIGATOR/src")

from navigator.inference.calibrate import Calibrator

BASE = "/Users/kalo/NAVIGATOR"
ARGS = (f"{BASE}/configs", f"{BASE}/data/compounds",
        f"{BASE}/data/formulations/native.json",
        f"{BASE}/data/calibration_space/search_space.json",
        f"{BASE}/data/calibration_space/mechanism_constraints.json",
        f"{BASE}/data/validation/f_lit.json")


def test_current_feasible():
    with Calibrator(*ARGS, workers=2) as cal:
        (objs, viol, aux), = cal.evaluate_many([cal.current_x])
    assert viol == 0, aux


def test_bad_mechanism_rejected():
    with Calibrator(*ARGS, workers=2) as cal:
        x = list(cal.current_x)
        # frac_large_pore au max (1e-3): tobra sort de sa bande mecanistique
        x[5] = -3.0
        (objs, viol, aux), = cal.evaluate_many([x])
    assert viol > 0, aux


if __name__ == "__main__":
    test_current_feasible()
    test_bad_mechanism_rejected()
    print("calibration OK")
