"""Tests unitaires socle."""

import sys
sys.path.insert(0, "/Users/kalo/NAVIGATOR/src")

from navigator.chemistry.descriptors import bohley_radius_nm
from navigator.chemistry.ionization import ionization_at_ph
from navigator.core.schemas import Compound
from navigator.gi.permeation.passive import goldman_factor, renkin_factor


def test_bohley():
    r = bohley_radius_nm(467.52)
    assert 0.5 < r < 1.0, r


def test_goldman_cation():
    m = goldman_factor(5.0, -60.0)
    assert m > 5.0, m  # cations fortement favorisés


def test_goldman_neutral():
    assert abs(goldman_factor(0.0) - 1.0) < 1e-9


def test_renkin_excluded():
    assert renkin_factor(1.0, 0.6) == 0.0


def test_ionization_tobra():
    c = Compound(id="t", mw=467.0, logp=-5.0, cs_mg_l=1000.0,
                 pkas=[{"type": "base", "pka": 8.0}])
    z, f_u, _ = ionization_at_ph(c, 6.5)
    assert z > 0.9, z
    assert f_u < 0.05, f_u


if __name__ == "__main__":
    test_bohley(); test_goldman_cation(); test_goldman_neutral()
    test_renkin_excluded(); test_ionization_tobra()
    print("unit OK")
