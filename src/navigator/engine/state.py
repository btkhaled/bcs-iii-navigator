"""État ODE: layout 5N+12."""

from __future__ import annotations

from dataclasses import dataclass


N_CUM = 9  # trans, para, endo, oct, pept1, ent, sglt, oatp, pgp
CUM_NAMES = ["trans", "para", "endo", "oct", "pept1", "ent", "sglt", "oatp", "pgp"]


@dataclass
class Layout:
    n: int

    @property
    def i_diss(self): return (0, self.n)
    @property
    def i_solid(self): return (self.n, 2 * self.n)
    @property
    def i_mucus(self): return (2 * self.n, 3 * self.n)
    @property
    def i_wall(self): return (3 * self.n, 4 * self.n)
    @property
    def i_teer(self): return (4 * self.n, 5 * self.n)
    @property
    def i_central(self): return 5 * self.n
    @property
    def i_elim(self): return 5 * self.n + 1
    @property
    def i_fecal(self): return 5 * self.n + 2
    @property
    def i_abs(self): return 5 * self.n + 3
    @property
    def i_cum(self): return (5 * self.n + 4, 5 * self.n + 4 + N_CUM)

    @property
    def size(self): return 5 * self.n + 4 + N_CUM


def y0_dose(layout: Layout, dose_mg: float, teer0_list):
    import numpy as np

    y = np.zeros(layout.size)
    # dose solide dans premier segment
    a, b = layout.i_solid
    y[a] = dose_mg
    a, b = layout.i_teer
    for i, t0 in enumerate(teer0_list):
        y[a + i] = t0
    return y
