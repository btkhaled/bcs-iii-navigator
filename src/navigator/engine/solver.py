"""Solveur LSODA/Radau/BDF."""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from ..core.errors import SolverError


def solve(dydt, y0, t_end_h: float, method: str = "LSODA", rtol: float = 1e-6, atol: float = 1e-9, n_eval: int = 400):
    m = {"LSODA": "LSODA", "RADAU": "Radau", "BDF": "BDF"}.get(method.upper(), "LSODA")
    t_eval = np.linspace(0, t_end_h, max(20, n_eval))
    try:
        sol = solve_ivp(dydt, (0, t_end_h), y0, method=m, t_eval=t_eval, rtol=rtol, atol=atol)
    except Exception as e:
        raise SolverError(str(e)) from e
    if not sol.success:
        raise SolverError(sol.message)
    # positivité douce (anti-bruit numérique)
    y = np.maximum(sol.y, 0.0)
    return sol.t, y
