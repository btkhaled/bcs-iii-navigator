"""Runner: compound + formulation -> SimResult."""

from __future__ import annotations

import numpy as np

from ..chemistry.ionization import ionization_at_ph
from ..core.provenance import run_id
from ..core.schemas import Compound, Formulation, SimResult, ResultBreakdown
from ..core.settings import Settings
from ..core.trace import Trace
from ..systemic.pk import pk_metrics
from .odes import make_dydt
from .solver import solve
from .state import Layout, y0_dose


def run(compound: Compound, formulation: Formulation, settings: Settings) -> SimResult:
    trace = Trace()
    layout = Layout(len(settings.physiology.segments))
    teer0 = [s.teer0 for s in settings.physiology.segments]
    y0 = y0_dose(layout, formulation.dose_mg, teer0)

    dydt, info = make_dydt(compound, formulation, settings)
    t, y = solve(dydt, y0, settings.numerics.t_end_h, settings.numerics.method,
                 settings.numerics.rtol, settings.numerics.atol, settings.numerics.n_eval)

    n = layout.n
    A_central = y[layout.i_central, :]
    auc, cmax, _ = pk_metrics(t, A_central, settings.systemic.vd_l)

    A_diss_f = float(np.sum(y[0:n, -1]))
    A_solid_f = float(np.sum(y[n:2 * n, -1]))
    A_muc_f = float(np.sum(y[2 * n:3 * n, -1]))
    A_wall_f = float(np.sum(y[3 * n:4 * n, -1]))
    A_fecal = float(y[layout.i_fecal, -1])
    A_central_f = float(y[layout.i_central, -1])
    A_elim = float(y[layout.i_elim, -1])
    A_abs = float(y[layout.i_abs, -1])
    gi_rest = A_diss_f + A_solid_f + A_muc_f + A_wall_f
    f_oral = max(0.0, min(1.0, A_abs / max(formulation.dose_mg, 1e-12)))
    mass_bal = abs(formulation.dose_mg - (gi_rest + A_central_f + A_elim + A_fecal))

    base = layout.i_cum[0]
    cum = [float(y[base + i, -1]) for i in range(9)]
    tot = sum(cum[0:8])
    if tot <= 0:
        fracs = [0.0] * 8
    else:
        fracs = [c / tot for c in cum[:8]]
    names = ["trans", "para", "endo", "oct", "pept1", "ent", "sglt", "oatp"]
    dom = names[int(np.argmax(fracs))] if tot > 0 else "none"
    bd = ResultBreakdown(trans=fracs[0], para=fracs[1], endo=fracs[2], oct=fracs[3],
                         pept1=fracs[4], ent=fracs[5], sglt=fracs[6], oatp=fracs[7], pgp=float(cum[8]))

    z, f_u, _ = ionization_at_ph(compound, settings.physiology.segments[1].ph if n > 1 else 6.5)
    rid = run_id(compound.model_dump(), formulation.model_dump(), settings.physiology.model_dump())
    warnings = []
    if mass_bal > max(settings.numerics.mass_tol_mg * max(formulation.dose_mg, 1.0), 0.5):
        warnings.append(f"mass_balance élevé: {mass_bal:.4f} mg")
    trace.log("F_oral", "A_abs/dose", {"A_abs": A_abs, "dose": formulation.dose_mg}, f_oral, "-")
    trace.log("Fg", "1/(1+fu*CLint*Tgut/Vgut)", {"fg": info["fg"]}, info["fg"], "-")
    trace.log("Fh", "Qh/(Qh+fu*CLint)", {"fh": info["fh"]}, info["fh"], "-")

    return SimResult(f_oral=f_oral, f_g=info["fg"], f_h=info["fh"], auc_mg_h_l=auc,
                     cmax_mg_l=cmax, breakdown=bd, dominant=dom, z=z,
                     f_unionized=f_u, mass_balance_mg=mass_bal, run_id=rid,
                     trace=trace.to_list(), warnings=warnings)
