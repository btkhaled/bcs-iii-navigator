"""Système ODE: dydt avec conservation masse par construction + limiteurs."""

from __future__ import annotations

import numpy as np

from ..core.schemas import Compound, Formulation
from ..core.settings import Settings
from ..gi.wall import fg_gut
from ..systemic.first_pass import fh_liver
from .fluxes import segment_fluxes
from .state import Layout

K_CAP_H = 200.0  # cap flux sortant: J <= A * K_CAP (anti-overshoot)


def _cap(j: float, a: float) -> float:
    if j <= 0 or a <= 0:
        return max(j, 0.0) if a > 0 else 0.0
    return min(j, a * K_CAP_H)


def make_dydt(compound: Compound, formulation: Formulation, settings: Settings):
    layout = Layout(len(settings.physiology.segments))
    fg = fg_gut(compound.fu_plasma, compound.clint_gut_ul_min_mg,
                settings.first_pass.tgut_h, settings.first_pass.vgut_l)
    fh = fh_liver(compound.fu_plasma, compound.clint_liver_ul_min_mg, settings.first_pass.qh_l_h)
    k_el = settings.systemic.cl_l_h / max(settings.systemic.vd_l, 1e-9)
    tau_on = settings.formulation_fx.tau_on_h
    tau_rec = settings.formulation_fx.tau_rec_h

    def dydt(t, y):
        n = layout.n
        dyd = np.zeros_like(y)
        J_to_blood_all = 0.0

        for j in range(n):
            A_diss = max(float(y[j]), 0.0)
            A_solid = max(float(y[n + j]), 0.0)
            A_muc = max(float(y[2 * n + j]), 0.0)
            A_wall = max(float(y[3 * n + j]), 0.0)
            teer = max(float(y[4 * n + j]), 1.0)
            fl, _ = segment_fluxes(j, t, A_diss, A_solid, A_muc, A_wall, teer,
                                   compound, formulation, settings, fg)
            # limiteurs par compartiment source
            J_trans = _cap(fl["J_trans"], A_diss)
            J_para = _cap(fl["J_para"], A_diss)
            J_endo_up = _cap(fl["J_endo"], A_diss)
            # transcytose: 1% rejoint le sang, 99% dégradé lysosomal
            f_tc = settings.permeation.endo_transcytosis
            J_endo = J_endo_up * f_tc
            J_endo_deg = J_endo_up * (1 - f_tc)
            J_oct = _cap(max(fl.get("J_oct", 0.0), 0.0), A_diss)
            J_pept1 = _cap(max(fl.get("J_pept1", 0.0), 0.0), A_diss)
            J_ent = fl.get("J_ent", 0.0)  # bidirectionnel, petit
            J_ent_pos = _cap(max(J_ent, 0.0), A_diss)
            J_sglt = _cap(max(fl.get("J_sglt", 0.0), 0.0), A_diss)
            J_oatp = _cap(max(fl.get("J_oatp", 0.0), 0.0), A_diss)
            J_lm = _cap(fl["J_lm"], A_diss)
            J_ml = _cap(fl["J_ml"], A_muc)
            J_muc_clear = _cap(fl["J_muc_clear"], A_muc)
            J_pgp = _cap(fl["J_pgp"], A_wall)
            J_bcrp = _cap(fl["J_bcrp"], A_wall)
            J_tr_diss = _cap(fl["J_tr_diss"], A_diss)
            J_tr_solid = _cap(fl["J_tr_solid"], A_solid)
            J_diss = _cap(fl["J_diss"], A_solid)
            J_prec = _cap(fl["J_prec"], A_diss)
            J_to_blood = _cap(fl["J_to_blood"], A_wall)
            J_wall_met = _cap(fl["J_wall_met"], A_wall)
            # renormalise si somme sorties lumen > disponible
            out_lumen = (J_trans + J_para + J_endo_up + J_oct + J_pept1 + J_ent_pos + J_sglt + J_oatp
                         + J_lm + J_tr_diss + J_prec)
            avail_lumen = A_diss * K_CAP_H + J_diss + J_ml + J_pgp + J_bcrp
            if out_lumen > avail_lumen and out_lumen > 0:
                s = avail_lumen / out_lumen
                J_trans *= s; J_para *= s; J_endo *= s; J_oct *= s; J_pept1 *= s
                J_ent_pos *= s; J_sglt *= s; J_oatp *= s; J_lm *= s; J_tr_diss *= s; J_prec *= s

            dyd[j] += (J_diss - J_prec - J_trans - J_para - J_endo_up - J_oct - J_pept1
                       - J_ent_pos - J_sglt - J_oatp - J_lm + J_ml + J_pgp + J_bcrp - J_tr_diss)
            dyd[n + j] += -J_diss + J_prec - J_tr_solid
            dyd[2 * n + j] += J_lm - J_ml - J_muc_clear
            J_in_wall = J_trans + J_para + J_endo + J_oct + J_pept1 + J_ent_pos + J_sglt + J_oatp
            dyd[3 * n + j] += J_in_wall - J_pgp - J_bcrp - J_to_blood - J_wall_met
            tgt = fl["teer_tgt"]
            tau = tau_on if tgt < teer else tau_rec
            dyd[4 * n + j] += (tgt - teer) / max(tau, 1e-9)

            if j + 1 < n:
                dyd[j + 1] += J_tr_diss
                dyd[n + j + 1] += J_tr_solid
            else:
                dyd[layout.i_fecal] += J_tr_diss + J_tr_solid
            dyd[layout.i_fecal] += J_muc_clear
            dyd[layout.i_elim] += J_wall_met + J_endo_deg

            base = layout.i_cum[0]
            dyd[base + 0] += J_trans
            dyd[base + 1] += J_para
            dyd[base + 2] += J_endo
            dyd[base + 3] += J_oct
            dyd[base + 4] += J_pept1
            dyd[base + 5] += J_ent_pos
            dyd[base + 6] += J_sglt
            dyd[base + 7] += J_oatp
            dyd[base + 8] += J_pgp

            J_to_blood_all += J_to_blood

        dyd[layout.i_central] += J_to_blood_all * fh - k_el * max(float(y[layout.i_central]), 0.0)
        dyd[layout.i_elim] += J_to_blood_all * (1 - fh) + k_el * max(float(y[layout.i_central]), 0.0)
        dyd[layout.i_abs] += J_to_blood_all * fh
        return dyd

    return dydt, {"fg": fg, "fh": fh, "k_el": k_el, "layout": layout}
