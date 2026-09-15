"""Efflux P-gp / BCRP depuis l'entérocyte vers la lumière."""

from __future__ import annotations

from .carriers import transporter_flux_mg_h


def pgp_efflux_mg_h(compound, z, c_cell_mg_l, area_cm2, vmax, km, expr=1.0, vmax_scale=1.0, inhibitor=0.0):
    # efflux = MM sur concentration cellulaire
    from ...core.units import pmol_cm2_min_to_mg_h, mg_l_to_um
    from .carriers import SUBSTRATE_FNS, mm_rate

    if not SUBSTRATE_FNS["PGP"](compound, z):
        return 0.0
    v = mm_rate(vmax * expr * vmax_scale, km, mg_l_to_um(max(c_cell_mg_l, 0.0), compound.mw))
    j = pmol_cm2_min_to_mg_h(v, area_cm2, compound.mw)
    if inhibitor > 0:
        j *= (1.0 - 0.9 * min(1.0, inhibitor))
    return max(0.0, j)


def bcrp_efflux_mg_h(compound, z, c_cell_mg_l, area_cm2, vmax, km, expr=1.0, vmax_scale=1.0):
    from ...core.units import pmol_cm2_min_to_mg_h, mg_l_to_um
    from .carriers import SUBSTRATE_FNS, mm_rate

    if not SUBSTRATE_FNS["BCRP"](compound, z):
        return 0.0
    v = mm_rate(vmax * expr * vmax_scale, km, mg_l_to_um(max(c_cell_mg_l, 0.0), compound.mw))
    return max(0.0, pmol_cm2_min_to_mg_h(v, area_cm2, compound.mw))
