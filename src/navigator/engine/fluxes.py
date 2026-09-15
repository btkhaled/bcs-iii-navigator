"""Flux mécanistiques par segment. Aucune constante en dur."""

from __future__ import annotations

import math
from typing import Dict, Tuple

from ..chemistry.descriptors import bohley_radius_nm
from ..chemistry.ionization import ionization_at_ph
from ..chemistry.permeability import pint0_cm_s
from ..chemistry.solubility import cs_at_ph
from ..core.schemas import Compound, Formulation
from ..core.settings import Settings
from ..core.units import cm_s_mg_l_cm2_to_mg_h
from ..formulation.composition import evaluate as eval_form
from ..formulation.hip import m_hip_dynamic
from ..formulation.release import effective_release
from ..formulation.sedds import sedds_solubilization
from ..gi.dissolution import noyes_whitney_dMdt_mg_h
from ..gi.mucus import f_mucus
from ..gi.permeation.carriers import transporter_flux_mg_h
from ..gi.permeation.efflux import bcrp_efflux_mg_h, pgp_efflux_mg_h
from ..gi.permeation.endocytosis import endo_flux_mg_h
from ..gi.permeation.enhancers import bliss_total, hill_drop, teer_multiplier, teer_target
from ..gi.permeation.passive import goldman_factor, renkin_factor, series_permeability, two_pore_factor, uwl_permeability_cm_s
from ..gi.precipitation import precipitation_rate_mg_h
from ..gi.transporters import vmax_km
from ..gi.wall import fg_gut


K_MUC_H = 2.0
K_MUC_CLEAR_H = 0.3
K_WALL_H = 3.0


def segment_fluxes(
    j: int,
    t_h: float,
    A_diss: float,
    A_solid: float,
    A_mucus: float,
    A_wall: float,
    teer: float,
    compound: Compound,
    formulation: Formulation,
    settings: Settings,
    fg: float,
) -> Tuple[Dict[str, float], Dict[str, float]]:
    seg = settings.physiology.segments[j]
    fx_cfg = settings.formulation_fx
    perm = settings.permeation

    V = seg.volume_l
    V_muc = V * 0.1
    V_wall = V * 0.1
    C_diss = max(A_diss, 0.0) / V
    C_muc = max(A_mucus, 0.0) / max(V_muc, 1e-9)
    C_wall = max(A_wall, 0.0) / max(V_wall, 1e-9)

    z, f_union, _ = ionization_at_ph(compound, seg.ph)
    cs = cs_at_ph(compound, seg.ph, bile_mm=seg.bile_mm + formulation.bile_mm)
    cs *= sedds_solubilization(formulation.lipid_frac, seg.bile_mm)

    f = eval_form(formulation)
    m_hip = m_hip_dynamic(formulation.m_hip, seg.ph, C_diss, t_h, formulation.hip_counterion)
    f_hip = min(0.95, (m_hip - 1.0) / (m_hip + 4.0)) if m_hip > 1 else 0.0
    f_trans_frac = f_union + (1 - f_union) * f_hip * 0.5
    C_free = C_diss * formulation.f_free * (1.0 - seg.fecal_bound)

    rel = effective_release(t_h, formulation.k_release_h, seg.ph, formulation.enteric, formulation.enteric_ph_threshold)

    radius_nm = bohley_radius_nm(compound.mw)
    # particule formulation pour dissolution (pas rayon moléculaire)
    diss_cfg = settings.dissolution
    part_um = diss_cfg.get("particle_radius_um", 25.0)
    if formulation.particle_nm > 0:
        part_nm = formulation.particle_nm
    else:
        part_nm = part_um * 1000.0
    h_um = diss_cfg.get("diffusion_layer_um", 30.0)
    rho = diss_cfg.get("density_mg_cm3", 1300.0)
    # dissolution sur solide disponible
    J_diss_raw = noyes_whitney_dMdt_mg_h(A_solid, C_diss, cs, part_nm, h_um, rho, compound.mw) * max(rel, 0.05)
    # cap anti-raideur: pas plus que le solide disponible / dt_min
    J_diss = min(J_diss_raw, max(A_solid, 0.0) * 50.0)
    J_prec = precipitation_rate_mg_h(C_diss, cs, V)

    # mucus corrigé de-novo (voie B)
    Fmuc_raw = f_mucus(formulation.particle_nm, seg.mucus_thickness_um, formulation.zeta_mv, fx_cfg.d_crit_nm)
    Fmuc = min(1.0, Fmuc_raw * formulation.denovo_mucus)

    # --- transcellulaire (échelle propre t_scale, distincte du para) ---
    p_int0 = pint0_cm_s(compound.logp, perm.qspr_a, perm.qspr_b) * perm.t_scale
    pip_factor = 1.0 + fx_cfg.pip_factor_per_mm * formulation.pip_mm
    p_mem = p_int0 * m_hip * f.lipid_factor * f.snac_factor * pip_factor * max(f_trans_frac, 1e-6)
    p_mem *= goldman_factor(z, perm.membrane_potential_mv) if perm.goldman else 1.0
    p_mem *= formulation.denovo_p_trans
    p_uwl = uwl_permeability_cm_s(compound.mw, perm.uwl_thickness_um)
    p_trans = series_permeability(p_mem, p_uwl)
    J_trans = cm_s_mg_l_cm2_to_mg_h(p_trans, C_free, seg.area_cm2) * Fmuc * rel

    # --- paracellulaire (2 pores: 0.6nm abondant + 5nm rare) ---
    phi = two_pore_factor(radius_nm, perm.pore_nm, perm.pore_large_nm, perm.frac_large_pore)
    m_teer = teer_multiplier(seg.teer0, teer, fx_cfg.teer_max_drop)
    m_chi = 1.0 + formulation.chitosan * fx_cfg.chitosan_teer_factor
    p_para = perm.p_base_cm_s * perm.p_scale * phi * m_teer * m_chi
    p_para *= formulation.denovo_p_para
    # sélectivité cationique jonctions: racine du facteur Goldman
    g = goldman_factor(z, perm.membrane_potential_mv) if perm.goldman else 1.0
    p_para *= math.sqrt(max(g, 1e-6))
    J_para = cm_s_mg_l_cm2_to_mg_h(p_para, C_free, seg.area_cm2) * Fmuc

    # --- endo ---
    J_endo = endo_flux_mg_h(C_muc, seg.area_cm2, compound.mw, perm.vmax_scale) * Fmuc

    # --- carriers influx ---
    J_carr: Dict[str, float] = {}
    for name in ("OCT3", "PEPT1", "ENT", "SGLT1", "OATP2B1"):
        vmax, km = vmax_km(settings, name)
        expr = seg.expr.get(name, 1.0)
        bidir = (name == "ENT")
        key = {"OCT3": "oct", "PEPT1": "pept1", "ENT": "ent", "SGLT1": "sglt", "OATP2B1": "oatp"}[name]
        J_carr[key] = transporter_flux_mg_h(
            name, compound, z, C_free, C_wall, seg.area_cm2, vmax, km,
            expr=expr, vmax_scale=perm.vmax_scale, bidirectional=bidir,
        ) * Fmuc * rel

    # --- efflux depuis wall ---
    vmax_pgp, km_pgp = vmax_km(settings, "PGP")
    vmax_bcrp, km_bcrp = vmax_km(settings, "BCRP")
    pgp_inhib = max(formulation.pgp_inhibitor, formulation.tpgs_mm * fx_cfg.tpgs_pgp_factor / 10.0)
    pgp_inhib = min(1.0, pgp_inhib)
    J_pgp = pgp_efflux_mg_h(compound, z, C_wall, seg.area_cm2, vmax_pgp, km_pgp,
                            expr=seg.expr.get("PGP", 1.0), vmax_scale=perm.vmax_scale,
                            inhibitor=pgp_inhib)
    J_bcrp = bcrp_efflux_mg_h(compound, z, C_wall, seg.area_cm2, vmax_bcrp, km_bcrp,
                              expr=seg.expr.get("BCRP", 1.0), vmax_scale=perm.vmax_scale)
    J_eff = J_pgp + J_bcrp

    # --- mucus exchange ---
    J_lm = K_MUC_H * max(A_diss, 0.0)
    J_ml = K_MUC_H * max(A_mucus, 0.0)
    J_muc_clear = K_MUC_CLEAR_H * max(A_mucus, 0.0)

    # --- wall to blood ---
    J_wall_out = K_WALL_H * max(A_wall, 0.0)
    J_to_blood = J_wall_out * fg
    J_wall_met = J_wall_out * (1 - fg)

    # --- transit ---
    k_tr = 1.0 / max(seg.transit_h, 1e-9)
    J_tr_diss = k_tr * max(A_diss, 0.0)
    J_tr_solid = k_tr * max(A_solid, 0.0)

    # --- TEER target ---
    drops = []
    enh_map = {"C10": formulation.c10_mm, "EDTA": formulation.edta_mm, "C8": formulation.c8_mm, "C12": formulation.c12_mm, "BILE": seg.bile_mm + formulation.bile_mm,
               "Labrasol": formulation.labrasol_mm, "NaCDC": formulation.nacdc_mm}
    for ename, conc in enh_map.items():
        h = fx_cfg.enhancers.get(ename)
        if h is None or conc <= 0:
            continue
        drops.append(hill_drop(conc * max(rel, 0.2), h.ec50_mm, h.dmax, h.n))
    e_tot = bliss_total(drops)
    teer_tgt = teer_target(seg.teer0, e_tot)

    fluxes = {
        "J_diss": J_diss, "J_prec": J_prec, "J_trans": max(J_trans, 0.0), "J_para": max(J_para, 0.0),
        "J_endo": max(J_endo, 0.0), **{f"J_{k}": max(v, 0.0) if k != "ent" else v for k, v in J_carr.items()},
        "J_pgp": J_pgp, "J_bcrp": J_bcrp, "J_eff": J_eff,
        "J_lm": J_lm, "J_ml": J_ml, "J_muc_clear": J_muc_clear,
        "J_to_blood": J_to_blood, "J_wall_met": J_wall_met,
        "J_tr_diss": J_tr_diss, "J_tr_solid": J_tr_solid,
        "teer_tgt": teer_tgt, "C_diss": C_diss, "C_wall": C_wall, "z": z,
        "f_union": f_union, "Fmuc": Fmuc, "rel": rel, "cs": cs,
    }
    return fluxes, {"teer_tgt": teer_tgt}
