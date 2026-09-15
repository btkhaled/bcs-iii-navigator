"""Transporteurs: Michaelis-Menten + grilles substrat + expression segmentaire."""

from __future__ import annotations

import math

from ...core.registry import TRANSPORTER_REGISTRY
from ...core.schemas import Compound
from ...core.units import pmol_cm2_min_to_mg_h, mg_l_to_um


def mm_rate(vmax: float, km: float, c: float) -> float:
    if vmax <= 0 or km <= 0 or c <= 0:
        return 0.0
    if not math.isfinite(vmax) or not math.isfinite(km) or not math.isfinite(c):
        return 0.0
    return vmax * c / (km + c)


# --- grilles substrat (mécanistiques, depuis la biologie des transporteurs) ---

def is_oct_substrate(c: Compound, z: float) -> bool:
    if not (0.2 <= z <= 2.0):
        return False
    if not (100 <= c.mw <= 350):
        return False
    if c.logp >= 2.0:
        return False
    return True


def is_pept1_substrate(c: Compound) -> bool:
    if c.is_biguanide:
        return False
    if not c.has_peptide_bond:
        return False
    if not c.has_alpha_amino:
        # pharmacophore PEPT1: amine alpha N-terminale + C=O + acide C-terminal
        return False
    if c.hbd < 3:
        return False
    return 150 <= c.mw <= 700


def is_ent_substrate(c: Compound) -> bool:
    if not c.has_azote_base:
        return False
    return 200 <= c.mw <= 350


def is_sglt_substrate(c: Compound) -> bool:
    if not c.is_polyol:
        return False
    if c.logp >= 0:
        return False
    if c.hbd < 3:
        return False
    return 150 <= c.mw <= 350


def is_oatp_substrate(c: Compound, z: float) -> bool:
    # anion ET lipophile (porte ET)
    if z >= -0.2:
        return False
    if c.logp < 1.0:
        return False
    return True


def is_pgp_substrate(c: Compound) -> bool:
    # amples/cationiques-amphiphiles uniquement : exclut polycations hydrophiles (aminosides)
    if c.logp < 0.0:
        return False
    return (c.n_o_count >= 8) and (c.mw > 400)


def is_bcrp_substrate(c: Compound, z: float) -> bool:
    return z < -0.2 and c.mw > 300


SUBSTRATE_FNS = {
    "OCT3": lambda c, z: is_oct_substrate(c, z),
    "PEPT1": lambda c, z: is_pept1_substrate(c),
    "ENT": lambda c, z: is_ent_substrate(c),
    "SGLT1": lambda c, z: is_sglt_substrate(c),
    "OATP2B1": lambda c, z: is_oatp_substrate(c, z),
    "PGP": lambda c, z: is_pgp_substrate(c),
    "BCRP": lambda c, z: is_bcrp_substrate(c, z),
}


def transporter_flux_mg_h(
    name: str,
    compound: Compound,
    z: float,
    c_lumen_mg_l: float,
    c_cell_mg_l: float,
    area_cm2: float,
    vmax_pmol_cm2_min: float,
    km_um: float,
    expr: float = 1.0,
    vmax_scale: float = 1.0,
    bidirectional: bool = False,
    pgp_inhibitor: float = 0.0,
) -> float:
    """Flux net influx (positif = vers entérocyte). P-gp/BCRP gérés en efflux séparé."""
    fn = SUBSTRATE_FNS.get(name)
    if fn is not None and not fn(compound, z):
        return 0.0
    # overrides molécule
    vmax = compound.vmax_override.get(name, vmax_pmol_cm2_min) * expr * vmax_scale
    km = compound.km_override.get(name, km_um)
    if vmax <= 0:
        return 0.0
    c_um = mg_l_to_um(max(c_lumen_mg_l, 0.0), compound.mw)
    v = mm_rate(vmax, km, c_um)  # pmol/cm2/min
    j = pmol_cm2_min_to_mg_h(v, area_cm2, compound.mw)
    if bidirectional:
        c_cell_um = mg_l_to_um(max(c_cell_mg_l, 0.0), compound.mw)
        v_back = mm_rate(vmax, km, c_cell_um)
        j -= pmol_cm2_min_to_mg_h(v_back, area_cm2, compound.mw)
    if name in ("PGP", "BCRP") and pgp_inhibitor > 0:
        j *= (1.0 - 0.9 * min(1.0, pgp_inhibitor))
    return max(0.0, j) if not bidirectional else j
