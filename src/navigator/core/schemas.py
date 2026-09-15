"""Schémas pydantic stricts. Tout JSON est validé ici, extra=forbid."""

from __future__ import annotations

from typing import Dict, List, Literal, Optional

from pydantic import BaseModel, Field


class Strict(BaseModel):
    model_config = {"extra": "forbid"}


class PkaGroup(Strict):
    type: Literal["acid", "base"]
    pka: float


class Compound(Strict):
    id: str
    mw: float = Field(gt=0)
    logp: float = Field(default=0.0)
    tpsa: float = Field(default=0.0, ge=0)
    hbd: int = Field(default=0, ge=0)
    hba: int = Field(default=0, ge=0)
    rotb: int = Field(default=0, ge=0)
    aromatic_rings: int = Field(default=0, ge=0)
    frac_sp3: float = Field(default=0.0, ge=0, le=1)
    n_o_count: int = Field(default=0, ge=0)
    formal_charge: float = Field(default=0.0)
    pkas: List[PkaGroup] = Field(default_factory=list)
    cs_mg_l: float = Field(gt=0, description="solubilité intrinsèque mg/L")
    has_peptide_bond: bool = False
    has_alpha_amino: bool = False
    is_biguanide: bool = False
    has_azote_base: bool = False
    is_polyol: bool = False
    fu_plasma: float = Field(default=0.9, gt=0, le=1)
    clint_liver_ul_min_mg: float = Field(default=0.0, ge=0)
    clint_gut_ul_min_mg: float = Field(default=0.0, ge=0)
    vmax_override: Dict[str, float] = Field(default_factory=dict)
    km_override: Dict[str, float] = Field(default_factory=dict)


class Formulation(Strict):
    id: str = "native"
    dose_mg: float = Field(gt=0)
    f_free: float = Field(default=1.0, gt=0, le=1)
    # HIP
    m_hip: float = Field(default=1.0, ge=1.0)
    hip_counterion: Optional[str] = None
    # SEDDS / NP
    lipid_frac: float = Field(default=0.0, ge=0, le=1)
    particle_nm: float = Field(default=0.0, ge=0)
    zeta_mv: float = Field(default=0.0)
    # enhancers (mM luminal cible)
    c10_mm: float = Field(default=0.0, ge=0)
    snac_mm: float = Field(default=0.0, ge=0)
    edta_mm: float = Field(default=0.0, ge=0)
    c8_mm: float = Field(default=0.0, ge=0)
    c12_mm: float = Field(default=0.0, ge=0)
    bile_mm: float = Field(default=0.0, ge=0)
    chitosan: float = Field(default=0.0, ge=0, le=1)
    polymer: Optional[str] = None
    # voie A: nouveaux excipients réels
    labrasol_mm: float = Field(default=0.0, ge=0)
    tpgs_mm: float = Field(default=0.0, ge=0)
    pip_mm: float = Field(default=0.0, ge=0)
    nacdc_mm: float = Field(default=0.0, ge=0)
    # voie B: excipient théorique de-novo (spec minimale viable)
    denovo_p_trans: float = Field(default=1.0, ge=1.0, le=10.0)
    denovo_p_para: float = Field(default=1.0, ge=1.0, le=10.0)
    denovo_mucus: float = Field(default=1.0, ge=1.0, le=5.0)
    # release
    k_release_h: float = Field(default=10.0, gt=0)
    enteric: float = Field(default=0.0, ge=0, le=1)
    enteric_ph_threshold: float = Field(default=5.5)
    pgp_inhibitor: float = Field(default=0.0, ge=0, le=1)


class SegmentPhys(Strict):
    name: str
    ph: float
    transit_h: float = Field(gt=0)
    volume_l: float = Field(gt=0)
    area_cm2: float = Field(gt=0)
    radius_cm: float = Field(gt=0)
    teer0: float = Field(gt=0)
    mucus_thickness_um: float = Field(default=100.0, ge=0)
    bile_mm: float = Field(default=3.0, ge=0)
    fecal_bound: float = Field(default=0.0, ge=0, le=0.99)
    # expression relative transporteurs (multiplicateur de Vmax)
    expr: Dict[str, float] = Field(default_factory=dict)


class Physiology(Strict):
    segments: List[SegmentPhys]
    gastric_emptying_h: float = Field(default=0.25, gt=0)
    colon_transit_h: float = Field(default=18.0, gt=0)


class TransporterDef(Strict):
    vmax_pmol_cm2_min: float = Field(ge=0)
    km_um: float = Field(gt=0)


class PermeationConfig(Strict):
    pore_nm: float = Field(gt=0)
    pore_large_nm: float = Field(default=5.0, gt=0)
    frac_large_pore: float = Field(default=0.001, ge=0, le=1)
    p_base_cm_s: float = Field(gt=0)
    qspr_a: float
    qspr_b: float
    uwl_thickness_um: float = Field(default=30.0, ge=0)
    membrane_potential_mv: float = Field(default=-60.0)
    goldman: bool = True
    vmax_scale: float = Field(default=1.0, gt=0)
    p_scale: float = Field(default=1.0, gt=0)
    t_scale: float = Field(default=100.0, gt=0)
    endo_transcytosis: float = Field(default=0.01, ge=0, le=1)


class FirstPassConfig(Strict):
    qh_l_h: float = Field(gt=0)
    vgut_l: float = Field(gt=0)
    tgut_h: float = Field(gt=0)


class SystemicConfig(Strict):
    vd_l: float = Field(gt=0)
    cl_l_h: float = Field(ge=0)


class NumericsConfig(Strict):
    method: str = "LSODA"
    t_end_h: float = Field(default=96.0, gt=0)
    rtol: float = Field(default=1e-6, gt=0)
    atol: float = Field(default=1e-9, gt=0)
    mass_tol_mg: float = Field(default=1e-6, gt=0)
    n_eval: int = Field(default=400, gt=10)


class ExcipientHill(Strict):
    ec50_mm: float = Field(gt=0)
    dmax: float = Field(ge=0, le=1)
    n: float = Field(default=1.5, gt=0)


class FormulationEffects(Strict):
    enhancers: Dict[str, ExcipientHill]
    snac_factor_per_mm: float = Field(default=0.02, ge=0)
    pip_factor_per_mm: float = Field(default=0.015, ge=0)
    tpgs_pgp_factor: float = Field(default=0.9, ge=0, le=1)
    lipid_factor_per_frac: float = Field(default=1.0, ge=0)
    chitosan_teer_factor: float = Field(default=2.0, ge=0)
    teer_max_drop: float = Field(default=10.0, gt=0)
    tau_on_h: float = Field(default=0.125, gt=0)
    tau_rec_h: float = Field(default=1.0, gt=0)
    d_crit_nm: float = Field(default=150.0, gt=0)


class ResultBreakdown(Strict):
    trans: float = 0.0
    para: float = 0.0
    endo: float = 0.0
    oct: float = 0.0
    sglt: float = 0.0
    ent: float = 0.0
    pept1: float = 0.0
    oatp: float = 0.0
    pgp: float = 0.0


class SimResult(Strict):
    f_oral: float
    f_g: float
    f_h: float
    auc_mg_h_l: float
    cmax_mg_l: float
    breakdown: ResultBreakdown
    dominant: str
    z: float
    f_unionized: float
    mass_balance_mg: float
    run_id: str
    trace: list = Field(default_factory=list)
    warnings: list = Field(default_factory=list)
