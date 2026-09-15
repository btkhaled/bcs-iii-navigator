"""Conversions d'unités centralisées.

Aucune conversion en dur ailleurs dans le code.
Tous les facteurs viennent d'ici (constantes physiques pures).
"""

# --- constantes physiques (CODATA, pas des paramètres métier) ---
AVOGADRO = 6.02214076e23
GAS_R = 8.314462618  # J / mol / K
FARADAY = 96485.33212  # C / mol
BODY_TEMP_K = 310.15  # 37 C

# --- conversions ---
# P[cm/s] * C[mg/L] * A[cm2] -> mg/h
# 1 cm/s * 1 cm2 = 1 cm3/s = 1 mL/s ; 1 mg/L = 1 ug/mL
# mg/h = P_cm_s * C_mg_L * A_cm2 * (1 L / 1000 mL) * 3600 s/h
#      = P * C * A * 3.6
CM_S_MG_L_CM2_TO_MG_H = 3.6

# Vmax[pmol/cm2/min] * A[cm2] -> mg/h pour une MW donnée :
# pmol/min * 1e-9 mmol/pmol... en pratique :
# mg/h = Vmax * A * 1e-9 (mmol/pmol->mol?) * MW(g/mol) * 1e3(mg/g) * 60(min/h) * 1e-12?
# Forme explicite : pmol = 1e-12 mol ; mg = MW * mol * 1e3
# mg/h = Vmax[pmol/cm2/min] * A * 1e-12 * MW * 1e3 * 60
PMOL_CM2_MIN_TO_MG_H_PER_MW = 1e-12 * 1e3 * 60.0  # x MW x A x Vmax

# uM -> mg/L : mg/L = uM * MW / 1000
UM_TO_MG_L_PER_MW = 1.0 / 1000.0

# mg/L -> uM : uM = mg/L * 1000 / MW
MG_L_TO_UM_PER_MW = 1000.0


def pmol_cm2_min_to_mg_h(vmax: float, area_cm2: float, mw: float) -> float:
    return vmax * area_cm2 * PMOL_CM2_MIN_TO_MG_H_PER_MW * mw


def um_to_mg_l(um: float, mw: float) -> float:
    return um * mw * UM_TO_MG_L_PER_MW


def mg_l_to_um(mg_l: float, mw: float) -> float:
    if mw <= 0:
        raise ValueError("MW must be > 0")
    return mg_l * MG_L_TO_UM_PER_MW / mw


def cm_s_mg_l_cm2_to_mg_h(p_cm_s: float, c_mg_l: float, a_cm2: float) -> float:
    return p_cm_s * c_mg_l * a_cm2 * CM_S_MG_L_CM2_TO_MG_H
