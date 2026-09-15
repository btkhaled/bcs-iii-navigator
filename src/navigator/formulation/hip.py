"""HIP dynamics: M_HIP dépend pH / concentration / temps (mécanistique simplifié)."""

from __future__ import annotations


def m_hip_dynamic(m_hip_base: float, ph: float, c_mg_l: float, t_h: float, counterion: str | None = None) -> float:
    # dissociation à pH extrêmes + dilution
    # facteur pH: optimal 5.5-7.0
    if ph < 4.5 or ph > 8.0:
        ph_f = 0.6
    elif 5.5 <= ph <= 7.0:
        ph_f = 1.0
    else:
        ph_f = 0.85
    # dilution: si très dilué, dissociation partielle
    dil_f = 1.0 if c_mg_l > 10 else 0.7 + 0.3 * (c_mg_l / 10.0)
    m = 1.0 + (m_hip_base - 1.0) * ph_f * dil_f
    return max(1.0, min(32.0, m))
