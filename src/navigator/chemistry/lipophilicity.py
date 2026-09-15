"""logD(pH) apparent."""

from __future__ import annotations

import math

from ..core.schemas import Compound


def logd_at_ph(compound: Compound, ph: float) -> float:
    logp = compound.logp
    # correction globale: somme des contributions (approximation sites indépendants)
    corr = 0.0
    for g in compound.pkas:
        if g.type == "acid":
            corr += math.log10(1.0 + 10.0 ** (ph - g.pka))
        else:
            corr += math.log10(1.0 + 10.0 ** (g.pka - ph))
    return logp - corr
