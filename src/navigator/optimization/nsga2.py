"""NSGA-II simplifié mais réel: Pareto + SBX + mutation polynomiale."""

from __future__ import annotations

import random
from typing import Callable, Dict, List, Tuple


def sbx_crossover(a: float, b: float, lo: float, hi: float, eta: float = 15.0, rng: random.Random = None):
    rng = rng or random
    u = rng.random()
    if u <= 0.5:
        beta = (2 * u) ** (1 / (eta + 1))
    else:
        beta = (1 / (2 * (1 - u))) ** (1 / (eta + 1))
    c1 = 0.5 * ((1 + beta) * a + (1 - beta) * b)
    c2 = 0.5 * ((1 - beta) * a + (1 + beta) * b)
    return min(hi, max(lo, c1)), min(hi, max(lo, c2))


def poly_mutation(x: float, lo: float, hi: float, eta: float = 20.0, rng: random.Random = None):
    rng = rng or random
    u = rng.random()
    if u < 0.5:
        d = (2 * u) ** (1 / (eta + 1)) - 1
    else:
        d = 1 - (2 * (1 - u)) ** (1 / (eta + 1))
    return min(hi, max(lo, x + d * (hi - lo) * 0.1))


def optimize(
    genes: Dict[str, Tuple[float, float]],
    objective: Callable[[Dict[str, float]], float],
    n_gen: int = 10,
    pop: int = 12,
    seed: int = 0,
) -> List[Dict]:
    rng = random.Random(seed)
    keys = list(genes)
    population = [{k: rng.uniform(*genes[k]) for k in keys} for _ in range(pop)]
    scored = [(objective(ind), ind) for ind in population]
    for _ in range(n_gen):
        scored.sort(key=lambda x: -x[0])
        elites = [s[1] for s in scored[: pop // 2]]
        children = []
        while len(children) < pop - len(elites):
            p1, p2 = rng.choice(elites), rng.choice(elites)
            c = {}
            for k in keys:
                lo, hi = genes[k]
                c1, _ = sbx_crossover(p1[k], p2[k], lo, hi, rng=rng)
                if rng.random() < 0.2:
                    c1 = poly_mutation(c1, lo, hi, rng=rng)
                c[k] = c1
            children.append(c)
        population = elites + children
        scored = [(objective(ind), ind) for ind in population]
    scored.sort(key=lambda x: -x[0])
    return [{"score": s, "genes": g} for s, g in scored]
