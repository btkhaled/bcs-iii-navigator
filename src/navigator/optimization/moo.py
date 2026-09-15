"""Vrai NSGA-II contraint (Deb 2002): tri non-dominé + crowding + tournoi.

`nsga2.py` (GA élitiste mono-objectif pour formulations) est inchangé.
Ce module sert la calibration multi-objectifs sous contraintes.
"""

from __future__ import annotations

import random
from typing import Callable, List, Tuple

from .nsga2 import poly_mutation, sbx_crossover

BIG = 1e6


def _dominates(a: List[float], b: List[float]) -> bool:
    return all(x <= y for x, y in zip(a, b)) and any(x < y for x, y in zip(a, b))


def _cdom(ai: Tuple[List[float], int], aj: Tuple[List[float], int]) -> bool:
    (oa, va), (ob, vb) = ai, aj
    if va == 0 and vb == 0:
        return _dominates(oa, ob)
    if va == 0:
        return True
    if vb == 0:
        return False
    return va < vb


def non_dominated_sort(items: List[Tuple[List[float], int]]) -> List[List[int]]:
    n = len(items)
    succ = [[] for _ in range(n)]
    ndom = [0] * n
    fronts: List[List[int]] = [[]]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if _cdom(items[i], items[j]):
                succ[i].append(j)
            elif _cdom(items[j], items[i]):
                ndom[i] += 1
        if ndom[i] == 0:
            fronts[0].append(i)
    k = 0
    while k < len(fronts) and fronts[k]:
        nxt = []
        for i in fronts[k]:
            for j in succ[i]:
                ndom[j] -= 1
                if ndom[j] == 0:
                    nxt.append(j)
        k += 1
        fronts.append(nxt)
    return [f for f in fronts if f]


def crowding(objs: List[List[float]], front: List[int]) -> dict:
    dist = {i: 0.0 for i in front}
    if len(front) <= 2:
        for i in front:
            dist[i] = float("inf")
        return dist
    m = len(objs[0])
    for k in range(m):
        vals = sorted(front, key=lambda i: objs[i][k])
        dist[vals[0]] = float("inf")
        dist[vals[-1]] = float("inf")
        rng = objs[vals[-1]][k] - objs[vals[0]][k]
        if rng <= 0:
            continue
        for a, b, c in zip(vals, vals[1:], vals[2:]):
            dist[b] += (objs[c][k] - objs[a][k]) / rng
    return dist


def optimize_moo(
    evaluate_many: Callable[[List[List[float]]], List[Tuple[List[float], int, dict]]],
    bounds: List[Tuple[float, float]],
    pop: int = 40,
    gen: int = 40,
    seed: int = 0,
    seed_with: List[List[float]] | None = None,
    cx_eta: float = 15.0,
    mut_eta: float = 20.0,
    mut_p: float = 0.2,
    log_every: int = 10,
    log_fn: Callable[[int, dict], None] | None = None,
) -> dict:
    """Espace normalisé [lo,hi] par gène. evaluate_many(xs) -> [(objs, viol, aux)]."""
    rng = random.Random(seed)

    def rand_ind():
        return [rng.uniform(lo, hi) for lo, hi in bounds]

    # LHS seeding for diversity: stratify the uniform sampling
    if pop >= 8:
        xs = []
        for i in range(pop):
            xs.append([(lo + (rng.random() + i) % pop / pop * (hi - lo)) for lo, hi in bounds])
        xs = xs[:pop]
    else:
        xs = [rand_ind() for _ in range(pop)]
    if seed_with:
        for i, w in enumerate(seed_with[:pop]):
            xs[i] = list(w)
    res = evaluate_many(xs)
    pop_state = [{"x": x, "objs": list(o), "viol": v, "aux": a} for x, (o, v, a) in zip(xs, res)]

    def tournament(pool):
        a, b = rng.choice(pool), rng.choice(pool)
        ia = (a["objs"], a["viol"])
        ib = (b["objs"], b["viol"])
        if _cdom(ia, ib) and not _cdom(ib, ia):
            return a
        if _cdom(ib, ia) and not _cdom(ia, ib):
            return b
        # même front -> crowding
        if a.get("crowd", 0) >= b.get("crowd", 0):
            return a
        return b

    for g in range(gen):
        items = [(p["objs"], p["viol"]) for p in pop_state]
        fronts = non_dominated_sort(items)
        objs_only = [p["objs"] for p in pop_state]
        for f in fronts:
            cd = crowding(objs_only, f)
            for i in f:
                pop_state[i]["crowd"] = cd[i]
                pop_state[i]["rank"] = fronts.index(f)
        feas = [p for p in pop_state if p["viol"] == 0]
        cur_best = min(feas, key=lambda p: p["objs"][0])["objs"][0] if feas else 1e9
        if g == 0:
            best_seen = cur_best
            stall = 0
        elif cur_best < best_seen - 1e-4:
            best_seen = cur_best
            stall = 0
        else:
            stall += 1
        if stall >= 12:
            if log_fn is not None:
                log_fn(g, {"best_rmse": cur_best, "viol": 0, "n_feas": len(feas), "early_stop": True})
            break
        # variation
        children = []
        while len(children) < pop:
            p1, p2 = tournament(pop_state), tournament(pop_state)
            c1, c2 = [], []
            for k, (lo, hi) in enumerate(bounds):
                v1, v2 = sbx_crossover(p1["x"][k], p2["x"][k], lo, hi, eta=cx_eta, rng=rng)
                c1.append(v1)
                c2.append(v2)
            for c in (c1, c2):
                for k, (lo, hi) in enumerate(bounds):
                    if rng.random() < mut_p:
                        c[k] = poly_mutation(c[k], lo, hi, eta=mut_eta, rng=rng)
                children.append(c)
                if len(children) >= pop:
                    break
        res_c = evaluate_many(children)
        combined = pop_state + [
            {"x": x, "objs": list(o), "viol": v, "aux": a} for x, (o, v, a) in zip(children, res_c)
        ]
        items = [(p["objs"], p["viol"]) for p in combined]
        fronts = non_dominated_sort(items)
        objs_only = [p["objs"] for p in combined]
        nxt = []
        for f in fronts:
            if len(nxt) + len(f) <= pop:
                nxt.extend(f)
            else:
                cd = crowding(objs_only, f)
                f_sorted = sorted(f, key=lambda i: -cd[i])
                nxt.extend(f_sorted[: pop - len(nxt)])
                break
        pop_state = [combined[i] for i in nxt]
        if log_fn is not None and (g % log_every == 0 or g == gen - 1):
            feas2 = [p for p in pop_state if p["viol"] == 0]
            best = min(feas2, key=lambda p: p["objs"][0]) if feas2 else min(pop_state, key=lambda p: (p["viol"], p["objs"][0]))
            log_fn(g, {"best_rmse": best["objs"][0], "viol": best["viol"], "n_feas": len(feas2)})

    items = [(p["objs"], p["viol"]) for p in pop_state]
    fronts = non_dominated_sort(items)
    objs_only = [p["objs"] for p in pop_state]
    for ri, f in enumerate(fronts):
        cd = crowding(objs_only, f)
        for i in f:
            pop_state[i]["rank"] = ri
            pop_state[i]["crowd"] = cd[i]
    return {"population": pop_state, "fronts": fronts}
