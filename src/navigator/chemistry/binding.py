"""Fraction libre / liaison."""

from __future__ import annotations


def free_fraction(fu_plasma: float) -> float:
    return max(1e-4, min(1.0, fu_plasma))
