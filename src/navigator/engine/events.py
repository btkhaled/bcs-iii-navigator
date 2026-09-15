"""Events: doses / repas (hook V1: dose unique t0, extension multi-dose)."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List


@dataclass
class DosingEvent:
    t_h: float
    amount_mg: float
    segment: int = 0


@dataclass
class Scenario:
    fed: bool = False
    events: List[DosingEvent] = field(default_factory=list)
