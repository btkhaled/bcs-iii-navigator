"""Audit trail: chaque calcul = (formule, entrées, sortie, unité)."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List


@dataclass
class TraceEntry:
    name: str
    formula: str
    inputs: Dict[str, Any]
    output: Any
    unit: str = ""


@dataclass
class Trace:
    entries: List[TraceEntry] = field(default_factory=list)

    def log(self, name: str, formula: str, inputs: Dict[str, Any], output: Any, unit: str = ""):
        self.entries.append(TraceEntry(name, formula, inputs, output, unit))

    def to_list(self) -> List[Dict[str, Any]]:
        return [
            {
                "name": e.name,
                "formula": e.formula,
                "inputs": e.inputs,
                "output": e.output,
                "unit": e.unit,
            }
            for e in self.entries
        ]
