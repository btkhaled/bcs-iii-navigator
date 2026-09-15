"""Use-cases applicatifs."""

from __future__ import annotations

from pathlib import Path

from ..core.loader import load_compound, load_formulation
from ..core.settings import Settings
from ..engine.runner import run


def predict(compound_path: str, formulation_path: str, config_dir: str):
    s = Settings(config_dir)
    c = load_compound(compound_path)
    f = load_formulation(formulation_path)
    return run(c, f, s)
