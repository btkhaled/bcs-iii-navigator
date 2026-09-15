"""Loader JSON unique avec validation."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .errors import ConfigError
from .schemas import Compound, Formulation


def load_json(path: str | Path) -> Any:
    p = Path(path)
    if not p.exists():
        raise ConfigError(f"JSON introuvable: {p}")
    return json.loads(p.read_text(encoding="utf-8"))


def load_compound(path: str | Path) -> Compound:
    return Compound(**load_json(path))


def load_formulation(path: str | Path) -> Formulation:
    return Formulation(**load_json(path))


def write_json(path: str | Path, obj: Any) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, indent=2, sort_keys=True), encoding="utf-8")
