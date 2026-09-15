"""Provenance: hash reproductible code+configs+inputs -> run_id."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


def canonical(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str)


def sha256_str(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def run_id(*parts: Any) -> str:
    h = hashlib.sha256()
    for p in parts:
        if isinstance(p, (dict, list)):
            h.update(canonical(p).encode())
        elif isinstance(p, Path):
            h.update(p.read_bytes() if p.exists() else str(p).encode())
        else:
            h.update(str(p).encode())
    return h.hexdigest()[:16]
