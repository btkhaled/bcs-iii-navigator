"""Batch: N runs."""

from __future__ import annotations

from .predict import predict


def batch(pairs, config_dir: str):
    return [predict(c, f, config_dir) for c, f in pairs]
