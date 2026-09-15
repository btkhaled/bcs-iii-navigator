"""Settings: charge tous les configs/*.json, validés."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

from .errors import ConfigError
from .schemas import (
    FirstPassConfig,
    FormulationEffects,
    NumericsConfig,
    PermeationConfig,
    Physiology,
    SystemicConfig,
)
from .loader import load_json


class Settings:
    def __init__(self, config_dir: str | Path):
        self.config_dir = Path(config_dir)
        try:
            self.physiology = Physiology(**load_json(self.config_dir / "gi_physiology.json"))
            self.permeation = PermeationConfig(**load_json(self.config_dir / "permeation.json"))
            self.first_pass = FirstPassConfig(**load_json(self.config_dir / "first_pass.json"))
            self.systemic = SystemicConfig(**load_json(self.config_dir / "systemic.json"))
            self.numerics = NumericsConfig(**load_json(self.config_dir / "numerics.json"))
            self.formulation_fx = FormulationEffects(
                **load_json(self.config_dir / "formulation_effects.json")
            )
            self.transporters: Dict[str, Dict[str, Any]] = load_json(
                self.config_dir / "transporters.json"
            )
            self.calibration: Dict[str, Any] = load_json(
                self.config_dir / "calibration.json"
            )
            self.dissolution: Dict[str, Any] = load_json(
                self.config_dir / "dissolution.json"
            )
        except FileNotFoundError as e:
            raise ConfigError(str(e)) from e
