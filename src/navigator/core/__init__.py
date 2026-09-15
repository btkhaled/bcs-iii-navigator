"""Core: infra logicielle, aucune science ici."""

from .errors import (
    ConfigError,
    DomainError,
    SolverError,
    ConvergenceError,
    ValidationError,
)

__all__ = [
    "ConfigError",
    "DomainError",
    "SolverError",
    "ConvergenceError",
    "ValidationError",
]
