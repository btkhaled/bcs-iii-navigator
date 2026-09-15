"""Hiérarchie d'exceptions métier."""


class NavigatorError(Exception):
    """Base."""


class ConfigError(NavigatorError):
    """Config JSON manquante / invalide."""


class DomainError(NavigatorError):
    """Valeur physiquement impossible."""


class SolverError(NavigatorError):
    """Échec d'intégration ODE."""


class ConvergenceError(SolverError):
    """Non-convergence / bilan masse violé."""


class ValidationError(NavigatorError):
    """Validation schéma / contrat."""
