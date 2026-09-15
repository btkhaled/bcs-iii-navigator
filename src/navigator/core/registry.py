"""Registry générique pour plugins (transporteurs, excipients, dissolution)."""

from __future__ import annotations

from typing import Callable, Dict, TypeVar

T = TypeVar("T")


class Registry:
    def __init__(self, name: str):
        self.name = name
        self._items: Dict[str, object] = {}

    def register(self, key: str):
        def deco(obj: T) -> T:
            if key in self._items:
                raise KeyError(f"[{self.name}] '{key}' déjà enregistré")
            self._items[key] = obj
            return obj

        return deco

    def get(self, key: str):
        try:
            return self._items[key]
        except KeyError as e:
            raise KeyError(f"[{self.name}] inconnu: '{key}'. Connus: {sorted(self._items)}") from e

    def keys(self):
        return sorted(self._items)

    def __contains__(self, key: str) -> bool:
        return key in self._items


DISSOLUTION_REGISTRY = Registry("dissolution")
TRANSPORTER_REGISTRY = Registry("transporter")
EXCIPIENT_REGISTRY = Registry("excipient")
KP_REGISTRY = Registry("kp")
