"""Compatibility imports for the independent lifecycle audit."""

from .audit import lifecycle as _implementation


def __getattr__(name):
    return getattr(_implementation, name)
