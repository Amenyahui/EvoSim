"""Compatibility entry point; implementation lives in plots.diagnostics."""

from .plots import diagnostics as _implementation


def __getattr__(name):
    return getattr(_implementation, name)


if __name__ == "__main__":
    _implementation.main()
