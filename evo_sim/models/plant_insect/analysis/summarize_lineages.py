"""Compatibility entry point; implementation lives in summaries.lineages."""

from .summaries import lineages as _implementation


def __getattr__(name):
    return getattr(_implementation, name)


if __name__ == "__main__":
    _implementation.main()
