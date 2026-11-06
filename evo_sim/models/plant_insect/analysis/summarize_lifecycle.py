"""Compatibility entry point; implementation lives in summaries.lifecycle."""

from .summaries import lifecycle as _implementation


def __getattr__(name):
    return getattr(_implementation, name)


if __name__ == "__main__":
    _implementation.main()
