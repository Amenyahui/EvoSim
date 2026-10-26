"""State paths and source identity for configured substeps."""

import hashlib
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def get_by_path(state, path):
    for key in path.split("/"):
        state = state[key]
    return state


def set_by_path(state, path, value):
    keys = path.split("/")
    for key in keys[:-1]:
        state = state[key]
    if keys[-1] not in state:
        raise KeyError(f"Undeclared state output: {path}")
    state[keys[-1]] = value


def source_hashes():
    package = PROJECT_ROOT / "evo_sim"
    return {
        p.relative_to(package).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(package.rglob("*"))
        if p.suffix in {".py", ".slim", ".eidos"}
    }
