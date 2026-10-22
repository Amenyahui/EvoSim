"""Read safe YAML; experiment files may extend a shared model configuration."""

from copy import deepcopy
from pathlib import Path

import yaml


def read_config(path, _parents=()):
    path = Path(path).resolve()
    if path in _parents:
        raise ValueError("Cyclic configuration inheritance")
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise ValueError("Configuration must be a mapping")
    parent = config.pop("extends", None)
    if parent is None:
        return config
    base = read_config(path.parent / parent, (*_parents, path))
    for key, value in config.items():
        if key in ("simulation_metadata", "parameters"):
            base[key] = {**base.get(key, {}), **value}
        else:
            base[key] = deepcopy(value)
    return base
