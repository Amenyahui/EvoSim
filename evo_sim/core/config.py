"""Read safe YAML; experiment files may extend a shared model configuration."""

from copy import deepcopy
from pathlib import Path

import yaml


def merge_config(base, overrides):
    """Merge nested mappings without aliases; lists/scalars (including null) replace."""
    if not isinstance(base, dict) or not isinstance(overrides, dict):
        raise ValueError("Configuration and overrides must be mappings")
    result = deepcopy(base)
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = merge_config(result[key], value)
        else:
            result[key] = deepcopy(value)
    return result


def read_config(path, _parents=(), *, input_loader=None):
    path = Path(path).resolve()
    if path in _parents:
        raise ValueError("Cyclic configuration inheritance")
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise ValueError("Configuration must be a mapping")
    parent = config.pop("extends", None)
    base = (
        read_config(path.parent / parent, (*_parents, path), input_loader=input_loader)
        if parent is not None
        else {}
    )
    files = config.pop("input_files", None)
    if files is not None:
        if input_loader is None:
            raise ValueError("CSV inputs require the model configuration reader")
        parameters = merge_config(base.get("parameters", {}), config.get("parameters", {}))
        base["parameters"], sources = input_loader(files, path.parent, parameters)
        base["input_sources"] = {**base.get("input_sources", {}), **sources}
    for key, value in config.items():
        if key in ("simulation_metadata", "parameters", "batch"):
            base[key] = merge_config(base.get(key, {}), value)
        else:
            base[key] = deepcopy(value)
    return base
