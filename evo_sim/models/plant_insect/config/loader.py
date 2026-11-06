"""Resolve model YAML inheritance and its declaring-file-relative CSV inputs."""

from os import PathLike

from evo_sim.core.config import read_config as read_yaml_config

from .input_loader import load_inputs


def read_config(path: str | PathLike) -> dict:
    """Resolve YAML inheritance and CSV paths relative to their declaring YAML."""
    return read_yaml_config(path, input_loader=load_inputs)
