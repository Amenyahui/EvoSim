import copy
from pathlib import Path

from evo_sim.models.plant_insect import Config, get_runner


def base_runner():
    return get_runner(
        Path(__file__).resolve().parents[1]
        / "evo_sim/models/plant_insect/configs/current_base.yaml"
    )


def config(**changes):
    values = copy.deepcopy(base_runner().config["parameters"])
    values.update(generations=2, plant_sizes=[20, 20], insect_sizes=[6])
    values.update(changes)
    return Config(**values)


def raises(function, error=ValueError):
    try:
        function()
    except error:
        return
    raise AssertionError(f"expected {error.__name__}")
