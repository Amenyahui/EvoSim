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
    c = Config(**values)
    # Inheritance fixtures use a pulse each round and ample local microsites;
    # perennial cadence and finite establishment have separate explicit checks.
    c.plant_demography.update(
        reproduction_interval=1, microsites_per_cell=max(1, sum(c.plant_sizes)) * 4
    )
    return c


def raises(function, error=ValueError):
    try:
        function()
    except error:
        return
    errors = error if isinstance(error, tuple) else (error,)
    raise AssertionError("expected " + " or ".join(item.__name__ for item in errors))
