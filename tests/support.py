import copy

from evo_sim.models.plant_insect import Config, get_runner


def config(**changes):
    values = copy.deepcopy(get_runner().config["parameters"])
    values.update(generations=2, plant_sizes=[20, 20], insect_sizes=[6])
    values.update(changes)
    return Config(**values)


def raises(function, error=ValueError):
    try:
        function()
    except error:
        return
    raise AssertionError(f"expected {error.__name__}")
