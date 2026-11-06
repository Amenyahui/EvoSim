"""Plant–herbivore coevolution model."""

from evo_sim.models.plant_insect.config.schema import Config

from .simulator import get_runner, simulate

__all__ = ["Config", "get_runner", "simulate"]
