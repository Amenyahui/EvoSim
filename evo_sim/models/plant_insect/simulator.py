"""Plant–insect model registration and episode lifecycle."""

from copy import deepcopy
from pathlib import Path

from evo_sim.core.registry import Registry
from evo_sim.core.runner import Runner

from .helpers.config import config_payload
from .helpers.demography import ComputationalLimit
from .helpers.inputs import read_config
from .helpers.runtime import fail, finalize, initialize_community
from .substeps.demography import CollectOffspring, SampleMaternallyProducedOffspring
from .substeps.ecology import ChooseHosts, ObserveCues, SettleFeeding
from .substeps.lifecycle import AdvanceGeneration
from .substeps.reproduction import CrossFertilizeAndDisperse

MODEL_DIR = Path(__file__).resolve().parent


def get_registry():
    registry = Registry()
    for function, name, category in (
        (initialize_community, "initialize_community", "initialization"),
        (AdvanceGeneration, "advance_generation", "transition"),
        (ObserveCues, "observe_cues", "observation"),
        (ChooseHosts, "choose_hosts", "policy"),
        (SettleFeeding, "settle_feeding", "transition"),
        (SampleMaternallyProducedOffspring, "sample_maternal_offspring", "policy"),
        (CollectOffspring, "collect_offspring", "transition"),
        (CrossFertilizeAndDisperse, "cross_fertilize_and_disperse", "transition"),
    ):
        registry.register(function, name, category)
    return registry


class PlantInsectRunner(Runner):
    def init(self, initial=None):
        super().init(initial)
        self.result = self.state["environment"]["runtime"].result
        self.state_trajectory = self.result["snapshots"]
        if self.state["terminated"]:
            finalize(self.state)
        return self.state

    def step(self, num_steps=None):
        if self.state is None:
            raise RuntimeError("Call runner.init() before step()")
        if num_steps is not None and (type(num_steps) is not int or num_steps < 0):
            raise ValueError("num_steps must be a nonnegative integer")
        total = self.config["simulation_metadata"]["num_steps_per_episode"]
        remaining = max(0, total - self.state["current_step"])
        try:
            super().step(remaining if num_steps is None else min(num_steps, remaining))
            # The terminal generation records adults/plants and reference observations only.
            if self.state["current_step"] >= total and not self.state["terminated"]:
                super().step(1)
        except (ComputationalLimit, MemoryError, ArithmeticError, ValueError) as error:
            fail(self.state, error)
        if self.state["terminated"]:
            finalize(self.state)
        return self.state


def get_runner(config_path=None, *, config=None):
    if config is not None and config_path is not None:
        raise ValueError("Provide config or config_path, not both")
    if config is None:
        config = read_config(config_path or MODEL_DIR / "configs/sexual_qtl_development.yaml")
    if config["parameters"].get("model_version") == "0.1.6":
        steps = list(config["substeps"].values())
        transitions = list(steps[0].get("transition", {}).values()) if len(steps) == 1 else []
        if (
            len(steps) != 1
            or steps[0]["active_agents"]
            or steps[0].get("observation")
            or steps[0].get("policy")
            or len(transitions) != 1
            or transitions[0]["generator"] != "advance_generation"
        ):
            raise ValueError("v0.1.6 requires the adult/larval advance_generation lifecycle")
    return PlantInsectRunner(config, get_registry())


def simulate(c, initial=None, *, runtime_config=None):
    """Convenience API; every generation executes the configured substeps."""
    config = (
        deepcopy(runtime_config)
        if runtime_config is not None
        else read_config(MODEL_DIR / "config.yaml")
    )
    config["parameters"] = config_payload(c)
    config["parameters"].pop("generations", None)
    config["simulation_metadata"]["num_steps_per_episode"] = c.generations
    runner = get_runner(config=config)
    runner.init(initial)
    runner.step()
    return runner.result
