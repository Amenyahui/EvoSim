"""Configuration-driven Observation -> Policy -> Transition scheduler."""

from copy import deepcopy

from .controller import Controller
from .initializer import Initializer


class Runner:
    def __init__(self, config, registry):
        self.config = deepcopy(config)
        self.registry = registry
        metadata = self.config["simulation_metadata"]
        if metadata.get("device", "cpu") != "cpu" or metadata.get("calibration", False):
            raise ValueError("This runtime supports NumPy CPU simulation without calibration")
        if metadata["num_substeps_per_step"] != len(self.config["substeps"]):
            raise ValueError("num_substeps_per_step must match configured substeps")
        if (
            type(metadata["num_steps_per_episode"]) is not int
            or metadata["num_steps_per_episode"] < 0
        ):
            raise ValueError("num_steps_per_episode must be a nonnegative integer")
        agents = self.config["state"]["agents"]
        for definition in self.config["substeps"].values():
            if any(agent not in agents for agent in definition["active_agents"]):
                raise ValueError("Substep names an undeclared agent type")
        self.initializer = Initializer(self.config, registry)
        self.controller = Controller(self.config, registry)
        self.state = None

    def init(self, initial=None):
        self.state = self.initializer.initialize(initial)
        return self.state

    def reset(self):
        return self.init()

    def step(self, num_steps=None):
        if self.state is None:
            raise RuntimeError("Call runner.init() before step()")
        if num_steps is None:
            num_steps = self.config["simulation_metadata"]["num_steps_per_episode"]
        if type(num_steps) is not int or num_steps < 0:
            raise ValueError("num_steps must be a nonnegative integer")
        for _ in range(num_steps):
            if self.state["terminated"]:
                break
            for substep, definition in self.config["substeps"].items():
                self.state["current_substep"] = substep
                action = {}
                for agent in definition["active_agents"]:
                    observation = self.controller.observe(self.state, agent)
                    action[agent] = self.controller.act(self.state, observation, agent)
                self.state = self.controller.progress(self.state, action)
                if self.state["terminated"]:
                    break
            if not self.state["terminated"]:
                self.state["current_step"] += 1
        return self.state

    def forward(self):
        self.init()
        return self.step()
