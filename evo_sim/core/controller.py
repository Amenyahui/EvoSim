"""Build configured components and apply their declared state updates."""

from .helpers import get_by_path, set_by_path


class Controller:
    def __init__(self, config, registry):
        self.config = config
        self.functions = {}
        for step, definition in config["substeps"].items():
            self.functions[step] = {}
            for category in ("observation", "policy", "transition"):
                groups = definition.get(category) or {}
                if category == "transition":
                    groups = {"shared": groups}
                self.functions[step][category] = {}
                for agent, functions in groups.items():
                    self.functions[step][category][agent] = [
                        registry.get(node["generator"], category)(
                            config,
                            node.get("input_variables", {}),
                            node.get("output_variables", []),
                            node.get("arguments"),
                        )
                        for node in (functions or {}).values()
                    ]

    def observe(self, state, agent_type):
        output = {}
        for function in self.functions[state["current_substep"]]["observation"].get(agent_type, []):
            output.update(function(state))
        return output

    def act(self, state, observation, agent_type):
        output = {}
        for function in self.functions[state["current_substep"]]["policy"].get(agent_type, []):
            output.update(function(state, observation))
        return output

    def progress(self, state, action):
        for function in self.functions[state["current_substep"]]["transition"].get("shared", []):
            values = function(state, action)
            if not set(values).issubset(function.output_variables):
                raise ValueError("Transition returned undeclared outputs")
            # Check every destination before committing any returned state values.
            for name in values:
                get_by_path(state, function.input_variables[name])
            # Failed reproduction leaves the current generation's agents intact.
            for name, value in values.items():
                set_by_path(state, function.input_variables[name], value)
        return state
