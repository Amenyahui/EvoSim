"""Explicit per-model registration, following AgentTorch's simulator.get_registry()."""


class Registry:
    def __init__(self):
        self.helpers = {
            key: {} for key in ("initialization", "observation", "policy", "transition")
        }

    def register(self, obj_source, name, key):
        if key not in self.helpers:
            raise ValueError(f"Unknown registration category: {key}")
        if name in self.helpers[key]:
            raise ValueError(f"Duplicate {key}: {name}")
        self.helpers[key][name] = obj_source

    def get(self, name, key):
        try:
            return self.helpers[key][name]
        except KeyError as error:
            raise ValueError(f"Unregistered {key}: {name}") from error
