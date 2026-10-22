"""Initialize a correlated population through its registered model initializer."""


class Initializer:
    def __init__(self, config, registry):
        self.config = config
        self.registry = registry
        self.state = None

    def initialize(self, initial=None):
        node = self.config["state"]["initialization_function"]
        function = self.registry.get(node["generator"], "initialization")
        self.state = function(self.config, initial=initial)
        return self.state
