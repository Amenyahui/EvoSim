"""AgentTorch-style forward contracts, using the model's NumPy state arrays."""

from abc import ABC, abstractmethod

from .helpers import get_by_path


class Substep(ABC):
    def __init__(self, config, input_variables, output_variables, arguments=None):
        self.config = config
        self.input_variables = input_variables
        self.output_variables = output_variables
        arguments = arguments or {}
        if arguments.get("learnable"):
            raise ValueError("This NumPy model does not implement learnable parameters")
        self.fixed_args = arguments.get("fixed", {})

    def input(self, state, name):
        return get_by_path(state, self.input_variables[name])

    def __call__(self, *args, **kwargs):
        return self.forward(*args, **kwargs)


class SubstepObservation(Substep):
    @abstractmethod
    def forward(self, state):
        raise NotImplementedError


class SubstepAction(Substep):
    @abstractmethod
    def forward(self, state, observation):
        raise NotImplementedError


class SubstepTransition(Substep):
    @abstractmethod
    def forward(self, state, action):
        raise NotImplementedError
