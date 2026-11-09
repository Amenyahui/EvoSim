"""Recruitment stage; preserve v0.1.12 draw and execution order."""

from dataclasses import dataclass

import numpy as np

from evo_sim.models.plant_insect.biology.agents import CommunityState

from .larval_development import LarvalDevelopment
from .offspring_genetics import OffspringGenetics
from .plant_dynamics import PlantDynamics
from .plant_reproduction import PlantReproduction


@dataclass
class Recruitment:
    """Values produced by this stage, before transition commit."""

    birth_state: CommunityState
    recruitment_state: CommunityState
    next_state: CommunityState


def recruit(
    current: CommunityState,
    *,
    offspring: OffspringGenetics,
    larvae: LarvalDevelopment,
    plant_update: PlantDynamics,
    plant_births: PlantReproduction,
) -> Recruitment:
    """Assemble retained plants and recruited insects without publishing state."""
    birth_state = CommunityState.from_arrays(dict(plant_births.pvalues, **offspring.hvalues))
    recruited_hvalues = {key: value[larvae.survived] for key, value in offspring.hvalues.items()}
    recruited_hvalues["insect_cells"] = larvae.adult_cells
    recruitment_state = CommunityState.from_arrays(dict(plant_births.pvalues, **recruited_hvalues))
    old_values = current.arrays()
    next_values = dict(recruited_hvalues)
    for key, value in plant_births.pvalues.items():
        preserved = (
            plant_update.recovered[plant_update.psurvived]
            if key == "plant_resource"
            else old_values[key][plant_update.psurvived]
        )
        next_values[key] = np.concatenate((preserved, value))
    next_state = CommunityState.from_arrays(next_values)
    next_state.validate()
    return Recruitment(
        birth_state=birth_state,
        recruitment_state=recruitment_state,
        next_state=next_state,
    )
