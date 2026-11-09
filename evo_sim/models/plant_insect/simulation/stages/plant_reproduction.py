"""Plant reproduction stage; preserve v0.1.12 draw and execution order."""

from dataclasses import dataclass

import numpy as np

from evo_sim.models.plant_insect.biology.agents import CommunityState
from evo_sim.models.plant_insect.records.contributions import contributions

from .adult_reproduction import draw_fathers
from .larval_development import LarvalDevelopment
from .offspring_genetics import OffspringGenetics
from .plant_dynamics import PlantDynamics


@dataclass
class PlantReproduction:
    """Values produced by this stage, before transition commit."""

    pmothers: np.ndarray
    pfathers: np.ndarray
    pvalues: dict
    pids: np.ndarray
    precords: list[dict]
    next_p: int
    pcells: np.ndarray


def reproduce_plants(
    runtime,
    current: CommunityState,
    t: int,
    *,
    offspring: OffspringGenetics,
    larvae: LarvalDevelopment,
    plant_update: PlantDynamics,
) -> PlantReproduction:
    """Draw fathers and inherit only established plant offspring through SLiM."""
    c = runtime.c
    pl = current.plants.species
    plants = current.plants
    pmothers = plant_update.propagule_mothers[plant_update.winners]
    pfathers = draw_fathers(plants, pmothers, plant_update.pcandidates, runtime, 0)
    runtime.stage = "plant_meiosis_mutation"
    pvalues, child_pl, pids, precords, next_p = offspring.genetic.fertilize(
        plants.genotypes, pl, plants.ids, pmothers, pfathers, "plant", runtime.next_ids[0], t + 1
    )
    pcells = plant_update.proposal_cells[plant_update.winners]
    pvalues.update(
        plant_species=child_pl,
        plant_ids=pids,
        plant_cells=pcells,
        plant_resource=c.plant_demography["newborn_resource_fraction"] * larvae.maximum[pmothers],
        plant_microsite=plant_update.proposal_sites[plant_update.winners],
    )
    contributions(plant_update.prows, pmothers, pfathers, np.ones(len(pids), bool))
    return PlantReproduction(
        pmothers=pmothers,
        pfathers=pfathers,
        pvalues=pvalues,
        pids=pids,
        precords=precords,
        next_p=next_p,
        pcells=pcells,
    )
