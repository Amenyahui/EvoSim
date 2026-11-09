"""Larval development stage; preserve v0.1.12 draw and execution order."""

from dataclasses import dataclass

import numpy as np

from evo_sim.models.plant_insect.biology import spatial
from evo_sim.models.plant_insect.biology.agents import CommunityState
from evo_sim.models.plant_insect.biology.insects import feeding
from evo_sim.models.plant_insect.records.contributions import contributions
from evo_sim.models.plant_insect.simulation.runtime import resource_vector

from .adult_reproduction import AdultReproduction
from .offspring_genetics import OffspringGenetics
from .oviposition import Oviposition


@dataclass
class LarvalDevelopment:
    """Values produced by this stage, before transition commit."""

    maximum: np.ndarray
    capacity: np.ndarray
    dmax: np.ndarray
    e: dict
    survived: np.ndarray
    sex: np.ndarray
    origins: np.ndarray
    adult_cells: np.ndarray


def develop_larvae(
    runtime,
    current: CommunityState,
    t: int,
    *,
    adults: AdultReproduction,
    oviposition: Oviposition,
    offspring: OffspringGenetics,
) -> LarvalDevelopment:
    """Allocate shared resources and draw larval survival, sex and dispersal."""
    c = runtime.c
    pl = current.plants.species
    plants = current.plants
    maximum = resource_vector(c.resource_capacity, c.plant_sizes, pl, "P")
    capacity = plants.resource
    dmax = resource_vector(c.max_feeding_demand, c.insect_sizes, offspring.child_hl, "H")
    runtime.stage = "larval_feeding_and_survival"
    e = feeding.settle_feeding(
        plants.defense,
        offspring.hvalues["I"],
        pl,
        offspring.child_hl,
        oviposition.V,
        oviposition.host,
        c,
        capacity,
        dmax,
        chemical=plants.chemical_traits,
        module_traits=offspring.hvalues["insect_module_traits"],
    )
    survived = np.zeros(len(offspring.hids), bool)
    for s, rng in enumerate(runtime.rng[2][1]):
        selected = offspring.child_hl == s
        survived[selected] = rng.random(selected.sum()) < e["survival_probability"][selected]
    sex = np.empty(len(offspring.hids), dtype=np.int8)
    for s, rng in enumerate(runtime.rng[7][1]):
        selected = offspring.child_hl == s
        sex[selected] = (
            rng.random(selected.sum())
            >= c.sexual_reproduction["insect_offspring_female_probability"]
        )
    origins = plants.cells[oviposition.host]
    adult_cells = spatial.disperse(
        origins[survived],
        c.spatial["insect_dispersal_radius"],
        c.spatial["grid_shape"],
        runtime.spatial_rng[1],
    )
    offspring.hvalues.update(
        insect_species=offspring.child_hl,
        insect_ids=offspring.hids,
        insect_cells=origins.copy(),
        insect_sex=sex,
        insect_birth_host_id=plants.ids[oviposition.host].copy(),
        insect_birth_host_species=pl[oviposition.host].copy(),
        insect_adult_cohort=np.full(len(offspring.hids), t + 1, dtype=int),
    )
    contributions(adults.hrows, oviposition.mothers, oviposition.fathers, survived)
    return LarvalDevelopment(
        maximum=maximum,
        capacity=capacity,
        dmax=dmax,
        e=e,
        survived=survived,
        sex=sex,
        origins=origins,
        adult_cells=adult_cells,
    )
