"""Plant dynamics stage; preserve v0.1.12 draw and execution order."""

from dataclasses import dataclass

import numpy as np

from evo_sim.core.runner import ComputationalLimit
from evo_sim.models.plant_insect.biology import spatial
from evo_sim.models.plant_insect.biology.agents import CommunityState
from evo_sim.models.plant_insect.biology.plants.demography import (
    leaf_recovery_rate,
    plant_resource_update,
)

from .adult_reproduction import sample_parents
from .larval_development import LarvalDevelopment
from .offspring_genetics import OffspringGenetics


@dataclass
class PlantDynamics:
    """Values produced by this stage, before transition commit."""

    post: np.ndarray
    defoliation: np.ndarray
    S: np.ndarray
    recovered: np.ndarray
    recovery_rate: np.ndarray
    pcounts: np.ndarray
    pcandidates: list[np.ndarray]
    prows: list[dict]
    psurvived: np.ndarray
    npropagules: int
    propagule_mothers: np.ndarray
    proposal_cells: np.ndarray
    winners: np.ndarray
    proposal_sites: np.ndarray


def update_plants(
    runtime,
    current: CommunityState,
    t: int,
    *,
    offspring: OffspringGenetics,
    larvae: LarvalDevelopment,
) -> PlantDynamics:
    """Update plant resources/survival, then propose and establish propagules."""
    c = runtime.c
    plants = current.plants
    runtime.stage = "plant_survival_and_reproduction"
    post, defoliation, S, recovered = plant_resource_update(
        larvae.capacity,
        larvae.e["X"],
        larvae.maximum,
        c.plant_demography,
        plants.chemical_traits[:, 2],
    )
    pcounts, pcandidates, prows, psurvived = sample_parents(
        plants, "plant", S, larvae.e["fp"] * recovered / larvae.maximum, runtime, t
    )
    npropagules = sum(map(int, pcounts))
    if npropagules + len(offspring.hids) > c.max_individuals_safety:
        raise ComputationalLimit(
            "propagule and larval allocation exceeds safety; no transition committed"
        )
    propagule_mothers = np.repeat(np.arange(len(plants)), pcounts)
    proposal_cells = spatial.disperse(
        plants.cells[propagule_mothers],
        c.spatial["plant_dispersal_radius"],
        c.spatial["grid_shape"],
        runtime.spatial_rng[0],
    )
    winners, proposal_sites = spatial.establish_propagules(
        proposal_cells, plants.microsite[psurvived], c, runtime.rng[9][0][0]
    )
    nplants = len(winners)
    retained = int(psurvived.sum())
    if retained + nplants + int(larvae.survived.sum()) > c.max_individuals_safety:
        raise ComputationalLimit(
            "retained and recruited population exceeds safety; no transition committed"
        )
    snapshot_values = 2 * (
        len(offspring.genetic.loci["plant"]) * (retained + nplants)
        + len(offspring.genetic.loci["insect"]) * int(larvae.survived.sum())
    )
    offspring.genetic.preflight(
        [nplants, 0],
        runtime.stored_genotype_values
        + offspring.hvalues["insect_genotypes"].size
        + snapshot_values,
    )
    return PlantDynamics(
        post=post,
        defoliation=defoliation,
        S=S,
        recovered=recovered,
        recovery_rate=leaf_recovery_rate(plants.chemical_traits[:, 2], c.plant_demography),
        pcounts=pcounts,
        pcandidates=pcandidates,
        prows=prows,
        psurvived=psurvived,
        npropagules=npropagules,
        propagule_mothers=propagule_mothers,
        proposal_cells=proposal_cells,
        winners=winners,
        proposal_sites=proposal_sites,
    )
