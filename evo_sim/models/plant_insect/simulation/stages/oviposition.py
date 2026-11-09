"""Oviposition stage; preserve v0.1.12 draw and execution order."""

from dataclasses import dataclass

import numpy as np

from evo_sim.models.plant_insect.biology import spatial
from evo_sim.models.plant_insect.biology.agents import CommunityState
from evo_sim.models.plant_insect.biology.insects import host_selection
from evo_sim.models.plant_insect.biology.plants import chemistry

from .adult_reproduction import AdultReproduction, draw_fathers


@dataclass
class Oviposition:
    """Values produced by this stage, before transition commit."""

    egg_mothers: np.ndarray
    V: np.ndarray
    choice: dict
    accepted: np.ndarray
    mothers: np.ndarray
    host: np.ndarray
    fathers: np.ndarray


def oviposit(runtime, current: CommunityState, *, adults: AdultReproduction) -> Oviposition:
    """Choose a host per potential egg, then draw accepted eggs' fathers."""
    c = runtime.c
    pl = current.plants.species
    hl = current.insects.species
    insects = current.insects
    plants = current.plants
    runtime.stage = "maternal_oviposition"
    egg_mothers = np.repeat(np.arange(len(insects)), adults.counts)
    V = chemistry.plant_voc(plants.defense, plants.voc_trait, c, runtime.B, plants.chemical_traits)
    V_cue = V[runtime.rng[8][0][0].permutation(len(V))] if c.shuffle_voc else V.copy()
    choice = host_selection.choose_hosts(
        insects.module_traits[egg_mothers],
        pl,
        hl[egg_mothers],
        V_cue,
        c,
        runtime.w,
        runtime.rng[1][0][0],
        runtime.rng[1][1],
        candidate_pools=spatial.candidate_pools(plants.cells, insects.cells[egg_mothers], c),
        PSM=plants.defense,
    )
    accepted = choice["host"] >= 0
    mothers, host = egg_mothers[accepted], choice["host"][accepted]
    fathers = draw_fathers(insects, mothers, adults.candidates, runtime, 1)
    return Oviposition(
        egg_mothers=egg_mothers,
        V=V,
        choice=choice,
        accepted=accepted,
        mothers=mothers,
        host=host,
        fathers=fathers,
    )
