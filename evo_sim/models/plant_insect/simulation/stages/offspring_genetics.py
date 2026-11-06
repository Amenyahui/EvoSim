"""Offspring genetics stage; preserve v0.1.11 draw and execution order."""

import copy
from dataclasses import dataclass

import numpy as np

from evo_sim.models.plant_insect.biology.agents import CommunityState
from evo_sim.models.plant_insect.biology.genetics.qtl import GeneticState

from .oviposition import Oviposition


@dataclass
class OffspringGenetics:
    """Values produced by this stage, before transition commit."""

    genetic: GeneticState
    hvalues: dict
    child_hl: np.ndarray
    hids: np.ndarray
    hrecords: list[dict]
    next_h: int


def inherit_offspring(
    runtime, current: CommunityState, t: int, *, oviposition: Oviposition
) -> OffspringGenetics:
    """Inherit accepted larvae through SLiM using a disposable genetic state."""
    hl = current.insects.species
    insects = current.insects
    runtime.genetic.preflight([0, len(oviposition.mothers)], runtime.stored_genotype_values)
    genetic = copy.copy(runtime.genetic)
    genetic.alleles = runtime.genetic.alleles.copy()
    genetic.events = runtime.genetic.events.copy()
    genetic.transmissions = runtime.genetic.transmissions.copy()
    genetic.rng = copy.deepcopy(runtime.genetic.rng)
    runtime.stage = "insect_meiosis_mutation"
    hvalues, child_hl, hids, hrecords, next_h = genetic.fertilize(
        insects.genotypes,
        hl,
        insects.ids,
        oviposition.mothers,
        oviposition.fathers,
        "insect",
        runtime.next_ids[1],
        t + 1,
    )
    return OffspringGenetics(
        genetic=genetic,
        hvalues=hvalues,
        child_hl=child_hl,
        hids=hids,
        hrecords=hrecords,
        next_h=next_h,
    )
