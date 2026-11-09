"""Adult reproduction stage; preserve v0.1.12 draw and execution order."""

import json
from dataclasses import dataclass

import numpy as np

from evo_sim.core.runner import ComputationalLimit
from evo_sim.models.plant_insect.biology import traits as modules
from evo_sim.models.plant_insect.biology.agents import CommunityState
from evo_sim.models.plant_insect.biology.insects.demography import (
    maternal_moments,
    mating_candidates,
)


def sample_parents(pop, kind, S, F, runtime, generation):
    """Plants survive this window; current adult insects have already survived it."""
    k = int(kind == "insect")
    c = runtime.c
    maternal = np.ones(len(pop), bool) if k == 0 else pop.sex == 0
    donor = np.ones(len(pop), bool) if k == 0 else pop.sex == 1
    if k == 0:
        indices = np.arange(len(pop))
        potential = [indices[indices != i] for i in indices]
    else:
        potential = mating_candidates(
            pop.species,
            pop.cells,
            donor,
            maternal,
            c.sexual_reproduction["insect_mate_radius"],
            c.spatial["grid_shape"],
        )
    mu = (runtime.lp if k == 0 else runtime.li)[pop.species] * F
    if not np.isfinite(mu).all():
        raise ComputationalLimit("nonfinite maternal offspring mean")
    q, expected, variance = maternal_moments(S, mu, potential, maternal)
    survived = np.ones(len(pop), bool)
    if k == 0:
        for s, rng in enumerate(runtime.rng[2][0]):
            mask = pop.species == s
            survived[mask] = rng.random(mask.sum()) < S[mask]
    candidates = [pool[survived[pool]] for pool in potential]
    eligible = maternal & survived & np.array([len(pool) > 0 for pool in candidates], bool)
    counts = np.zeros(len(pop), dtype=np.int64)
    for s, rng in enumerate(runtime.rng[3][k]):
        mask = eligible & (pop.species == s)
        try:
            counts[mask] = rng.poisson(mu[mask])
        except ValueError as error:
            raise ComputationalLimit("offspring mean exceeds Poisson sampler capacity") from error
    rows = []
    for i in range(len(pop)):
        rows.append(
            dict(
                generation=generation,
                reproduction_round=generation,
                kind=kind,
                species=f"{'P' if k == 0 else 'H'}{pop.species[i] + 1}",
                species_id=int(pop.species[i]),
                individual_id=int(pop.ids[i]),
                sex_or_reproductive_role="hermaphrodite"
                if k == 0
                else ("female" if pop.sex[i] == 0 else "male"),
                is_maternal_role=bool(maternal[i]),
                S_i=float(S[i]),
                F_i=float(F[i]),
                W_i=float(S[i] * F[i]),
                lambda_s=float((runtime.lp if k == 0 else runtime.li)[pop.species[i]]),
                conditional_mu=float(mu[i]),
                mate_probability=float(q[i]),
                expected_maternal_offspring=float(expected[i]) if k else None,
                maternal_offspring_variance=float(variance[i]) if k else None,
                expected_propagules=float(expected[i]) if not k else None,
                propagule_variance=float(variance[i]) if not k else None,
                propagule_count=int(counts[i]) if not k else None,
                reproduction_active=True,
                survived_to_mating_stage=bool(survived[i]),
                eligible_as_mother=bool(maternal[i] and survived[i]),
                eligible_as_donor=bool(donor[i] and survived[i]),
                candidate_mate_count=len(candidates[i]),
                has_mate=bool(len(candidates[i])),
                potential_egg_count=int(counts[i]) if k else None,
                offspring_count_as_mother=0,
                offspring_count_sired=0,
                sampling_applied=True,
                contribution_committed=False,
                reason="not_maternal_role"
                if not maternal[i]
                else "did_not_survive"
                if not survived[i]
                else "no_compatible_mate"
                if not len(candidates[i])
                else "eligible",
            )
        )
    for i, row in enumerate(rows):
        pool = candidates[i]
        probabilities = np.full(len(pool), 1 / len(pool)) if len(pool) else np.empty(0)
        row.update(
            mating_row=int(pop.cells[i, 0]),
            mating_col=int(pop.cells[i, 1]),
            candidate_mate_ids=json.dumps(pop.ids[pool].tolist()),
            candidate_mate_probabilities=json.dumps(probabilities.tolist()),
        )
    return counts, candidates, rows, survived


def draw_fathers(pop, mothers, candidates, runtime, k):
    return np.array(
        [runtime.rng[5][k][pop.species[m]].choice(candidates[m]) for m in mothers], dtype=int
    )


@dataclass
class AdultReproduction:
    """Values produced by this stage, before transition commit."""

    counts: np.ndarray
    candidates: list[np.ndarray]
    hrows: list[dict]
    potential_count: int


def prepare_adults(runtime, current: CommunityState, t: int) -> AdultReproduction:
    """Sample adult maternal egg counts without repeating past larval survival."""
    c = runtime.c
    insects = current.insects
    runtime.stage = "adult_mating_and_potential_eggs"
    fi, processing_cost, sensory_cost = modules.maternal_fecundity(insects.module_traits, c)
    counts, candidates, hrows, _ = sample_parents(
        insects, "insect", np.ones(len(insects)), fi, runtime, t
    )
    for i, row in enumerate(hrows):
        row.update(
            processing_expression_cost=float(processing_cost[i]),
            sensory_expression_cost=float(sensory_cost[i]),
        )
    potential_count = sum(map(int, counts))
    if potential_count > c.max_individuals_safety:
        raise ComputationalLimit(
            f"{potential_count} potential eggs exceed safety; no transition committed"
        )
    return AdultReproduction(
        counts=counts,
        candidates=candidates,
        hrows=hrows,
        potential_count=potential_count,
    )
