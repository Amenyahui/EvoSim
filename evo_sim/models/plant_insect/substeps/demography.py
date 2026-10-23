import numpy as np

from evo_sim.core.substep import SubstepAction, SubstepTransition

from ..helpers.demography import ComputationalLimit, maternal_moments, mating_candidates
from ..helpers.runtime import community_state


class SampleMaternallyProducedOffspring(SubstepAction):
    def forward(self, state, observation):
        runtime = self.input(state, "runtime")
        community = community_state(self.input(state, "agents"))
        runtime.stage = "demography"
        k = self.fixed_args["taxon"]
        kind = "plant" if k == 0 else "insect"
        pop = community.plants if k == 0 else community.insects
        labels, ids = pop.species, pop.ids
        S = runtime.e["S_" + kind]
        F = runtime.e["fp" if k == 0 else "fi"]
        lam = (runtime.lp if k == 0 else runtime.li)[labels]
        phi = runtime.plant_phi if k == 0 else np.ones(len(labels))
        mu = lam * F * phi
        if not np.isfinite(mu).all():
            raise ComputationalLimit("nonfinite offspring mean")
        valid = np.ones(len(labels), bool) if k == 0 else runtime.e["host"] >= 0
        cells = pop.cells.copy()
        if k == 1:
            cells[valid] = community.plants.cells[runtime.e["host"][valid]]
        mother = valid & (np.ones(len(labels), bool) if k == 0 else pop.sex == 0)
        donor = valid & (F > 0) & (np.ones(len(labels), bool) if k == 0 else pop.sex == 1)
        radius = runtime.c.sexual_reproduction["pollen_radius" if k == 0 else "insect_mate_radius"]
        potential = mating_candidates(
            labels, cells, donor, mother, radius, runtime.c.spatial["grid_shape"]
        )
        q, expected, variance = maternal_moments(S, mu, potential, mother)
        survived = np.zeros(len(labels), bool)
        births = np.zeros(len(labels), dtype=np.int64)
        candidates = [np.empty(0, dtype=int) for _ in labels]
        if runtime.update:
            for s, rng in enumerate(runtime.rng[2][k]):
                mask = labels == s
                survived[mask] = rng.random(mask.sum()) < S[mask]
            candidates = [pool[survived[pool]] for pool in potential]
            eligible = mother & survived & np.array([len(pool) > 0 for pool in candidates], bool)
            for s, rng in enumerate(runtime.rng[3][k]):
                mask = eligible & (labels == s)
                try:
                    births[mask] = rng.poisson(mu[mask])
                except ValueError as error:
                    raise ComputationalLimit(
                        "offspring mean exceeds Poisson sampler capacity"
                    ) from error
        rows = []
        for i in range(len(labels)):
            rows.append(
                dict(
                    generation=state["current_step"],
                    kind=kind,
                    species=f"{'P' if k == 0 else 'H'}{labels[i] + 1}",
                    species_id=int(labels[i]),
                    individual_id=int(ids[i]),
                    sex_or_reproductive_role="seed_and_pollen"
                    if k == 0
                    else ("female" if pop.sex[i] == 0 else "male"),
                    S_i=S[i],
                    F_i=F[i],
                    W_i=S[i] * F[i],
                    phi_i=phi[i],
                    lambda_s=lam[i],
                    conditional_mu=mu[i],
                    mate_probability=q[i],
                    expected_maternal_offspring=expected[i],
                    maternal_offspring_variance=variance[i],
                    survived_to_mating_stage=bool(survived[i]) if runtime.update else None,
                    eligible_as_mother=bool(mother[i] and survived[i]) if runtime.update else None,
                    eligible_as_donor=bool(donor[i] and survived[i]) if runtime.update else None,
                    candidate_mate_count=len(candidates[i]) if runtime.update else None,
                    has_mate=bool(len(candidates[i])) if runtime.update else None,
                    offspring_count_as_mother=int(births[i]) if runtime.update else None,
                    offspring_count_sired=0 if runtime.update else None,
                    sampling_applied=runtime.update,
                    reason=(
                        "not_maternal_role"
                        if not mother[i]
                        else "did_not_survive"
                        if not survived[i]
                        else "no_compatible_mate"
                        if not len(candidates[i])
                        else "eligible"
                    )
                    if runtime.update
                    else "observation_only",
                )
            )
        runtime.result["parents"].extend(rows)
        return {
            self.output_variables[0]: dict(
                survived_to_reproduce=survived,
                offspring_count=births,
                candidates=candidates,
                cells=cells,
                parent_rows=rows,
            )
        }


class CollectOffspring(SubstepTransition):
    def forward(self, state, action):
        runtime = self.input(state, "runtime")
        runtime.sampled = [action[agent]["offspring"] for agent in ("plants", "insects")]
        return {}
