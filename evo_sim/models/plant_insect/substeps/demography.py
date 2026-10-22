import numpy as np

from evo_sim.core.substep import SubstepAction, SubstepTransition

from ..helpers.demography import draw_offspring
from ..helpers.runtime import community_state


class SampleOffspring(SubstepAction):
    def forward(self, state, observation):
        runtime = self.input(state, "runtime")
        community = community_state(self.input(state, "agents"))
        runtime.stage = "demography"
        k = self.fixed_args["taxon"]
        kind = "plant" if k == 0 else "insect"
        pl, hl = community.plants.species, community.insects.species
        labels = pl if k == 0 else hl
        S = runtime.e["S_plant"] if k == 0 else runtime.e["S_insect"]
        F = runtime.e["fp"] if k == 0 else runtime.e["fi"]
        lam = (runtime.lp if k == 0 else runtime.li)[labels]
        factors = runtime.plant_phi if k == 0 else np.ones(len(labels))
        values = {
            name: np.zeros(
                len(labels),
                dtype=bool
                if name == "survived_to_reproduce"
                else np.int64
                if name == "offspring_count"
                else float,
            )
            for name in (
                "survived_to_reproduce",
                "offspring_count",
                "conditional_offspring_mean_mu",
                "expected_offspring",
                "offspring_variance",
            )
        }
        with np.errstate(over="raise", invalid="raise"):
            values["conditional_offspring_mean_mu"] = lam * F * factors
            values["expected_offspring"] = S * values["conditional_offspring_mean_mu"]
            values["offspring_variance"] = (
                values["expected_offspring"]
                + S * (1 - S) * values["conditional_offspring_mean_mu"] ** 2
            )
        if runtime.update:
            for s in range(len(runtime.c.plant_sizes if k == 0 else runtime.c.insect_sizes)):
                mask = labels == s
                drawn = draw_offspring(
                    S[mask],
                    F[mask],
                    lam[mask],
                    factors[mask],
                    runtime.rng[2][k][s],
                    runtime.rng[3][k][s],
                )
                for name, value in drawn.items():
                    values[name][mask] = value
        ids = community.plants.ids if k == 0 else community.insects.ids
        prefix = "P" if k == 0 else "H"
        for i in range(len(labels)):
            runtime.result["parents"].append(
                dict(
                    generation=state["current_step"],
                    kind=kind,
                    species=f"{prefix}{labels[i] + 1}",
                    species_id=int(labels[i]),
                    individual_id=int(ids[i]),
                    parent_genotype_id=f"{prefix}:{ids[i]}",
                    S_i=S[i],
                    F_i=F[i],
                    W_i=S[i] * F[i],
                    phi_i=factors[i],
                    lambda_s=lam[i],
                    conditional_offspring_mean_mu=values["conditional_offspring_mean_mu"][i],
                    expected_offspring=values["expected_offspring"][i],
                    offspring_variance=values["offspring_variance"][i],
                    survived_to_reproduce=bool(values["survived_to_reproduce"][i])
                    if runtime.update
                    else None,
                    offspring_count=int(values["offspring_count"][i]) if runtime.update else None,
                    sampling_applied=runtime.update,
                    survival_stream=f"phase2/{kind}/species{labels[i]}",
                    offspring_stream=f"phase3/{kind}/species{labels[i]}",
                )
            )
        return {self.output_variables[0]: values}


class CollectOffspring(SubstepTransition):
    def forward(self, state, action):
        runtime = self.input(state, "runtime")
        runtime.sampled = [action[agent]["offspring"] for agent in ("plants", "insects")]
        return {}
