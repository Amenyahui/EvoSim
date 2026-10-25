import copy
import json

import numpy as np

from evo_sim.core.substep import SubstepAction, SubstepObservation, SubstepTransition

from ..helpers import ecology as shared
from ..helpers import spatial, voc
from ..helpers.demography import density_phi
from ..helpers.observations import standardized
from ..helpers.runtime import community_state, mean_or_na, resource_vector


class ObserveCues(SubstepObservation):
    def forward(self, state):
        runtime = self.input(state, "runtime")
        community = community_state(self.input(state, "agents"))
        runtime.stage = "ecology"
        snapshot = community.snapshot(state["current_step"])
        runtime.result["snapshots"].append(snapshot)
        runtime.stored_genotype_values += (
            community.plants.genotypes.size + community.insects.genotypes.size
        )
        runtime.result["genetic_summary"].extend(
            runtime.genetic.summary(
                community,
                state["current_step"],
                (len(runtime.c.plant_sizes), len(runtime.c.insect_sizes)),
            )
        )
        runtime.result["haplotype_summary"].extend(
            runtime.genetic.haplotype_summary(
                community,
                state["current_step"],
                (len(runtime.c.plant_sizes), len(runtime.c.insect_sizes)),
            )
        )
        D = community.plants.defense
        V_emit = voc.map_voc(
            D, community.plants.voc_trait, runtime.B, runtime.c.rho_voc, runtime.c.voc_bounds
        )
        V_cue = (
            V_emit[runtime.rng[1][0][0].permutation(len(D))]
            if runtime.c.shuffle_voc
            else V_emit.copy()
        )
        return {self.output_variables[0]: V_emit, self.output_variables[1]: V_cue}


class ChooseHosts(SubstepAction):
    def forward(self, state, observation):
        runtime = self.input(state, "runtime")
        community = community_state(self.input(state, "agents"))
        options = {}
        if runtime.c.spatial is not None:
            options = dict(
                candidate_pools=spatial.candidate_pools(
                    community.plants.cells, community.insects.cells, runtime.c
                ),
            )
        choice = shared.choose_hosts(
            community.insects.preference,
            community.plants.species,
            community.insects.species,
            observation["V_cue"],
            runtime.c,
            runtime.w,
            runtime.rng[1][0][0],
            runtime.rng[1][1],
            **options,
        )
        return {self.output_variables[0]: dict(choice=choice, **observation)}


class SettleFeeding(SubstepTransition):
    def forward(self, state, action):
        runtime = self.input(state, "runtime")
        community = community_state(self.input(state, "agents"))
        chosen = action["insects"]["host_choice"]
        c = runtime.c
        pl, hl = community.plants.species, community.insects.species
        capacity = resource_vector(c.resource_capacity, c.plant_sizes, pl, "P")
        dmax = resource_vector(c.max_feeding_demand, c.insect_sizes, hl, "H")
        e = shared.settle_feeding(
            community.plants.defense,
            community.insects.resistance,
            pl,
            hl,
            chosen["V_emit"],
            chosen["choice"]["host"],
            c,
            capacity,
            dmax,
        )
        e["wi"] = e["g"] * e["fi"]
        e.update(
            chosen["choice"],
            V_emit=chosen["V_emit"],
            V_cue=chosen["V_cue"],
            R=capacity,
            dmax=dmax,
            S_plant=c.sP0 + (1 - c.sP0) * (1 - e["h"]),
            S_insect=e["g"],
        )
        runtime.e = e
        snapshot = runtime.result["snapshots"][-1]
        snapshot.update(
            R=runtime.e["R"].copy(),
            dmax=runtime.e["dmax"].copy(),
            V_emit=runtime.e["V_emit"].copy(),
            V_cue=runtime.e["V_cue"].copy(),
        )
        pl, hl = (community.plants.species, community.insects.species)
        runtime.counts = [
            np.bincount(labels, minlength=len(sizes))
            for labels, sizes in ((pl, runtime.c.plant_sizes), (hl, runtime.c.insect_sizes))
        ]
        runtime.stage = "plant_density"
        runtime.phi = density_phi(runtime.alpha, runtime.counts[0])
        runtime.plant_phi = runtime.phi[pl]
        if runtime.c.spatial is not None:
            runtime.plant_phi = spatial.local_competition(
                community.plants.cells, pl, runtime.alpha, runtime.c
            )
            runtime.phi = np.array(
                [mean_or_na(runtime.plant_phi[pl == s]) for s in range(len(runtime.c.plant_sizes))]
            )
        runtime.update = state["current_step"] < runtime.c.generations and len(pl) + len(hl) > 0
        runtime.reproduction_checkpoint = (
            (copy.deepcopy(runtime.rng), copy.deepcopy(runtime.spatial_rng))
            if runtime.update
            else None
        )
        runtime.start = len(runtime.result["populations"])
        for kind, labels, sizes, fitness, gain, damage, amount, phis in [
            (
                "plant",
                pl,
                runtime.c.plant_sizes,
                runtime.e["wp"],
                None,
                runtime.e["h"],
                runtime.e["X"],
                runtime.phi,
            ),
            (
                "insect",
                hl,
                runtime.c.insect_sizes,
                runtime.e["wi"],
                runtime.e["g"],
                None,
                runtime.e["x"],
                np.ones(len(runtime.c.insect_sizes)),
            ),
        ]:
            prefix = "P" if kind == "plant" else "H"
            for s in range(len(sizes)):
                mask = labels == s
                name = f"{prefix}{s + 1}"
                if not mask.any() and name not in runtime.extinct:
                    runtime.extinct[name] = state["current_step"]
                runtime.result["populations"].append(
                    dict(
                        generation=state["current_step"],
                        kind=kind,
                        species=name,
                        N_start=int(mask.sum()),
                        N_female=int(((community.insects.sex == 0) & mask).sum())
                        if kind == "insect"
                        else None,
                        N_male=int(((community.insects.sex == 1) & mask).sum())
                        if kind == "insect"
                        else None,
                        N_survived_to_reproduce=None,
                        N_offspring=None,
                        N_next=None,
                        update_applied=False,
                        extinct_at_generation=runtime.extinct.get(name),
                        zero_fitness_count=int((fitness[mask] == 0).sum()),
                        mean_g=mean_or_na(gain[mask]) if gain is not None else np.nan,
                        mean_h=mean_or_na(damage[mask]) if damage is not None else np.nan,
                        total_resource_transfer=float(amount[mask].sum()),
                        density_phi=phis[s],
                        population_status="extinct" if not mask.any() else "extant",
                    )
                )
        for i, p in enumerate(runtime.e["host"]):
            row = dict(
                generation=state["current_step"],
                seed=runtime.c.seed,
                insect_id=int(community.insects.ids[i]),
                insect_row=i,
                species=f"H{hl[i] + 1}",
                species_id=int(hl[i]),
                host_row=int(p),
                host_id=int(community.plants.ids[p]) if p >= 0 else -1,
                plant_species=f"P{pl[p] + 1}" if p >= 0 else "",
                accepted=bool(p >= 0),
                reason=runtime.e["reason"][i],
                candidate_count=int(runtime.e["candidate_count"][i]),
                choice_order=int(runtime.e["choice_order"][i]),
                q=runtime.e["q"][i],
                d=runtime.e["d"][i],
                x=runtime.e["x"][i],
                g=runtime.e["g"][i],
                dmax=runtime.e["dmax"][i],
                R=runtime.e["R"][p] if p >= 0 else 0.0,
                gamma=runtime.e["gamma"][p] if p >= 0 else 1.0,
                host_h=runtime.e["h"][p] if p >= 0 else 0.0,
                I_sum=community.insects.resistance[i].sum(),
                W_insect=runtime.e["wi"][i],
            )
            if runtime.observe:
                for pidx in range(len(runtime.c.plant_sizes)):
                    row[f"P{pidx + 1}_candidate_share"] = runtime.e["candidate_share"][i, pidx]
                    row[f"P{pidx + 1}_conditional_probability"] = runtime.e["conditional_share"][
                        i, pidx
                    ]
            if runtime.c.spatial is not None:
                row.update(
                    search_row=int(community.insects.cells[i, 0]),
                    search_col=int(community.insects.cells[i, 1]),
                    host_cell_row=int(community.plants.cells[p, 0]) if p >= 0 else -1,
                    host_cell_col=int(community.plants.cells[p, 1]) if p >= 0 else -1,
                    host_assignment="local_choice",
                )
                row["candidate_plant_ids"] = json.dumps(
                    community.plants.ids[runtime.e["candidate_rows"][i]].tolist()
                )
            runtime.result["events"].append(row)
        for p in range(len(pl)):
            runtime.result["plants"].append(
                dict(
                    generation=state["current_step"],
                    plant_id=int(community.plants.ids[p]),
                    plant_row=p,
                    species=f"P{pl[p] + 1}",
                    species_id=int(pl[p]),
                    R=runtime.e["R"][p],
                    T=runtime.e["total_demand"][p],
                    gamma=runtime.e["gamma"][p],
                    X=runtime.e["X"][p],
                    h=runtime.e["h"][p],
                    herbivore_load=int(runtime.e["herbivore_load"][p]),
                    D_sum=community.plants.defense[p].sum(),
                    V_emit_mean=runtime.e["V_emit"][p].mean(),
                    defense_cost=runtime.e["defense_cost"][p],
                    voc_cost=runtime.e["voc_cost"][p],
                    W_plant=runtime.e["wp"][p],
                )
            )
            if runtime.c.spatial is not None:
                runtime.result["plants"][-1].update(
                    cell_row=int(community.plants.cells[p, 0]),
                    cell_col=int(community.plants.cells[p, 1]),
                    density_phi=runtime.plant_phi[p],
                )
        runtime.network_start = len(runtime.result["networks"])
        for h in range(len(runtime.c.insect_sizes)):
            total_resource = float(runtime.e["x"][hl == h].sum())
            for p in range(len(runtime.c.plant_sizes)):
                matched = (hl == h) & (runtime.e["host"] >= 0)
                matched &= np.array(
                    [pl[k] == p if k >= 0 else False for k in runtime.e["host"]], bool
                )
                actual = float(runtime.e["x"][matched].sum())
                size = int(runtime.counts[1][h])
                runtime.result["networks"].append(
                    dict(
                        generation=state["current_step"],
                        insect_species=f"H{h + 1}",
                        plant_species=f"P{p + 1}",
                        N_accept=int(matched.sum()),
                        Q_actual=actual,
                        L_per_insect=actual / size if size else np.nan,
                        resource_share=actual / total_resource if total_resource else np.nan,
                        maternal_offspring=None,
                        paternal_offspring=None,
                        genetic_contribution=None,
                        maternal_share=None,
                        paternal_share=None,
                        genetic_contribution_share=None,
                        contribution_committed=False,
                    )
                )
                if runtime.observe:
                    runtime.result["networks"][-1]["candidate_exposure"] = int(
                        runtime.e["candidate_exposure"][h, pl == p].sum()
                    )
        if runtime.observe:
            rows, probability, q = standardized(
                community, runtime.c, runtime.B, runtime.w, state["current_step"]
            )
            runtime.result["diagnostics"].extend(rows)
            snapshot.update(
                reference_p=probability,
                reference_q=q,
                reference_plant_species=np.array(
                    [f"P{p + 1}" for p in range(len(runtime.c.plant_sizes))]
                ),
            )
        return {}
