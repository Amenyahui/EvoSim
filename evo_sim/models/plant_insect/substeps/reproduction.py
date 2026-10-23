import copy

import numpy as np

from evo_sim.core.substep import SubstepTransition

from ..helpers import spatial
from ..helpers.agents import CommunityState
from ..helpers.demography import ComputationalLimit
from ..helpers.runtime import agent_state, community_state


class CrossFertilizeAndDisperse(SubstepTransition):
    def forward(self, state, action):
        runtime = self.input(state, "runtime")
        community = community_state(self.input(state, "agents"))
        pl, hl = community.plants.species, community.insects.species
        if not len(pl) + len(hl):
            runtime.status = dict(
                status="all_species_extinct",
                generation=state["current_step"],
                message="no individuals; absorbing state, no immigration",
            )
            return {"terminated": True}
        if not runtime.update:
            runtime.status = dict(
                status="completed",
                generation=state["current_step"],
                message="requested number of generation updates completed",
            )
            return {"terminated": True}
        for row in runtime.result["populations"][runtime.start :]:
            k = 0 if row["kind"] == "plant" else 1
            labels = pl if k == 0 else hl
            s = int(row["species"][1:]) - 1
            row["N_survived_to_reproduce"] = int(
                runtime.sampled[k]["survived_to_reproduce"][labels == s].sum()
            )
            row["N_offspring"] = sum(
                int(x) for x in runtime.sampled[k]["offspring_count"][labels == s]
            )
        total = sum(int(x) for group in runtime.sampled for x in group["offspring_count"])
        runtime.stage = "offspring_allocation"
        if total > runtime.c.max_individuals_safety:
            raise ComputationalLimit(
                f"sampled {total} offspring exceeds safety {runtime.c.max_individuals_safety}; no offspring truncated or committed"
            )
        runtime.genetic.preflight(
            [sum(int(x) for x in group["offspring_count"]) for group in runtime.sampled],
            runtime.stored_genotype_values,
        )
        # ponytail: copy catalog lists, sharing immutable Alleles; persistent catalogs if scale demands it.
        runtime.stage = "meiosis_mutation"
        genetic = copy.copy(runtime.genetic)
        genetic.alleles = runtime.genetic.alleles.copy()
        genetic.events = runtime.genetic.events.copy()
        genetic.transmissions = runtime.genetic.transmissions.copy()
        genetic.rng = copy.deepcopy(runtime.genetic.rng)
        new_state, new_records, next_ids, sired_counts = {}, [], [], []
        for k, kind in enumerate(("plant", "insect")):
            pop = community.plants if k == 0 else community.insects
            sampled = runtime.sampled[k]
            mothers = np.repeat(np.arange(len(pop)), sampled["offspring_count"])
            fathers = np.empty(len(mothers), dtype=int)
            for i, mother in enumerate(mothers):
                fathers[i] = runtime.rng[5][k][pop.species[mother]].choice(
                    sampled["candidates"][mother]
                )
            children, labels, ids, records, next_id = genetic.fertilize(
                pop.genotypes,
                pop.species,
                pop.ids,
                mothers,
                fathers,
                kind,
                runtime.next_ids[k],
                state["current_step"] + 1,
            )
            next_ids.append(next_id)
            sired_counts.append(np.bincount(fathers, minlength=len(pop)))
            new_state.update(children)
            new_state[kind + "_species"] = labels
            new_state[kind + "_ids"] = ids
            cells = spatial.offspring_cells(
                community,
                runtime.e["host"],
                sampled["offspring_count"],
                kind,
                runtime.c,
                runtime.spatial_rng[k],
            )
            new_state[kind + "_cells"] = cells
            sex = np.full(len(ids), -1, dtype=np.int8)
            if k == 1:
                probability = runtime.c.sexual_reproduction["insect_offspring_female_probability"]
                for s, rng in enumerate(runtime.rng[7][k]):
                    mask = labels == s
                    sex[mask] = (rng.random(mask.sum()) >= probability).astype(np.int8)
                new_state["insect_sex"] = sex
            for i, record in enumerate(records):
                record.update(
                    cell_row=int(cells[i, 0]),
                    cell_col=int(cells[i, 1]),
                    maternal_cell=sampled["cells"][mothers[i]].tolist(),
                    paternal_cell=sampled["cells"][fathers[i]].tolist(),
                    offspring_cell=cells[i].tolist(),
                    offspring_sex=int(sex[i]) if k == 1 else None,
                )
            new_records.extend(records)
        runtime.stage = "offspring_validation"
        community_next = CommunityState.from_arrays(new_state)
        community_next.validate()
        # Prepare every derived record and returned agent mapping before publishing anything.
        committed_records = runtime.result["offspring"] + new_records
        parent_groups = [[dict(row) for row in group["parent_rows"]] for group in runtime.sampled]
        for rows, counts in zip(parent_groups, sired_counts):
            for row, count in zip(rows, counts):
                row["offspring_count_sired"] = int(count)
        population_rows = [dict(row) for row in runtime.result["populations"][runtime.start :]]
        extinct = runtime.extinct.copy()
        for row in population_rows:
            k = 0 if row["kind"] == "plant" else 1
            parents = [p for p in parent_groups[k] if p["species"] == row["species"]]
            row["N_mothers_with_mates"] = sum(
                p["eligible_as_mother"] and p["has_mate"] for p in parents
            )
            row["N_seed_or_female_parents"] = sum(
                p["offspring_count_as_mother"] > 0 for p in parents
            )
            row["N_pollen_or_male_parents"] = sum(p["offspring_count_sired"] > 0 for p in parents)
            row["N_next"] = row["N_offspring"]
            row["update_applied"] = True
            if row["N_next"] == 0 and row["species"] not in extinct:
                extinct[row["species"]] = state["current_step"] + 1
                row["extinct_at_generation"] = state["current_step"] + 1
        parents_start = len(runtime.result["parents"]) - len(pl) - len(hl)
        committed_parents = (
            runtime.result["parents"][:parents_start] + parent_groups[0] + parent_groups[1]
        )
        committed_populations = runtime.result["populations"][: runtime.start] + population_rows
        output = {self.output_variables[0]: agent_state(community_next)}
        runtime.genetic = genetic
        runtime.next_ids = next_ids
        runtime.extinct = extinct
        runtime.result["offspring"] = committed_records
        runtime.result["parents"] = committed_parents
        runtime.result["populations"] = committed_populations
        runtime.reproduction_checkpoint = None
        return output
