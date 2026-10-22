from evo_sim.core.substep import SubstepTransition

from ..helpers import spatial
from ..helpers.agents import CommunityState
from ..helpers.demography import ComputationalLimit, inherit
from ..helpers.runtime import agent_state, community_state


class InheritAndDisperse(SubstepTransition):
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
        new_state = {}
        new_records = []
        if runtime.genetic is not None:
            runtime.genetic.preflight(
                [sum(int(x) for x in group["offspring_count"]) for group in runtime.sampled],
                runtime.stored_genotype_values,
            )
        for k, (kind, trait_names) in enumerate(
            (("plant", ("D", "plant_voc_independent_trait")), ("insect", ("I", "A")))
        ):
            traits = {name: community.arrays()[name] for name in trait_names}
            if runtime.genetic is not None:
                children, labels, ids, records, runtime.next_ids[k] = runtime.genetic.inherit(
                    (community.plants if k == 0 else community.insects).genotypes,
                    pl if k == 0 else hl,
                    community.plants.ids if k == 0 else community.insects.ids,
                    runtime.sampled[k]["offspring_count"],
                    kind,
                    runtime.next_ids[k],
                    state["current_step"] + 1,
                )
            else:
                children, labels, ids, records, runtime.next_ids[k] = inherit(
                    traits,
                    pl if k == 0 else hl,
                    community.plants.ids if k == 0 else community.insects.ids,
                    runtime.sampled[k]["offspring_count"],
                    runtime.c,
                    kind,
                    runtime.rng[4][k],
                    runtime.next_ids[k],
                    state["current_step"] + 1,
                )
            new_state.update(children)
            new_state[kind + "_species"] = labels
            new_state[kind + "_ids"] = ids
            if runtime.c.spatial is not None:
                cells = spatial.offspring_cells(
                    community,
                    runtime.e["host"],
                    runtime.sampled[k]["offspring_count"],
                    kind,
                    runtime.c,
                    runtime.spatial_rng[k],
                )
                new_state[kind + "_cells"] = cells
                for record, cell in zip(records, cells):
                    record.update(cell_row=int(cell[0]), cell_col=int(cell[1]))
            new_records.extend(records)
        for row in runtime.result["populations"][runtime.start :]:
            row["N_next"] = row["N_offspring"]
            row["update_applied"] = True
            if row["N_next"] == 0 and row["species"] not in runtime.extinct:
                runtime.extinct[row["species"]] = state["current_step"] + 1
                row["extinct_at_generation"] = state["current_step"] + 1
        runtime.result["offspring"].extend(new_records)
        community = CommunityState.from_arrays(new_state)
        return {self.output_variables[0]: agent_state(community)}
