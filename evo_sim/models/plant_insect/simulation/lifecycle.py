"""Current adult reproduction -> larval resource window -> adult recruitment.

One disposable transition owns the coupled plant/insect update and publishes its
stage ledgers only after both taxa and the next snapshot fit the storage budget.
"""

import copy

from evo_sim.core.substep import SubstepTransition
from evo_sim.models.plant_insect.biology.agents import CommunityState
from evo_sim.models.plant_insect.records.events import generation_records
from evo_sim.models.plant_insect.records.observations import observe
from evo_sim.models.plant_insect.simulation.runtime import (
    agent_state,
    community_state,
)

from .stages.adult_reproduction import prepare_adults
from .stages.larval_development import develop_larvae
from .stages.offspring_genetics import inherit_offspring
from .stages.oviposition import oviposit
from .stages.plant_dynamics import update_plants
from .stages.plant_reproduction import reproduce_plants
from .stages.recruitment import recruit


def transition(runtime, current: CommunityState, t: int):
    """Execute existing stage order; return disposable state and unpublished ledgers."""
    adults = prepare_adults(runtime, current, t)
    oviposition = oviposit(runtime, current, adults=adults)
    offspring = inherit_offspring(runtime, current, t, oviposition=oviposition)
    larvae = develop_larvae(
        runtime, current, t, adults=adults, oviposition=oviposition, offspring=offspring
    )
    plant_update = update_plants(runtime, current, t, offspring=offspring, larvae=larvae)
    plant_births = reproduce_plants(
        runtime, current, t, offspring=offspring, larvae=larvae, plant_update=plant_update
    )
    recruitment = recruit(
        current,
        offspring=offspring,
        larvae=larvae,
        plant_update=plant_update,
        plant_births=plant_births,
    )
    records, extinct = generation_records(
        runtime,
        current,
        t,
        adults=adults,
        oviposition=oviposition,
        offspring=offspring,
        larvae=larvae,
        plant_update=plant_update,
        plant_births=plant_births,
        recruitment=recruitment,
    )
    return (
        recruitment.next_state,
        offspring.genetic,
        [plant_births.next_p, offspring.next_h],
        records,
        extinct,
    )


class AdvanceGeneration(SubstepTransition):
    def forward(self, state, action):
        runtime = self.input(state, "runtime")
        current = community_state(self.input(state, "agents"))
        t = state["current_step"]
        observe(runtime, current, t)
        if t >= runtime.c.generations or not len(current.plants) + len(current.insects):
            runtime.status = dict(
                status="completed" if t >= runtime.c.generations else "all_species_extinct",
                generation=t,
                message="adult/plant snapshot; no larval window executed",
            )
            return {"terminated": True}
        runtime.reproduction_checkpoint = (
            copy.deepcopy(runtime.rng),
            copy.deepcopy(runtime.spatial_rng),
        )
        next_state, genetic, next_ids, records, extinct = transition(runtime, current, t)
        # Prepare allocation-prone derived values before publishing the transition.
        committed = {name: runtime.result[name] + rows for name, rows in records.items()}
        stored_values = runtime.stored_genotype_values + sum(
            2 * len(genetic.loci[r["kind"]]) for r in records["offspring"]
        )
        output = {"agents": agent_state(next_state)}
        runtime.result.update(committed)
        runtime.stored_genotype_values = stored_values
        runtime.genetic, runtime.next_ids, runtime.extinct = genetic, next_ids, extinct
        runtime.reproduction_checkpoint = None
        return output
