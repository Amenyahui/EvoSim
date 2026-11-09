"""Unpublished generation event ledgers and genetic contribution records."""

import json
from typing import TYPE_CHECKING

import numpy as np

from evo_sim.models.plant_insect.biology.agents import CommunityState
from evo_sim.models.plant_insect.records.contributions import allele_contributions

if TYPE_CHECKING:
    from evo_sim.models.plant_insect.simulation.stages.adult_reproduction import AdultReproduction
    from evo_sim.models.plant_insect.simulation.stages.larval_development import LarvalDevelopment
    from evo_sim.models.plant_insect.simulation.stages.offspring_genetics import OffspringGenetics
    from evo_sim.models.plant_insect.simulation.stages.oviposition import Oviposition
    from evo_sim.models.plant_insect.simulation.stages.plant_dynamics import PlantDynamics
    from evo_sim.models.plant_insect.simulation.stages.plant_reproduction import PlantReproduction
    from evo_sim.models.plant_insect.simulation.stages.recruitment import Recruitment


def host_name(label):
    return "unknown" if label < 0 else f"P{label + 1}"


def source_flows(parents, mothers, fathers, destinations, labels, survived, generation, sizes):
    """Nuclear parent-role accounting, not an effective gene-flow estimator."""
    rows = []
    backgrounds = range(-1, sizes[0])
    for s in range(sizes[1]):
        for stage, selected in (("birth", labels == s), ("recruitment", (labels == s) & survived)):
            for source in backgrounds:
                for destination in backgrounds:
                    at_destination = selected & (destinations == destination)
                    count = (
                        sum(
                            np.count_nonzero(at_destination & (parents[role] == source))
                            for role in (mothers, fathers)
                        )
                        / 2
                    )
                    rows.append(
                        dict(
                            reproduction_round=generation,
                            adult_cohort=generation + 1,
                            insect_species=f"H{s + 1}",
                            stage=stage,
                            source_host_species=host_name(source),
                            destination_host_species=host_name(destination),
                            genetic_contribution=count,
                        )
                    )
    return rows


def generation_records(
    runtime,
    current: CommunityState,
    t: int,
    *,
    adults: "AdultReproduction",
    oviposition: "Oviposition",
    offspring: "OffspringGenetics",
    larvae: "LarvalDevelopment",
    plant_update: "PlantDynamics",
    plant_births: "PlantReproduction",
    recruitment: "Recruitment",
) -> tuple[dict, dict]:
    """Format recorded stage outcomes and contributions before atomic publication."""
    c = runtime.c
    pl, hl = current.plants.species, current.insects.species
    insects, plants = current.insects, current.plants
    new = {
        key: []
        for key in (
            "offspring",
            "parents",
            "populations",
            "events",
            "plants",
            "networks",
            "oviposition",
            "larvae",
            "host_flows",
            "lagged_networks",
            "allele_contributions",
            "plant_proposals",
        )
    }
    propagule_indices = (
        np.concatenate([np.arange(n) for n in plant_update.pcounts])
        if plant_update.npropagules
        else np.empty(0, int)
    )
    established_ids = dict(zip(plant_update.winners.tolist(), plant_births.pids.tolist()))
    for j, mother in enumerate(plant_update.propagule_mothers):
        new["plant_proposals"].append(
            dict(
                event_id=f"{t}:{int(plants.ids[mother])}:{int(propagule_indices[j])}",
                reproduction_round=t,
                offspring_generation=t + 1,
                maternal_id=int(plants.ids[mother]),
                propagule_index=int(propagule_indices[j]),
                cell_row=int(plant_update.proposal_cells[j, 0]),
                cell_col=int(plant_update.proposal_cells[j, 1]),
                established=bool(plant_update.proposal_sites[j] >= 0),
                microsite=int(plant_update.proposal_sites[j])
                if plant_update.proposal_sites[j] >= 0
                else None,
                offspring_id=int(established_ids[j]) if j in established_ids else None,
            )
        )
    egg_indices = (
        np.concatenate([np.arange(n) for n in adults.counts])
        if adults.potential_count
        else np.empty(0, int)
    )
    egg_children = dict(zip(np.flatnonzero(oviposition.accepted), offspring.hids))
    for j, m in enumerate(oviposition.egg_mothers):
        p = oviposition.choice["host"][j]
        new["oviposition"].append(
            dict(
                event_id=f"{t}:{int(insects.ids[m])}:{int(egg_indices[j])}",
                offspring_id=int(egg_children[j]) if j in egg_children else None,
                reproduction_round=t,
                adult_cohort=t + 1,
                maternal_id=int(insects.ids[m]),
                egg_index=int(egg_indices[j]),
                host_id=int(plants.ids[p]) if p >= 0 else -1,
                plant_species=host_name(pl[p]) if p >= 0 else "unknown",
                accepted=bool(p >= 0),
                reason=str(oviposition.choice["reason"][j]),
                candidate_count=int(oviposition.choice["candidate_count"][j]),
                candidate_plant_ids=json.dumps(
                    plants.ids[oviposition.choice["candidate_rows"][j]].tolist()
                ),
                choice_order=int(oviposition.choice["choice_order"][j]),
                contact_plant_ids=json.dumps(
                    plants.ids[oviposition.choice["contact_rows"][j]].tolist()
                ),
                contact_probabilities=json.dumps(oviposition.choice["contact_probabilities"][j]),
                contact_accepted=json.dumps(oviposition.choice["contact_accepted"][j]),
            )
        )
    accepted_indices = np.flatnonzero(oviposition.accepted)
    for i, record in enumerate(offspring.hrecords):
        p, m, f = oviposition.host[i], oviposition.mothers[i], oviposition.fathers[i]
        egg = accepted_indices[i]
        known = insects.birth_host_species[m] >= 0 and insects.birth_host_species[f] >= 0
        different = bool(insects.birth_host_species[m] != insects.birth_host_species[f])
        event_id = f"{t}:{int(insects.ids[m])}:{int(egg_indices[egg])}"
        record.update(
            oviposition_event_id=event_id,
            egg_index=int(egg_indices[egg]),
            reproduction_round=t,
            adult_cohort=t + 1,
            birth_host_id=int(plants.ids[p]),
            birth_host_species=host_name(pl[p]),
            maternal_birth_host_species=host_name(insects.birth_host_species[m]),
            paternal_birth_host_species=host_name(insects.birth_host_species[f]),
            cross_background_pair=different if known else None,
            pair_background_status="unknown" if not known else "different" if different else "same",
            survived_to_adult=bool(larvae.survived[i]),
            offspring_sex=int(larvae.sex[i]),
            maternal_cell=insects.cells[m].tolist(),
            paternal_cell=insects.cells[f].tolist(),
            offspring_cell=larvae.origins[i].tolist(),
            cell_row=int(larvae.origins[i, 0]),
            cell_col=int(larvae.origins[i, 1]),
            genotype=json.dumps(offspring.hvalues["insect_genotypes"][i].tolist()),
        )
        event = dict(
            oviposition_event_id=event_id,
            egg_index=int(egg_indices[egg]),
            generation=t,
            reproduction_round=t,
            adult_cohort=t + 1,
            seed=c.seed,
            insect_id=int(offspring.hids[i]),
            insect_row=i,
            maternal_id=record["maternal_id"],
            paternal_id=record["paternal_id"],
            species=f"H{offspring.child_hl[i] + 1}",
            species_id=int(offspring.child_hl[i]),
            host_row=int(p),
            host_id=int(plants.ids[p]),
            plant_species=host_name(pl[p]),
            birth_host_id=int(plants.ids[p]),
            birth_host_species=host_name(pl[p]),
            accepted=True,
            reason="accepted",
            candidate_count=int(oviposition.choice["candidate_count"][accepted_indices[i]]),
            choice_order=int(oviposition.choice["choice_order"][accepted_indices[i]]),
            q=float(larvae.e["q"][i]),
            d=float(larvae.e["d"][i]),
            x=float(larvae.e["x"][i]),
            g=float(larvae.e["g"][i]),
            dmax=float(larvae.dmax[i]),
            R=float(larvae.capacity[p]),
            gamma=float(larvae.e["gamma"][p]),
            host_h=float(larvae.e["h"][p]),
            I_sum=float(offspring.hvalues["I"][i].sum()),
            W_insect=float(larvae.e["survival_probability"][i]),
            survived_to_adult=bool(larvae.survived[i]),
            insect_life_stage="larva",
            same_species_others=int(larvae.e["same_species_others"][i]),
            other_species_others=int(larvae.e["other_species_others"][i]),
        )
        for name in ("dose", "processing_capacity", "processed", "residual", "antinutrient"):
            for axis in range(2):
                event[f"{name}_{axis + 1}"] = float(larvae.e[name][i, axis])
        for name in ("nutrition_survival", "toxicity_survival", "antinutrition_survival"):
            event[name] = float(larvae.e[name][i])
        new["events"].append(event)
        new["larvae"].append(dict(event))
    for i, record in enumerate(plant_births.precords):
        record.update(
            reproduction_round=t,
            adult_cohort=t + 1,
            survived_to_adult=True,
            offspring_cell=plant_births.pcells[i].tolist(),
            cell_row=int(plant_births.pcells[i, 0]),
            cell_col=int(plant_births.pcells[i, 1]),
            maternal_cell=plants.cells[plant_births.pmothers[i]].tolist(),
            paternal_cell=plants.cells[plant_births.pfathers[i]].tolist(),
            genotype=json.dumps(plant_births.pvalues["plant_genotypes"][i].tolist()),
        )
    new["offspring"] = plant_births.precords + offspring.hrecords
    new["parents"] = plant_update.prows + adults.hrows
    for i in range(len(plants)):
        new["plants"].append(
            dict(
                generation=t,
                plant_id=int(plants.ids[i]),
                plant_row=i,
                species=host_name(pl[i]),
                species_id=int(pl[i]),
                R=float(larvae.capacity[i]),
                T=float(larvae.e["total_demand"][i]),
                gamma=float(larvae.e["gamma"][i]),
                X=float(larvae.e["X"][i]),
                h=float(larvae.e["h"][i]),
                herbivore_load=int(larvae.e["herbivore_load"][i]),
                D_sum=float(plants.defense[i].sum()),
                V_emit_mean=float(oviposition.V[i].mean()),
                defense_cost=float(larvae.e["defense_cost"][i]),
                voc_cost=float(larvae.e["voc_cost"][i]),
                W_plant=float(plant_update.prows[i]["W_i"]),
                resource_capacity=float(larvae.maximum[i]),
                resource_after_feeding=float(plant_update.post[i]),
                tolerance_Z=float(plants.chemical_traits[i, 2]),
                recovery_rate=float(plant_update.recovery_rate[i]),
                resource_after_recovery=float(plant_update.recovered[i])
                if plant_update.psurvived[i]
                else None,
                defoliation=float(plant_update.defoliation[i]),
                survival_probability=float(plant_update.S[i]),
                survived=bool(plant_update.psurvived[i]),
                resource_recovery=float(plant_update.recovered[i] - plant_update.post[i])
                if plant_update.psurvived[i]
                else None,
                resource_next=float(plant_update.recovered[i])
                if plant_update.psurvived[i]
                else None,
                microsite=int(plants.microsite[i]),
                reproduction_active=True,
                cell_row=int(plants.cells[i, 0]),
                cell_col=int(plants.cells[i, 1]),
            )
        )
    extinct = runtime.extinct.copy()
    for k, (kind, before, after, rows, live) in enumerate(
        (
            (
                "plant",
                plants,
                recruitment.next_state.plants,
                plant_update.prows,
                plant_update.psurvived,
            ),
            (
                "insect",
                insects,
                recruitment.next_state.insects,
                adults.hrows,
                np.ones(len(insects), bool),
            ),
        )
    ):
        sizes = c.plant_sizes if k == 0 else c.insect_sizes
        for s in range(len(sizes)):
            par = [r for r in rows if r["species_id"] == s]
            start = int(np.count_nonzero(before.species == s))
            nxt = int(np.count_nonzero(after.species == s))
            name = f"{'P' if k == 0 else 'H'}{s + 1}"
            if not start:
                extinct.setdefault(name, t)
            elif not nxt:
                extinct.setdefault(name, t + 1)
            new["populations"].append(
                dict(
                    generation=t,
                    kind=kind,
                    species=name,
                    N_start=start,
                    N_survived_to_reproduce=int(live[before.species == s].sum()),
                    N_offspring=sum(r["offspring_count_as_mother"] for r in par),
                    N_next=nxt,
                    N_retained=int(live[before.species == s].sum()) if not k else 0,
                    N_propagules=sum(int(r["propagule_count"]) for r in par) if not k else None,
                    N_potential_eggs=sum(int(r["potential_egg_count"]) for r in par) if k else None,
                    update_applied=True,
                    extinct_at_generation=extinct.get(name),
                    zero_fitness_count=sum(r["W_i"] == 0 for r in par),
                    population_status="extant" if start else "extinct",
                    N_female=int(np.count_nonzero((before.species == s) & (insects.sex == 0)))
                    if k
                    else None,
                    N_male=int(np.count_nonzero((before.species == s) & (insects.sex == 1)))
                    if k
                    else None,
                    N_mothers_with_mates=sum(
                        r["eligible_as_mother"] and r["has_mate"] for r in par
                    ),
                    N_seed_or_female_parents=sum(r["offspring_count_as_mother"] > 0 for r in par),
                    N_pollen_or_male_parents=sum(r["offspring_count_sired"] > 0 for r in par),
                    total_resource_transfer=float(larvae.e["X"][pl == s].sum())
                    if not k
                    else float(larvae.e["x"][offspring.child_hl == s].sum()),
                )
            )
    for h in range(len(c.insect_sizes)):
        total_resource = float(larvae.e["x"][offspring.child_hl == h].sum())
        for p in range(len(c.plant_sizes)):
            selected = (offspring.child_hl == h) & (pl[oviposition.host] == p)
            resource = float(larvae.e["x"][selected].sum())
            new["networks"].append(
                dict(
                    generation=t,
                    reproduction_round=t,
                    adult_cohort=t + 1,
                    insect_species=f"H{h + 1}",
                    plant_species=host_name(p),
                    N_accept=int(selected.sum()),
                    Q_actual=resource,
                    resource_share=resource / total_resource if total_resource else np.nan,
                    L_per_insect=resource / int(selected.sum()) if selected.any() else np.nan,
                    successful_recruits=int(np.count_nonzero(selected & larvae.survived)),
                    birth_count=int(selected.sum()),
                    contribution_committed=True,
                    edge_semantics="larval_resource_and_recruitment",
                )
            )
        for p in range(-1, len(c.plant_sizes)):
            selected = (hl == h) & (insects.birth_host_species == p)
            par = [adults.hrows[i] for i in np.flatnonzero(selected)]
            new["lagged_networks"].append(
                dict(
                    reproduction_round=t,
                    adult_cohort=t + 1,
                    insect_species=f"H{h + 1}",
                    plant_species=host_name(p),
                    maternal_offspring=sum(r["offspring_count_as_mother"] for r in par),
                    paternal_offspring=sum(r["offspring_count_sired"] for r in par),
                    genetic_contribution=sum(r["realized_genetic_contribution"] for r in par),
                    recruitment_contribution=sum(
                        r["realized_recruitment_contribution"] for r in par
                    ),
                    edge_semantics="parent_birth_host_lagged_contribution",
                )
            )
    new["host_flows"] = source_flows(
        insects.birth_host_species,
        oviposition.mothers,
        oviposition.fathers,
        pl[oviposition.host],
        offspring.child_hl,
        larvae.survived,
        t,
        (len(c.plant_sizes), len(c.insect_sizes)),
    )
    for stage, after, parent_groups in (
        ("birth", recruitment.birth_state, [plant_update.prows, adults.hrows]),
        (
            "recruitment",
            recruitment.recruitment_state,
            [
                [
                    dict(r, realized_genetic_contribution=r["realized_recruitment_contribution"])
                    for r in rows
                ]
                for rows in (plant_update.prows, adults.hrows)
            ],
        ),
    ):
        rows = allele_contributions(
            current,
            after,
            parent_groups,
            offspring.genetic,
            t,
            (len(c.plant_sizes), len(c.insect_sizes)),
            plant_births.precords,
        )
        new["allele_contributions"].extend(dict(row, stage=stage) for row in rows)
    return new, extinct
