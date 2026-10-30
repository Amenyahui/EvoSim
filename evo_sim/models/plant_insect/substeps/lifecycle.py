"""v0.1.6 adult reproduction -> larval resource window -> adult recruitment.

One disposable transition owns the coupled plant/insect update. Biological kernels
and SLiM transport are shared with the previous model; only event order changes.
"""

import copy
import json

import numpy as np

from evo_sim.core.substep import SubstepTransition

from ..helpers import chemistry, modules, spatial
from ..helpers import ecology as feeding
from ..helpers.agents import CommunityState
from ..helpers.contributions import allele_contributions, commit_parent_contributions
from ..helpers.demography import ComputationalLimit, maternal_moments, mating_candidates
from ..helpers.observations import standardized
from ..helpers.runtime import agent_state, community_state, resource_vector


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


def sample_parents(pop, kind, S, F, phi, runtime, generation):
    """Plants survive this window; current adult insects have already survived it."""
    k = int(kind == "insect")
    c = runtime.c
    maternal = np.ones(len(pop), bool) if k == 0 else pop.sex == 0
    donor = np.ones(len(pop), bool) if k == 0 else pop.sex == 1
    radius = c.sexual_reproduction["pollen_radius" if k == 0 else "insect_mate_radius"]
    potential = mating_candidates(
        pop.species, pop.cells, donor, maternal, radius, c.spatial["grid_shape"]
    )
    mu = (runtime.lp if k == 0 else runtime.li)[pop.species] * F * phi
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
                sex_or_reproductive_role="seed_and_pollen"
                if k == 0
                else ("female" if pop.sex[i] == 0 else "male"),
                is_maternal_role=bool(maternal[i]),
                S_i=float(S[i]),
                F_i=float(F[i]),
                W_i=float(S[i] * F[i]),
                phi_i=float(phi[i]),
                lambda_s=float((runtime.lp if k == 0 else runtime.li)[pop.species[i]]),
                conditional_mu=float(mu[i]),
                mate_probability=float(q[i]),
                expected_maternal_offspring=float(expected[i]),
                maternal_offspring_variance=float(variance[i]),
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
    if modules.enabled(c):
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


def contributions(rows, mothers, fathers, survived):
    bm = np.bincount(mothers, minlength=len(rows))
    bf = np.bincount(fathers, minlength=len(rows))
    rm = np.bincount(mothers[survived], minlength=len(rows))
    rf = np.bincount(fathers[survived], minlength=len(rows))
    for i, row in enumerate(rows):
        row.update(
            offspring_count_as_mother=int(bm[i]),
            offspring_count_sired=int(bf[i]),
            recruited_count_as_mother=int(rm[i]),
            recruited_count_sired=int(rf[i]),
            realized_recruitment_contribution=float((rm[i] + rf[i]) / 2),
        )
    commit_parent_contributions([rows])


def observe(runtime, current, generation):
    snapshot = current.snapshot(generation)
    snapshot["V_emit"] = chemistry.plant_voc(
        current.plants.defense,
        current.plants.voc_trait,
        runtime.c,
        runtime.B,
        current.plants.chemical_traits,
    )
    snapshot["R"] = resource_vector(
        runtime.c.resource_capacity, runtime.c.plant_sizes, current.plants.species, "P"
    )
    runtime.result["snapshots"].append(snapshot)
    runtime.stored_genotype_values += current.plants.genotypes.size + current.insects.genotypes.size
    sizes = (len(runtime.c.plant_sizes), len(runtime.c.insect_sizes))
    runtime.result["genetic_summary"].extend(runtime.genetic.summary(current, generation, sizes))
    runtime.result["haplotype_summary"].extend(
        runtime.genetic.haplotype_summary(current, generation, sizes)
    )
    if runtime.observe:
        rows, probability, q = standardized(current, runtime.c, runtime.B, runtime.w, generation)
        runtime.result["diagnostics"].extend(rows)
        snapshot.update(
            reference_p=probability,
            reference_q=q,
            reference_plant_species=np.array(
                runtime.c.observation_reference_traits["reference_ids"]
                if modules.enabled(runtime.c)
                else [f"P{p + 1}" for p in range(sizes[0])]
            ),
        )
    for h in range(sizes[1]):
        for p in range(-1, sizes[0]):
            count = int(
                np.count_nonzero(
                    (current.insects.species == h) & (current.insects.birth_host_species == p)
                )
            )
            runtime.result["host_groups"].append(
                dict(
                    generation=generation,
                    insect_species=f"H{h + 1}",
                    birth_host_species=host_name(p),
                    N=count,
                    status="present" if count else "absent",
                )
            )
    return snapshot


def transition(runtime, current, t):
    c = runtime.c
    pl, hl = current.plants.species, current.insects.species
    insects, plants = current.insects, current.plants
    modern = modules.enabled(c)
    # Adults retain the birth-host-plus-dispersal cells when the local mate pool is built.
    runtime.stage = "adult_mating_and_potential_eggs"
    if modern:
        fi, processing_cost, sensory_cost = modules.maternal_fecundity(insects.module_traits, c)
    else:
        fi = np.maximum(0, 1 - c.insect_cost * insects.resistance.sum(axis=1))
    counts, candidates, hrows, _ = sample_parents(
        insects, "insect", np.ones(len(insects)), fi, np.ones(len(insects)), runtime, t
    )
    if modern:
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
    runtime.stage = "maternal_oviposition"
    egg_mothers = np.repeat(np.arange(len(insects)), counts)
    V = chemistry.plant_voc(plants.defense, plants.voc_trait, c, runtime.B, plants.chemical_traits)
    V_cue = V[runtime.rng[1][0][0].permutation(len(V))] if c.shuffle_voc else V.copy()
    choice = feeding.choose_hosts(
        (insects.module_traits if modern else insects.preference)[egg_mothers],
        pl,
        hl[egg_mothers],
        V_cue,
        c,
        runtime.w,
        runtime.rng[1][0][0],
        runtime.rng[1][1],
        candidate_pools=spatial.candidate_pools(plants.cells, insects.cells[egg_mothers], c),
    )
    accepted = choice["host"] >= 0
    mothers, host = egg_mothers[accepted], choice["host"][accepted]
    fathers = draw_fathers(insects, mothers, candidates, runtime, 1)
    runtime.genetic.preflight([0, len(mothers)], runtime.stored_genotype_values)
    genetic = copy.copy(runtime.genetic)
    genetic.alleles = runtime.genetic.alleles.copy()
    genetic.events = runtime.genetic.events.copy()
    genetic.transmissions = runtime.genetic.transmissions.copy()
    genetic.rng = copy.deepcopy(runtime.genetic.rng)
    runtime.stage = "insect_meiosis_mutation"
    hvalues, child_hl, hids, hrecords, next_h = genetic.fertilize(
        insects.genotypes, hl, insects.ids, mothers, fathers, "insect", runtime.next_ids[1], t + 1
    )
    capacity = resource_vector(c.resource_capacity, c.plant_sizes, pl, "P")
    dmax = resource_vector(c.max_feeding_demand, c.insect_sizes, child_hl, "H")
    runtime.stage = "larval_feeding_and_survival"
    e = feeding.settle_feeding(
        plants.defense,
        hvalues["I"],
        pl,
        child_hl,
        V,
        host,
        c,
        capacity,
        dmax,
        chemical=plants.chemical_traits,
        module_traits=hvalues.get("insect_module_traits"),
    )
    survived = np.zeros(len(hids), bool)
    for s, rng in enumerate(runtime.rng[2][1]):
        selected = child_hl == s
        survived[selected] = rng.random(selected.sum()) < e["g"][selected]
    sex = np.empty(len(hids), dtype=np.int8)
    for s, rng in enumerate(runtime.rng[7][1]):
        selected = child_hl == s
        sex[selected] = (
            rng.random(selected.sum())
            >= c.sexual_reproduction["insect_offspring_female_probability"]
        )
    origins = plants.cells[host]
    adult_cells = spatial.disperse(
        origins[survived],
        c.spatial["insect_dispersal_radius"],
        c.spatial["grid_shape"],
        runtime.spatial_rng[1],
    )
    hvalues.update(
        insect_species=child_hl,
        insect_ids=hids,
        insect_cells=origins.copy(),
        insect_sex=sex,
        insect_birth_host_id=plants.ids[host].copy(),
        insect_birth_host_species=pl[host].copy(),
        insect_adult_cohort=np.full(len(hids), t + 1, dtype=int),
    )
    contributions(hrows, mothers, fathers, survived)
    runtime.stage = "plant_survival_and_reproduction"
    phi = spatial.local_competition(plants.cells, pl, runtime.alpha, c)
    S = c.sP0 + (1 - c.sP0) * (1 - e["h"])
    pcounts, pcandidates, prows, psurvived = sample_parents(
        plants, "plant", S, e["fp"], phi, runtime, t
    )
    nplants = sum(map(int, pcounts))
    if nplants + len(hids) > c.max_individuals_safety:
        raise ComputationalLimit("birth allocation exceeds safety; no transition committed")
    genetic.preflight(
        [nplants, 0], runtime.stored_genotype_values + hvalues["insect_genotypes"].size
    )
    pmothers = np.repeat(np.arange(len(plants)), pcounts)
    pfathers = draw_fathers(plants, pmothers, pcandidates, runtime, 0)
    runtime.stage = "plant_meiosis_mutation"
    pvalues, child_pl, pids, precords, next_p = genetic.fertilize(
        plants.genotypes, pl, plants.ids, pmothers, pfathers, "plant", runtime.next_ids[0], t + 1
    )
    pcells = spatial.disperse(
        plants.cells[pmothers],
        c.spatial["plant_dispersal_radius"],
        c.spatial["grid_shape"],
        runtime.spatial_rng[0],
    )
    pvalues.update(plant_species=child_pl, plant_ids=pids, plant_cells=pcells)
    contributions(prows, pmothers, pfathers, np.ones(len(pids), bool))
    birth_values = dict(pvalues, **hvalues)
    birth_state = CommunityState.from_arrays(birth_values)
    next_values = dict(pvalues, **{key: value[survived] for key, value in hvalues.items()})
    next_values["insect_cells"] = adult_cells
    next_state = CommunityState.from_arrays(next_values)
    next_state.validate()
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
        )
    }
    egg_indices = (
        np.concatenate([np.arange(n) for n in counts]) if potential_count else np.empty(0, int)
    )
    egg_children = dict(zip(np.flatnonzero(accepted), hids))
    for j, m in enumerate(egg_mothers):
        p = choice["host"][j]
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
                reason=str(choice["reason"][j]),
                candidate_count=int(choice["candidate_count"][j]),
                candidate_plant_ids=json.dumps(plants.ids[choice["candidate_rows"][j]].tolist()),
                choice_order=int(choice["choice_order"][j]),
            )
        )
    accepted_indices = np.flatnonzero(accepted)
    for i, record in enumerate(hrecords):
        p, m, f = host[i], mothers[i], fathers[i]
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
            survived_to_adult=bool(survived[i]),
            offspring_sex=int(sex[i]),
            maternal_cell=insects.cells[m].tolist(),
            paternal_cell=insects.cells[f].tolist(),
            offspring_cell=origins[i].tolist(),
            cell_row=int(origins[i, 0]),
            cell_col=int(origins[i, 1]),
            genotype=json.dumps(hvalues["insect_genotypes"][i].tolist()),
        )
        event = dict(
            oviposition_event_id=event_id,
            egg_index=int(egg_indices[egg]),
            generation=t,
            reproduction_round=t,
            adult_cohort=t + 1,
            seed=c.seed,
            insect_id=int(hids[i]),
            insect_row=i,
            maternal_id=record["maternal_id"],
            paternal_id=record["paternal_id"],
            species=f"H{child_hl[i] + 1}",
            species_id=int(child_hl[i]),
            host_row=int(p),
            host_id=int(plants.ids[p]),
            plant_species=host_name(pl[p]),
            birth_host_id=int(plants.ids[p]),
            birth_host_species=host_name(pl[p]),
            accepted=True,
            reason="accepted",
            candidate_count=int(choice["candidate_count"][accepted_indices[i]]),
            choice_order=int(choice["choice_order"][accepted_indices[i]]),
            q=float(e["q"][i]),
            d=float(e["d"][i]),
            x=float(e["x"][i]),
            g=float(e["g"][i]),
            dmax=float(dmax[i]),
            R=float(capacity[p]),
            gamma=float(e["gamma"][p]),
            host_h=float(e["h"][p]),
            I_sum=float(hvalues["I"][i].sum()),
            W_insect=float(e["g"][i]),
            survived_to_adult=bool(survived[i]),
            insect_life_stage="larva",
            same_species_others=int(e["same_species_others"][i]),
            other_species_others=int(e["other_species_others"][i]),
        )
        new["events"].append(event)
        new["larvae"].append(dict(event))
    for i, record in enumerate(precords):
        record.update(
            reproduction_round=t,
            adult_cohort=t + 1,
            survived_to_adult=True,
            offspring_cell=pcells[i].tolist(),
            cell_row=int(pcells[i, 0]),
            cell_col=int(pcells[i, 1]),
            maternal_cell=plants.cells[pmothers[i]].tolist(),
            paternal_cell=plants.cells[pfathers[i]].tolist(),
            genotype=json.dumps(pvalues["plant_genotypes"][i].tolist()),
        )
    new["offspring"] = precords + hrecords
    new["parents"] = prows + hrows
    for i in range(len(plants)):
        new["plants"].append(
            dict(
                generation=t,
                plant_id=int(plants.ids[i]),
                plant_row=i,
                species=host_name(pl[i]),
                species_id=int(pl[i]),
                R=float(capacity[i]),
                T=float(e["total_demand"][i]),
                gamma=float(e["gamma"][i]),
                X=float(e["X"][i]),
                h=float(e["h"][i]),
                herbivore_load=int(e["herbivore_load"][i]),
                D_sum=float(plants.defense[i].sum()),
                V_emit_mean=float(V[i].mean()),
                defense_cost=float(e["defense_cost"][i]),
                voc_cost=float(e["voc_cost"][i]),
                W_plant=float(e["wp"][i]),
                density_phi=float(phi[i]),
                cell_row=int(plants.cells[i, 0]),
                cell_col=int(plants.cells[i, 1]),
            )
        )
    extinct = runtime.extinct.copy()
    for k, (kind, before, after, rows, live) in enumerate(
        (
            ("plant", plants, next_state.plants, prows, psurvived),
            ("insect", insects, next_state.insects, hrows, np.ones(len(insects), bool)),
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
                    total_resource_transfer=float(e["X"][pl == s].sum())
                    if not k
                    else float(e["x"][child_hl == s].sum()),
                    density_phi=float(phi[pl == s].mean())
                    if not k and start
                    else 1.0
                    if k
                    else np.nan,
                )
            )
    for h in range(len(c.insect_sizes)):
        total_resource = float(e["x"][child_hl == h].sum())
        for p in range(len(c.plant_sizes)):
            selected = (child_hl == h) & (pl[host] == p)
            resource = float(e["x"][selected].sum())
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
                    successful_recruits=int(np.count_nonzero(selected & survived)),
                    birth_count=int(selected.sum()),
                    contribution_committed=True,
                    edge_semantics="larval_resource_and_recruitment",
                )
            )
        for p in range(-1, len(c.plant_sizes)):
            selected = (hl == h) & (insects.birth_host_species == p)
            par = [hrows[i] for i in np.flatnonzero(selected)]
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
        mothers,
        fathers,
        pl[host],
        child_hl,
        survived,
        t,
        (len(c.plant_sizes), len(c.insect_sizes)),
    )
    for stage, after, parent_groups in (
        ("birth", birth_state, [prows, hrows]),
        (
            "recruitment",
            next_state,
            [
                [
                    dict(r, realized_genetic_contribution=r["realized_recruitment_contribution"])
                    for r in rows
                ]
                for rows in (prows, hrows)
            ],
        ),
    ):
        rows = allele_contributions(
            current, after, parent_groups, genetic, t, (len(c.plant_sizes), len(c.insect_sizes))
        )
        new["allele_contributions"].extend(dict(row, stage=stage) for row in rows)
    return next_state, genetic, [next_p, next_h], new, extinct


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
