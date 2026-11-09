"""Independent reconstruction of saved adult/larval life-history ledgers.

Only configuration validation and generic assertion helpers are shared with the
simulator. Equations, event joins and random draws below use saved raw records.
"""

import json
from functools import partial

import numpy as np

from .lifecycle_checks.events import (
    check_rows,
    close,
    distance,
    host_name,
    read_records,
    unique_keys,
)
from .lifecycle_checks.feeding import (
    chemical_outcomes,
    contact_probability,
    module_choice,
    module_profiles,
)
from .lifecycle_checks.genetics import audit_children, read_genetics
from .lifecycle_checks.population import audit_parents, audit_snapshot_summaries, check_population
from .lifecycle_checks.resources import resource_shares, vector


def audit_run_lifecycle(path, meta):
    """Audit all births, including dead larvae, separately from adult recruits."""
    c, names, tables, states, status = read_records(path, meta)
    if not states:
        assert status == "computational_limit"
        assert all(tables[name].empty for name in names if name != "allele_catalog")
        assert meta.get("stored_genotype_values", 0) == 0
        return dict(
            status=status,
            states=0,
            transitions=0,
            max_overrun=0.0,
            max_balance_error=0.0,
            genetic_audit=True,
            contribution_audit=True,
            contribution_transitions=0,
            alleles_checked=len(tables["allele_catalog"]),
            mutations_checked=0,
            offspring_rows=0,
        )
    spec, loci, catalog, mutations, mutation_lookup, phenotype = read_genetics(
        path, meta, c, tables
    )
    rng, spatial_rng, shape, first, microsites = check_population(
        path, c, states, loci, catalog, phenotype, spec
    )
    tx = tables["gamete_transmissions"]
    retained_values = len(tx) + sum(
        state[kind + "_genotypes"].size for state in states for kind in ("plant", "insect")
    )
    if "stored_genotype_values" in meta:
        assert retained_values == meta["stored_genotype_values"]
        assert retained_values <= spec["max_genotype_values_safety"]
    unique_keys(
        tx, ("generation", "kind", "offspring_id", "role", "copy", "locus_id"), "transmissions"
    )
    tx_lookup = {
        (r.generation, r.kind, r.offspring_id, r.role, r.copy, r.locus_id): r
        for r in tx.itertuples()
    }
    seen_tx, seen_mutations = set(), set()
    calls = {(r["generation"], r["kind"]): r for r in meta["slim_calls"]}
    assert len(calls) == len(meta["slim_calls"])
    seen_calls = set()
    expected = {
        name: {}
        for name in names
        if name not in ("allele_catalog", "mutations", "gamete_transmissions")
    }
    keys = dict(
        events=("generation", "insect_id"),
        larvae=("generation", "insect_id"),
        plants=("generation", "plant_id"),
        plant_proposals=("reproduction_round", "maternal_id", "propagule_index"),
        parents=("generation", "kind", "individual_id"),
        offspring=("generation", "kind", "offspring_id"),
        populations=("generation", "kind", "species"),
        networks=("generation", "insect_species", "plant_species"),
        oviposition=("reproduction_round", "maternal_id", "egg_index"),
        host_flows=(
            "reproduction_round",
            "insect_species",
            "stage",
            "source_host_species",
            "destination_host_species",
        ),
        host_groups=("generation", "insect_species", "birth_host_species"),
        lagged_networks=("reproduction_round", "insect_species", "plant_species"),
        allele_contributions=("generation", "kind", "species", "stage", "locus_id", "allele_id"),
    )

    def add(name, key, **values):
        assert key not in expected[name], ("duplicate expected key", name, key)
        expected[name][key] = values

    transitions = 0
    first_extinct = {}
    next_id = [len(first["plant_ids"]), len(first["insect_ids"])]
    for t, state in enumerate(states):
        pl, hl = state["plant_species"], state["insect_species"]
        pids, hids = state["plant_ids"], state["insect_ids"]
        pcells, hcells = state["plant_cells"], state["insect_cells"]
        for h in range(len(c.insect_sizes)):
            for p in range(-1, len(c.plant_sizes)):
                n = int(np.count_nonzero((hl == h) & (state["insect_birth_host_species"] == p)))
                add(
                    "host_groups",
                    (t, f"H{h + 1}", host_name(p)),
                    N=n,
                    status="present" if n else "absent",
                )
        if t + 1 == len(states):
            continue
        transitions += 1
        D = state["D"]
        chemical = state["plant_chemical_traits"]
        V = chemical[:, 3:5]
        processing_cost = c.insect_cost * state["insect_module_traits"][:, :2].sum(axis=1)
        sensory_cost = c.functional_modules["sensory_cost"] * state["insect_module_traits"][
            :, 4:6
        ].sum(axis=1)
        fi = np.maximum(0, 1 - processing_cost - sensory_cost)
        cue = V[rng[8, 0, 0].permutation(len(V))] if c.shuffle_voc else V
        current_rows, birth_genotypes, child_rows, live_masks = {}, {}, {}, {}

        parents = partial(
            audit_parents,
            hcells=hcells,
            state=state,
            c=c,
            shape=shape,
            rng=rng,
            add=add,
            t=t,
            expected=expected,
            processing_cost=processing_cost,
            sensory_cost=sensory_cost,
        )

        hcounts, hpools, _ = parents("insect", np.ones(len(hids)), fi, np.ones(len(hids)))
        egg_mothers = np.repeat(np.arange(len(hids)), hcounts)
        choice_rng = rng[1, 0, 0]
        hosts = np.full(len(egg_mothers), -1, int)
        ranks = choice_rng.permutation(len(egg_mothers))
        egg_indices = (
            np.concatenate([np.arange(n) for n in hcounts])
            if len(egg_mothers)
            else np.empty(0, int)
        )
        for rank, j in enumerate(ranks):
            m = egg_mothers[j]
            available = np.flatnonzero(
                distance(hcells[m], pcells, shape) <= c.spatial["search_radius"]
            )
            stream = rng[1, 1, int(hl[m])]
            candidates = stream.choice(
                available, min(c.candidate_size, len(available)), replace=False
            )
            reason = "no_candidate"
            visited, acceptance_probabilities, contact_results = [], [], []
            if len(candidates):
                remaining = candidates.copy()
                while len(remaining):
                    probabilities = module_choice(
                        cue[remaining], state["insect_module_traits"][m], c
                    )
                    selected = stream.choice(len(probabilities), p=probabilities)
                    if selected == len(remaining):
                        break
                    plant = int(remaining[selected])
                    probability = float(
                        contact_probability(D[[plant]], state["insect_module_traits"][m], c)[0]
                    )
                    contact_ok = bool(stream.random() < probability)
                    visited.append(plant)
                    acceptance_probabilities.append(probability)
                    contact_results.append(contact_ok)
                    if contact_ok:
                        hosts[j] = plant
                        break
                    remaining = np.delete(remaining, selected)
                reason = "accepted" if hosts[j] >= 0 else "rejected"
            p = hosts[j]
            add(
                "oviposition",
                (t, int(hids[m]), int(egg_indices[j])),
                event_id=f"{t}:{int(hids[m])}:{int(egg_indices[j])}",
                offspring_id=None,
                adult_cohort=t + 1,
                host_id=int(pids[p]) if p >= 0 else -1,
                plant_species=host_name(pl[p]) if p >= 0 else "unknown",
                accepted=bool(p >= 0),
                reason=reason,
                candidate_count=len(candidates),
                choice_order=rank,
                candidate_plant_ids=json.dumps(pids[candidates].tolist()),
                contact_plant_ids=json.dumps(pids[np.asarray(visited, int)].tolist()),
                contact_accepted=json.dumps(contact_results),
            )
            recorded = tables["oviposition"]
            recorded = recorded[
                (recorded.reproduction_round == t)
                & (recorded.maternal_id == hids[m])
                & (recorded.egg_index == egg_indices[j])
            ]
            assert len(recorded) == 1
            close(
                json.loads(recorded.iloc[0].contact_probabilities),
                acceptance_probabilities,
                "contact probabilities",
            )
        accepted = hosts >= 0
        for i, j in enumerate(np.flatnonzero(accepted)):
            key = (t, int(hids[egg_mothers[j]]), int(egg_indices[j]))
            expected["oviposition"][key]["offspring_id"] = next_id[1] + i
        mothers, hosts = egg_mothers[accepted], hosts[accepted]
        fathers = np.array([rng[5, 1, int(hl[m])].choice(hpools[m]) for m in mothers], int)

        children = partial(
            audit_children,
            tables=tables,
            t=t,
            next_id=next_id,
            loci=loci,
            phenotype=phenotype,
            catalog=catalog,
            calls=calls,
            seen_calls=seen_calls,
            c=c,
            meta=meta,
            add=add,
            tx_lookup=tx_lookup,
            seen_tx=seen_tx,
            mutation_lookup=mutation_lookup,
            spec=spec,
            seen_mutations=seen_mutations,
            expected=expected,
            birth_genotypes=birth_genotypes,
            child_rows=child_rows,
            state=state,
        )

        hg, hchildren = children("insect", mothers, fathers)
        child_hl = hl[mothers]
        child_I = module_profiles(phenotype(hg, "insect"), c)[0]
        R = state["plant_resource"]
        K = vector(c.resource_capacity, pl, "P")
        dmax = vector(c.max_feeding_demand, child_hl, "H")
        demand = dmax.copy()
        total, gamma, x, X = resource_shares(hosts, demand, R)
        h = np.divide(X, R, out=np.zeros_like(X, dtype=float), where=R > 0)
        g = x / dmax
        exposure = chemical_outcomes(D[hosts], child_I, x, dmax, c)
        q = chemical_outcomes(D[hosts], child_I, dmax, dmax, c)["survival_probability"]
        live = np.zeros(len(hchildren), bool)
        sex = np.empty(len(hchildren), np.int8)
        for s in range(len(c.insect_sizes)):
            selected = child_hl == s
            live[selected] = (
                rng[2, 1, s].random(selected.sum()) < exposure["survival_probability"][selected]
            )
            sex[selected] = (
                rng[7, 1, s].random(selected.sum())
                >= c.sexual_reproduction["insect_offspring_female_probability"]
            )
        adult_cells = (
            pcells[hosts[live]]
            + spatial_rng[1].integers(
                -c.spatial["insect_dispersal_radius"],
                c.spatial["insect_dispersal_radius"] + 1,
                size=(live.sum(), 2),
            )
        ) % shape
        origins = pcells[hosts]
        accepted_eggs = np.flatnonzero(accepted)
        for i, row in enumerate(hchildren.itertuples()):
            p, m, f = hosts[i], mothers[i], fathers[i]
            a, b = state["insect_birth_host_species"][[m, f]]
            values = dict(
                birth_host_id=int(pids[p]),
                birth_host_species=host_name(pl[p]),
                maternal_birth_host_species=host_name(a),
                paternal_birth_host_species=host_name(b),
                cross_background_pair=None if min(a, b) < 0 else bool(a != b),
                pair_background_status="unknown"
                if min(a, b) < 0
                else "same"
                if a == b
                else "different",
                survived_to_adult=bool(live[i]),
                offspring_sex=int(sex[i]),
            )
            expected["offspring"][t + 1, "insect", row.offspring_id].update(values)
            for field, value in (
                ("maternal_cell", hcells[m]),
                ("paternal_cell", hcells[f]),
                ("offspring_cell", origins[i]),
            ):
                close(
                    json.loads(getattr(row, field)),
                    value,
                    ("offspring location", t, row.offspring_id, field),
                )
            egg = accepted_eggs[i]
            ov = tables["oviposition"]
            ov = ov[
                (ov.reproduction_round == t)
                & (ov.maternal_id == row.maternal_id)
                & (ov.egg_index == egg_indices[egg])
            ]
            assert len(ov) == 1
            assert row.oviposition_event_id == ov.iloc[0].event_id
            assert ov.iloc[0].offspring_id == row.offspring_id
            expected["offspring"][t + 1, "insect", row.offspring_id]["egg_index"] = int(
                egg_indices[egg]
            )
            same = int(np.count_nonzero((hosts == p) & (child_hl == child_hl[i]))) - 1
            other = int(np.count_nonzero(hosts == p)) - same - 1
            event = dict(
                insect_row=i,
                seed=c.seed,
                accepted=True,
                reason="accepted",
                candidate_count=expected["oviposition"][t, row.maternal_id, int(egg_indices[egg])][
                    "candidate_count"
                ],
                choice_order=expected["oviposition"][t, row.maternal_id, int(egg_indices[egg])][
                    "choice_order"
                ],
                reproduction_round=t,
                adult_cohort=t + 1,
                species=f"H{child_hl[i] + 1}",
                species_id=int(child_hl[i]),
                maternal_id=row.maternal_id,
                paternal_id=row.paternal_id,
                host_row=int(p),
                host_id=int(pids[p]),
                plant_species=host_name(pl[p]),
                birth_host_id=int(pids[p]),
                birth_host_species=host_name(pl[p]),
                q=q[i],
                d=demand[i],
                x=x[i],
                g=g[i],
                dmax=dmax[i],
                R=R[p],
                gamma=gamma[p],
                host_h=h[p],
                I_sum=child_I[i].sum(),
                W_insect=exposure["survival_probability"][i],
                survived_to_adult=bool(live[i]),
                insect_life_stage="larva",
                same_species_others=same,
                other_species_others=other,
                oviposition_event_id=row.oviposition_event_id,
                egg_index=int(egg_indices[egg]),
            )
            for name in ("dose", "processing_capacity", "processed", "residual", "antinutrient"):
                for axis in range(2):
                    event[f"{name}_{axis + 1}"] = exposure[name][i, axis]
            for name in ("nutrition_survival", "toxicity_survival", "antinutrition_survival"):
                event[name] = exposure[name][i]
            add("events", (t, row.offspring_id), **event)
            add("larvae", (t, row.offspring_id), **event)
        live_masks["insect"] = live
        fp = np.maximum(
            0, 1 - c.plant_cost * D.sum(axis=1) - c.voc_cost * chemical[:, 3:5].sum(axis=1)
        )
        post = np.maximum(0, R - X)
        defoliation = np.divide(X, R, out=np.zeros_like(X, dtype=float), where=R > 0)
        spec_plant = c.plant_demography
        Sp = np.exp(
            -spec_plant["background_mortality_rate"]
            - spec_plant["defoliation_mortality_rate"] * defoliation
        )
        recovery_rate = (
            spec_plant["resource_recovery_rate"]
            + spec_plant["tolerance_recovery_rate"] * chemical[:, 2]
        )
        replenished = K - (K - post) * np.exp(-recovery_rate)
        pcounts, ppools, psurvived = parents("plant", Sp, fp * replenished / K, np.ones(len(pids)))
        proposal_mothers = np.repeat(np.arange(len(pids)), pcounts)
        proposal_indices = (
            np.concatenate([np.arange(n) for n in pcounts])
            if len(proposal_mothers)
            else np.empty(0, int)
        )
        proposal_cells = (
            pcells[proposal_mothers]
            + spatial_rng[0].integers(
                -c.spatial["plant_dispersal_radius"],
                c.spatial["plant_dispersal_radius"] + 1,
                size=(len(proposal_mothers), 2),
            )
        ) % shape
        occupied = set(map(int, state["plant_microsite"][psurvived]))
        proposal_sites = np.full(len(proposal_mothers), -1, int)
        for i in rng[9, 0, 0].permutation(len(proposal_mothers)):
            cell = proposal_cells[i]
            start = int(cell[0] * shape[1] + cell[1]) * microsites
            vacancy = next(
                (slot for slot in range(start, start + microsites) if slot not in occupied), None
            )
            if vacancy is not None:
                proposal_sites[i] = vacancy
                occupied.add(vacancy)
        established = np.flatnonzero(proposal_sites >= 0)
        pm = proposal_mothers[established]
        pf = np.array([rng[5, 0, int(pl[m])].choice(ppools[m]) for m in pm], int)
        _, pchildren = children("plant", pm, pf)
        child_pcells = proposal_cells[established]
        child_sites = proposal_sites[established]
        child_ids = dict(zip(established, pchildren.offspring_id))
        for i, mother in enumerate(proposal_mothers):
            add(
                "plant_proposals",
                (t, int(pids[mother]), int(proposal_indices[i])),
                event_id=f"{t}:{int(pids[mother])}:{int(proposal_indices[i])}",
                offspring_generation=t + 1,
                cell_row=int(proposal_cells[i, 0]),
                cell_col=int(proposal_cells[i, 1]),
                established=bool(proposal_sites[i] >= 0),
                microsite=int(proposal_sites[i]) if proposal_sites[i] >= 0 else None,
                offspring_id=int(child_ids[i]) if i in child_ids else None,
            )
        recovery = replenished - post
        live_masks["plant"] = np.ones(len(pchildren), bool)
        for i, ident in enumerate(pids):
            add(
                "plants",
                (t, int(ident)),
                species=host_name(pl[i]),
                R=R[i],
                T=total[i],
                gamma=gamma[i],
                X=X[i],
                h=h[i],
                herbivore_load=int(np.count_nonzero(hosts == i)),
                D_sum=D[i].sum(),
                V_emit_mean=V[i].mean(),
                defense_cost=c.plant_cost * D[i].sum(),
                voc_cost=c.voc_cost * chemical[i, 3:5].sum(),
                W_plant=Sp[i] * fp[i] * replenished[i] / K[i],
                resource_capacity=K[i],
                resource_after_feeding=post[i],
                tolerance_Z=chemical[i, 2],
                recovery_rate=recovery_rate[i],
                resource_after_recovery=replenished[i] if psurvived[i] else None,
                defoliation=defoliation[i],
                survival_probability=Sp[i],
                survived=bool(psurvived[i]),
                resource_recovery=recovery[i] if psurvived[i] else None,
                resource_next=replenished[i] if psurvived[i] else None,
                microsite=int(state["plant_microsite"][i]),
                reproduction_active=True,
                cell_row=int(pcells[i, 0]),
                cell_col=int(pcells[i, 1]),
            )
        for kind, mm, ff, survivor, live_parents in (
            ("insect", mothers, fathers, live, np.ones(len(hids), bool)),
            ("plant", pm, pf, np.ones(len(pm), bool), psurvived),
        ):
            k = int(kind == "insect")
            ids, labels = state[kind + "_ids"], state[kind + "_species"]
            bm = np.bincount(mm, minlength=len(ids))
            bf = np.bincount(ff, minlength=len(ids))
            rm, rf = (
                np.bincount(mm[survivor], minlength=len(ids)),
                np.bincount(ff[survivor], minlength=len(ids)),
            )
            weights, recruited = (bm + bf) / 2, (rm + rf) / 2
            current_rows[kind] = (bm, bf, weights, rm, rf, recruited)
            for i, ident in enumerate(ids):
                expected["parents"][t, kind, int(ident)].update(
                    offspring_count_as_mother=int(bm[i]),
                    offspring_count_sired=int(bf[i]),
                    recruited_count_as_mother=int(rm[i]),
                    recruited_count_sired=int(rf[i]),
                    realized_genetic_contribution=weights[i],
                    realized_recruitment_contribution=recruited[i],
                )
            nxt = states[t + 1]
            rows = child_rows[kind]
            retained = np.flatnonzero(psurvived) if k == 0 else np.empty(0, int)
            for suffix, newborn_values in (
                ("_ids", rows.offspring_id.to_numpy()[survivor]),
                ("_species", labels[mm][survivor]),
                ("_genotypes", birth_genotypes[kind][survivor]),
                ("_cells", adult_cells if k else child_pcells),
            ):
                wanted = (
                    newborn_values
                    if k
                    else np.concatenate((state[kind + suffix][retained], newborn_values))
                )
                np.testing.assert_array_equal(nxt[kind + suffix], wanted)
            if k == 0:
                np.testing.assert_array_equal(
                    nxt["plant_microsite"],
                    np.concatenate((state["plant_microsite"][retained], child_sites)),
                )
                close(
                    nxt["plant_resource"],
                    np.concatenate(
                        (
                            replenished[retained],
                            vector(c.resource_capacity, labels[mm], "P")
                            * spec_plant["newborn_resource_fraction"],
                        )
                    ),
                    "next plant stock",
                )
                close(
                    nxt["plant_chemical_traits"][: len(retained)],
                    chemical[retained],
                    "retained chemical traits",
                )
            if k:
                for name, values in (
                    ("insect_sex", sex),
                    ("insect_birth_host_id", pids[hosts]),
                    ("insect_birth_host_species", pl[hosts]),
                    ("insect_adult_cohort", np.full(len(rows), t + 1)),
                ):
                    np.testing.assert_array_equal(nxt[name], values[survivor])
            else:
                for i, row in enumerate(rows.itertuples()):
                    expected["offspring"][t + 1, kind, row.offspring_id]["survived_to_adult"] = True
                    for field, value in (
                        ("maternal_cell", pcells[mm[i]]),
                        ("paternal_cell", pcells[ff[i]]),
                        ("offspring_cell", child_pcells[i]),
                    ):
                        close(
                            json.loads(getattr(row, field)),
                            value,
                            ("plant child location", row.offspring_id, field),
                        )
            sizes = c.insect_sizes if k else c.plant_sizes
            for s in range(len(sizes)):
                name = f"{'H' if k else 'P'}{s + 1}"
                selected = labels == s
                child_selected = labels[mm] == s
                n, born, recruited_next = (
                    int(selected.sum()),
                    int(child_selected.sum()),
                    int(np.count_nonzero(child_selected & survivor)),
                )
                if k:
                    assert np.isclose(weights[selected].sum(), born) and np.isclose(
                        recruited[selected].sum(), recruited_next
                    )
                n_retained = int(live_parents[selected].sum()) if k == 0 else 0
                n_next = recruited_next + n_retained
                if n == 0:
                    first_extinct.setdefault(name, t)
                elif n_next == 0:
                    first_extinct.setdefault(name, t + 1)
                add(
                    "populations",
                    (t, kind, name),
                    N_start=n,
                    N_offspring=born,
                    N_next=n_next,
                    N_survived_to_reproduce=int(live_parents[selected].sum()),
                    update_applied=True,
                    population_status="extant" if n else "extinct",
                    extinct_at_generation=first_extinct.get(name),
                    total_resource_transfer=float(x[child_hl == s].sum())
                    if k
                    else float(X[pl == s].sum()),
                    N_retained=n_retained,
                    N_propagules=int(pcounts[selected].sum()) if k == 0 else None,
                )
                if k:
                    expected["populations"][t, kind, name]["N_potential_eggs"] = int(
                        hcounts[selected].sum()
                    )
                for stage, chosen, ww in (
                    ("birth", np.ones(len(rows), bool), weights),
                    ("recruitment", survivor, recruited),
                ):
                    before = state[kind + "_genotypes"][selected]
                    after = birth_genotypes[kind][child_selected & chosen]
                    count = len(after)
                    if not k:
                        children_for_label = child_selected & chosen
                        ww = (
                            np.bincount(mm[children_for_label], minlength=len(ids))
                            + np.bincount(ff[children_for_label], minlength=len(ids))
                        ) / 2
                    for l, locus in enumerate(loci[kind]):
                        alleles = np.union1d(before[:, :, l], after[:, :, l])
                        if not k:
                            alleles = np.union1d(alleles, state[kind + "_genotypes"][ww > 0, :, l])
                        for allele in alleles:
                            z = (before[:, :, l] == allele).sum(axis=1) / 2
                            p0 = float(z.mean()) if n else np.nan
                            if k:
                                pw = float(ww[selected] @ z / count) if count else np.nan
                            else:
                                all_z = (state[kind + "_genotypes"][:, :, l] == allele).sum(
                                    axis=1
                                ) / 2
                                pw = float(ww @ all_z / count) if count else np.nan
                            p1 = float(np.mean(after[:, :, l] == allele)) if count else np.nan
                            add(
                                "allele_contributions",
                                (t, kind, name, stage, locus["locus_id"], int(allele)),
                                offspring_generation=t + 1,
                                N_parent=n,
                                N_offspring=count,
                                p_parent=p0,
                                p_parent_weighted=pw,
                                p_offspring=p1,
                                delta_parent_contribution=pw - p0,
                                delta_transmission=p1 - pw,
                            )
        for s in range(len(c.insect_sizes)):
            species = f"H{s + 1}"
            bm, bf, weights, rm, rf, rw = current_rows["insect"]
            total_resource = x[child_hl == s].sum()
            for p in range(len(c.plant_sizes)):
                selected = (child_hl == s) & (pl[hosts] == p)
                resource = float(x[selected].sum())
                add(
                    "networks",
                    (t, species, host_name(p)),
                    reproduction_round=t,
                    adult_cohort=t + 1,
                    N_accept=int(selected.sum()),
                    birth_count=int(selected.sum()),
                    Q_actual=resource,
                    resource_share=resource / total_resource if total_resource else np.nan,
                    L_per_insect=resource / selected.sum() if selected.any() else np.nan,
                    successful_recruits=int(np.count_nonzero(selected & live)),
                    contribution_committed=True,
                )
            for p in range(-1, len(c.plant_sizes)):
                selected = (hl == s) & (state["insect_birth_host_species"] == p)
                add(
                    "lagged_networks",
                    (t, species, host_name(p)),
                    adult_cohort=t + 1,
                    maternal_offspring=int(bm[selected].sum()),
                    paternal_offspring=int(bf[selected].sum()),
                    genetic_contribution=weights[selected].sum(),
                    recruitment_contribution=rw[selected].sum(),
                )
                for destination in range(-1, len(c.plant_sizes)):
                    for stage, chosen in (
                        ("birth", np.ones(len(live), bool)),
                        ("recruitment", live),
                    ):
                        child_selected = (child_hl == s) & (pl[hosts] == destination) & chosen
                        amount = (
                            np.count_nonzero(
                                child_selected & (state["insect_birth_host_species"][mothers] == p)
                            )
                            + np.count_nonzero(
                                child_selected & (state["insect_birth_host_species"][fathers] == p)
                            )
                        ) / 2
                        add(
                            "host_flows",
                            (t, species, stage, host_name(p), host_name(destination)),
                            adult_cohort=t + 1,
                            genetic_contribution=amount,
                        )
    audit_snapshot_summaries(path, c, states, catalog, phenotype)
    assert seen_tx == set(tx_lookup), ("orphan/missing transmission", set(tx_lookup) ^ seen_tx)
    assert seen_mutations == set(mutation_lookup), "orphan/missing mutation"
    assert seen_calls == set(calls), "orphan/missing SLiM call"
    for name in expected:
        check_rows(tables[name], keys[name], expected[name], name)
    assert states[-1]["generation"] == meta["status"]["generation"]
    if status == "completed":
        assert transitions == c.generations and len(states) == c.generations + 1
    elif status == "all_species_extinct":
        assert len(states[-1]["plant_ids"]) + len(states[-1]["insect_ids"]) == 0
    else:
        assert status in ("computational_limit", "backend_error", "arithmetic_error")
        assert meta["status"].get("message")
    return dict(
        status=status,
        states=len(states),
        transitions=transitions,
        event_rows=len(tables["events"]),
        plant_rows=len(tables["plants"]),
        offspring_rows=len(tables["offspring"]),
        max_overrun=0.0,
        max_balance_error=0.0,
        genetic_audit=True,
        contribution_audit=True,
        contribution_transitions=transitions,
        allele_contribution_rows=len(tables["allele_contributions"]),
        alleles_checked=len(catalog),
        mutations_checked=len(mutations),
        stored_genotype_values=retained_values,
    )
