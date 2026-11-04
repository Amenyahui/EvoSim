"""Read-only matched mating opportunities and observed QTL-copy ancestry."""

import json
from collections import Counter, defaultdict

import numpy as np

UNKNOWN = "unknown"
UNASSIGNED = "unassigned"
WEIGHTS = ("conditional_mu", "potential_egg_count", "offspring_count_as_mother")


def host_groups(state):
    values = state.get("insect_birth_host_species")
    if values is None:
        return None
    return {
        int(i): f"P{int(h) + 1}" if h >= 0 else UNKNOWN for i, h in zip(state["insect_ids"], values)
    }


def known(label):
    return label not in (UNKNOWN, UNASSIGNED)


def flag(value):
    if value in (True, "True", 1):
        return True
    if value in (False, "False", 0):
        return False
    raise ValueError("invalid maternal eligibility flag")


def opportunities(states, parents, mappings):
    required = {
        "generation",
        "kind",
        "individual_id",
        "eligible_as_mother",
        "candidate_mate_ids",
        "candidate_mate_probabilities",
        *WEIGHTS,
    }
    rows = []
    for t in range(max(states)):
        state = states[t]
        host = host_groups(state)
        if not required <= set(parents) or "insect_sex" not in state or host is None:
            rows.append(dict(reproduction_round=t, status="not_saved"))
            continue
        ids = state["insect_ids"]
        pools = dict(zip(map(int, ids), map(int, state["insect_species"])))
        males = [int(i) for i, sex in zip(ids, state["insect_sex"]) if sex == 1]
        ledger = parents[(parents.generation == t) & (parents.kind == "insect")]
        if ledger.individual_id.duplicated().any() or set(ledger.individual_id) != set(ids):
            raise ValueError("parent ledger does not match snapshot")
        females = {int(i) for i, sex in zip(ids, state["insect_sex"]) if sex == 0}
        for system, labels in (("birth_host", host), ("candidate_cluster", mappings[t, "insect"])):
            accum = defaultdict(lambda: defaultdict(float))
            for mother in ledger.to_dict("records"):
                if not flag(mother["eligible_as_mother"]):
                    continue
                ident = int(mother["individual_id"])
                if ident not in females:
                    raise ValueError("eligible insect mother is not female")
                local = json.loads(mother["candidate_mate_ids"])
                probs = np.asarray(json.loads(mother["candidate_mate_probabilities"]), dtype=float)
                global_pool = [i for i in males if pools[i] == pools[ident]]
                if (
                    len(local) != len(probs)
                    or len(set(local)) != len(local)
                    or not np.isfinite(probs).all()
                    or any(i not in global_pool for i in local)
                    or (len(local) and not np.allclose(probs, 1 / len(local)))
                ):
                    raise ValueError("invalid local male pool or sampling probabilities")
                source = labels.get(ident, UNKNOWN)
                unknown_local = sum(not known(labels.get(i, UNKNOWN)) for i in local)
                unknown_global = sum(not known(labels.get(i, UNKNOWN)) for i in global_pool)
                complete = bool(local) and known(source) and unknown_global == 0
                targets = sorted(
                    {
                        labels[i]
                        for i in pools
                        if pools[i] == pools[ident] and known(labels.get(i, UNKNOWN))
                    }
                    | ({source} if known(source) else set())
                )
                # An all-cross row permits more than two groups without conflating direction.
                targets.append("__all_cross__")
                for weight_name in WEIGHTS:
                    weight = float(mother[weight_name])
                    if not np.isfinite(weight) or weight < 0:
                        raise ValueError("invalid maternal weight")
                    for target in targets:
                        key = (system, pools[ident], source, target, weight_name)
                        a = accum[key]
                        a["mothers"] += 1
                        a["weight"] += weight
                        a["no_mate_mothers"] += not local
                        a["no_mate_weight"] += weight if not local else 0
                        a["unknown_mother_count"] += not known(source)
                        a["unknown_mother_weight"] += weight if not known(source) else 0
                        a["unknown_local_edges"] += unknown_local
                        a["local_edges"] += len(local)
                        a["unknown_global_edges"] += unknown_global
                        a["global_edges"] += len(global_pool)
                        if not complete:
                            a["excluded_mothers"] += 1
                            a["excluded_weight"] += weight
                            continue
                        a["matched_mothers"] += 1
                        a["matched_weight"] += weight

                        def match(i):
                            return (
                                labels[i] != source
                                if target == "__all_cross__"
                                else labels[i] == target
                            )

                        a["local_sum"] += weight * sum(match(i) for i in local) / len(local)
                        a["global_sum"] += (
                            weight * sum(match(i) for i in global_pool) / len(global_pool)
                        )
            if not accum:
                rows.append(
                    dict(reproduction_round=t, group_system=system, status="no_eligible_mothers")
                )
            for (system, pool, source, target, weight_name), a in sorted(accum.items()):
                w = a["matched_weight"]
                p_local = a["local_sum"] / w if w else np.nan
                p_global = a["global_sum"] / w if w else np.nan
                rows.append(
                    dict(
                        reproduction_round=t,
                        group_system=system,
                        compatibility_pool=pool,
                        mother_group=source,
                        target_group=target,
                        weight_type=weight_name,
                        **{
                            k: a[k]
                            for k in (
                                "mothers",
                                "weight",
                                "no_mate_mothers",
                                "no_mate_weight",
                                "unknown_mother_count",
                                "unknown_mother_weight",
                                "unknown_local_edges",
                                "local_edges",
                                "unknown_global_edges",
                                "global_edges",
                                "excluded_mothers",
                                "excluded_weight",
                                "matched_mothers",
                                "matched_weight",
                            )
                        },
                        p_local=p_local,
                        p_global=p_global,
                        ri_relative=1 - p_local / p_global
                        if p_global > 0 and target != source
                        else np.nan,
                        status="available",
                    )
                )
    return rows


def ancestry(states, offspring, mappings, reference, loci, transmissions):
    if type(reference) is not int or reference not in states:
        raise ValueError("ancestry reference generation must be an available snapshot")
    loci = loci[loci.organism == "insect"].locus_id.tolist()
    if not loci or len(set(loci)) != len(loci):
        raise ValueError("insect QTL locus list is missing or duplicated")
    required = {
        "generation",
        "kind",
        "offspring_id",
        "role",
        "parent_id",
        "locus_id",
        "inherited_allele_id",
        "child_allele_id",
        "transmitted_homolog_before_mutation",
    }
    if not required <= set(transmissions) or not {"genotype", "birth_host_species"} <= set(
        offspring
    ):
        raise ValueError(
            "QTL ancestry requires complete transmission and offspring genotype ledgers"
        )
    tx = transmissions[transmissions.kind == "insect"]
    rows = []
    ref = states[reference]
    host = host_groups(ref)
    if host is None:
        raise ValueError("QTL ancestry requires saved reference birth hosts")
    systems = {"birth_host": host, "candidate_cluster": mappings[reference, "insect"]}
    for system, sources in systems.items():
        source_labels = set(sources.values()) | {UNKNOWN, UNASSIGNED}
        current = {
            int(i): np.full((2, len(loci)), sources.get(int(i), UNKNOWN), dtype=object)
            for i in ref["insect_ids"]
        }

        def emit(t, stage, labels, arrays):
            counts = defaultdict(Counter)
            for ident, array in arrays.items():
                target = labels.get(ident, UNKNOWN)
                counts[target].update(array.ravel())
            if not counts:
                rows.append(
                    dict(
                        reference_generation=reference,
                        generation=t,
                        k=t - reference,
                        source_system=system,
                        target_system="birth_host" if stage == "birth" else "candidate_cluster",
                        target_group=UNKNOWN,
                        target_lineage_id=UNASSIGNED,
                        stage=stage,
                        qtl_copies=0,
                        total_qtl_copies=0,
                        fraction=np.nan,
                        unknown_fraction=np.nan,
                        status="empty_population",
                    )
                )
            for target, counts_by_source in sorted(counts.items()):
                total = sum(counts_by_source.values())
                for source in sorted(source_labels):
                    rows.append(
                        dict(
                            reference_generation=reference,
                            generation=t,
                            k=t - reference,
                            source_system=system,
                            source_group=source,
                            target_system="birth_host" if stage == "birth" else "candidate_cluster",
                            target_group=target,
                            target_lineage_id=UNASSIGNED,
                            stage=stage,
                            qtl_copies=counts_by_source[source],
                            total_qtl_copies=total,
                            fraction=counts_by_source[source] / total,
                            unknown_fraction=sum(counts_by_source[s] for s in (UNKNOWN, UNASSIGNED))
                            / total,
                            status="available",
                        )
                    )

        emit(reference, "reference", mappings[reference, "insect"], current)
        for t in range(reference, max(states)):
            state = states[t]
            parent_genotypes = {
                int(i): g for i, g in zip(state["insect_ids"], state["insect_genotypes"])
            }
            births = offspring[(offspring.kind == "insect") & (offspring.reproduction_round == t)]
            records = tx[tx.generation == t + 1]
            if len(records) != len(births) * 2 * len(loci):
                raise ValueError("incomplete QTL transmission records")
            keyed = {}
            for record in records.to_dict("records"):
                key = (int(record["offspring_id"]), record["role"], record["locus_id"])
                if key in keyed:
                    raise ValueError("duplicate QTL transmission")
                keyed[key] = record
            born = {}
            birth_hosts = {}
            for child in births.to_dict("records"):
                ident = int(child["offspring_id"])
                genotype = np.asarray(json.loads(child["genotype"]))
                if genotype.shape != (2, len(loci)):
                    raise ValueError("invalid offspring QTL genotype")
                array = np.empty((2, len(loci)), dtype=object)
                for h, role in enumerate(("maternal", "paternal")):
                    parent = int(child[role + "_id"])
                    if parent not in current:
                        raise ValueError("missing QTL ancestry parent")
                    for l, locus in enumerate(loci):
                        record = keyed.pop((ident, role, locus), None)
                        if record is None:
                            raise ValueError("missing QTL transmission")
                        homolog = int(record["transmitted_homolog_before_mutation"])
                        if (
                            homolog not in (0, 1)
                            or int(record["parent_id"]) != parent
                            or int(record["inherited_allele_id"])
                            != parent_genotypes[parent][homolog, l]
                            or int(record["child_allele_id"]) != genotype[h, l]
                        ):
                            raise ValueError("inconsistent QTL transmission")
                        array[h, l] = current[parent][homolog, l]
                born[ident] = array
                birth_hosts[ident] = child.get("birth_host_species", UNKNOWN)
            if keyed:
                raise ValueError("unexpected QTL transmission")
            emit(t + 1, "birth", birth_hosts, born)
            next_ids = list(map(int, states[t + 1]["insect_ids"]))
            current = {i: born[i] for i in next_ids}
            child_genotypes = {
                int(row["offspring_id"]): np.asarray(json.loads(row["genotype"]))
                for row in births.to_dict("records")
            }
            if any(
                not np.array_equal(child_genotypes[i], g)
                for i, g in zip(next_ids, states[t + 1]["insect_genotypes"])
            ):
                raise ValueError("recruited QTL genotype disagrees with snapshot")
            emit(t + 1, "recruitment", mappings[t + 1, "insect"], current)
    return rows
