"""Read-only operational QTL clusters, ancestry events and stage-explicit event networks.

Run with --input <batch> --design <observation.yaml>. Outputs are confined to the
batch's data/summary/lineages directory. These clusters are not confirmed species.
"""

import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from evo_sim.core.io import atomic_csv

from .common import Batch, audit_file, batch_parser, read_settings, run_path, sha256, write_json

UNASSIGNED = "unassigned"
UNKNOWN = "unknown"
KINDS = ("plant", "insect")
TABLE_COLUMNS = {
    "memberships": "generation kind individual_id compatibility_pool cluster_id lineage_id component_size continuity_generations persistent_candidate",
    "clusters": "generation kind cluster_id compatibility_pool size lineage_id continuity_generations persistent_candidate",
    "transitions": "kind source_generation target_generation source_cluster target_cluster recruited_genetic_contribution",
    "lineage_events": "generation kind event source_clusters target_cluster reason event_scope",
    "populations": "generation kind compatibility_pool individual_count detected_clusters pool_state extinction_observed",
    "mating": "reproduction_round offspring_cohort parent_generation kind stage offspring_weighted_pairs classified_parent_pairs cross_cluster_parent_pairs unknown_parent_pairs unclassified_parent_pairs cross_fraction_of_classified opportunity_status",
    "mating_opportunities": "reproduction_round kind eligible_mothers classified_mothers candidate_edges cross_cluster_candidate_edges classified_probability_mass cross_probability_mass probability_weighted_cross_fraction opportunity_status",
    "networks": "reproduction_round insect_generation plant_generation edge_semantics insect_cluster plant_cluster actual_consumption birth_count recruitment_count",
    "network_balances": "reproduction_round birth_count recruitment_count actual_consumption maternal_network_births maternal_network_recruits maternal_network_consumption recruit_network_recruits",
}


def read_design(design):
    """Validate the three observation-only thresholds; no implicit defaults."""
    values = (
        dict(design)
        if isinstance(design, dict)
        else yaml.safe_load(Path(design).read_text(encoding="utf-8"))
    )
    keys = {"distance_threshold", "minimum_cluster_size", "persistence_generations"}
    if not isinstance(values, dict) or set(values) != keys:
        raise ValueError(f"lineage design requires exactly {sorted(keys)}")
    distance = values["distance_threshold"]
    if type(distance) not in (float, int) or not np.isfinite(distance) or not 0 <= distance <= 1:
        raise ValueError("distance_threshold must be finite in [0, 1]")
    for key in keys - {"distance_threshold"}:
        if type(values[key]) is not int or values[key] < (
            2 if key == "minimum_cluster_size" else 1
        ):
            raise ValueError(f"invalid {key}")
    return values


def cluster_genotypes(ids, pools, genotypes, *, kind, generation, design):
    """Threshold graph components using unphased allele-copy mismatch in [0, 1]."""
    ids, pools, g = np.asarray(ids), np.asarray(pools), np.asarray(genotypes)
    if (
        ids.ndim != 1
        or ids.dtype.kind not in "iu"
        or np.any(ids < 0)
        or pools.shape != ids.shape
        or pools.dtype.kind not in "iu"
        or np.any(pools < 0)
        or len(set(ids.tolist())) != len(ids)
        or g.ndim != 3
        or g.shape[:2] != (len(ids), 2)
        or not g.shape[2]
        or g.dtype.kind not in "iu"
        or np.any(g < 0)
    ):
        raise ValueError("invalid snapshot identities or diploid QTL genotypes")
    root = np.arange(len(ids))

    def find(i):
        while root[i] != i:
            root[i] = root[root[i]]
            i = root[i]
        return i

    # ponytail: O(N²L) time, O(NL) working memory; partition/index if large batches need it.
    for i in range(len(ids)):
        candidates = np.flatnonzero((pools == pools[i]) & (np.arange(len(ids)) > i))
        a, b = g[i], g[candidates]
        direct = (b[:, 0] == a[0]).astype(int) + (b[:, 1] == a[1])
        swapped = (b[:, 0] == a[1]).astype(int) + (b[:, 1] == a[0])
        distance = 1 - np.maximum(direct, swapped).sum(axis=1) / (2 * g.shape[2])
        for j in candidates[distance <= design["distance_threshold"]]:
            root[find(int(j))] = find(i)
    components = defaultdict(list)
    for i in range(len(ids)):
        components[find(i)].append(i)
    assignment, groups, sizes = {}, {}, {}
    for members in sorted(components.values(), key=lambda rows: min(ids[rows])):
        label = f"{kind}:g{generation:06d}:c{int(min(ids[members]))}"
        if len(members) < design["minimum_cluster_size"]:
            label = UNASSIGNED
        else:
            groups[label] = dict(
                compatibility_pool=int(pools[members[0]]), size=len(members), members=members
            )
        for i in members:
            assignment[int(ids[i])] = label
            sizes[int(ids[i])] = len(members)
    return assignment, groups, sizes


def _flag(value):
    if value in (True, "True", 1):
        return True
    if value in (False, "False", 0):
        return False
    raise ValueError(f"invalid survived_to_adult value {value!r}")


def _source(mapping, ident):
    return mapping.get(int(ident), UNKNOWN) if pd.notna(ident) else UNKNOWN


def _summarize_run(folder, design):
    """Return tables from one saved run without calling simulation functions."""
    metadata = json.loads((folder / "metadata.json").read_text(encoding="utf-8"))
    if metadata["version"] != "0.1.7" or metadata.get("genetics_version") not in {
        "slim-qtl-3",
        "slim-qtl-4",
    }:
        raise ValueError("lineages require v0.1.7 / slim-qtl-3 or slim-qtl-4 saved ledgers")
    hashes = {}

    def checked(name):
        path = folder / name
        digest = sha256(path)
        if metadata["output_sha256"].get(name) != digest:
            raise ValueError(f"saved input hash mismatch: {name}")
        hashes[name] = digest
        return path

    states = {}
    for file in sorted((folder / "snapshots").glob("g*.npz")):
        with np.load(checked(file.relative_to(folder).as_posix()), allow_pickle=False) as saved:
            state = {key: saved[key] for key in saved.files}
        t = int(state["generation"])
        if t in states:
            raise ValueError("duplicate snapshot generation")
        states[t] = state
    if not states or sorted(states) != list(range(max(states) + 1)):
        raise ValueError("lineage tracking requires consecutive snapshots starting at generation 0")
    offspring = pd.read_csv(checked("offspring.csv"), float_precision="round_trip")
    events = pd.read_csv(checked("events.csv"), float_precision="round_trip")
    parents = pd.read_csv(checked("parents.csv"), float_precision="round_trip")
    required = {
        "kind",
        "offspring_id",
        "maternal_id",
        "paternal_id",
        "reproduction_round",
        "adult_cohort",
        "survived_to_adult",
        "birth_host_id",
    }
    if not required <= set(offspring):
        raise ValueError(f"offspring ledger missing {sorted(required - set(offspring))}")
    if offspring.duplicated(["kind", "offspring_id"]).any():
        raise ValueError("duplicate offspring identity")
    if not offspring.empty and (
        not offspring.adult_cohort.eq(offspring.reproduction_round + 1).all()
        or not offspring.reproduction_round.isin(range(max(states))).all()
        or not offspring.kind.isin(KINDS).all()
    ):
        raise ValueError("invalid offspring time or kind")
    offspring = offspring.copy()
    offspring["recruited"] = offspring.survived_to_adult.map(_flag)
    if not {"generation", "insect_id", "host_id", "x"} <= set(events):
        raise ValueError("events require generation, insect_id, host_id and x")
    if events.duplicated(["generation", "insect_id"]).any():
        raise ValueError("duplicate larval event")
    if not events.generation.isin(range(max(states))).all():
        raise ValueError("event outside committed rounds")
    if not np.isfinite(events.x).all() or (events.x < 0).any():
        raise ValueError("invalid resource transfer")

    tables = {name: [] for name in TABLE_COLUMNS}
    mappings, groups, continuity = {}, {}, {}
    all_pools = {
        kind: sorted({int(p) for s in states.values() for p in s[kind + "_species"]})
        for kind in KINDS
    }
    for t in sorted(states):
        state = states[t]
        for kind in KINDS:
            mapping, current, component_sizes = cluster_genotypes(
                state[kind + "_ids"],
                state[kind + "_species"],
                state[kind + "_genotypes"],
                kind=kind,
                generation=t,
                design=design,
            )
            mappings[t, kind], groups[t, kind] = mapping, current
            incoming, outgoing, weights = defaultdict(set), defaultdict(set), Counter()
            if t:
                births = offspring[(offspring.kind == kind) & (offspring.adult_cohort == t)]
                recruits = births[births.recruited]
                if set(recruits.offspring_id) != set(mapping):
                    raise ValueError(f"{kind} recruits do not match snapshot {t}")
                previous = mappings[t - 1, kind]
                for child in recruits.itertuples(index=False):
                    target = mapping[int(child.offspring_id)]
                    for ident in (child.maternal_id, child.paternal_id):
                        source = _source(previous, ident)
                        weights[source, target] += 0.5
                        if source in groups[t - 1, kind] and target in current:
                            incoming[target].add(source)
                            outgoing[source].add(target)
                for (source, target), weight in sorted(weights.items()):
                    tables["transitions"].append(
                        dict(
                            kind=kind,
                            source_generation=t - 1,
                            target_generation=t,
                            source_cluster=source,
                            target_cluster=target,
                            recruited_genetic_contribution=weight,
                        )
                    )
                for old, details in groups[t - 1, kind].items():
                    if not outgoing[old]:
                        pool_alive = np.any(
                            state[kind + "_species"] == details["compatibility_pool"]
                        )
                        represented = any(source == old for source, _ in weights)
                        tables["lineage_events"].append(
                            dict(
                                generation=t,
                                kind=kind,
                                event="termination",
                                source_clusters=json.dumps([old]),
                                target_cluster=UNASSIGNED,
                                reason="pool_extinct"
                                if not pool_alive
                                else (
                                    "observation_not_detected"
                                    if represented
                                    else "no_recruited_descendants"
                                ),
                                event_scope="operational_component",
                            )
                        )
            for label, details in current.items():
                sources = sorted(incoming[label])
                split = any(len(outgoing[source]) > 1 for source in sources)
                merge = len(sources) > 1
                event = (
                    "split_merge"
                    if split and merge
                    else "merge"
                    if merge
                    else ("split" if split else "continue" if sources else "origin")
                )
                if event == "continue":
                    lineage, age = continuity[sources[0]]
                    age += 1
                else:
                    lineage, age = label, 1
                continuity[label] = lineage, age
                details.update(lineage_id=lineage, continuity_generations=age)
                tables["clusters"].append(
                    dict(
                        generation=t,
                        kind=kind,
                        cluster_id=label,
                        compatibility_pool=details["compatibility_pool"],
                        size=details["size"],
                        lineage_id=lineage,
                        continuity_generations=age,
                        persistent_candidate=age >= design["persistence_generations"],
                    )
                )
                tables["lineage_events"].append(
                    dict(
                        generation=t,
                        kind=kind,
                        event=event,
                        source_clusters=json.dumps(sources),
                        target_cluster=label,
                        reason="",
                        event_scope="operational_component",
                    )
                )
            for ident, pool in zip(state[kind + "_ids"], state[kind + "_species"]):
                label = mapping[int(ident)]
                lineage, age = continuity.get(label, (UNASSIGNED, 0))
                tables["memberships"].append(
                    dict(
                        generation=t,
                        kind=kind,
                        individual_id=int(ident),
                        compatibility_pool=int(pool),
                        cluster_id=label,
                        lineage_id=lineage,
                        component_size=component_sizes[int(ident)],
                        continuity_generations=age,
                        persistent_candidate=age >= design["persistence_generations"],
                    )
                )
            for pool in all_pools[kind]:
                n = int(np.count_nonzero(state[kind + "_species"] == pool))
                tables["populations"].append(
                    dict(
                        generation=t,
                        kind=kind,
                        compatibility_pool=pool,
                        individual_count=n,
                        detected_clusters=sum(
                            d["compatibility_pool"] == pool for d in current.values()
                        ),
                        pool_state="present" if n else "empty_at_snapshot",
                        extinction_observed=bool(
                            t and n == 0 and np.any(states[t - 1][kind + "_species"] == pool)
                        ),
                    )
                )

    for t in range(max(states)):
        for kind in KINDS:
            born = offspring[(offspring.kind == kind) & (offspring.reproduction_round == t)]
            mapping = mappings[t, kind]
            opportunity = dict(
                reproduction_round=t,
                kind=kind,
                eligible_mothers=0,
                classified_mothers=0,
                candidate_edges=0,
                cross_cluster_candidate_edges=0,
                classified_probability_mass=0.0,
                cross_probability_mass=0.0,
            )
            pool_fields = {
                "generation",
                "kind",
                "individual_id",
                "eligible_as_mother",
                "candidate_mate_ids",
                "candidate_mate_probabilities",
            }
            available = pool_fields <= set(parents)
            if available:
                ledger = parents[(parents.generation == t) & (parents.kind == kind)]
                for mother in ledger.to_dict("records"):
                    if not _flag(mother["eligible_as_mother"]):
                        continue
                    candidates = json.loads(mother["candidate_mate_ids"])
                    probabilities = np.asarray(json.loads(mother["candidate_mate_probabilities"]))
                    if (
                        len(candidates) != len(probabilities)
                        or len(set(candidates)) != len(candidates)
                        or not np.isfinite(probabilities).all()
                        or (probabilities < 0).any()
                        or (len(candidates) and not np.isclose(probabilities.sum(), 1))
                    ):
                        raise ValueError("invalid saved candidate mate distribution")
                    if not candidates:
                        continue
                    opportunity["eligible_mothers"] += 1
                    source = _source(mapping, mother["individual_id"])
                    if source in (UNASSIGNED, UNKNOWN):
                        continue
                    opportunity["classified_mothers"] += 1
                    for father, probability in zip(candidates, probabilities):
                        target = _source(mapping, father)
                        if target not in (UNASSIGNED, UNKNOWN):
                            opportunity["candidate_edges"] += 1
                            opportunity["classified_probability_mass"] += probability
                            if source != target:
                                opportunity["cross_cluster_candidate_edges"] += 1
                                opportunity["cross_probability_mass"] += probability
            mass = opportunity["classified_probability_mass"]
            opportunity.update(
                probability_weighted_cross_fraction=opportunity["cross_probability_mass"] / mass
                if mass
                else np.nan,
                opportunity_status="saved_maternal_pool" if available else "not_saved",
            )
            if not available:
                for key in set(opportunity) - {"kind", "reproduction_round", "opportunity_status"}:
                    opportunity[key] = np.nan
            tables["mating_opportunities"].append(opportunity)
            for stage in ("birth", "recruitment"):
                sample = born if stage == "birth" else born[born.recruited]
                classified, cross, unknown = 0, 0, 0
                for child in sample.itertuples(index=False):
                    m, f = (
                        _source(mapping, ident) for ident in (child.maternal_id, child.paternal_id)
                    )
                    unknown += UNKNOWN in (m, f)
                    if m not in (UNASSIGNED, UNKNOWN) and f not in (UNASSIGNED, UNKNOWN):
                        classified += 1
                        cross += m != f
                tables["mating"].append(
                    dict(
                        reproduction_round=t,
                        offspring_cohort=t + 1,
                        parent_generation=t,
                        kind=kind,
                        stage=stage,
                        offspring_weighted_pairs=len(sample),
                        classified_parent_pairs=classified,
                        cross_cluster_parent_pairs=cross,
                        unknown_parent_pairs=unknown,
                        unclassified_parent_pairs=len(sample) - classified,
                        cross_fraction_of_classified=cross / classified if classified else np.nan,
                        opportunity_status="separate_maternal_pool_table"
                        if available
                        else "not_saved",
                    )
                )
        born = offspring[(offspring.kind == "insect") & (offspring.reproduction_round == t)]
        ev = events[events.generation == t]
        if set(ev.insect_id) != set(born.offspring_id):
            raise ValueError("larval events do not match all births including dead larvae")
        lookup = born.set_index("offspring_id")
        edges = defaultdict(lambda: np.zeros(3))
        for event in ev.itertuples(index=False):
            child = lookup.loc[event.insect_id]
            if int(event.host_id) != int(child.birth_host_id):
                raise ValueError("larval event / birth host mismatch")
            plant = _source(mappings[t, "plant"], event.host_id)
            maternal = _source(mappings[t, "insect"], child.maternal_id)
            edges["maternal_source", t, maternal, plant] += (
                float(event.x),
                1,
                int(child.recruited),
            )
            if child.recruited:
                recruit = _source(mappings[t + 1, "insect"], event.insect_id)
                edges["recruit_cluster", t + 1, recruit, plant] += (0, 0, 1)
        for (semantics, source_t, insect, plant), weight in sorted(edges.items()):
            tables["networks"].append(
                dict(
                    reproduction_round=t,
                    insect_generation=source_t,
                    plant_generation=t,
                    edge_semantics=semantics,
                    insect_cluster=insect,
                    plant_cluster=plant,
                    actual_consumption=weight[0] if semantics == "maternal_source" else np.nan,
                    birth_count=int(weight[1]) if semantics == "maternal_source" else np.nan,
                    recruitment_count=int(weight[2]),
                )
            )
        maternal_weights = sum(
            (weight for key, weight in edges.items() if key[0] == "maternal_source"),
            start=np.zeros(3),
        )
        recruit_count = sum(
            weight[2] for key, weight in edges.items() if key[0] == "recruit_cluster"
        )
        expected = [float(ev.x.sum()), len(born), int(born.recruited.sum())]
        if not np.allclose(maternal_weights, expected) or recruit_count != expected[2]:
            raise ValueError("operational network ledger imbalance")
        tables["network_balances"].append(
            dict(
                reproduction_round=t,
                birth_count=expected[1],
                recruitment_count=expected[2],
                actual_consumption=expected[0],
                maternal_network_births=int(maternal_weights[1]),
                maternal_network_recruits=int(maternal_weights[2]),
                maternal_network_consumption=float(maternal_weights[0]),
                recruit_network_recruits=int(recruit_count),
            )
        )
    return tables, hashes


def summarize(batch, design):
    """Summarize an audited v0.1.7 batch; return output tables for further inspection."""
    batch = batch if isinstance(batch, Batch) else Batch.open(Path(batch), allow_lifecycle=True)
    if read_settings(batch.path).get("version") != "0.1.7":
        raise ValueError("lineages require a v0.1.7 batch")
    settings = read_design(design)
    collected, inputs, processing = defaultdict(list), {}, []
    for run in batch.status.itertuples(index=False):
        folder = run_path(batch.path, run.case, run.replicate)
        if run.status == "invalid_parameter":
            processing.append(
                dict(case=run.case, replicate=run.replicate, status="invalid_parameter")
            )
            continue
        tables, hashes = _summarize_run(folder, settings)
        for name, rows in tables.items():
            collected[name].extend(
                dict(case=run.case, replicate=run.replicate, **row) for row in rows
            )
        for name, digest in hashes.items():
            inputs[(folder / name).relative_to(batch.path).as_posix()] = digest
        inputs[(folder / "metadata.json").relative_to(batch.path).as_posix()] = sha256(
            folder / "metadata.json"
        )
        processing.append(dict(case=run.case, replicate=run.replicate, status="summarized"))
    output = batch.summary_dir / "lineages"
    output.mkdir(parents=True, exist_ok=True)
    frames = {
        name: pd.DataFrame(collected[name], columns=["case", "replicate", *columns.split()])
        for name, columns in TABLE_COLUMNS.items()
    }
    frames["run_processing"] = pd.DataFrame(processing)
    for name, frame in frames.items():
        atomic_csv(frame, output / f"{name}.csv")
    inputs["settings.json"] = sha256(batch.path / "settings.json")
    inputs["data/run_status.csv"] = sha256(batch.path / "data/run_status.csv")
    inputs[audit_file(batch.path).relative_to(batch.path).as_posix()] = sha256(
        audit_file(batch.path)
    )
    write_json(
        output / "metadata.json",
        dict(
            observation_schema="operational-qtl-lineages-1",
            analysis_source_sha256=sha256(Path(__file__)),
            design=settings,
            design_sha256=sha256(design) if not isinstance(design, dict) else None,
            input_sha256=inputs,
            biological_species_confirmed=False,
            distance="1 - locuswise maximum matching allele copies / (2 * number of loci)",
            opportunity_definition=(
                "eligible_mothers have a saved nonempty pool; candidate_edges count classified "
                "mother/father pairs; probability masses sum saved weights over those pairs"
            ),
            limitations=[
                "single-link components can merge through one bridging genotype",
                "QTL identity is a selected limited panel, not whole-genome reproductive isolation",
                "persistence is consecutive unbranched ancestry; split/merge reset the counter",
                "thresholds and persistence are observation choices, not biological species criteria",
                "actual mating weights offspring; saved opportunity pools weight eligible mothers equally",
                "maternal_source networks describe maternal ancestry, not larval cluster membership",
                "recruit_cluster networks link t+1 recruits to host plant clusters at t",
                "pool extinction is distinct from loss of a detected operational cluster",
            ],
        ),
    )
    return frames


def main():
    parser = batch_parser(__doc__)
    parser.add_argument("--design", type=Path, required=True)
    args = parser.parse_args()
    summarize(args.input, args.design)


if __name__ == "__main__":
    main()
