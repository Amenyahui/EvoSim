"""Independent reconstruction from saved genetic ledgers; no runtime G->P calls."""

import hashlib
import json
from itertools import combinations

import numpy as np
import pandas as pd


def audit_genetics(path, c, meta):
    spec = c.genetics
    mapping = spec["phenotype_map"]
    lower, upper = mapping["bounds"]
    assert meta["genetics_version"] == spec["schema_version"]
    slim = spec["schema_version"] == "slim-qtl-1"
    calls = {}
    if slim:
        for call in meta["slim_calls"]:
            key = (call["generation"], call["kind"])
            assert key not in calls
            seeds = np.random.SeedSequence(
                [c.seed, 140, 4, ("plant", "insect").index(call["kind"]), call["generation"]]
            ).generate_state(2)
            assert call["meiosis_seed"] == int(seeds[0]) + 1
            assert call["mutation_seed"] == int(seeds[1]) + 1
            assert "SLiM version 5.2" in call["version"]
            assert len(call["binary_sha256"]) == len(call["script_sha256"]) == 64
            assert (
                call["script_sha256"]
                == meta["source_hash"]["models/plant_insect/helpers/slim_reproduction.slim"]
            )
            calls[key] = call
    assert meta["ecology_rule_version"] == c.model_version
    encoded = json.dumps(spec, sort_keys=True, separators=(",", ":")).encode()
    assert hashlib.sha256(encoded).hexdigest() == meta["genetic_architecture_sha256"]
    loci = pd.read_csv(path / "genetic_loci.csv", float_precision="round_trip")
    assert len(loci) == len(spec["loci"])
    for actual, expected in zip(loci.to_dict("records"), spec["loci"]):
        for key, value in expected.items():
            assert (json.loads(actual[key]) if isinstance(value, list) else actual[key]) == value
    alleles = pd.read_csv(path / "allele_catalog.csv", float_precision="round_trip")
    assert alleles.allele_id.tolist() == list(range(len(alleles)))
    effects = np.array([json.loads(x) for x in alleles.effect])
    assert np.isfinite(effects).all()
    keys = alleles.locus_id.tolist()
    founder_count = len(spec["alleles"])
    for i, a in enumerate(spec["alleles"]):
        assert keys[i] == a["locus_id"] and pd.isna(alleles.parent_allele_id.iloc[i])
        assert alleles.birth_generation.iloc[i] == 0
        np.testing.assert_array_equal(effects[i], a["effect"])
    events = pd.read_csv(path / "mutations.csv", float_precision="round_trip")
    assert sorted(events.allele_id.tolist()) == list(range(founder_count, len(alleles)))
    assert not events.duplicated(["generation", "kind", "offspring_id", "locus_id", "copy"]).any()
    by_locus = {x["locus_id"]: x for x in spec["loci"]}
    for row in events.itertuples():
        i, parent = int(row.allele_id), int(row.parent_allele_id)
        locus = by_locus[row.locus_id]
        assert 0 <= parent < i and keys[i] == keys[parent] == row.locus_id
        assert row.kind == locus["organism"] and spec["mutation_enabled"][row.kind]
        assert locus["mutation_rate"] > 0 and row.copy in (0, 1)
        assert row.generation >= 1 and row.generation == alleles.birth_generation.iloc[i]
        assert row.generation > alleles.birth_generation.iloc[parent]
        assert row.parent_allele_id == alleles.parent_allele_id.iloc[i]
        delta = np.asarray(json.loads(row.delta))
        np.testing.assert_allclose(effects[i], effects[parent] + delta, rtol=1e-12, atol=1e-14)
        deterministic = np.asarray(locus["mutation_sd"]) == 0
        np.testing.assert_array_equal(
            delta[deterministic], np.asarray(locus["mutation_mean"])[deterministic]
        )
        if slim:
            seed = calls[row.generation, row.kind]["mutation_seed"]
            assert row.mutation_stream == f"SLiM/seed{seed}/{row.kind}/generation{row.generation}"
        else:
            assert (
                row.mutation_stream
                == f"seed140/phase1/{row.kind}/species{int(row.species[1:]) - 1}"
            )
    for i, key in enumerate(keys):
        assert key in by_locus
        other = list(set(range(effects.shape[1])) - set(by_locus[key]["targets"]))
        assert np.all(effects[i, other] == 0)
    transmissions = pd.read_csv(path / "gamete_transmissions.csv", float_precision="round_trip")
    transmission_groups = transmissions.groupby(["generation", "kind"], sort=False).indices
    assert sum(map(len, transmission_groups.values())) == len(transmissions)
    seen_transmissions = 0
    children = pd.read_csv(path / "offspring.csv", float_precision="round_trip", low_memory=False)
    summary = pd.read_csv(path / "genetic_summary.csv", float_precision="round_trip")
    summary_groups = summary.groupby(["generation", "species"], sort=False)
    haplotypes = pd.read_csv(path / "haplotype_summary.csv", float_precision="round_trip")
    expected_haplotypes = []
    expected_summary_rows = 0
    seen_events = 0
    previous = None
    with (
        np.load(path / "genotypes.npz") as exported_g,
        np.load(path / "phenotypes.npz") as exported_p,
    ):
        expected_g, expected_p = set(), set()
        for file in sorted((path / "snapshots").glob("*.npz")):
            with np.load(file) as loaded:
                state = {key: loaded[key] for key in loaded.files}
            t = int(state["generation"])
            for k, (kind, sizes, names) in enumerate(
                (
                    ("plant", c.plant_sizes, ("D", "plant_voc_independent_trait")),
                    ("insect", c.insect_sizes, ("I", "A")),
                )
            ):
                g = state[kind + "_genotypes"]
                ids, labels = state[kind + "_ids"], state[kind + "_species"]
                local = [x for x in spec["loci"] if x["organism"] == kind]
                assert g.shape == (len(ids), 2, len(local)) and g.dtype.kind in "iu"
                assert ((g >= 0) & (g < len(alleles))).all()
                raw = np.broadcast_to(
                    np.asarray(spec["background"][kind], float), (len(ids), effects.shape[1])
                ).copy()
                for l, locus in enumerate(local):
                    for h in range(2):
                        assert all(keys[int(a)] == locus["locus_id"] for a in g[:, h, l])
                        assert np.all(alleles.birth_generation.to_numpy()[g[:, h, l]] <= t)
                        raw += effects[g[:, h, l]] * mapping["effect_scale"]
                derived = np.clip(raw, lower, upper)
                for name, values in zip(names, (derived[:, : c.n], derived[:, c.n :])):
                    np.testing.assert_allclose(values, state[name], rtol=1e-12, atol=1e-14)
                    key = f"g{t:06d}_{name}"
                    expected_p.add(key)
                    np.testing.assert_array_equal(exported_p[key], state[name])
                for name in (kind + "_genotypes", kind + "_ids", kind + "_species"):
                    key = f"g{t:06d}_{name}"
                    expected_g.add(key)
                    np.testing.assert_array_equal(exported_g[key], state[name])
                prefix = "P" if kind == "plant" else "H"
                if t == 0:
                    for s, size in enumerate(sizes):
                        rng = np.random.default_rng(np.random.SeedSequence([c.seed, 140, 0, k, s]))
                        for l, locus in enumerate(local):
                            freq = spec["founder_frequencies"][f"{prefix}{s + 1}"][
                                locus["locus_id"]
                            ]
                            p = np.array([x["frequency"] for x in freq], float)
                            expected = rng.choice(
                                [x["allele_id"] for x in freq], size=(size, 2), p=p / p.sum()
                            )
                            np.testing.assert_array_equal(g[labels == s, :, l], expected)
                else:
                    rows = children[
                        (children.generation == t) & (children.kind == kind)
                    ].sort_values("offspring_id")
                    assert rows.offspring_id.tolist() == ids.tolist()
                    lookup = {int(x): i for i, x in enumerate(previous[kind + "_ids"])}
                    maternal = np.array([lookup[int(x)] for x in rows.maternal_id], dtype=int)
                    paternal = np.array([lookup[int(x)] for x in rows.paternal_id], dtype=int)
                    assert np.all(maternal != paternal)
                    for indices in (maternal, paternal):
                        np.testing.assert_array_equal(previous[kind + "_species"][indices], labels)
                    selected = events[(events.generation == t) & (events.kind == kind)]
                    seen_events += len(selected)
                    tx = transmissions.iloc[transmission_groups.get((t, kind), [])]
                    assert len(tx) == len(ids) * 2 * len(local), (
                        str(path),
                        t,
                        kind,
                        len(tx),
                        len(ids),
                        len(local),
                    )
                    if slim and len(ids):
                        assert calls[t, kind]["offspring_count"] == len(ids)
                    assert not tx.duplicated(["offspring_id", "role", "locus_id"]).any()
                    seen_transmissions += len(tx)
                    child_lookup = {int(x): i for i, x in enumerate(ids)}
                    locus_lookup = {x["locus_id"]: i for i, x in enumerate(local)}
                    event_lookup = {
                        (int(r.offspring_id), int(r.copy), r.locus_id): r
                        for r in selected.itertuples()
                    }
                    ancestry = {}
                    for row in tx.itertuples():
                        assert row.role in ("maternal", "paternal")
                        h = 0 if row.role == "maternal" else 1
                        i, l = child_lookup[row.offspring_id], locus_lookup[row.locus_id]
                        parent = (maternal if h == 0 else paternal)[i]
                        assert row.parent_id == previous[kind + "_ids"][parent]
                        source = int(row.transmitted_homolog_before_mutation)
                        assert source in (0, 1)
                        if slim:
                            assert pd.isna(row.crossover_count)
                        else:
                            assert row.crossover_count >= 0 and row.crossover_count == int(
                                row.crossover_count
                            )
                        assert (
                            row.chromosome_id == local[l]["chromosome_id"]
                            and row.position_cM == local[l]["position_cM"]
                        )
                        assert row.switch_before_locus in (True, False)
                        if not slim:
                            assert bool(row.switch_before_locus) == bool(row.crossover_count % 2)
                        old = previous[kind + "_genotypes"][parent, source, l]
                        assert old == row.inherited_allele_id and row.child_allele_id == g[i, h, l]
                        event = event_lookup.get((int(row.offspring_id), h, row.locus_id))
                        if event is None:
                            assert pd.isna(row.mutation_event_id) and old == row.child_allele_id
                        else:
                            assert (
                                row.mutation_event_id
                                == event.event_id
                                == event.allele_id
                                == row.child_allele_id
                            )
                            assert event.parent_allele_id == old and event.role == row.role
                        ancestry.setdefault((i, h, row.chromosome_id), []).append(row)
                    for group in ancestry.values():
                        group.sort(key=lambda r: r.position_cM)
                        assert not group[0].switch_before_locus
                        if not slim:
                            assert group[0].crossover_count == 0
                        for prev, curr in zip(group, group[1:]):
                            assert curr.transmitted_homolog_before_mutation == (
                                int(prev.transmitted_homolog_before_mutation)
                                ^ int(curr.switch_before_locus)
                            )
                    for i, row in enumerate(rows.itertuples()):
                        for h, role in enumerate(("maternal", "paternal")):
                            assert getattr(row, "mutations_" + role) == sum(
                                key[0] == row.offspring_id and key[1] == h for key in event_lookup
                            )
                    for locus in local:
                        if spec["mutation_enabled"][kind] and locus["mutation_rate"] == 1:
                            assert (selected.locus_id == locus["locus_id"]).sum() == 2 * len(ids)
                groups = {}
                for l, locus in enumerate(local):
                    groups.setdefault(locus["chromosome_id"], []).append(l)
                for indices in groups.values():
                    indices.sort(key=lambda l: local[l]["position_cM"])
                    for left, right in combinations(indices, 2):
                        for species in range(len(sizes)):
                            copies = g[labels == species][:, :, [left, right]].reshape(-1, 2)
                            pairs, counts = np.unique(copies, axis=0, return_counts=True)
                            for pair, count in zip(pairs, counts):
                                expected_haplotypes.append(
                                    dict(
                                        generation=t,
                                        kind=kind,
                                        species=f"{prefix}{species + 1}",
                                        locus_left=local[left]["locus_id"],
                                        locus_right=local[right]["locus_id"],
                                        allele_left=int(pair[0]),
                                        allele_right=int(pair[1]),
                                        frequency=float(count / len(copies)),
                                    )
                                )
                for s in range(len(sizes)):
                    mask = labels == s
                    n = int(mask.sum())
                    ss = summary_groups.get_group((t, f"{prefix}{s + 1}"))
                    assert ss.N.eq(n).all() and ss.kind.eq(kind).all()
                    expected_count = 0
                    offset = 0
                    for name in names:
                        for j in range(state[name].shape[1]):
                            rows = ss[(ss.record_type == "trait") & (ss.trait == f"{name}_{j + 1}")]
                            assert len(rows) == 1, (str(path), t, kind, s, name, j, len(rows))
                            expected_count += 1
                            np.testing.assert_allclose(
                                rows.raw_mean,
                                raw[mask, offset + j].mean() if n else np.nan,
                                equal_nan=True,
                            )
                            x = state[name][mask, j]
                            unbounded = raw[mask, offset + j]
                            expected = (
                                [
                                    x.mean(),
                                    x.var(),
                                    np.mean((unbounded < lower) | (unbounded > upper)),
                                ]
                                if n
                                else [np.nan] * 3
                            )
                            np.testing.assert_allclose(
                                rows[["mean", "variance", "clipped_fraction"]].iloc[0],
                                expected,
                                rtol=1e-12,
                                atol=1e-14,
                                equal_nan=True,
                            )
                        offset += state[name].shape[1]
                    for l, locus in enumerate(local):
                        rows = ss[(ss.record_type == "allele") & (ss.locus_id == locus["locus_id"])]
                        active = alleles[
                            (alleles.locus_id == locus["locus_id"])
                            & (alleles.birth_generation <= t)
                        ].allele_id.tolist()
                        assert rows.allele_id.tolist() == active
                        expected_count += len(active)
                        copies = g[mask, :, l]
                        f = (
                            np.array([np.count_nonzero(copies == a) / (2 * n) for a in active])
                            if n
                            else np.full(len(active), np.nan)
                        )
                        np.testing.assert_allclose(
                            rows.frequency, f, rtol=1e-12, atol=1e-14, equal_nan=True
                        )
                        he = 1 - f @ f if n else np.nan
                        ho = np.mean(copies[:, 0] != copies[:, 1]) if n else np.nan
                        np.testing.assert_allclose(
                            rows.gene_diversity, he, rtol=1e-12, atol=1e-14, equal_nan=True
                        )
                        np.testing.assert_allclose(
                            rows.observed_heterozygosity, ho, rtol=1e-12, atol=1e-14, equal_nan=True
                        )
                        m = events[
                            (events.generation == t)
                            & (events.species == f"{prefix}{s + 1}")
                            & (events.locus_id == locus["locus_id"])
                        ]
                        assert rows.mutation_count.eq(len(m)).all()
                    assert len(ss) == expected_count
                    expected_summary_rows += expected_count
            previous = state
        assert set(exported_g.files) == expected_g and set(exported_p.files) == expected_p
    assert seen_events == len(events) and expected_summary_rows == len(summary)
    assert seen_transmissions == len(transmissions)
    if slim:
        expected_calls = set(zip(children.generation, children.kind))
        assert set(calls) == expected_calls
    if expected_haplotypes:
        pd.testing.assert_frame_equal(
            haplotypes, pd.DataFrame(expected_haplotypes), check_exact=False, rtol=1e-12, atol=1e-14
        )
    else:
        assert haplotypes.empty
    return dict(genetic_audit=True, alleles_checked=len(alleles), mutations_checked=len(events))


def genetic_limit(c, states, populations):
    """Independently check a reported genetic allocation stop."""
    if c.genetics is None or c.genetics.get("mode") != "sexual_qtl":
        return False
    g = c.genetics
    n_loci = [sum(x["organism"] == kind for x in g["loci"]) for kind in ("plant", "insect")]
    if not states:
        return (
            sum(2 * l * sum(n) for l, n in zip(n_loci, (c.plant_sizes, c.insect_sizes)))
            > g["max_genotype_values_safety"]
        )
    stored = sum(s[k + "_genotypes"].size for s in states for k in ("plant", "insect"))
    counts = [
        int(populations.loc[populations.kind == k, "N_offspring"].sum())
        for k in ("plant", "insect")
    ]
    new_values = sum(2 * l * n for l, n in zip(n_loci, counts))
    # The runtime reserves the all-eligible-copies-mutate upper bound before allocating either taxon.
    # Include extinct alleles; the current genotypes alone do not give the catalog size.
    existing = populations.attrs["allele_count"]
    reserve = sum(
        2 * n * sum(x["organism"] == k and x["mutation_rate"] > 0 for x in g["loci"])
        for k, n in zip(("plant", "insect"), counts)
        if g["mutation_enabled"][k]
    )
    return (
        stored + new_values > g["max_genotype_values_safety"]
        or existing + reserve > g["max_alleles_safety"]
    )
