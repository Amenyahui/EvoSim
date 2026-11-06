"""Independent genetics checks; equations preserved from v0.1.11 audit."""

import hashlib
import json
from functools import partial

import numpy as np
import pandas as pd

from .events import (
    close,
    compare,
    flag,
    unique_keys,
)


def audit_phenotype(g, kind, clip=True, *, spec, loci, catalog, effects):
    raw = np.broadcast_to(
        np.asarray(spec["background"][kind], float), (len(g), len(spec["background"][kind]))
    ).copy()
    assert g.shape == (len(g), 2, len(loci[kind])) and g.dtype.kind in "iu"
    assert ((g >= 0) & (g < len(catalog))).all()
    for j, locus in enumerate(loci[kind]):
        for role in (0, 1):
            assert all(catalog.iloc[int(a)].locus_id == locus["locus_id"] for a in g[:, role, j])
            raw += (
                np.array([effects[int(a)] for a in g[:, role, j]]).reshape(len(g), raw.shape[1])
                * spec["phenotype_map"]["effect_scale"]
            )
    return np.clip(raw, *spec["phenotype_map"]["bounds"]) if clip else raw


def audit_children(
    kind,
    mothers,
    fathers,
    *,
    tables,
    t,
    next_id,
    loci,
    phenotype,
    catalog,
    calls,
    seen_calls,
    c,
    meta,
    add,
    tx_lookup,
    seen_tx,
    mutation_lookup,
    spec,
    seen_mutations,
    expected,
    birth_genotypes,
    child_rows,
    state,
):
    k = int(kind == "insect")
    ids, labels, g = (state[kind + x] for x in ("_ids", "_species", "_genotypes"))
    rows = tables["offspring"]
    rows = rows[(rows.generation == t + 1) & (rows.kind == kind)].sort_values("offspring_id")
    np.testing.assert_array_equal(
        rows.offspring_id, np.arange(next_id[k], next_id[k] + len(mothers))
    )
    next_id[k] += len(mothers)
    np.testing.assert_array_equal(rows.maternal_id, ids[mothers])
    np.testing.assert_array_equal(rows.paternal_id, ids[fathers])
    decoded = [json.loads(v) for v in rows.genotype]
    if decoded:
        assert np.asarray(decoded).dtype.kind in "iu", "noninteger newborn genotype"
    child_g = np.asarray(decoded, int).reshape(len(rows), 2, len(loci[kind]))
    phenotype(child_g, kind)
    assert (catalog.birth_generation.to_numpy()[child_g] <= t + 1).all()
    if len(rows):
        call = calls[t + 1, kind]
        seen_calls.add((t + 1, kind))
        seeds = np.random.SeedSequence([c.seed, 140, 4, k, t + 1]).generate_state(2)
        assert call["offspring_count"] == len(rows)
        assert (
            call["meiosis_seed"] == int(seeds[0]) + 1 and call["mutation_seed"] == int(seeds[1]) + 1
        )
        assert "SLiM version 5.2" in call["version"]
        assert (
            call["script_sha256"]
            == meta["source_hash"][
                "models/plant_insect/biology/genetics/slim_reproduction.slim"
                if "models/plant_insect/biology/genetics/slim_reproduction.slim"
                in meta["source_hash"]
                else "models/plant_insect/helpers/slim_reproduction.slim"
            ]
        )
    for i, row in enumerate(rows.itertuples()):
        values = dict(
            species=f"{'P' if k == 0 else 'H'}{labels[mothers[i]] + 1}",
            reproduction_round=t,
            adult_cohort=t + 1,
            maternal_id=int(ids[mothers[i]]),
            paternal_id=int(ids[fathers[i]]),
        )
        add("offspring", (t + 1, kind, row.offspring_id), **values)
        ancestry = {}
        mutation_counts = [0, 0]
        for role, parent in enumerate((mothers[i], fathers[i])):
            for l, locus in enumerate(loci[kind]):
                key = (
                    t + 1,
                    kind,
                    row.offspring_id,
                    "maternal" if role == 0 else "paternal",
                    role,
                    locus["locus_id"],
                )
                transmission = tx_lookup[key]
                seen_tx.add(key)
                source = int(transmission.transmitted_homolog_before_mutation)
                assert source in (0, 1) and transmission.parent_id == ids[parent]
                assert transmission.chromosome_id == locus["chromosome_id"]
                compare(transmission.position_cM, locus["position_cM"], key)
                assert pd.isna(transmission.crossover_count)
                old, new = int(g[parent, source, l]), int(child_g[i, role, l])
                assert (
                    transmission.inherited_allele_id == old and transmission.child_allele_id == new
                )
                ancestry[role, locus["locus_id"]] = (
                    source,
                    flag(transmission.switch_before_locus),
                )
                mutation_key = (t + 1, kind, row.offspring_id, locus["locus_id"], role)
                mutation = mutation_lookup.get(mutation_key)
                if old == new:
                    assert mutation is None and pd.isna(transmission.mutation_event_id)
                    assert (
                        not k or not spec["mutation_enabled"][kind] or locus["mutation_rate"] != 1
                    )
                else:
                    assert (
                        mutation is not None
                        and mutation.parent_allele_id == old
                        and mutation.allele_id == new
                    )
                    assert transmission.mutation_event_id == new
                    assert mutation.species == row.species
                    assert mutation.event_id == new and mutation.role == key[3]
                    assert mutation.chromosome_id == locus["chromosome_id"]
                    assert (
                        mutation.mutation_stream
                        == f"SLiM/seed{calls[t + 1, kind]['mutation_seed']}/{kind}/generation{t + 1}"
                    )
                    seen_mutations.add(mutation_key)
                    mutation_counts[role] += 1
        for role in (0, 1):
            for chromosome in {l["chromosome_id"] for l in loci[kind]}:
                ordered = sorted(
                    (l for l in loci[kind] if l["chromosome_id"] == chromosome),
                    key=lambda l: l["position_cM"],
                )
                previous = None
                for locus in ordered:
                    source, switched = ancestry[role, locus["locus_id"]]
                    assert switched == (previous is not None and previous != source)
                    previous = source
        expected["offspring"][t + 1, kind, row.offspring_id].update(
            mutations_maternal=mutation_counts[0], mutations_paternal=mutation_counts[1]
        )
    birth_genotypes[kind], child_rows[kind] = child_g, rows
    return child_g, rows


def read_genetics(path, meta, c, tables):
    """Validate the allele catalog and bind an independently reconstructed phenotype map."""
    spec = c.genetics
    architecture = json.dumps(spec, sort_keys=True, separators=(",", ":")).encode()
    assert hashlib.sha256(architecture).hexdigest() == meta["genetic_architecture_sha256"]
    loci = {
        kind: [l for l in spec["loci"] if l["organism"] == kind] for kind in ("plant", "insect")
    }
    locus_lookup = {l["locus_id"]: l for l in spec["loci"]}
    saved_loci = pd.read_csv(path / "genetic_loci.csv", float_precision="round_trip")
    assert len(saved_loci) == len(spec["loci"])
    for actual, wanted in zip(saved_loci.to_dict("records"), spec["loci"]):
        for key, value in wanted.items():
            assert (json.loads(actual[key]) if isinstance(value, list) else actual[key]) == value
    catalog = tables["allele_catalog"]
    assert catalog.allele_id.tolist() == list(range(len(catalog)))
    effects = [np.asarray(json.loads(v), float) for v in catalog.effect]
    assert all(np.isfinite(e).all() for e in effects)
    for i, row in enumerate(spec["alleles"]):
        assert catalog.iloc[i].locus_id == row["locus_id"]
        assert pd.isna(catalog.iloc[i].parent_allele_id)
        assert catalog.iloc[i].birth_generation == 0
        np.testing.assert_array_equal(effects[i], row["effect"])
    for i, row in enumerate(catalog.itertuples()):
        locus = locus_lookup[row.locus_id]
        assert len(effects[i]) == len(spec["background"][locus["organism"]])
        other = list(set(range(len(effects[i]))) - set(locus["targets"]))
        assert np.all(effects[i][other] == 0)
    mutations = tables["mutations"]
    unique_keys(mutations, ("generation", "kind", "offspring_id", "locus_id", "copy"), "mutations")
    assert sorted(mutations.allele_id.tolist()) == list(range(len(spec["alleles"]), len(catalog)))
    mutation_lookup = {}
    for row in mutations.itertuples():
        a, old = int(row.allele_id), int(row.parent_allele_id)
        locus = locus_lookup[row.locus_id]
        assert (
            0 <= old < a and catalog.iloc[a].locus_id == catalog.iloc[old].locus_id == row.locus_id
        )
        assert row.kind == locus["organism"] and spec["mutation_enabled"][row.kind]
        assert locus["mutation_rate"] > 0 and row.copy in (0, 1)
        assert catalog.iloc[a].parent_allele_id == old
        assert (
            row.generation == catalog.iloc[a].birth_generation > catalog.iloc[old].birth_generation
        )
        delta = np.asarray(json.loads(row.delta), float)
        close(effects[a], effects[old] + delta, ("mutation effect", a))
        fixed = np.asarray(locus["mutation_sd"]) == 0
        close(delta[fixed], np.asarray(locus["mutation_mean"])[fixed], ("fixed mutation", a))
        mutation_lookup[row.generation, row.kind, row.offspring_id, row.locus_id, row.copy] = row

    phenotype = partial(
        audit_phenotype,
        spec=spec,
        loci=loci,
        catalog=catalog,
        effects=effects,
    )
    return spec, loci, catalog, mutations, mutation_lookup, phenotype
