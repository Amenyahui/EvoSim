"""Independent reconstruction of saved adult/larval life-history ledgers.

Only configuration validation and generic assertion helpers are shared with the
simulator. Equations, event joins and random draws below use saved raw records.
"""

import hashlib
import json
from itertools import combinations

import numpy as np
import pandas as pd

from .audit_contributions import close, unique_keys
from .helpers.config import Config, config_payload


def flag(value):
    assert value in (True, False, "True", "False"), ("invalid boolean", value)
    return value in (True, "True")


def host_name(value):
    return "unknown" if value < 0 else f"P{value + 1}"


def distance(origin, targets, shape):
    delta = np.abs(np.asarray(targets) - origin)
    return np.minimum(delta, shape - delta).max(axis=-1)


def compare(actual, expected, context):
    if expected is None or (isinstance(expected, float) and np.isnan(expected)):
        assert pd.isna(actual), (context, actual, expected)
    elif isinstance(expected, (bool, np.bool_)):
        assert flag(actual) == expected, (context, actual, expected)
    elif isinstance(expected, str):
        assert actual == expected, (context, actual, expected)
    else:
        close(actual, expected, context)


def check_rows(frame, keys, expected, name):
    """Require exact key coverage; check every declared expected value."""
    unique_keys(frame, keys, name)
    indexed = {tuple(row[k] for k in keys): row for row in frame.to_dict("records")}
    assert set(indexed) == set(expected), (name, "missing/extra keys", set(expected) ^ set(indexed))
    for key, values in expected.items():
        for column, value in values.items():
            assert column in frame.columns, (name, "missing column", column)
            compare(indexed[key][column], value, (name, key, column))


def vector(value, labels, prefix):
    return (
        np.array([value[f"{prefix}{s + 1}"] for s in labels], float)
        if isinstance(value, dict)
        else np.full(len(labels), value, float)
    )


def audit_snapshot_summaries(path, c, states, catalog, phenotype):
    """Snapshot exports and reference assays have no larval-stage assumptions."""
    summary = pd.read_csv(path / "genetic_summary.csv", float_precision="round_trip")
    haplotypes = pd.read_csv(path / "haplotype_summary.csv", float_precision="round_trip")
    mutations = pd.read_csv(path / "mutations.csv", float_precision="round_trip")
    hap_expected, expected_count = {}, 0
    with np.load(path / "phenotypes.npz") as exported:
        export_keys = set()
        for state in states:
            t = int(state["generation"])
            for kind, sizes in (("plant", c.plant_sizes), ("insect", c.insect_sizes)):
                g, labels = state[kind + "_genotypes"], state[kind + "_species"]
                raw = phenotype(g, kind, clip=False)
                local = [l for l in c.genetics["loci"] if l["organism"] == kind]
                traits = ("plant_chemical_traits",) if kind == "plant" else ("I", "A")
                for name in ("D",) + traits if kind == "plant" else traits:
                    key = f"g{t:06d}_{name}"
                    export_keys.add(key)
                    np.testing.assert_array_equal(exported[key], state[name])
                for s in range(len(sizes)):
                    species = f"{'P' if kind == 'plant' else 'H'}{s + 1}"
                    selected = labels == s
                    n = int(selected.sum())
                    own = summary[
                        (summary.generation == t)
                        & (summary.kind == kind)
                        & (summary.species == species)
                    ]
                    assert own.N.eq(n).all()
                    offset = 0
                    for name in traits:
                        for j in range(state[name].shape[1]):
                            rows = own[
                                (own.record_type == "trait") & (own.trait == f"{name}_{j + 1}")
                            ]
                            assert len(rows) == 1
                            expected_count += 1
                            x, unbounded = state[name][selected, j], raw[selected, offset + j]
                            lower, upper = c.genetics["phenotype_map"]["bounds"]
                            for column, value in dict(
                                mean=x.mean() if n else np.nan,
                                variance=x.var() if n else np.nan,
                                raw_mean=unbounded.mean() if n else np.nan,
                                clipped_fraction=np.mean((unbounded < lower) | (unbounded > upper))
                                if n
                                else np.nan,
                            ).items():
                                compare(
                                    rows.iloc[0][column],
                                    value,
                                    ("summary", t, species, name, j, column),
                                )
                        offset += state[name].shape[1]
                    for l, locus in enumerate(local):
                        rows = own[
                            (own.record_type == "allele") & (own.locus_id == locus["locus_id"])
                        ]
                        active = catalog[
                            (catalog.locus_id == locus["locus_id"])
                            & (catalog.birth_generation <= t)
                        ].allele_id.tolist()
                        assert rows.allele_id.tolist() == active
                        expected_count += len(rows)
                        copies = g[selected, :, l]
                        frequencies = (
                            np.array([np.count_nonzero(copies == a) / (2 * n) for a in active])
                            if n
                            else np.full(len(active), np.nan)
                        )
                        close(rows.frequency, frequencies, ("summary frequencies", t, species, l))
                        close(
                            rows.observed_heterozygosity,
                            np.mean(copies[:, 0] != copies[:, 1]) if n else np.nan,
                            "heterozygosity",
                        )
                        close(
                            rows.gene_diversity,
                            1 - frequencies @ frequencies if n else np.nan,
                            "gene diversity",
                        )
                        m = mutations[
                            (mutations.generation == t)
                            & (mutations.species == species)
                            & (mutations.locus_id == locus["locus_id"])
                        ]
                        assert rows.mutation_count.eq(len(m)).all()
                    for chromosome in {l["chromosome_id"] for l in local}:
                        indices = sorted(
                            (i for i, l in enumerate(local) if l["chromosome_id"] == chromosome),
                            key=lambda i: local[i]["position_cM"],
                        )
                        for left, right in combinations(indices, 2):
                            copies = g[selected][:, :, [left, right]].reshape(-1, 2)
                            pairs, counts = np.unique(copies, axis=0, return_counts=True)
                            for pair, count in zip(pairs, counts):
                                hap_expected[
                                    t,
                                    kind,
                                    species,
                                    local[left]["locus_id"],
                                    local[right]["locus_id"],
                                    int(pair[0]),
                                    int(pair[1]),
                                ] = dict(frequency=float(count / len(copies)))
        assert set(exported.files) == export_keys
    assert expected_count == len(summary)
    check_rows(
        haplotypes,
        (
            "generation",
            "kind",
            "species",
            "locus_left",
            "locus_right",
            "allele_left",
            "allele_right",
        ),
        hap_expected,
        "haplotypes",
    )
    if not c.record_observations:
        return
    diagnostics = pd.read_csv(path / "diagnostics.csv", float_precision="round_trip")
    expected = {}
    ref = np.asarray(c.observation_reference_traits["plant_chemical_traits"], float)
    D = ref[:, :2]
    V = np.column_stack((c.eta_shared * ref[:, 1], ref[:, 2]))
    weights = np.asarray(c.weights) / sum(c.weights)
    for state in states:
        t = int(state["generation"])
        A, I = state["A"], state["I"]
        d = ((A[:, None, :] - V[None, :, :]) ** 2) @ weights
        if c.olfaction_mode == "neutral" or c.tau == 0:
            probability = np.full(d.shape, 1 / len(D))
        else:
            logits = -((d - d.min(axis=1, keepdims=True)) / c.sigma) / c.sigma / 2
            mass = np.exp(logits * c.tau)
            probability = mass / mass.sum(axis=1, keepdims=True)
        ratios = np.divide(
            I[:, None, :],
            D[None, :, :],
            out=np.full((len(I), *D.shape), np.inf),
            where=D[None, :, :] > 0,
        )
        q = np.minimum(1, ratios.min(axis=2))
        close(state["reference_p"], probability, ("reference choice", t))
        close(state["reference_q"], q, ("reference performance", t))
        assert state["reference_plant_species"].tolist() == [f"P{p + 1}" for p in range(len(D))]
        for s in range(len(c.insect_sizes)):
            selected = state["insect_species"] == s
            n = int(selected.sum())
            for p in range(len(D)):
                a, b = probability[selected, p], q[selected, p]
                joint = a * b
                expected[t, f"H{s + 1}", f"P{p + 1}"] = dict(
                    insect_count=n,
                    panel="equal_ancestral_species",
                    preference=a.mean() if n else np.nan,
                    q_unlimited=b.mean() if n else np.nan,
                    joint_pq=joint.mean() if n else np.nan,
                    cov_pq=joint.mean() - a.mean() * b.mean() if n else np.nan,
                    variance_p=a.var() if n else np.nan,
                    variance_q=b.var() if n else np.nan,
                    variance_pq=joint.var() if n else np.nan,
                )
    check_rows(
        diagnostics,
        ("generation", "insect_species", "reference_plant_species"),
        expected,
        "diagnostics",
    )


def audit_run_v016(path, meta):
    """Audit all births, including dead larvae, separately from adult recruits."""
    assert __debug__, "scientific audits require assertions"
    c = Config(**meta["config"])
    c.validate()
    assert c.model_version == meta["version"] == "0.1.6"
    assert meta["genetics_version"] == "slim-qtl-2"
    assert meta["ecology_rule_version"] == "0.1.6"
    payload = json.dumps(config_payload(c), sort_keys=True, separators=(",", ":")).encode()
    assert hashlib.sha256(payload).hexdigest() == meta["config_sha256"]
    assert json.loads((path / "status.json").read_text(encoding="utf-8")) == meta["status"]
    for name, digest in meta["output_sha256"].items():
        assert hashlib.sha256((path / name).read_bytes()).hexdigest() == digest, name
    names = (
        "events",
        "plants",
        "parents",
        "offspring",
        "populations",
        "networks",
        "oviposition",
        "larvae",
        "host_flows",
        "host_groups",
        "lagged_networks",
        "allele_contributions",
        "gamete_transmissions",
        "mutations",
        "allele_catalog",
    )
    tables = {
        name: pd.read_csv(path / f"{name}.csv", float_precision="round_trip", low_memory=False)
        for name in names
    }
    assert {f"{name}.csv" for name in names} <= set(meta["output_sha256"])
    assert {
        "genetic_loci.csv",
        "genetic_summary.csv",
        "haplotype_summary.csv",
        "genotypes.npz",
        "phenotypes.npz",
    } <= set(meta["output_sha256"])
    if c.record_observations:
        assert "diagnostics.csv" in meta["output_sha256"]
    files = sorted((path / "snapshots").glob("*.npz"))
    assert {file.relative_to(path).as_posix() for file in files} == {
        name for name in meta["output_sha256"] if name.startswith("snapshots/")
    }
    states = []
    for file in files:
        with np.load(file) as loaded:
            states.append(dict(loaded))
    assert [int(s["generation"]) for s in states] == list(range(len(states)))
    status = meta["status"]["status"]
    if not states:
        assert status == "computational_limit"
        assert all(tables[name].empty for name in names if name != "allele_catalog")
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

    def phenotype(g, kind, clip=True):
        raw = np.broadcast_to(
            np.asarray(spec["background"][kind], float), (len(g), len(spec["background"][kind]))
        ).copy()
        assert g.shape == (len(g), 2, len(loci[kind])) and g.dtype.kind in "iu"
        assert ((g >= 0) & (g < len(catalog))).all()
        for j, locus in enumerate(loci[kind]):
            for role in (0, 1):
                assert all(
                    catalog.iloc[int(a)].locus_id == locus["locus_id"] for a in g[:, role, j]
                )
                raw += (
                    np.array([effects[int(a)] for a in g[:, role, j]]).reshape(len(g), raw.shape[1])
                    * spec["phenotype_map"]["effect_scale"]
                )
        return np.clip(raw, *spec["phenotype_map"]["bounds"]) if clip else raw

    rng = {
        (phase, k, s): np.random.default_rng(np.random.SeedSequence([c.seed, 130, phase, k, s]))
        for phase in range(8)
        for k, sizes in enumerate((c.plant_sizes, c.insect_sizes))
        for s in range(len(sizes))
    }
    spatial_rng = [
        np.random.default_rng(np.random.SeedSequence([c.seed, 150, 1, k])) for k in range(2)
    ]
    shape = np.asarray(c.spatial["grid_shape"])
    first = states[0]
    for k, kind in enumerate(("plant", "insect")):
        sizes = c.plant_sizes if k == 0 else c.insect_sizes
        ids, labels, g = (first[kind + x] for x in ("_ids", "_species", "_genotypes"))
        np.testing.assert_array_equal(ids, np.arange(sum(sizes)))
        np.testing.assert_array_equal(labels, np.repeat(np.arange(len(sizes)), sizes))
        given = c.spatial[f"initial_{kind}_cells"]
        expected_cells = (
            np.random.default_rng(
                np.random.SeedSequence([c.seed, 150, 0 if k == 0 else 2])
            ).integers(0, shape, size=(len(ids), 2))
            if (given is None)
            else np.asarray(given).reshape(-1, 2)
        )
        np.testing.assert_array_equal(first[kind + "_cells"], expected_cells)
        for s, size in enumerate(sizes):
            founder_rng = np.random.default_rng(np.random.SeedSequence([c.seed, 140, 0, k, s]))
            for l, locus in enumerate(loci[kind]):
                freq = spec["founder_frequencies"][f"{'P' if k == 0 else 'H'}{s + 1}"][
                    locus["locus_id"]
                ]
                p = np.array([v["frequency"] for v in freq])
                expected = founder_rng.choice(
                    [v["allele_id"] for v in freq], (size, 2), p=p / p.sum()
                )
                np.testing.assert_array_equal(g[labels == s, :, l], expected)
    for key, value in (
        ("insect_birth_host_id", -1),
        ("insect_birth_host_species", -1),
        ("insect_adult_cohort", 0),
    ):
        np.testing.assert_array_equal(first[key], np.full(len(first["insect_ids"]), value))
    for s in range(len(c.insect_sizes)):
        mask = first["insect_species"] == s
        expected = (
            rng[6, 1, s].random(mask.sum())
            >= c.sexual_reproduction["insect_founder_female_probability"]
        )
        np.testing.assert_array_equal(first["insect_sex"][mask], expected.astype(np.int8))
    for state in states:
        for kind in ("plant", "insect"):
            ids, labels, cells = (state[kind + x] for x in ("_ids", "_species", "_cells"))
            assert len(ids) == len(np.unique(ids)) == len(labels) == len(cells)
            assert ((cells >= 0) & (cells < shape)).all()
            values = phenotype(state[kind + "_genotypes"], kind)
            if kind == "plant":
                close(
                    state["plant_chemical_traits"], values, (state["generation"], kind, "phenotype")
                )
                close(state["D"], values[:, :2], "plant defense")
            else:
                close(state["I"], values[:, : c.n], "insect resistance")
                close(state["A"], values[:, c.n :], "insect preference")
            assert (
                catalog.birth_generation.to_numpy()[state[kind + "_genotypes"]]
                <= int(state["generation"])
            ).all()
    with np.load(path / "genotypes.npz") as exported:
        expected_keys = {
            f"g{t:06d}_{kind}{suffix}"
            for t in range(len(states))
            for kind in ("plant", "insect")
            for suffix in ("_genotypes", "_ids", "_species")
        }
        assert set(exported.files) == expected_keys
        for t, state in enumerate(states):
            for kind in ("plant", "insect"):
                for suffix in ("_genotypes", "_ids", "_species"):
                    np.testing.assert_array_equal(
                        exported[f"g{t:06d}_{kind}{suffix}"], state[kind + suffix]
                    )
    tx = tables["gamete_transmissions"]
    unique_keys(tx, ("generation", "kind", "offspring_id", "role", "locus_id"), "transmissions")
    tx_lookup = {
        (r.generation, r.kind, r.offspring_id, r.role, r.locus_id): r for r in tx.itertuples()
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
        V = np.column_stack((c.eta_shared * chemical[:, 1], chemical[:, 2]))
        cue = V[rng[1, 0, 0].permutation(len(V))] if c.shuffle_voc else V
        fi = np.maximum(0, 1 - c.insect_cost * state["I"].sum(axis=1))
        current_rows, birth_genotypes, child_rows, live_masks = {}, {}, {}, {}

        def parents(kind, S, F, phi):
            k = int(kind == "insect")
            labels, ids, cells = (state[kind + x] for x in ("_species", "_ids", "_cells"))
            maternal = np.ones(len(ids), bool) if k == 0 else state["insect_sex"] == 0
            donor = np.ones(len(ids), bool) if k == 0 else state["insect_sex"] == 1
            radius = c.sexual_reproduction["pollen_radius" if k == 0 else "insect_mate_radius"]
            potential = []
            for i in range(len(ids)):
                mask = donor & (labels == labels[i]) & (distance(cells[i], cells, shape) <= radius)
                mask[i] = False
                potential.append(np.flatnonzero(mask) if maternal[i] else np.empty(0, int))
            survived = np.ones(len(ids), bool)
            if k == 0:
                for s in range(len(c.plant_sizes)):
                    selected = labels == s
                    survived[selected] = rng[2, 0, s].random(selected.sum()) < S[selected]
            pools = [pool[survived[pool]] for pool in potential]
            lam = vector(
                c.plant_offspring_lambda if k == 0 else c.insect_offspring_lambda,
                labels,
                "P" if k == 0 else "H",
            )
            mu = lam * F * phi
            counts = np.zeros(len(ids), int)
            eligible = maternal & survived & np.array([len(pool) > 0 for pool in pools], bool)
            for s in range(len(c.plant_sizes if k == 0 else c.insect_sizes)):
                selected = eligible & (labels == s)
                counts[selected] = rng[3, k, s].poisson(mu[selected])
            for i, ident in enumerate(ids):
                q = 1 - np.prod(1 - S[potential[i]])
                prob = S[i] * q
                add(
                    "parents",
                    (t, kind, int(ident)),
                    species=f"{'P' if k == 0 else 'H'}{labels[i] + 1}",
                    species_id=int(labels[i]),
                    is_maternal_role=bool(maternal[i]),
                    S_i=S[i],
                    F_i=F[i],
                    W_i=S[i] * F[i],
                    phi_i=phi[i],
                    lambda_s=lam[i],
                    conditional_mu=mu[i],
                    mate_probability=q,
                    expected_maternal_offspring=prob * mu[i] if maternal[i] else None,
                    maternal_offspring_variance=prob * mu[i] + prob * (1 - prob) * mu[i] ** 2
                    if maternal[i]
                    else None,
                    survived_to_mating_stage=bool(survived[i]),
                    eligible_as_mother=bool(maternal[i] and survived[i]),
                    eligible_as_donor=bool(donor[i] and survived[i]),
                    candidate_mate_count=len(pools[i]),
                    has_mate=bool(len(pools[i])),
                    sampling_applied=True,
                    contribution_committed=True,
                )
                if k:
                    expected["parents"][t, kind, int(ident)]["potential_egg_count"] = counts[i]
            return counts, pools, survived

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
            if len(candidates):
                weights = np.asarray(c.weights) / np.sum(c.weights)
                distances = (cue[candidates] - state["A"][m]) ** 2 @ weights
                if c.choice_mode == "forced":
                    relative = (
                        np.zeros(len(candidates))
                        if c.olfaction_mode == "neutral" or c.tau == 0
                        else -((distances - distances.min()) / c.sigma) / c.sigma / 2
                    )
                    logits = (relative - relative.max()) * c.tau
                    values = np.exp(logits)
                    probabilities = np.append(values / values.sum(), 0.0)
                else:
                    scores = (
                        np.full(len(candidates), c.b)
                        if c.olfaction_mode == "neutral"
                        else c.b - (distances / c.sigma) / c.sigma / 2
                    )
                    logits = (
                        np.append(scores, c.outside_utility / c.tau)
                        if c.tau >= 1
                        else np.append(
                            scores * c.tau if c.tau else np.zeros(len(scores)), c.outside_utility
                        )
                    )
                    logits = (logits - logits.max()) * (c.tau if c.tau >= 1 else 1)
                    values = np.exp(logits)
                    probabilities = values / values.sum()
                    probabilities /= probabilities.sum()
                selected = stream.choice(len(probabilities), p=probabilities)
                if selected < len(candidates):
                    hosts[j] = candidates[selected]
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
            )
        accepted = hosts >= 0
        for i, j in enumerate(np.flatnonzero(accepted)):
            key = (t, int(hids[egg_mothers[j]]), int(egg_indices[j]))
            expected["oviposition"][key]["offspring_id"] = next_id[1] + i
        mothers, hosts = egg_mothers[accepted], hosts[accepted]
        fathers = np.array([rng[5, 1, int(hl[m])].choice(hpools[m]) for m in mothers], int)

        def children(kind, mothers, fathers):
            k = int(kind == "insect")
            ids, labels, g = (state[kind + x] for x in ("_ids", "_species", "_genotypes"))
            rows = tables["offspring"]
            rows = rows[(rows.generation == t + 1) & (rows.kind == kind)].sort_values(
                "offspring_id"
            )
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
                    call["meiosis_seed"] == int(seeds[0]) + 1
                    and call["mutation_seed"] == int(seeds[1]) + 1
                )
                assert "SLiM version 5.2" in call["version"]
                assert (
                    call["script_sha256"]
                    == meta["source_hash"]["models/plant_insect/helpers/slim_reproduction.slim"]
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
                            transmission.inherited_allele_id == old
                            and transmission.child_allele_id == new
                        )
                        ancestry[role, locus["locus_id"]] = (
                            source,
                            flag(transmission.switch_before_locus),
                        )
                        mutation_key = (t + 1, kind, row.offspring_id, locus["locus_id"], role)
                        mutation = mutation_lookup.get(mutation_key)
                        if old == new:
                            assert mutation is None and pd.isna(transmission.mutation_event_id)
                            assert not spec["mutation_enabled"][kind] or locus["mutation_rate"] != 1
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

        hg, hchildren = children("insect", mothers, fathers)
        child_hl = hl[mothers]
        child_I = phenotype(hg, "insect")[:, : c.n]
        R = vector(c.resource_capacity, pl, "P")
        dmax = vector(c.max_feeding_demand, child_hl, "H")
        ratios = np.divide(child_I, D[hosts], out=np.full_like(child_I, np.inf), where=D[hosts] > 0)
        q = np.minimum(1, ratios.min(axis=1))
        demand = q * dmax
        total = np.bincount(hosts, weights=demand, minlength=len(pids))
        gamma = np.ones(len(pids))
        limited = total > R
        gamma[limited] = R[limited] / total[limited]
        x = demand * gamma[hosts]
        X = np.bincount(hosts, weights=x, minlength=len(pids))
        h = np.divide(X, R, out=np.zeros_like(X, dtype=float), where=R > 0)
        g = x / dmax
        live = np.zeros(len(hchildren), bool)
        sex = np.empty(len(hchildren), np.int8)
        for s in range(len(c.insect_sizes)):
            selected = child_hl == s
            live[selected] = rng[2, 1, s].random(selected.sum()) < g[selected]
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
                W_insect=g[i],
                survived_to_adult=bool(live[i]),
                insect_life_stage="larva",
                same_species_others=same,
                other_species_others=other,
                oviposition_event_id=row.oviposition_event_id,
                egg_index=int(egg_indices[egg]),
            )
            add("events", (t, row.offspring_id), **event)
            add("larvae", (t, row.offspring_id), **event)
        live_masks["insect"] = live
        fp = np.maximum(0, 1 - c.plant_cost * D.sum(axis=1) - c.voc_cost * chemical[:, 2])
        Sp = c.sP0 + (1 - c.sP0) * (1 - h)
        alpha = np.array(
            [
                [
                    c.plant_competition_alpha[f"P{s + 1}"][f"P{r + 1}"]
                    for r in range(len(c.plant_sizes))
                ]
                for s in range(len(c.plant_sizes))
            ]
        )
        phi = np.array(
            [
                1
                / (
                    1
                    + alpha[pl[i]]
                    @ np.bincount(
                        pl[distance(cell, pcells, shape) <= c.spatial["competition_radius"]],
                        minlength=len(c.plant_sizes),
                    )
                )
                for i, cell in enumerate(pcells)
            ]
        )
        pcounts, ppools, psurvived = parents("plant", Sp, fp, phi)
        pm = np.repeat(np.arange(len(pids)), pcounts)
        pf = np.array([rng[5, 0, int(pl[m])].choice(ppools[m]) for m in pm], int)
        _, pchildren = children("plant", pm, pf)
        child_pcells = (
            pcells[pm]
            + spatial_rng[0].integers(
                -c.spatial["plant_dispersal_radius"],
                c.spatial["plant_dispersal_radius"] + 1,
                size=(len(pm), 2),
            )
        ) % shape
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
                voc_cost=c.voc_cost * chemical[i, 2],
                W_plant=Sp[i] * fp[i],
                density_phi=phi[i],
                cell_row=int(pcells[i, 0]),
                cell_col=int(pcells[i, 1]),
            )
        for kind, mm, ff, survivor, live_parents in (
            ("insect", mothers, fathers, live, np.ones(len(hids), bool)),
            ("plant", pm, pf, np.ones(len(pm), bool), psurvived),
        ):
            k = int(kind == "insect")
            ids, labels = state[kind + "_ids"], state[kind + "_species"]
            bm, bf = np.bincount(mm, minlength=len(ids)), np.bincount(ff, minlength=len(ids))
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
            np.testing.assert_array_equal(
                nxt[kind + "_ids"], rows.offspring_id.to_numpy()[survivor]
            )
            np.testing.assert_array_equal(nxt[kind + "_species"], labels[mm][survivor])
            np.testing.assert_array_equal(nxt[kind + "_genotypes"], birth_genotypes[kind][survivor])
            np.testing.assert_array_equal(nxt[kind + "_cells"], adult_cells if k else child_pcells)
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
                n, born, n_next = (
                    int(selected.sum()),
                    int(child_selected.sum()),
                    int(np.count_nonzero(child_selected & survivor)),
                )
                assert np.isclose(weights[selected].sum(), born) and np.isclose(
                    recruited[selected].sum(), n_next
                )
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
                    density_phi=1.0 if k else float(phi[selected].mean()) if n else np.nan,
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
                    for l, locus in enumerate(loci[kind]):
                        for allele in np.union1d(before[:, :, l], after[:, :, l]):
                            z = (before[:, :, l] == allele).sum(axis=1) / 2
                            p0 = float(z.mean()) if n else np.nan
                            pw = float(ww[selected] @ z / count) if count else np.nan
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
    )
