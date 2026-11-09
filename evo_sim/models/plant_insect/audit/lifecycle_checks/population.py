"""Independent population checks; equations preserved from v0.1.12 audit."""

import json
from itertools import combinations

import numpy as np
import pandas as pd

from .events import (
    check_rows,
    close,
    compare,
    distance,
)
from .feeding import chemical_outcomes, contact_probability, module_choice, module_profiles
from .resources import vector


def audit_snapshot_summaries(path, c, states, catalog, phenotype):
    """Snapshot exports and reference assays have no larval-stage assumptions."""
    summary = pd.read_csv(path / "genetic_summary.csv", float_precision="round_trip")
    summary_groups = summary.groupby(["generation", "kind", "species"], sort=False)
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
                traits = (
                    ("plant_chemical_traits",) if kind == "plant" else ("insect_module_traits",)
                )
                export_traits = (("D",) if kind == "plant" else ("I", "A")) + traits
                for name in export_traits:
                    key = f"g{t:06d}_{name}"
                    export_keys.add(key)
                    np.testing.assert_array_equal(exported[key], state[name])
                for s in range(len(sizes)):
                    species = f"{'P' if kind == 'plant' else 'H'}{s + 1}"
                    selected = labels == s
                    n = int(selected.sum())
                    own = summary_groups.get_group((t, kind, species))
                    assert own.N.eq(n).all()
                    offset = 0
                    for name in traits:
                        trait_values = (
                            state[name][:, None] if state[name].ndim == 1 else state[name]
                        )
                        for j in range(trait_values.shape[1]):
                            rows = own[
                                (own.record_type == "trait") & (own.trait == f"{name}_{j + 1}")
                            ]
                            assert len(rows) == 1, (
                                "trait summary key",
                                t,
                                kind,
                                species,
                                name,
                                j,
                                len(rows),
                            )
                            expected_count += 1
                            x, unbounded = trait_values[selected, j], raw[selected, offset + j]
                            lower, upper = c.genetics["phenotype_map"]["bounds"]
                            for column, value in dict(
                                mean=x.mean() if n else np.nan,
                                variance=x.var() if n else np.nan,
                                raw_mean=unbounded.mean() if n else np.nan,
                                clipped_fraction=(
                                    np.mean((unbounded < lower) | (unbounded > upper))
                                    if n
                                    else np.nan
                                ),
                            ).items():
                                compare(
                                    rows.iloc[0][column],
                                    value,
                                    ("summary", t, species, name, j, column),
                                )
                        offset += trait_values.shape[1]
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
    V = ref[:, 3:5]
    reference_ids = c.observation_reference_traits["reference_ids"]
    for state in states:
        t = int(state["generation"])
        I = state["I"]
        probability = np.array(
            [
                module_choice(V, x, c, forced=True)[:-1] * contact_probability(D, x, c)
                for x in state["insect_module_traits"]
            ]
        ).reshape(len(I), len(D))
        maximum = np.broadcast_to(
            vector(c.max_feeding_demand, state["insect_species"], "H")[:, None], (len(I), len(D))
        )
        q = chemical_outcomes(D[None, :, :], I[:, None, :], maximum, maximum, c)[
            "survival_probability"
        ]
        close(state["reference_p"], probability, ("reference choice", t))
        close(state["reference_q"], q, ("reference performance", t))
        assert state["reference_plant_species"].tolist() == reference_ids
        for s in range(len(c.insect_sizes)):
            selected = state["insect_species"] == s
            n = int(selected.sum())
            for p in range(len(D)):
                a, b = probability[selected, p], q[selected, p]
                joint = a * b
                expected[t, f"H{s + 1}", reference_ids[p]] = dict(
                    insect_count=n,
                    panel="fixed_chemical_reference",
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


def audit_parents(
    kind,
    S,
    F,
    phi,
    *,
    hcells,
    state,
    c,
    shape,
    rng,
    add,
    t,
    expected,
    processing_cost,
    sensory_cost,
):
    k = int(kind == "insect")
    labels, ids, cells = (state[kind + x] for x in ("_species", "_ids", "_cells"))
    if k:
        cells = hcells
    maternal = np.ones(len(ids), bool) if k == 0 else state["insect_sex"] == 0
    donor = np.ones(len(ids), bool) if k == 0 else state["insect_sex"] == 1
    radius = c.sexual_reproduction["insect_mate_radius"] if k else 0
    active = True
    potential = []
    for i in range(len(ids)):
        mask = (
            donor & (labels == labels[i]) & (distance(cells[i], cells, shape) <= radius)
            if k
            else donor.copy()
        )
        mask[i] = False
        potential.append(np.flatnonzero(mask) if active and maternal[i] else np.empty(0, int))
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
    mu = lam * F * int(active)
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
            lambda_s=lam[i],
            conditional_mu=mu[i],
            mate_probability=q,
            expected_maternal_offspring=prob * mu[i] if k and maternal[i] else None,
            maternal_offspring_variance=prob * mu[i] + prob * (1 - prob) * mu[i] ** 2
            if k and maternal[i]
            else None,
            survived_to_mating_stage=bool(survived[i]),
            eligible_as_mother=bool(active and maternal[i] and survived[i]),
            eligible_as_donor=bool(active and donor[i] and survived[i]),
            candidate_mate_count=len(pools[i]),
            has_mate=bool(len(pools[i])),
            reproduction_active=active,
            propagule_count=int(counts[i]) if k == 0 else None,
            expected_propagules=prob * mu[i] if k == 0 else None,
            propagule_variance=prob * mu[i] + prob * (1 - prob) * mu[i] ** 2 if k == 0 else None,
            sampling_applied=True,
            contribution_committed=True,
        )
        if k:
            expected["parents"][t, kind, int(ident)]["potential_egg_count"] = counts[i]
    for i, ident in enumerate(ids):
        pool = pools[i]
        probabilities = np.full(len(pool), 1 / len(pool)) if len(pool) else np.empty(0)
        expected["parents"][t, kind, int(ident)].update(
            mating_row=int(cells[i, 0]),
            mating_col=int(cells[i, 1]),
            candidate_mate_ids=json.dumps(ids[pool].tolist()),
            candidate_mate_probabilities=json.dumps(probabilities.tolist()),
        )
        if k:
            expected["parents"][t, kind, int(ident)].update(
                processing_expression_cost=processing_cost[i],
                sensory_expression_cost=sensory_cost[i],
            )
    return counts, pools, survived


def check_population(path, c, states, loci, catalog, phenotype, spec):
    """Reconstruct founders and validate every saved population/genotype snapshot."""
    rng = {
        (phase, k, s): np.random.default_rng(np.random.SeedSequence([c.seed, 130, phase, k, s]))
        for phase in range(10)
        for k, sizes in enumerate((c.plant_sizes, c.insect_sizes))
        for s in range(len(sizes))
    }
    spatial_rng = [
        np.random.default_rng(np.random.SeedSequence([c.seed, 150, 1, k])) for k in range(2)
    ]
    shape = np.asarray(c.spatial["grid_shape"])
    first = states[0]
    microsites = c.plant_demography["microsites_per_cell"]
    total_sites = int(np.prod(shape)) * microsites
    for k, kind in enumerate(("plant", "insect")):
        sizes = c.plant_sizes if k == 0 else c.insect_sizes
        ids, labels, g = (first[kind + x] for x in ("_ids", "_species", "_genotypes"))
        np.testing.assert_array_equal(ids, np.arange(sum(sizes)))
        np.testing.assert_array_equal(labels, np.repeat(np.arange(len(sizes)), sizes))
        given = c.spatial[f"initial_{kind}_cells"]
        initial_rng = np.random.default_rng(
            np.random.SeedSequence([c.seed, 150, 0 if k == 0 else 2])
        )
        if k == 0:
            if given is None:
                slots = initial_rng.choice(total_sites, len(ids), replace=False)
                cells = slots // microsites
                expected_cells = np.column_stack((cells // shape[1], cells % shape[1]))
            else:
                expected_cells = np.asarray(given, int).reshape(-1, 2)
                occupied, slots = set(), []
                for cell in expected_cells:
                    start = int(cell[0] * shape[1] + cell[1]) * microsites
                    slot = next(i for i in range(start, start + microsites) if i not in occupied)
                    occupied.add(slot)
                    slots.append(slot)
            np.testing.assert_array_equal(first["plant_microsite"], slots)
            close(
                first["plant_resource"], vector(c.resource_capacity, labels, "P"), "founder stock"
            )
        else:
            expected_cells = (
                initial_rng.integers(0, shape, size=(len(ids), 2))
                if given is None
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
        sites = state["plant_microsite"]
        stocks = state["plant_resource"]
        capacities = vector(c.resource_capacity, state["plant_species"], "P")
        assert sites.dtype.kind in "iu" and len(sites) == len(np.unique(sites))
        assert ((sites >= 0) & (sites < total_sites)).all()
        np.testing.assert_array_equal(
            sites // microsites,
            state["plant_cells"][:, 0] * shape[1] + state["plant_cells"][:, 1],
        )
        assert stocks.shape == capacities.shape and np.isfinite(stocks).all()
        assert ((stocks >= 0) & (stocks <= capacities)).all()
        close(state["R"], stocks, "snapshot stock")
        for kind in ("plant", "insect"):
            ids, labels, cells = (state[kind + x] for x in ("_ids", "_species", "_cells"))
            assert len(ids) == len(np.unique(ids)) == len(labels) == len(cells)
            assert ((cells >= 0) & (cells < shape)).all()
            values = phenotype(state[kind + "_genotypes"], kind)
            if kind == "plant":
                close(
                    state["plant_chemical_traits"],
                    values[:, :5],
                    (state["generation"], kind, "phenotype"),
                )
                close(state["D"], values[:, :2], "plant defense")
            else:
                close(state["insect_module_traits"], values, "insect modules")
                I, A = module_profiles(values, c)
                close(state["I"], I, "effective processing")
                close(state["A"], A, "sensory profile")
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
    return rng, spatial_rng, shape, first, microsites
