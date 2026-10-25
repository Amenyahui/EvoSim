"""Independent v0.1.2 ledger checks, using snapshots and child identities only.

Do not import production observation/contribution helpers here: their output is
what this module checks, not an independent source of expected values.
"""

import numpy as np
import pandas as pd

ALLELE_KEYS = ("generation", "offspring_generation", "kind", "species", "locus_id", "allele_id")
ALLELE_VALUES = (
    "N_parent",
    "N_offspring",
    "p_parent",
    "p_parent_weighted",
    "p_offspring",
    "delta_parent_contribution",
    "delta_transmission",
)
NETWORK_VALUES = (
    "maternal_offspring",
    "paternal_offspring",
    "genetic_contribution",
    "maternal_share",
    "paternal_share",
    "genetic_contribution_share",
)


def close(actual, expected, context):
    assert np.allclose(
        np.asarray(actual, float),
        np.asarray(expected, float),
        rtol=1e-11,
        atol=1e-14,
        equal_nan=True,
    ), context


def require_columns(frame, columns, name):
    assert set(columns) <= set(frame.columns), (
        name,
        "missing columns",
        set(columns) - set(frame.columns),
    )


def unique_keys(frame, columns, name):
    require_columns(frame, columns, name)
    assert frame[list(columns)].notna().all().all(), (name, "missing key")
    assert not frame.duplicated(list(columns)).any(), (name, "duplicate key")


def boolean(values, name):
    mapped = values.map({True: True, False: False, "True": True, "False": False})
    assert mapped.notna().all(), (name, "invalid boolean")
    return mapped.to_numpy(bool)


def audit_contributions(path, c, states, tables):
    """Check complete committed transitions, including extinction and empty ledgers."""
    allele = pd.read_csv(path / "allele_contributions.csv", float_precision="round_trip")
    require_columns(allele, ALLELE_KEYS + ALLELE_VALUES, "allele_contributions")
    unique_keys(allele, ALLELE_KEYS, "allele_contributions")
    parents, children = tables["parents"], tables["offspring"]
    populations, networks = tables["populations"], tables["networks"]
    require_columns(parents, ("realized_genetic_contribution", "contribution_committed"), "parents")
    require_columns(
        networks, ("resource_share", "contribution_committed") + NETWORK_VALUES, "networks"
    )
    for name, frame, keys in (
        ("parents", parents, ("generation", "kind", "individual_id")),
        ("offspring", children, ("generation", "kind", "offspring_id")),
        ("populations", populations, ("generation", "kind", "species")),
        ("networks", networks, ("generation", "insect_species", "plant_species")),
    ):
        unique_keys(frame, keys, name)
        assert set(frame.generation) <= {int(s["generation"]) for s in states}, (
            name,
            "unknown generation",
        )
    assert set(parents.kind) <= {"plant", "insect"}
    assert set(children.kind) <= {"plant", "insect"}
    expected_alleles = {}
    transitions = 0
    for index, state in enumerate(states):
        t = int(state["generation"])
        pop = populations[populations.generation == t]
        par = parents[parents.generation == t]
        net = networks[networks.generation == t]
        update = boolean(pop.update_applied, "populations.update_applied")
        committed = bool(len(update) and update.all())
        assert not update.any() or committed, (t, "partial commit")
        assert committed == (index + 1 < len(states)), (t, "commit/snapshot mismatch")
        assert (
            boolean(par.contribution_committed, "parents.contribution_committed") == committed
        ).all()
        assert (
            boolean(net.contribution_committed, "networks.contribution_committed") == committed
        ).all()
        if not committed:
            assert par.realized_genetic_contribution.isna().all(), (
                t,
                "uncommitted parent contribution",
            )
            assert net[list(NETWORK_VALUES)].isna().all().all(), (
                t,
                "uncommitted network contribution",
            )
        following = states[index + 1] if committed else None
        kid = children[children.generation == t + 1]
        if not committed:
            assert kid.empty
        else:
            transitions += 1
        weights_by_kind = {}
        for kind, sizes, prefix in (("plant", c.plant_sizes, "P"), ("insect", c.insect_sizes, "H")):
            ids, labels = state[kind + "_ids"], state[kind + "_species"]
            rows = par[par.kind == kind]
            assert set(rows.individual_id) <= set(ids), (t, kind, "unknown parent")
            if committed:
                assert set(rows.individual_id) == set(ids), (t, kind, "missing parent")
            lookup = {int(v): i for i, v in enumerate(ids)}
            for row in rows.itertuples():
                assert row.species == f"{prefix}{int(labels[lookup[row.individual_id]]) + 1}"
            if not committed:
                continue
            rows = rows.set_index("individual_id").reindex(ids)
            offspring = kid[kid.kind == kind]
            next_ids, next_labels = following[kind + "_ids"], following[kind + "_species"]
            assert set(offspring.offspring_id) == set(next_ids), (t, kind, "missing offspring")
            next_lookup = dict(zip(map(int, next_ids), map(int, next_labels)))
            maternal, paternal = np.zeros(len(ids)), np.zeros(len(ids))
            for child in offspring.itertuples():
                assert child.maternal_id in lookup and child.paternal_id in lookup
                m, f = lookup[child.maternal_id], lookup[child.paternal_id]
                assert m != f
                s = next_lookup[child.offspring_id]
                assert labels[m] == labels[f] == s
                assert child.species == f"{prefix}{s + 1}"
                maternal[m] += 1
                paternal[f] += 1
            weight = (maternal + paternal) / 2
            close(rows.offspring_count_as_mother, maternal, (t, kind, "mother count"))
            close(rows.offspring_count_sired, paternal, (t, kind, "sired count"))
            close(rows.realized_genetic_contribution, weight, (t, kind, "realized contribution"))
            weights_by_kind[kind] = (maternal, paternal, weight)
            loci = [l for l in c.genetics["loci"] if l["organism"] == kind]
            for s in range(len(sizes)):
                species = f"{prefix}{s + 1}"
                selected, next_selected = labels == s, next_labels == s
                n, n_next = int(selected.sum()), int(next_selected.sum())
                close(weight[selected].sum(), n_next, (t, species, "contribution conservation"))
                population = pop[(pop.kind == kind) & (pop.species == species)]
                assert len(population) == 1
                close(
                    population.iloc[0][["N_start", "N_offspring", "N_next"]],
                    [n, n_next, n_next],
                    (t, species, "population counts"),
                )
                parent_g = state[kind + "_genotypes"][selected]
                child_g = following[kind + "_genotypes"][next_selected]
                for l, locus in enumerate(loci):
                    active = np.union1d(parent_g[:, :, l], child_g[:, :, l])
                    for a in active:
                        z = (parent_g[:, :, l] == a).sum(axis=1) / 2
                        p = float(z.mean()) if n else np.nan
                        weighted = float(weight[selected] @ z / n_next) if n_next else np.nan
                        p_next = (
                            float(np.count_nonzero(child_g[:, :, l] == a) / (2 * n_next))
                            if n_next
                            else np.nan
                        )
                        key = (t, t + 1, kind, species, locus["locus_id"], int(a))
                        expected_alleles[key] = (
                            n,
                            n_next,
                            p,
                            weighted,
                            p_next,
                            weighted - p,
                            p_next - weighted,
                        )
        events = tables["events"][tables["events"].generation == t]
        for h in range(len(c.insect_sizes)):
            species = f"H{h + 1}"
            edges = net[net.insect_species == species]
            if net.empty:
                continue
            assert set(edges.plant_species) == {f"P{p + 1}" for p in range(len(c.plant_sizes))}
            selected = events[events.species == species]
            total_resource = float(selected.x.sum())
            mother_total = father_total = genetic_total = 0.0
            if committed:
                mask = state["insect_species"] == h
                m, f, w = weights_by_kind["insect"]
                mother_total, father_total, genetic_total = (
                    m[mask].sum(),
                    f[mask].sum(),
                    w[mask].sum(),
                )
            host_lookup = dict(zip(events.insect_id, events.plant_species))
            for edge in edges.itertuples():
                resource = float(
                    selected.loc[selected.plant_species == edge.plant_species, "x"].sum()
                )
                close(edge.Q_actual, resource, (t, species, edge.plant_species, "resource"))
                close(
                    edge.resource_share,
                    resource / total_resource if total_resource else np.nan,
                    (t, species, "resource share"),
                )
                if committed:
                    on_host = np.array(
                        [
                            label == h and host_lookup.get(int(i)) == edge.plant_species
                            for i, label in zip(state["insect_ids"], state["insect_species"])
                        ],
                        bool,
                    )
                    mother, father, genetic = m[on_host].sum(), f[on_host].sum(), w[on_host].sum()
                    expected = (
                        mother,
                        father,
                        genetic,
                        mother / mother_total if mother_total else np.nan,
                        father / father_total if father_total else np.nan,
                        genetic / genetic_total if genetic_total else np.nan,
                    )
                    close(
                        [getattr(edge, k) for k in NETWORK_VALUES],
                        expected,
                        (t, species, edge.plant_species, "host contribution"),
                    )
            if committed:
                close(
                    edges.genetic_contribution.sum(),
                    genetic_total,
                    (t, species, "host conservation"),
                )
    actual_keys = set(allele[list(ALLELE_KEYS)].itertuples(index=False, name=None))
    assert actual_keys == set(expected_alleles), (
        "allele_contributions",
        "missing/extra keys",
        set(expected_alleles) - actual_keys,
        actual_keys - set(expected_alleles),
    )
    for row in allele.itertuples(index=False):
        key = tuple(getattr(row, k) for k in ALLELE_KEYS)
        close(
            [getattr(row, k) for k in ALLELE_VALUES],
            expected_alleles[key],
            ("allele_contributions", key),
        )
    return dict(
        contribution_audit=True,
        contribution_transitions=transitions,
        allele_contribution_rows=len(allele),
    )
