"""Analytic examples distinguish transmission, mutation, survival, mixing and retention."""

import json

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect.analysis.summaries.allele_change import (
    COMPONENTS,
    decompose_transition,
)
from tests.support import raises


def fixture(kind, parents, labels, children, pairs, inherited, retained, recruited):
    """One-locus diploid records; children follow parents in each taxon's ID space."""
    before, after = {"generation": 0}, {"generation": 1}
    born, tx = [], []
    for own in ("plant", "insect"):
        count = len(parents) if own == kind else 0
        before[own + "_ids"] = np.arange(count)
        before[own + "_species"] = np.array(labels if count else [], dtype=int)
        before[own + "_genotypes"] = np.array(parents if count else [], dtype=int).reshape(
            count, 2, 1
        )
        if own != kind:
            after[own + "_ids"] = np.empty(0, int)
            after[own + "_species"] = np.empty(0, int)
            after[own + "_genotypes"] = np.empty((0, 2, 1), int)
    for i, (child, pair, gamete) in enumerate(zip(children, pairs, inherited, strict=True)):
        born.append(
            dict(
                generation=1,
                kind=kind,
                offspring_id=len(parents) + i,
                maternal_id=pair[0],
                paternal_id=pair[1],
                genotype=json.dumps([[child[0]], [child[1]]]),
            )
        )
        for role in (0, 1):
            tx.append(
                dict(
                    generation=1,
                    kind=kind,
                    offspring_id=len(parents) + i,
                    copy=role,
                    locus_id="L",
                    inherited_allele_id=gamete[role],
                    child_allele_id=child[role],
                    parent_id=pair[role],
                )
            )
    after[kind + "_ids"] = np.array(retained + [len(parents) + i for i in recruited], int)
    after[kind + "_species"] = np.array(
        [labels[i] for i in retained] + [labels[pairs[i][0]] for i in recruited], int
    )
    after[kind + "_genotypes"] = np.array(
        [parents[i] for i in retained] + [children[i] for i in recruited], int
    ).reshape(-1, 2, 1)
    births = pd.DataFrame(
        born, columns="generation kind offspring_id maternal_id paternal_id genotype".split()
    )
    transmissions = pd.DataFrame(
        tx,
        columns="generation kind offspring_id copy locus_id inherited_allele_id child_allele_id parent_id".split(),
    )
    return (
        before,
        after,
        births,
        transmissions,
        {"plant": ["L"], "insect": ["L"]},
        {"plant": 2, "insect": 1},
    )


def test_decomposition_separates_meiosis_mutation_and_family_survival():
    inputs = fixture(
        "insect",
        [[0, 1], [0, 0]],
        [0, 0],
        [[2, 0], [1, 0]],
        [[0, 1], [0, 1]],
        [[1, 0], [1, 0]],
        [],
        [0],
    )
    table = decompose_transition(*inputs).query("scope == 'all'").set_index("allele_id")
    old, new = table.loc[1], table.loc[2]
    np.testing.assert_allclose(
        old[
            ["p_parent", "p_parent_weighted", "p_gamete", "p_birth", "p_recruit", "p_next"]
        ].to_numpy(float),
        [0.25, 0.25, 0.5, 0.25, 0, 0],
    )
    np.testing.assert_allclose(
        old[["delta_meiosis", "delta_mutation", "delta_viability"]].to_numpy(float),
        [0.25, -0.25, -0.25],
    )
    assert new.delta_mutation == 0.25 and new.delta_viability == 0.25 and new.delta_total == 0.5
    assert np.allclose(table.closure_error, 0)
    # No mutation: within-family survival remains a distinct component.
    inputs = fixture(
        "insect",
        [[0, 1], [0, 0]],
        [0, 0],
        [[1, 0], [0, 0]],
        [[0, 1], [0, 1]],
        [[1, 0], [0, 0]],
        [],
        [0],
    )
    row = decompose_transition(*inputs).query("scope == 'all' and allele_id == 1").iloc[0]
    assert row.delta_meiosis == row.delta_mutation == 0 and row.delta_viability == 0.25


def test_plant_source_mixing_and_retention_are_not_reproductive_selection():
    inputs = fixture("plant", [[0, 0], [1, 1]], [0, 1], [[0, 1]], [[0, 1]], [[0, 1]], [0], [0])
    table = decompose_transition(*inputs)
    pooled = table.query("scope == 'all' and allele_id == 1").iloc[0]
    assert pooled.p_parent == 0.5 and pooled.p_next == 0.25
    assert pooled.delta_source_mixing == pooled.delta_parent_reweighting == 0
    assert pooled.delta_retention == -0.25
    maternal = table.query("species == 'P1' and allele_id == 1").iloc[0]
    assert maternal.p_parent == 0 and maternal.p_source_mixed == maternal.p_parent_weighted == 0.5
    assert maternal.delta_source_mixing == 0.25 and maternal.delta_parent_reweighting == 0
    assert maternal.p_next == 0.25 and maternal.N_next == 2
    extinct = table.query("species == 'P2'")
    assert extinct[COMPONENTS + ["p_next", "delta_total"]].isna().all().all()


def test_no_births_preserve_survival_selection_and_extinction_is_undefined():
    inputs = fixture("plant", [[0, 0], [1, 1]], [0, 0], [], [], [], [1], [])
    row = decompose_transition(*inputs).query("scope == 'all' and allele_id == 1").iloc[0]
    assert row.N_birth == row.N_recruit == 0 and row.N_retained == 1
    assert np.isnan(row.p_birth) and np.isnan(row.p_gamete)
    assert row.delta_retention == row.delta_total == 0.5 and row.p_next == 1
    np.testing.assert_array_equal(row[COMPONENTS[:5]].to_numpy(float), 0)
    inputs = fixture("plant", [[0, 0], [1, 1]], [0, 0], [], [], [], [], [])
    table = decompose_transition(*inputs).query("scope == 'all'")
    assert table[COMPONENTS + ["delta_total", "p_next"]].isna().all().all()


def test_decomposition_rejects_inconsistent_parent_child_and_snapshot_keys():
    inputs = fixture("insect", [[0, 1], [0, 0]], [0, 0], [[1, 0]], [[0, 1]], [[1, 0]], [], [0])
    for column, value in (("parent_id", 999), ("child_allele_id", 999), ("copy", 3)):
        broken = inputs[3].copy()
        broken.loc[0, column] = value
        raises(lambda: decompose_transition(inputs[0], inputs[1], inputs[2], broken, *inputs[4:]))
    after = dict(inputs[1], insect_genotypes=np.array([[[0], [0]]]))
    raises(lambda: decompose_transition(inputs[0], after, *inputs[2:]))
