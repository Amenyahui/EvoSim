"""Hand-check matched opportunities and transmitted QTL ancestry."""

import json

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect.analysis.lineage_observations import ancestry, opportunities
from tests.support import raises


def test_matched_opportunities_weights_unknown_and_no_mate():
    states = {
        0: dict(
            insect_ids=np.arange(6),
            insect_species=np.zeros(6, int),
            insect_sex=np.array([0, 0, 0, 1, 1, 1]),
            insect_birth_host_species=np.array([0, 0, 0, 0, 1, 1]),
        ),
        1: {},
    }
    mappings = {(0, "insect"): {0: "A", 1: "A", 2: "A", 3: "A", 4: "B", 5: "B"}}
    parents = []
    for ident, pool, mu, eggs, born in [(0, [3], 2, 2, 1), (1, [3, 4], 1, 1, 1), (2, [], 4, 0, 0)]:
        parents.append(
            dict(
                generation=0,
                kind="insect",
                individual_id=ident,
                eligible_as_mother=True,
                candidate_mate_ids=json.dumps(pool),
                candidate_mate_probabilities=json.dumps([1 / len(pool) for _ in pool]),
                conditional_mu=mu,
                potential_egg_count=eggs,
                offspring_count_as_mother=born,
            )
        )
    parents.extend(
        dict(generation=0, kind="insect", individual_id=i, eligible_as_mother=False)
        for i in (3, 4, 5)
    )
    table = pd.DataFrame(opportunities(states, pd.DataFrame(parents), mappings))
    row = table.query(
        "group_system == 'candidate_cluster' and mother_group == 'A' and target_group == 'B' and weight_type == 'conditional_mu'"
    ).iloc[0]
    assert np.isclose(row.p_local, 1 / 6) and np.isclose(row.p_global, 2 / 3)
    assert np.isclose(row.ri_relative, 0.75)
    assert row.no_mate_mothers == 1 and row.no_mate_weight == 4 and row.matched_weight == 3
    row = table.query(
        "group_system == 'candidate_cluster' and mother_group == 'A' and target_group == 'B' and weight_type == 'offspring_count_as_mother'"
    ).iloc[0]
    assert row.p_local == 0.25
    mappings[0, "insect"][5] = "unassigned"
    table = pd.DataFrame(opportunities(states, pd.DataFrame(parents), mappings))
    row = table.query(
        "group_system == 'candidate_cluster' and target_group == 'B' and weight_type == 'conditional_mu'"
    ).iloc[0]
    assert row.matched_mothers == 0 and row.unknown_global_edges == 3 and np.isnan(row.ri_relative)
    parents[0]["eligible_as_mother"] = "invalid"
    raises(lambda: opportunities(states, pd.DataFrame(parents), mappings))


def test_qtl_copy_ancestry_recombination_mutation_and_missing():
    states = {
        0: dict(
            insect_ids=np.array([0, 1]),
            insect_birth_host_species=np.array([0, 1]),
            insect_genotypes=np.array([[[0, 1], [2, 3]], [[4, 5], [6, 7]]]),
        ),
        1: dict(insect_ids=np.array([2]), insect_genotypes=np.array([[[0, 3], [6, 8]]])),
    }
    mappings = {(0, "insect"): {0: "A", 1: "B"}, (1, "insect"): {2: "C"}}
    children = pd.DataFrame(
        [
            dict(
                kind="insect",
                reproduction_round=0,
                offspring_id=2,
                maternal_id=0,
                paternal_id=1,
                birth_host_species="P1",
                genotype=json.dumps([[0, 3], [6, 8]]),
            )
        ]
    )
    loci = pd.DataFrame([dict(organism="insect", locus_id=x) for x in ("L1", "L2")])
    transmissions = pd.DataFrame(
        [
            dict(
                generation=1,
                kind="insect",
                offspring_id=2,
                role=role,
                parent_id=parent,
                locus_id=locus,
                transmitted_homolog_before_mutation=homolog,
                inherited_allele_id=old,
                child_allele_id=new,
            )
            for role, parent, locus, homolog, old, new in [
                ("maternal", 0, "L1", 0, 0, 0),
                ("maternal", 0, "L2", 1, 3, 3),
                ("paternal", 1, "L1", 1, 6, 6),
                ("paternal", 1, "L2", 1, 7, 8),
            ]
        ]
    )
    result = pd.DataFrame(ancestry(states, children, mappings, 0, loci, transmissions))
    birth = result.query(
        "source_system == 'candidate_cluster' and stage == 'birth' and source_group == 'A'"
    ).iloc[0]
    assert birth.target_group == "P1" and birth.fraction == 0.5 and birth.total_qtl_copies == 4
    baseline = result.query(
        "source_system == 'candidate_cluster' and stage == 'reference' and target_group == 'B' and source_group == 'A'"
    ).iloc[0]
    assert baseline.fraction == 0
    raises(lambda: ancestry(states, children, mappings, 0, loci, transmissions.iloc[:-1]))
    invalid = transmissions.copy()
    invalid.loc[0, "inherited_allele_id"] = 99
    raises(lambda: ancestry(states, children, mappings, 0, loci, invalid))


def test_two_generation_copy_paths_and_empty_population():
    states = {
        0: dict(
            insect_ids=np.array([0, 1]),
            insect_birth_host_species=np.array([0, -1]),
            insect_genotypes=np.array([[[0], [1]], [[2], [3]]]),
        ),
        1: dict(insect_ids=np.array([2, 3]), insect_genotypes=np.array([[[0], [2]], [[1], [3]]])),
        2: dict(insect_ids=np.array([4]), insect_genotypes=np.array([[[0], [1]]])),
        3: dict(insect_ids=np.array([], dtype=int), insect_genotypes=np.empty((0, 2, 1), int)),
    }
    mappings = {
        (0, "insect"): {0: "A", 1: "B"},
        (1, "insect"): {2: "C", 3: "C"},
        (2, "insect"): {4: "D"},
        (3, "insect"): {},
    }
    children = pd.DataFrame(
        [
            dict(
                kind="insect",
                reproduction_round=t,
                offspring_id=i,
                maternal_id=m,
                paternal_id=f,
                birth_host_species="P1",
                genotype=json.dumps(g),
            )
            for t, i, m, f, g in [
                (0, 2, 0, 1, [[0], [2]]),
                (0, 3, 0, 1, [[1], [3]]),
                (1, 4, 2, 3, [[0], [1]]),
            ]
        ]
    )
    loci = pd.DataFrame([dict(organism="insect", locus_id="L")])
    transmissions = pd.DataFrame(
        [
            dict(
                generation=t,
                kind="insect",
                offspring_id=i,
                role=role,
                parent_id=parent,
                locus_id="L",
                transmitted_homolog_before_mutation=h,
                inherited_allele_id=allele,
                child_allele_id=allele,
            )
            for t, i, role, parent, h, allele in [
                (1, 2, "maternal", 0, 0, 0),
                (1, 2, "paternal", 1, 0, 2),
                (1, 3, "maternal", 0, 1, 1),
                (1, 3, "paternal", 1, 1, 3),
                (2, 4, "maternal", 2, 0, 0),
                (2, 4, "paternal", 3, 0, 1),
            ]
        ]
    )
    result = pd.DataFrame(ancestry(states, children, mappings, 0, loci, transmissions))
    selected = result.query(
        "source_system == 'candidate_cluster' and source_group == 'A' and generation == 2 and stage == 'recruitment'"
    ).iloc[0]
    assert selected.fraction == 1  # Both parents are half A, but both transmitted their A copy.
    unknown = result.query(
        "source_system == 'birth_host' and generation == 1 and stage == 'recruitment'"
    ).iloc[0]
    assert unknown.unknown_fraction == 0.5
    empty = result.query("generation == 3 and stage == 'recruitment'")
    assert (empty.status == "empty_population").all() and empty.fraction.isna().all()


def test_opportunity_baselines_and_pool_boundaries():
    state = dict(
        insect_ids=np.arange(4),
        insect_species=np.array([0, 0, 0, 1]),
        insect_sex=np.array([0, 1, 1, 1]),
        insect_birth_host_species=np.array([0, 0, 1, 1]),
    )
    mappings = {(0, "insect"): {0: "A", 1: "A", 2: "B", 3: "B"}}
    parents = [
        dict(
            generation=0,
            kind="insect",
            individual_id=i,
            eligible_as_mother=i == 0,
            candidate_mate_ids=json.dumps([1, 2]),
            candidate_mate_probabilities=json.dumps([0.5, 0.5]),
            conditional_mu=1,
            potential_egg_count=0,
            offspring_count_as_mother=1,
        )
        for i in range(4)
    ]

    def calculate():
        frame = pd.DataFrame(opportunities({0: state, 1: {}}, pd.DataFrame(parents), mappings))
        return frame.query(
            "group_system == 'candidate_cluster' and target_group == '__all_cross__' and weight_type == 'conditional_mu'"
        ).iloc[0]

    assert calculate().ri_relative == 0  # The male in another pool is excluded.
    parents[0]["candidate_mate_ids"] = "[2]"
    parents[0]["candidate_mate_probabilities"] = "[1]"
    assert calculate().ri_relative == -1
    mappings[0, "insect"][2] = "A"
    assert calculate().p_global == 0 and np.isnan(calculate().ri_relative)
    parents[0]["conditional_mu"] = 0
    assert np.isnan(calculate().p_local)
    parents[0]["candidate_mate_ids"] = "[3]"
    raises(lambda: calculate())


def test_target_with_no_males_is_a_zero_baseline():
    state = dict(
        insect_ids=np.array([0, 1, 2]),
        insect_species=np.zeros(3, int),
        insect_sex=np.array([0, 0, 1]),
        insect_birth_host_species=np.array([0, 1, 0]),
    )
    mapping = {(0, "insect"): {0: "A", 1: "B", 2: "A"}}
    parents = pd.DataFrame(
        [
            dict(
                generation=0,
                kind="insect",
                individual_id=i,
                eligible_as_mother=i < 2,
                candidate_mate_ids="[2]",
                candidate_mate_probabilities="[1]",
                conditional_mu=1,
                potential_egg_count=1,
                offspring_count_as_mother=1,
            )
            for i in range(3)
        ]
    )
    rows = pd.DataFrame(opportunities({0: state, 1: {}}, parents, mapping))
    row = rows.query(
        "group_system == 'candidate_cluster' and mother_group == 'A' and target_group == 'B' and weight_type == 'conditional_mu'"
    ).iloc[0]
    assert row.p_local == row.p_global == 0 and np.isnan(row.ri_relative)
