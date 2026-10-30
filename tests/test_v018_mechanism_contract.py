"""Ecological probability, parental roles and native genetic replay contracts."""

import copy
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from evo_sim.models.plant_insect import get_runner, simulate
from evo_sim.models.plant_insect.audit import audit_run
from evo_sim.models.plant_insect.batch import save_result
from evo_sim.models.plant_insect.helpers import ecology, modules, traits
from evo_sim.models.plant_insect.substeps.lifecycle import draw_fathers, sample_parents
from tests.support import raises
from tests.test_v017_modules import CONFIG, module_config


def test_larval_survival_probability_boundaries_and_cost_separation():
    probability = ecology.larval_survival_probability
    np.testing.assert_array_equal(probability([0, 0.25, 1]), [0, 0.25, 1])
    assert probability([]).shape == (0,)
    for invalid in (-0.01, 1.01, np.nan, np.inf, -np.inf):
        raises(lambda: probability([invalid]))
    c = module_config()
    values = np.array([[0.6, 0.4, 0, 1, 0.5, 0.2, 0, 1]] * 2)
    I = modules.effective(values, c.functional_modules)[0]
    for resource, expected in ((0.0, 0.0), (2.0, 1.0)):
        result = ecology.settle_feeding(
            np.zeros((1, 2)),
            I,
            np.zeros(1, int),
            np.zeros(2, int),
            np.zeros((1, 2)),
            np.zeros(2, int),
            c,
            np.array([resource]),
            np.ones(2),
            chemical=np.zeros((1, 3)),
            module_traits=values,
        )
        np.testing.assert_array_equal(result["survival_probability"], [expected] * 2)
        np.testing.assert_array_equal(result["survival_probability"], result["x"])
        assert (result["fi"] < 1).all()
        if resource:
            assert (result["wi"] < result["survival_probability"]).all()


def test_competing_larvae_conserve_resources_and_fraction():
    host = np.array([0, 0, 1, 1, -1])
    q = np.array([1.0, 0.5, 1.0, 0.0, np.nan])
    capacity, dmax = np.array([1.0, 0.0]), np.array([2.0, 2.0, 1.0, 1.0, 1.0])
    e = ecology.allocate_resources(host, q, capacity, dmax)
    np.testing.assert_allclose(e["x"], [2 / 3, 1 / 3, 0, 0, 0])
    np.testing.assert_array_equal(e["g"], e["x"] / dmax)
    np.testing.assert_array_equal(ecology.larval_survival_probability(e["g"]), e["g"])
    assert (e["X"] <= capacity).all() and e["x"].sum() <= capacity.sum()
    rng = np.random.default_rng(18)
    for _ in range(30):
        host = rng.integers(0, 4, 100)
        q, capacity, dmax = rng.random(100), rng.random(4), rng.uniform(0.01, 2, 100)
        e = ecology.allocate_resources(host, q, capacity, dmax)
        assert (e["X"] <= capacity).all()
        np.testing.assert_array_equal(e["X"], np.bincount(host, weights=e["x"], minlength=4))
        np.testing.assert_array_equal(e["g"], e["x"] / dmax)
        assert (e["g"] <= q).all()


def test_limiting_component_boundaries_and_processing_monotonicity():
    np.testing.assert_array_equal(traits.consumption([[0, 0], [1, 0]], [[0, 0], [0, 1]]), [1, 0])
    np.testing.assert_array_equal(
        traits.consumption([[0, 1], [1, 1]], [[0, 0.5], [1, 0]]), [0.5, 0]
    )
    assert traits.consumption([[1e-320, 0]], [[1, 0]])[0] == 1
    assert traits.consumption(np.empty((0, 2)), np.empty((0, 2))).shape == (0,)
    for bad in (-1.0, np.nan, np.inf):
        raises(lambda: traits.consumption([[bad, 1]], [[1, 1]]))
        raises(lambda: traits.consumption([[1, 1]], [[bad, 1]]))
    raises(lambda: traits.consumption([], []))
    c = module_config()
    rng = np.random.default_rng(18)
    values, D = rng.random((100, 8)), rng.random((100, 2))
    enhanced = values.copy()
    enhanced[:, :2] = (enhanced[:, :2] + 1) / 2
    I = modules.effective(values, c.functional_modules)[0]
    more = modules.effective(enhanced, c.functional_modules)[0]
    assert (I.sum(axis=1) <= values[:, :2].sum(axis=1)).all()
    assert (more >= I).all()
    assert (traits.consumption(D, more) >= traits.consumption(D, I)).all()


def test_mating_roles_distance_compatibility_and_per_egg_fathers():
    class Population(SimpleNamespace):
        def __len__(self):
            return len(self.ids)

    pop = Population(
        ids=np.arange(6),
        species=np.array([0, 0, 0, 1, 0, 0]),
        sex=np.array([0, 1, 1, 1, 0, 1]),
        cells=np.array([[0, 0], [9, 9], [5, 5], [0, 0], [0, 0], [0, 1]]),
    )
    runtime = SimpleNamespace(
        c=module_config(),
        li=np.array([4.0, 4.0]),
        rng={
            3: {1: [np.random.default_rng(1), np.random.default_rng(2)]},
            5: {1: [np.random.default_rng(3), np.random.default_rng(4)]},
        },
    )
    counts, pools, rows, _ = sample_parents(
        pop, "insect", np.ones(6), np.ones(6), np.ones(6), runtime, 0
    )
    assert [p.tolist() for p in pools] == [[1, 5], [], [], [], [1, 5], []]
    assert not counts[pop.sex == 1].any()
    assert rows[0]["candidate_mate_count"] == 2
    mothers = np.array([0, 4, 0, 4, 0, 4])
    expected_rng = copy.deepcopy(runtime.rng[5][1][0])
    expected = np.array([expected_rng.choice(pools[m]) for m in mothers])
    np.testing.assert_array_equal(draw_fathers(pop, mothers, pools, runtime, 1), expected)
    assert runtime.rng[5][1][0].bit_generator.state == expected_rng.bit_generator.state
    pop.sex[:] = 0
    counts, pools, rows, _ = sample_parents(
        pop, "insect", np.ones(6), np.ones(6), np.ones(6), runtime, 0
    )
    assert not counts.any() and all(not len(p) for p in pools)
    assert all(r["reason"] == "no_compatible_mate" for r in rows)


def test_native_survival_rng_event_consistency_and_replay():
    runner = get_runner(CONFIG)
    runner.config["simulation_metadata"]["num_steps_per_episode"] = 2
    runner.init()
    runtime = runner.state["environment"]["runtime"]
    for generation in range(2):
        expected_rng = copy.deepcopy(runtime.rng[2][1][0])
        runner.step(1)
        rows = [r for r in runner.result["events"] if r["generation"] == generation]
        assert rows
        probabilities = np.array([r["W_insect"] for r in rows])
        np.testing.assert_array_equal(probabilities, [r["g"] for r in rows])
        np.testing.assert_array_equal(probabilities, [r["x"] / r["dmax"] for r in rows])
        np.testing.assert_array_equal(
            [r["survived_to_adult"] for r in rows], expected_rng.random(len(rows)) < probabilities
        )
        assert runtime.rng[2][1][0].bit_generator.state == expected_rng.bit_generator.state
    result = runner.result
    assert result["metadata"]["status"]["status"] == "completed"
    repeated = get_runner(CONFIG)
    repeated.config["simulation_metadata"]["num_steps_per_episode"] = 2
    repeated.init()
    repeated.step()
    for key in ("events", "offspring", "mutations", "gamete_transmissions"):
        assert result[key] == repeated.result[key]
    for first, second in zip(result["snapshots"], repeated.result["snapshots"], strict=True):
        for key in first:
            np.testing.assert_array_equal(first[key], second[key])
    with TemporaryDirectory() as temporary:
        path = Path(temporary) / "run"
        save_result(result, path)
        audit = audit_run(path)
        assert audit["genetic_audit"] and audit["contribution_audit"]
    # Check actual inherited copies against each recorded parent's pre-mutation homolog.
    locus_index = {}
    for row in result["genetic_loci"]:
        key = row["organism"]
        locus_index.setdefault(key, {})[row["locus_id"]] = len(locus_index.get(key, {}))
    for tx in result["gamete_transmissions"]:
        state = result["snapshots"][tx["generation"] - 1]
        kind = tx["kind"]
        parent = np.flatnonzero(state[kind + "_ids"] == tx["parent_id"])[0]
        allele = state[kind + "_genotypes"][
            parent, tx["transmitted_homolog_before_mutation"], locus_index[kind][tx["locus_id"]]
        ]
        assert allele == tx["inherited_allele_id"]


def test_lifecycle_consumes_explicit_probability_and_records_it():
    # Wiring test: a replacement kernel must drive both the draw and event field.
    c = module_config(choice_mode="forced")
    with patch.object(
        ecology, "larval_survival_probability", side_effect=lambda g: np.zeros_like(g)
    ):
        result = simulate(c)
    assert result["metadata"]["status"]["status"] == "completed"
    assert result["events"] and any(r["g"] > 0 for r in result["events"])
    assert all(r["W_insect"] == 0 and not r["survived_to_adult"] for r in result["events"])
    assert not len(result["snapshots"][-1]["insect_ids"])
    assert all(json.loads(r["genotype"]) for r in result["offspring"] if r["kind"] == "insect")
