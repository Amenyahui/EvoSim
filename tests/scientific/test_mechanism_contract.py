"""Ecological probability, parental roles and native genetic replay contracts."""

import copy
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from evo_sim.models.plant_insect import get_runner, simulate
from evo_sim.models.plant_insect.audit.run import audit_run
from evo_sim.models.plant_insect.batch import save_result
from evo_sim.models.plant_insect.biology.insects import feeding as _feeding
from evo_sim.models.plant_insect.simulation.stages.adult_reproduction import (
    draw_fathers,
    sample_parents,
)
from tests.scientific.test_functional_modules import CONFIG, module_config


def test_competing_larvae_conserve_resources_and_fraction():
    host = np.repeat(np.arange(3), [42, 61, 1000])
    capacity = np.array([30.0, 30.0, 300.0])
    for order in (np.arange(len(host)), np.arange(len(host))[::-1]):
        result = _feeding.allocate_resources(host[order], capacity, np.ones(len(host)))
        assert (result["X"] <= capacity).all()
        np.testing.assert_allclose(result["X"], capacity, rtol=1e-12)
        for plant in range(3):
            shares = result["x"][host[order] == plant]
            assert (shares == shares[0]).all()
    result = _feeding.allocate_resources(np.array([0, -1]), [1], [2, 2])
    np.testing.assert_array_equal(result["d"], [2, 0])
    np.testing.assert_array_equal(result["x"], [1, 0])


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
    counts, pools, rows, _ = sample_parents(pop, "insect", np.ones(6), np.ones(6), runtime, 0)
    assert [p.tolist() for p in pools] == [[1, 5], [], [], [], [1, 5], []]
    assert not counts[pop.sex == 1].any()
    assert rows[0]["candidate_mate_count"] == 2
    mothers = np.array([0, 4, 0, 4, 0, 4])
    expected_rng = copy.deepcopy(runtime.rng[5][1][0])
    expected = np.array([expected_rng.choice(pools[m]) for m in mothers])
    np.testing.assert_array_equal(draw_fathers(pop, mothers, pools, runtime, 1), expected)
    assert runtime.rng[5][1][0].bit_generator.state == expected_rng.bit_generator.state
    pop.sex[:] = 0
    counts, pools, rows, _ = sample_parents(pop, "insect", np.ones(6), np.ones(6), runtime, 0)
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
        np.testing.assert_allclose(
            probabilities,
            [
                r["nutrition_survival"] * r["toxicity_survival"] * r["antinutrition_survival"]
                for r in rows
            ],
        )
        assert all(r["g"] == r["x"] / r["dmax"] for r in rows)
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
    c = module_config(choice_mode="forced")
    original = _feeding.chemical_survival

    def no_survival(*args, **kwargs):
        result = original(*args, **kwargs)
        result["survival_probability"] = np.zeros_like(result["survival_probability"])
        return result

    with patch.object(_feeding, "chemical_survival", side_effect=no_survival):
        result = simulate(c)
    assert result["metadata"]["status"]["status"] == "completed"
    assert result["events"] and any(row["g"] > 0 for row in result["events"])
    assert all(row["W_insect"] == 0 and not row["survived_to_adult"] for row in result["events"])
    assert not len(result["snapshots"][-1]["insect_ids"])
