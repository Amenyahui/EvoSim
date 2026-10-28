"""Stage-causal counterexamples and native cross-generation v0.1.6 checks."""

import copy
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from evo_sim.models.plant_insect import get_runner, simulate
from evo_sim.models.plant_insect.helpers import ecology
from evo_sim.models.plant_insect.helpers.genetics import GeneticState
from evo_sim.models.plant_insect.substeps.lifecycle import source_flows
from tests.support import config, raises


def test_directed_background_sources_are_not_cross_pair_counts():
    # Same-source parents changing host, two sources entering a third host, unknown,
    # and a dead child all have different accounting consequences.
    rows = source_flows(
        np.array([0, 0, 1, -1]),
        np.array([0, 0, 3, 0]),
        np.array([1, 2, 0, 2]),
        np.array([1, 2, 0, 2]),
        np.zeros(4, int),
        np.array([True, True, True, False]),
        4,
        (3, 1),
    )
    lookup = {
        (r["stage"], r["source_host_species"], r["destination_host_species"]): r[
            "genetic_contribution"
        ]
        for r in rows
    }
    assert lookup["recruitment", "P1", "P2"] == 1
    assert lookup["recruitment", "P1", "P3"] == 0.5
    assert lookup["recruitment", "P2", "P3"] == 0.5
    assert lookup["recruitment", "unknown", "P1"] == 0.5
    assert lookup["birth", "P1", "P3"] == 1
    assert sum(r["genetic_contribution"] for r in rows if r["stage"] == "birth") == 4
    assert sum(r["genetic_contribution"] for r in rows if r["stage"] == "recruitment") == 3


def test_ten_births_three_recruits_keep_all_transmissions_and_ids():
    runner = get_runner()
    runner.config["simulation_metadata"]["num_steps_per_episode"] = 1
    runner.config["parameters"].update(
        plant_sizes=[1, 0],
        insect_sizes=[2],
        choice_mode="forced",
        resource_capacity={"P1": 3.0, "P2": 0.0},
        insect_cost=1.0,
    )
    runner.config["parameters"]["spatial"].update(
        initial_plant_cells=[[0, 0]],
        initial_insect_cells=[[0, 0], [0, 0]],
        insect_dispersal_radius=0,
    )
    runner.init()
    runtime = runner.state["environment"]["runtime"]
    plants, insects = runner.state["agents"]["plants"], runner.state["agents"]["insects"]
    plants["defense"][:] = 0
    plants["chemical_traits"][:] = 0
    insects["sex"][:] = [0, 1]
    insects["resistance"][:] = [[0, 0], [1, 1]]
    runtime.rng[3][1][0] = SimpleNamespace(poisson=lambda mu: np.full(len(mu), 10, int))
    runtime.rng[2][1][0] = SimpleNamespace(random=lambda n: np.linspace(0.05, 0.95, n))
    runner.step()
    result = runner.result
    assert result["metadata"]["status"]["status"] == "completed"
    children = [r for r in result["offspring"] if r["kind"] == "insect"]
    assert len(children) == 10 and sum(r["survived_to_adult"] for r in children) == 3
    assert result["snapshots"][1]["insect_ids"].tolist() == [2, 3, 4]
    assert all(
        r["pair_background_status"] == "unknown" and r["cross_background_pair"] is None
        for r in children
    )
    parents = [r for r in result["parents"] if r["kind"] == "insect"]
    assert parents[1]["F_i"] == 0 and parents[1]["eligible_as_donor"]
    assert [r["offspring_count_as_mother"] for r in parents] == [10, 0]
    assert [r["offspring_count_sired"] for r in parents] == [0, 10]
    assert sum(r["realized_genetic_contribution"] for r in parents) == 10
    assert sum(r["realized_recruitment_contribution"] for r in parents) == 3
    tx = [r for r in result["gamete_transmissions"] if r["kind"] == "insect"]
    assert len(tx) == 10 * 2 * 8
    assert {r["offspring_id"] for r in tx} == {r["offspring_id"] for r in children}
    assert len({r["event_id"] for r in result["oviposition"]}) == 10
    assert {r["event_id"] for r in result["oviposition"]} == {
        r["oviposition_event_id"] for r in children
    }
    for stage, count in (("birth", 10), ("recruitment", 3)):
        assert (
            sum(r["genetic_contribution"] for r in result["host_flows"] if r["stage"] == stage)
            == count
        )
    assert all(
        r["generation"] == 0 for r in result["events"] + result["networks"] + result["parents"]
    )


def test_native_two_rounds_use_maternal_A_child_I_and_no_adult_resampling():
    c = config(
        generations=2,
        insect_sizes=[20],
        choice_mode="forced",
        sP0=1,
        resource_capacity={"P1": 1000.0, "P2": 1000.0},
    )
    c.spatial.update(initial_plant_cells=[[0, 0]] * 40, initial_insect_cells=[[0, 0]] * 20)
    seen_A, seen_I = [], []
    choose, settle = ecology.choose_hosts, ecology.settle_feeding

    def choices(A, *args, **kwargs):
        seen_A.append(A.copy())
        return choose(A, *args, **kwargs)

    def feeding(D, I, *args, **kwargs):
        seen_I.append(I.copy())
        return settle(D, I, *args, **kwargs)

    with (
        patch.object(ecology, "choose_hosts", choices),
        patch.object(ecology, "settle_feeding", feeding),
    ):
        result = simulate(c)
    assert result["metadata"]["status"]["status"] == "completed"
    assert len(seen_A) == len(seen_I) == 2
    catalog = GeneticState(c)
    # Rebuild every born child's I from the exported actual inherited/mutated alleles.
    effects = [np.array(__import__("json").loads(r["effect"])) for r in result["allele_catalog"]]
    for t, snapshot in enumerate(result["snapshots"][:-1]):
        eggs = [r for r in result["oviposition"] if r["reproduction_round"] == t]
        lookup = {int(i): j for j, i in enumerate(snapshot["insect_ids"])}
        np.testing.assert_array_equal(
            seen_A[t], [snapshot["A"][lookup[r["maternal_id"]]] for r in eggs]
        )
        children = [
            r for r in result["offspring"] if r["kind"] == "insect" and r["reproduction_round"] == t
        ]
        derived = []
        for child in children:
            g = np.array(__import__("json").loads(child["genotype"]), int)
            raw = np.array(catalog.spec["background"]["insect"], float)
            for allele in g.ravel():
                raw += catalog.spec["phenotype_map"]["effect_scale"] * effects[allele]
            derived.append(np.clip(raw, *catalog.spec["phenotype_map"]["bounds"])[: c.n])
        np.testing.assert_allclose(seen_I[t], np.array(derived).reshape(-1, c.n))
        parents = [r for r in result["parents"] if r["generation"] == t and r["kind"] == "insect"]
        assert all(r["S_i"] == 1 and r["survived_to_mating_stage"] for r in parents)
    assert all(np.all(s["insect_birth_host_species"] == -1) for s in result["snapshots"][:1])


def test_zero_resource_no_mates_rejection_terminal_and_atomicity():
    allocation = ecology.allocate_resources(
        np.array([0]), np.array([1.0]), np.array([0.0]), np.array([1.0])
    )
    assert allocation["x"].tolist() == allocation["g"].tolist() == allocation["h"].tolist() == [0.0]
    empty = ecology.allocate_resources(
        np.array([], int), np.array([]), np.array([0.0]), np.array([])
    )
    assert empty["h"].tolist() == [0.0]
    raises(
        lambda: ecology.allocate_resources(
            np.array([0]), np.array([1.0]), np.array([-1.0]), np.array([1.0])
        )
    )
    c = config(
        generations=1,
        insect_sizes=[10],
        choice_mode="forced",
        resource_capacity={"P1": 0.0, "P2": 0.0},
    )
    c.spatial.update(initial_plant_cells=[[0, 0]] * 40, initial_insect_cells=[[0, 0]] * 10)
    result = simulate(c)
    assert result["larvae"] and not any(r["survived_to_adult"] for r in result["larvae"])
    assert all(r["x"] == 0 for r in result["larvae"])
    assert not len(result["snapshots"][-1]["insect_ids"])
    c.sexual_reproduction["insect_founder_female_probability"] = 1.0
    result = simulate(c)
    assert not result["oviposition"] and not result["larvae"]
    assert all(r["potential_egg_count"] == 0 for r in result["parents"] if r["kind"] == "insect")
    c.sexual_reproduction["insect_founder_female_probability"] = 0.5
    c.outside_utility = 1e6
    c.choice_mode = "outside"
    result = simulate(c)
    assert result["oviposition"] and not result["larvae"]
    assert all(not r["accepted"] for r in result["oviposition"])
    runner = get_runner()
    runner.config["simulation_metadata"]["num_steps_per_episode"] = 0
    runner.init()
    runtime = runner.state["environment"]["runtime"]
    streams = copy.deepcopy(runtime.rng)
    runner.step()
    assert not any(
        runner.result[k]
        for k in ("events", "parents", "larvae", "plants", "networks", "oviposition")
    )
    for phase in streams:
        for left, right in zip(streams[phase], runtime.rng[phase]):
            for a, b in zip(left, right):
                assert a.bit_generator.state == b.bit_generator.state
    runner = get_runner()
    runner.init()
    runtime = runner.state["environment"]["runtime"]
    streams = copy.deepcopy(runtime.rng)
    with patch(
        "evo_sim.models.plant_insect.substeps.lifecycle.source_flows",
        side_effect=ValueError("late ledger failure"),
    ):
        runner.step(1)
    assert runner.result["metadata"]["status"]["status"] == "arithmetic_error"
    assert not any(
        runner.result[k]
        for k in ("events", "parents", "larvae", "plants", "networks", "offspring", "host_flows")
    )
    for phase in streams:
        for left, right in zip(streams[phase], runtime.rng[phase]):
            for a, b in zip(left, right):
                assert a.bit_generator.state == b.bit_generator.state
