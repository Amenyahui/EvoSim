import numpy as np

from evo_sim.models.plant_insect import get_runner, simulate
from evo_sim.models.plant_insect.helpers.demography import maternal_moments, mating_candidates
from evo_sim.models.plant_insect.helpers.spatial import disperse, distance, local_competition
from evo_sim.models.plant_insect.substeps.demography import SampleMaternallyProducedOffspring
from tests.support import config


def test_periodic_mates_and_competition():
    cells = np.array([[0, 0], [9, 9], [0, 0], [5, 5]])
    labels = np.array([0, 0, 1, 0])
    candidates = mating_candidates(labels, cells, np.ones(4, bool), np.ones(4, bool), 1, [10, 10])
    assert [x.tolist() for x in candidates] == [[1], [0], [], []]
    assert distance([0, 0], [[9, 9], [5, 5]], [10, 10]).tolist() == [1, 5]
    c = config()
    c.spatial.update(grid_shape=[1, 1], competition_radius=1)
    alpha = np.array([[0.2, 0.1], [0.3, 0.4]])
    phi = local_competition(np.zeros((4, 2), int), np.array([0, 1, 0, 1]), alpha, c)
    np.testing.assert_allclose(phi, [1 / 1.6, 1 / 2.4, 1 / 1.6, 1 / 2.4])
    origin = np.array([[0, 0], [9, 9]])
    a = disperse(origin, 1, [10, 10], np.random.default_rng(8))
    b = disperse(origin, 1, [10, 10], np.random.default_rng(8))
    np.testing.assert_array_equal(a, b)


def test_maternal_sampling_moments_and_survival():
    runner = get_runner()
    runner.config["parameters"]["plant_sizes"] = [3, 0]
    runner.init()
    runtime = runner.state["environment"]["runtime"]
    runner.state["agents"]["plants"]["cells"][:] = 0
    runtime.update = True
    runtime.plant_phi = np.ones(3)
    runtime.e = dict(S_plant=np.array([0.6, 0.4, 0.5]), fp=np.ones(3))
    component = SampleMaternallyProducedOffspring(
        runner.config,
        {"runtime": "environment/runtime", "agents": "agents"},
        ["offspring"],
        {"fixed": {"taxon": 0}},
    )
    samples = []
    for _ in range(12000):
        runtime.result["parents"] = []
        values = component.forward(runner.state, {})["offspring"]
        samples.append(values["offspring_count"][0])
        assert np.all(values["offspring_count"][~values["survived_to_reproduce"]] == 0)
    q, mean, variance = maternal_moments(
        np.array([0.6, 0.4, 0.5]),
        np.full(3, 2.0),
        [np.array([1, 2]), np.array([0]), np.array([], int)],
        np.ones(3, bool),
    )
    np.testing.assert_allclose(q, [0.7, 0.6, 0])
    assert abs(np.mean(samples) - mean[0]) < 0.04
    assert abs(np.var(samples) - variance[0]) < 0.1
    print(
        "maternal moments: N=12000",
        np.mean(samples),
        np.var(samples),
        "expected",
        mean[0],
        variance[0],
    )
    runtime.e["S_plant"] = np.array([1.0, 0.0, 0.0])
    values = component.forward(runner.state, {})["offspring"]
    assert not values["offspring_count"].any()
    assert values["survived_to_reproduce"].tolist() == [True, False, False]
    runtime.e["S_plant"][:] = 1
    runtime.e["fp"][:] = 0
    assert not component.forward(runner.state, {})["offspring"]["offspring_count"].any()


def test_single_sex_singleton_and_extinction():
    for probability in (0.0, 1.0):
        c = config()
        c.sexual_reproduction["insect_founder_female_probability"] = probability
        r = simulate(c)
        assert all(
            x["offspring_count_as_mother"] in (0, None)
            for x in r["parents"]
            if x["kind"] == "insect"
        )
        assert r["snapshots"][1]["insect_genotypes"].shape == (0, 2, 8)
    r = simulate(config(plant_sizes=[1, 0], insect_sizes=[0]))
    assert not r["offspring"] and r["metadata"]["status"]["status"] == "all_species_extinct"
    r = simulate(config(plant_sizes=[0, 0], insect_sizes=[0]))
    assert r["metadata"]["status"]["generation"] == 0
    r = simulate(config(plant_sizes=[0, 0]))
    assert all(x["x"] == 0 for x in r["events"])
    assert r["snapshots"][-1]["insect_genotypes"].shape == (0, 2, 8)


def test_default_sexual_founders_can_reproduce():
    # Regression for the old 6-founder configuration with no compatible pair.
    runner = get_runner()
    runner.config["simulation_metadata"]["num_steps_per_episode"] = 1
    runner.init()
    runner.step()
    parents = [
        r for r in runner.result["parents"] if r["generation"] == 0 and r["kind"] == "insect"
    ]
    assert len(parents) == 40
    assert any(r["expected_maternal_offspring"] > 0 for r in parents)
    assert sum(r["offspring_count_as_mother"] for r in parents) > 0
