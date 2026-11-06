import numpy as np

from evo_sim.models.plant_insect import simulate
from evo_sim.models.plant_insect.biology.insects.demography import mating_candidates
from evo_sim.models.plant_insect.biology.spatial import disperse, distance
from evo_sim.models.plant_insect.simulation.runtime import community_state
from evo_sim.models.plant_insect.simulation.stages.adult_reproduction import sample_parents
from tests.support import base_runner, config


def test_periodic_mates_and_dispersal():
    cells = np.array([[0, 0], [9, 9], [0, 0], [5, 5]])
    labels = np.array([0, 0, 1, 0])
    candidates = mating_candidates(labels, cells, np.ones(4, bool), np.ones(4, bool), 1, [10, 10])
    assert [x.tolist() for x in candidates] == [[1], [0], [], []]
    assert distance([0, 0], [[9, 9], [5, 5]], [10, 10]).tolist() == [1, 5]
    origin = np.array([[0, 0], [9, 9]])
    a = disperse(origin, 1, [10, 10], np.random.default_rng(8))
    b = disperse(origin, 1, [10, 10], np.random.default_rng(8))
    np.testing.assert_array_equal(a, b)


def test_maternal_sampling_moments_and_survival():
    runner = base_runner()
    runner.config["parameters"]["plant_sizes"] = [2, 1]
    runner.init()
    runtime = runner.state["environment"]["runtime"]
    runner.state["agents"]["plants"]["cells"][:] = [[0, 0], [5, 5], [9, 9]]
    runtime.update = True
    runtime.e = dict(S_plant=np.array([0.6, 0.4, 0.5]), fp=np.ones(3))
    population = community_state(runner.state["agents"]).plants

    def sample():
        counts, candidates, rows, survived = sample_parents(
            population,
            "plant",
            runtime.e["S_plant"],
            runtime.e["fp"],
            runtime,
            0,
        )
        return dict(
            offspring_count=counts, survived_to_reproduce=survived, rows=rows, candidates=candidates
        )

    samples = []
    for _ in range(12000):
        runtime.result["parents"] = []
        values = sample()
        samples.append(values["offspring_count"][0])
        assert np.all(values["offspring_count"][~values["survived_to_reproduce"]] == 0)
    p = 0.6 * (1 - (1 - 0.4) * (1 - 0.5))
    mean = p * 2.0
    variance = mean + p * (1 - p) * 2.0**2
    assert np.isclose(values["rows"][0]["expected_propagules"], mean)
    assert np.isclose(values["rows"][0]["propagule_variance"], variance)
    assert abs(np.mean(samples) - mean) < 0.04
    assert abs(np.var(samples) - variance) < 0.1
    runtime.e["S_plant"] = np.array([1.0, 0.0, 0.0])
    values = sample()
    assert values["survived_to_reproduce"].tolist() == [True, False, False]
    assert not values["offspring_count"].any()
    runtime.e["S_plant"][:] = 1
    values = sample()
    assert values["candidates"][0].tolist() == [1, 2]
    runtime.e["fp"][:] = 0
    assert not sample()["offspring_count"].any()


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
        assert r["snapshots"][1]["insect_genotypes"].shape == (0, 2, 16)
    c = config(plant_sizes=[1, 0], insect_sizes=[0])
    c.plant_offspring_lambda["P1"] = 100.0
    r = simulate(c)
    assert not r["offspring"]  # A lone surviving plant cannot outcross.
    r = simulate(config(plant_sizes=[0, 0], insect_sizes=[0]))
    assert r["metadata"]["status"]["generation"] == 0
    r = simulate(config(plant_sizes=[0, 0]))
    assert all(x["x"] == 0 for x in r["events"])
    assert r["snapshots"][-1]["insect_genotypes"].shape == (0, 2, 16)


def test_default_sexual_founders_can_reproduce():
    # Regression for the old 6-founder configuration with no compatible pair.
    runner = base_runner()
    runner.config["simulation_metadata"]["num_steps_per_episode"] = 1
    runner.init()
    runner.step()
    parents = [
        r for r in runner.result["parents"] if r["generation"] == 0 and r["kind"] == "insect"
    ]
    assert len(parents) == 40
    assert any(r["expected_maternal_offspring"] > 0 for r in parents)
    assert sum(r["offspring_count_as_mother"] for r in parents) > 0
