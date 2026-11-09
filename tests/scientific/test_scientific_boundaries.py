"""Mechanism reference and the scientific limits of resource and choice closures."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np

from evo_sim.models.plant_insect import Config, get_runner, simulate
from evo_sim.models.plant_insect.biology import spatial
from evo_sim.models.plant_insect.biology import traits as modules
from evo_sim.models.plant_insect.biology.insects import sensory as _sensory

CONFIG = (
    Path(__file__).resolve().parents[2]
    / "evo_sim/models/plant_insect/configs/current_single_ancestor.yaml"
)


def reference_config():
    return Config(**dict(get_runner(CONFIG).config["parameters"], generations=0))


def test_completed_feeding_window_does_not_repeat_host_death_on_larvae():
    c = reference_config()
    c.generations, c.plant_sizes, c.insect_sizes = 1, [2], [20]
    c.resource_capacity, c.choice_mode = {"P1": 100.0}, "forced"
    c.plant_offspring_lambda = {"P1": c.plant_offspring_lambda["P1"]}
    c.parameter_provenance["plant_offspring_lambda"] = {
        "P1": c.parameter_provenance["plant_offspring_lambda"]["P1"]
    }
    c.genetics["founder_frequencies"].pop("P2", None)
    c.plant_demography.update(background_mortality_rate=1000.0, microsites_per_cell=2)
    c.spatial.update(initial_plant_cells=[[0, 0]] * 2, initial_insect_cells=[[0, 0]] * 20)
    result = simulate(c)
    assert result["metadata"]["status"]["status"] == "completed"
    dead_hosts = {
        parent["individual_id"]
        for parent in result["parents"]
        if parent["kind"] == "plant" and not parent["survived_to_mating_stage"]
    }
    assert len(dead_hosts) == 2
    assert any(
        larva["host_id"] in dead_hosts and larva["survived_to_adult"] for larva in result["events"]
    )
    diagnostic = result["metadata"]["founder_feeding_diagnostics"]
    assert diagnostic["pair_count"] == 40 and 0 <= diagnostic["survival_mean"] <= 1


def test_relative_avoidance_and_outside_logit_are_distinct():
    c = reference_config()
    sensory = np.array([0.6, 0.4, 0, 1, 0, 1, 0, 1, 0.5, 0.5])
    scores = modules.scores(np.array([[0, 0], [0, 1]]), sensory, c)
    assert 0 < scores[1] < scores[0]
    probabilities = _sensory.choice_probabilities(scores, c.tau, c.outside_utility, "outside")
    assert probabilities[1] < probabilities[0] and probabilities[-1] < 1 / 3
    assert _sensory.choice_probabilities(scores, 1, 10, "outside")[-1] > 0.99
    assert _sensory.choice_probabilities(scores, 0, 10, "outside")[-1] > 0.99
    for valence in ([1, 1], [-1, -1], [0, 0]):
        modules.validate(dict(c.functional_modules, sensory_valence=valence))
    c.b = 0
    assert modules.scores(np.array([[0, 1]]), sensory, c)[0] < 0


def test_periodic_displacement_residues_are_not_triangular():
    steps = np.arange(-6, 7)
    displacements = np.column_stack((steps, steps))
    rng = SimpleNamespace(integers=lambda low, high, size: displacements.copy())
    positions = spatial.disperse(np.zeros_like(displacements), 6, [10, 10], rng)
    counts = np.bincount(positions[:, 0], minlength=10)
    assert counts.max() - counts.min() == 1
    np.testing.assert_array_equal(counts, [1, 1, 1, 1, 2, 2, 2, 1, 1, 1])
