"""Mechanism reference and the scientific limits of resource and choice closures."""

from itertools import product
from math import comb
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from evo_sim.models.plant_insect import Config, get_runner, simulate
from evo_sim.models.plant_insect.biology import traits as modules
from evo_sim.models.plant_insect.helpers import ecology, spatial, traits, voc
from evo_sim.models.plant_insect.helpers.demography import plant_resource_update

CONFIG = (
    Path(__file__).resolve().parents[1]
    / "evo_sim/models/plant_insect/configs/current_single_ancestor.yaml"
)


def reference_config():
    return Config(**dict(get_runner(CONFIG).config["parameters"], generations=0))


def founder_reference_support():
    """Exact equal-frequency four-copy support; unrelated sensory axes integrate to one."""
    c = reference_config()
    scale = c.genetics["phenotype_map"]["effect_scale"]
    centers = []
    for axis in range(4):
        effects = [
            allele["effect"][axis]
            for allele in c.genetics["alleles"]
            if allele["locus_id"].startswith("insect_") and allele["effect"][axis] != 0
        ]
        assert len(effects) == 4
        centers.append(scale * 4 * np.mean(effects))
    levels = np.arange(-2, 3) * 0.025
    values = np.array(list(product(*(center + levels for center in centers))))
    insect_traits = np.column_stack((values, np.zeros((len(values), 4))))
    defense = np.array(list(product(0.35 + levels, 0.55 + levels)))
    frequency = np.array([comb(4, k) / 16 for k in range(5)])
    weights = (
        np.array([np.prod(frequency[list(indices)]) for indices in product(range(5), repeat=4)])[
            :, None
        ]
        * np.array([np.prod(frequency[list(indices)]) for indices in product(range(5), repeat=2)])[
            None, :
        ]
    )
    return c, np.array(centers), insect_traits, defense, weights


def test_current_founder_reference_crosses_both_defense_thresholds():
    c, centers, values, defense, weights = founder_reference_support()
    tuning = np.exp(
        -0.5 * ((centers[2:, None] - [0, 1]) / c.functional_modules["processing_width"]) ** 2
    )
    allocation = tuning / (c.functional_modules["processing_idle"] + tuning.sum(axis=1)[:, None])
    np.testing.assert_allclose(allocation.T @ centers[:2], [0.35, 0.55], atol=1e-14)
    response = modules.effective(values, c.functional_modules)[0]
    q = traits.consumption(defense[None, :, :], response[:, None, :])
    assert q.size == 15625 and 0.75 < q.min() < 0.76 and q.max() == 1
    assert np.count_nonzero(q < 1 - 1e-12) == 11828
    assert abs(weights[q < 1 - 1e-12].sum() - 0.7512261867523193) < 1e-12
    original = values.copy()
    original[:, :2] += 0.8 - centers[:2]
    original_response = modules.effective(original, c.functional_modules)[0]
    original_ratios = original_response[:, None, :] / defense[None, :, :]
    assert np.isclose(original_ratios.min(), 1.1135204602329019)
    np.testing.assert_array_equal(
        traits.consumption(defense[None, :, :], original_response[:, None, :]), np.ones_like(q)
    )
    for axis in range(2):
        stronger = defense.copy()
        stronger[:, axis] += 1e-4
        assert (traits.consumption(stronger[None, :, :], response[:, None, :]) < q - 1e-12).any()
    # The preserved four independent 0.5-frequency copies retain coordinate SD=.025.
    assert abs(np.sqrt(4 * 0.5**2 * 0.025**2) - 0.025) < 1e-15


def test_founder_diagnostic_reports_single_larva_resource_masking_and_empty_state():
    defense = np.array([[0.5, 0.5]])
    response = np.array([[0.25, 1.0], [1.0, 0.25]])
    before = defense.copy(), response.copy()
    ample = ecology.founder_feeding_diagnostics(defense, response, [10.0], [1.0, 1.0])
    assert ample["pair_count"] == 2 and ample["q_mean"] == 0.5
    assert ample["q_sensitive_fraction"] == [0.5, 0.5]
    assert ample["single_larva_g_sensitive_fraction"] == [0.5, 0.5]
    limited = ecology.founder_feeding_diagnostics(defense, response, [0.1], [1.0, 1.0])
    assert limited["q_sensitive_fraction"] == [0.5, 0.5]
    assert limited["single_larva_g_sensitive_fraction"] == [0.0, 0.0]
    assert limited["single_larva_resource_masked_q_sensitive_fraction"] == [0.5, 0.5]
    assert limited["single_larva_resource_limited_fraction"] == 1
    np.testing.assert_array_equal(defense, before[0])
    np.testing.assert_array_equal(response, before[1])
    saturated = ecology.founder_feeding_diagnostics(defense, [[1.0, 1.0]], [10.0], [1.0])
    assert saturated["q_unsaturated_fraction"] == 0
    assert saturated["q_sensitive_fraction"] == [0.0, 0.0]
    empty = ecology.founder_feeding_diagnostics(np.empty((0, 2)), response, [], [1.0, 1.0])
    assert empty["pair_count"] == 0 and empty["status"] == "no_initial_pairs"


def test_heterogeneous_shares_change_under_defense_and_expected_recruitment_is_not_a_cap():
    host, demand = np.zeros(2, int), np.ones(2)
    response = [[1.0, 1.0], [0.4, 1.0]]
    weak = ecology.allocate_resources(host, traits.consumption([0.3, 0.5], response), [1], demand)
    strong = ecology.allocate_resources(host, traits.consumption([1, 1], response), [1], demand)
    np.testing.assert_allclose(weak["g"], [0.5, 0.5])
    np.testing.assert_allclose(strong["g"], [5 / 7, 2 / 7])
    np.testing.assert_array_equal(weak["h"], strong["h"])
    # Independent Bernoulli trials can recruit both larvae despite expected sum=1.
    assert np.prod(strong["g"]) > 0 and np.isclose(strong["g"].sum(), 1)
    homogeneous = ecology.allocate_resources(host, [1, 1], [1], [1, 2])
    heterogeneous = ecology.allocate_resources(host, [1, 0.5], [1], [1, 2])
    assert np.isclose(homogeneous["g"].sum(), 2 / 3)
    assert np.isclose(heterogeneous["g"].sum(), 3 / 4)


def test_zero_stock_has_zero_herbivory_and_resource_mortality_still_applies():
    allocation = ecology.allocate_resources(np.array([0, 0]), [1, 0.5], [0], [1, 1])
    np.testing.assert_array_equal(allocation["h"], [0])
    np.testing.assert_array_equal(allocation["g"], [0, 0])
    c = reference_config()
    _, deficit, survival, _ = plant_resource_update([0], allocation["X"], [1], c.plant_demography)
    np.testing.assert_array_equal(deficit, [1])
    np.testing.assert_allclose(
        survival,
        np.exp(
            -sum(
                c.plant_demography[key]
                for key in ("background_mortality_rate", "defoliation_mortality_rate")
            )
        ),
    )


def test_completed_feeding_window_does_not_repeat_host_death_on_larvae():
    c = reference_config()
    c.generations, c.plant_sizes, c.insect_sizes = 1, [2], [20]
    c.resource_capacity, c.choice_mode = {"P1": 1.0}, "forced"
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
    assert diagnostic["pair_count"] == 40 and all(diagnostic["q_sensitive_fraction"])


def test_relative_avoidance_and_outside_logit_are_distinct():
    c = reference_config()
    sensory = np.array([0.6, 0.4, 0, 1, 0, 1, 0, 1])
    scores = modules.scores(np.array([[0, 0], [0, 1]]), sensory, c)
    assert 0 < scores[1] < scores[0]
    probabilities = voc.choice_probabilities(scores, c.tau, c.outside_utility, "outside")
    assert probabilities[1] < probabilities[0] and probabilities[-1] < 1 / 3
    assert voc.choice_probabilities(scores, 1, 10, "outside")[-1] > 0.99
    assert voc.choice_probabilities(scores, 0, 10, "outside")[-1] > 0.99
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
