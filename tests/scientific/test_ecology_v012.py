"""Analytic v0.1.12 scientific acceptance and independent saved-record checks."""

import copy
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from evo_sim.models.plant_insect import simulate
from evo_sim.models.plant_insect.audit.run import audit_run
from evo_sim.models.plant_insect.batch import save_result
from evo_sim.models.plant_insect.biology.genetics.qtl import GeneticState
from evo_sim.models.plant_insect.biology.insects.feeding import chemical_survival
from evo_sim.models.plant_insect.biology.insects.host_selection import (
    choose_hosts,
    contact_acceptance,
)
from evo_sim.models.plant_insect.biology.plants.demography import (
    leaf_recovery_rate,
    plant_resource_update,
)
from tests.scientific.test_independent_audit import reject_rehashed_tampering
from tests.support import config, raises


def test_psm_voc_independent_genetic_changes_and_optional_pleiotropy():
    c = config()
    c.validate()
    genetic = GeneticState(c)
    genome = genetic.founders([1, 0], "plant")
    baseline = genetic.phenotype(genome, "plant")[0]["plant_chemical_traits"]
    for axis in range(5):
        altered = copy.deepcopy(genetic)
        locus = next(
            index for index, entry in enumerate(altered.loci["plant"]) if entry["targets"] == [axis]
        )
        allele = int(genome[0, 0, locus])
        effect = list(altered.alleles[allele].effect)
        effect[axis] += 0.02
        altered.alleles[allele] = replace(altered.alleles[allele], effect=tuple(effect))
        changed = altered.phenotype(genome, "plant")[0]["plant_chemical_traits"]
        difference = changed - baseline
        assert difference[0, axis] > 0
        np.testing.assert_array_equal(np.delete(difference, axis, axis=1), 0)
    locus = c.genetics["loci"][0]
    locus["targets"].append(3)
    for allele in c.genetics["alleles"]:
        if allele["locus_id"] == locus["locus_id"]:
            allele["effect"][3] = -allele["effect"][0]
    c.validate()
    genotype = GeneticState(c)
    assert genotype.phenotype(genotype.founders([1, 0], "plant"), "plant")[1].shape == (1, 5)


def test_contact_weights_and_bounded_search_after_rejection():
    c = config(choice_mode="forced", plant_sizes=[2, 0])
    traits = np.full((1, 10), 0.5)
    psm = np.array([[1.0, 0.0], [0.25, 0.75]])
    traits[0, 8:10] = [1, 0]
    probability = contact_acceptance(psm, traits[0], c.contact_response)
    np.testing.assert_allclose(probability, 1 / (1 + np.exp(-np.array([1.0, -0.5]))))
    assert np.all((probability >= 0) & (probability <= 1))
    assert probability[0] > 0.5 > probability[1] and probability.sum() > 1
    for row, expected in zip(psm, probability, strict=True):
        np.testing.assert_array_equal(
            contact_acceptance(row[None, :], traits[0], c.contact_response), [expected]
        )
    traits[0, 8:10] = [0, 0]
    assert contact_acceptance(psm, traits[0], c.contact_response).sum() < 1

    class SearchStream:
        def choice(self, population, size=None, replace=True, p=None):
            return np.asarray(population)[:size] if size is not None else 0

        def random(self):
            return 0.5

    traits[0, 8:10] = [0, 1]
    result = choose_hosts(
        traits,
        np.zeros(2, int),
        np.zeros(1, int),
        np.zeros((2, 2)),
        c,
        None,
        np.random.default_rng(1),
        [SearchStream()],
        candidate_pools=[np.array([0, 1])],
        PSM=psm,
    )
    assert result["host"].tolist() == [1]
    assert result["contact_rows"] == [[0, 1]] and result["contact_accepted"] == [[False, True]]
    c.contact_response["intercept"] = -1000
    result = choose_hosts(
        traits,
        np.zeros(2, int),
        np.zeros(1, int),
        np.zeros((2, 2)),
        c,
        None,
        np.random.default_rng(1),
        [SearchStream()],
        candidate_pools=[np.array([0, 1])],
        PSM=psm,
    )
    assert result["host"].tolist() == [-1] and len(result["contact_rows"][0]) == 2


def test_oviposition_terminal_reasons_preserve_contact_history():
    class SearchStream:
        def __init__(self, selections):
            self.selections = iter(selections)

        def choice(self, population, size=None, replace=True, p=None):
            return np.asarray(population)[:size] if size is not None else next(self.selections)

        def random(self):
            return 0.5

    traits = np.full((1, 10), 0.5)
    traits[0, 8:10] = [0, 1]
    psm = np.array([[1.0, 0.0], [0.25, 0.75]])
    for selections, candidates, intercept, reason, visited in (
        ([1], [0, 1], 0, "accepted", [1]),
        ([0, 0], [0, 1], -1000, "rejected", [0, 1]),
        ([2], [0, 1], 0, "outside_option", []),
        ([0, 1], [0, 1], 0, "outside_option", [0]),
        ([], [], 0, "no_candidate", []),
    ):
        c = config(plant_sizes=[2, 0], choice_mode="outside")
        c.contact_response["intercept"] = intercept
        result = choose_hosts(
            traits,
            np.zeros(2, int),
            np.zeros(1, int),
            np.zeros((2, 2)),
            c,
            None,
            np.random.default_rng(1),
            [SearchStream(selections)],
            candidate_pools=[np.array(candidates, int)],
            PSM=psm,
        )
        assert result["reason"].tolist() == [reason]
        assert result["contact_rows"] == [visited]
        assert result["host_found"].tolist() == [reason == "accepted"]
        assert result["host"].tolist() == ([1] if reason == "accepted" else [-1])
        assert result["contact_accepted"] == [[reason == "accepted"] * len(visited)]


def test_processing_reduces_residual_but_not_intake_or_demand():
    spec = config().chemical_processing
    first = chemical_survival([[0.8, 0.4]], [[0.1, 0.1]], [0.5], [1.0], spec)
    second = chemical_survival([[0.8, 0.4]], [[1.0, 1.0]], [0.5], [1.0], spec)
    np.testing.assert_allclose(first["dose"], [[0.4, 0.2]])
    np.testing.assert_allclose(first["processed"], [[0.1, 0.1]])
    np.testing.assert_allclose(first["residual"], [[0.3, 0.1]])
    assert not second["residual"].any()
    assert second["survival_probability"][0] > first["survival_probability"][0]
    from evo_sim.models.plant_insect.biology.insects.feeding import settle_feeding

    c = config()
    values = np.full((2, 10), 0.5)
    results = []
    for processing in (np.zeros((2, 2)), np.ones((2, 2))):
        results.append(
            settle_feeding(
                np.array([[0.8, 0.4]]),
                processing,
                np.zeros(1, int),
                np.zeros(2, int),
                np.zeros((1, 2)),
                np.zeros(2, int),
                c,
                np.array([1.0]),
                np.ones(2),
                chemical=np.array([[0.8, 0.4, 0.5, 0, 0]]),
                module_traits=values,
            )
        )
    for key in ("d", "total_demand", "x", "g", "X"):
        np.testing.assert_array_equal(results[0][key], results[1][key])
    np.testing.assert_array_equal(results[0]["d"], [1, 1])
    assert (results[0]["residual"] > results[1]["residual"]).all()


def test_tolerance_recovers_faster_at_equal_damage_without_changing_mortality():
    spec = config().plant_demography
    post, damage, survival, recovery = plant_resource_update([8, 8], [4, 4], [10, 10], spec, [0, 1])
    np.testing.assert_array_equal(post, [4, 4])
    np.testing.assert_array_equal(damage, [0.5, 0.5])
    rate = leaf_recovery_rate([0, 1], damage, spec)
    np.testing.assert_allclose(
        rate, spec["resource_recovery_rate"] + spec["tolerance_recovery_rate"] * np.array([0, 0.5])
    )
    assert survival[0] == survival[1] and recovery[1] > recovery[0]
    assert np.all(post <= recovery) and np.all(recovery <= 10)


def test_tolerance_does_not_change_undamaged_or_later_baseline_recovery():
    spec = config().plant_demography
    for stock in (0, 2.5, 10):
        post, damage, survival, recovered = plant_resource_update(
            [stock, stock], [0, 0], [10, 10], spec, [0, 1]
        )
        np.testing.assert_array_equal(damage, [0, 0])
        np.testing.assert_array_equal(
            leaf_recovery_rate([0, 1], damage, spec), [spec["resource_recovery_rate"]] * 2
        )
        np.testing.assert_allclose(
            recovered, 10 - (10 - post) * np.exp(-spec["resource_recovery_rate"])
        )
        assert recovered[0] == recovered[1] and (recovered <= 10).all()
        np.testing.assert_allclose(survival, np.exp(-spec["background_mortality_rate"]))
    # After a damaged round, a later unbitten round has only the baseline rate.
    *_, first = plant_resource_update([8, 8], [4, 4], [10, 10], spec, [0, 1])
    _, damage, _, second = plant_resource_update(first, [0, 0], [10, 10], spec, [0, 1])
    np.testing.assert_allclose(
        (10 - second) / (10 - first), np.exp(-spec["resource_recovery_rate"])
    )
    np.testing.assert_array_equal(damage, [0, 0])


def test_recovery_rejects_invalid_damage_tolerance_and_resource_budgets():
    spec = config().plant_demography
    for z, h in (([0], [0, 1]), ([np.nan], [0]), ([0], [np.inf]), ([-1], [0]), ([0], [1.1])):
        raises(lambda z=z, h=h: leaf_recovery_rate(z, h, spec))
    for stock, consumed, capacity in (([11], [0], [10]), ([1], [2], [10]), ([-1], [0], [10])):
        raises(lambda: plant_resource_update(stock, consumed, capacity, spec, [0.5]))
    assert leaf_recovery_rate([], [], spec).shape == (0,)


def test_newborn_low_stock_without_feeding_is_not_defoliation():
    spec = config().plant_demography
    _, damage, survival, _ = plant_resource_update([2.5, 0], [0, 0], [10, 10], spec, [0.5, 0.5])
    np.testing.assert_array_equal(damage, [0, 0])
    np.testing.assert_allclose(survival, np.exp(-spec["background_mortality_rate"]))


def test_new_axes_native_transmission_recombination_and_mutation():
    c = config()
    for kind, axes in (("plant", {2, 3, 4}), ("insect", {8, 9})):
        c.genetics["mutation_enabled"][kind] = False
        genetic = GeneticState(c)
        genome = genetic.founders(
            [2] + [0] * (len(c.plant_sizes) - 1) if kind == "plant" else [2], kind
        )
        for index, locus in enumerate(genetic.loci[kind]):
            ids = [
                entry["allele_id"]
                for entry in c.genetics["founder_frequencies"]["P1" if kind == "plant" else "H1"][
                    locus["locus_id"]
                ]
            ]
            genome[:, 0, index] = ids[0]
            genome[:, 1, index] = ids[1]
        genetic.fertilize(
            genome,
            np.zeros(2, int),
            np.arange(2),
            np.zeros(128, int),
            np.ones(128, int),
            kind,
            2,
            1,
        )
        for axis in axes:
            own = {locus["locus_id"] for locus in genetic.loci[kind] if axis in locus["targets"]}
            assert any(
                row["switch_before_locus"]
                for row in genetic.transmissions
                if row["locus_id"] in own
            )
        c.genetics["mutation_enabled"][kind] = True
        for locus in c.genetics["loci"]:
            if locus["organism"] == kind:
                locus["mutation_rate"] = 1
        genetic = GeneticState(c)
        genetic.fertilize(
            genome, np.zeros(2, int), np.arange(2), np.zeros(4, int), np.ones(4, int), kind, 2, 1
        )
        for axis in axes:
            assert any(row["delta"][axis] != 0 for row in genetic.events)


def test_nutrition_toxicity_antinutrition_are_independent_controls():
    spec = copy.deepcopy(config().chemical_processing)
    spec.update(toxicity_weights=[0, 0], antinutrition_weights=[0, 0])
    baseline = chemical_survival([[1, 1]], [[0, 0]], [0.5], [1.0], spec)
    np.testing.assert_allclose(baseline["survival_probability"], [0.5])
    for field, factor in (
        ("toxicity_weights", np.exp(-0.5)),
        ("antinutrition_weights", np.exp(-0.5)),
    ):
        altered = dict(spec, **{field: [1, 0]})
        result = chemical_survival([[1, 1]], [[0, 0]], [0.5], [1.0], altered)
        np.testing.assert_allclose(result["survival_probability"], 0.5 * factor)
        np.testing.assert_array_equal(result["nutrition_survival"], baseline["nutrition_survival"])
    spec["nutrition_exponent"] = 0
    np.testing.assert_array_equal(
        chemical_survival([[1, 1]], [[0, 0]], [0.5], [1.0], spec)["survival_probability"], [1]
    )
    assert chemical_survival([[1, 1]], [[0, 0]], [0], [1.0], spec)["survival_probability"][0] == 0
    raises(lambda: chemical_survival([[1, 1]], [[0, 0]], [2], [1.0], spec))


def test_new_records_replay_independently_and_reject_rehashed_tampering():
    c = config(plant_sizes=[4, 4], insect_sizes=[12], generations=2, choice_mode="outside")
    c.spatial.update(initial_plant_cells=[[0, 0]] * 8, initial_insect_cells=[[0, 0]] * 12)
    result = simulate(c)
    assert result["metadata"]["status"]["status"] == "completed" and result["events"]
    assert any(row["reason"] == "outside_option" for row in result["oviposition"])
    for row in result["plants"]:
        np.testing.assert_allclose(row["defoliation"], row["X"] / row["R"] if row["R"] else 0)
        np.testing.assert_allclose(
            row["recovery_rate"],
            c.plant_demography["resource_recovery_rate"]
            + c.plant_demography["tolerance_recovery_rate"]
            * row["tolerance_Z"]
            * row["defoliation"],
        )
        np.testing.assert_allclose(
            row["W_plant"],
            row["survival_probability"]
            * max(0, 1 - row["defense_cost"] - row["voc_cost"])
            * (
                row["resource_capacity"]
                - (row["resource_capacity"] - row["resource_after_feeding"])
                * np.exp(-row["recovery_rate"])
            )
            / row["resource_capacity"],
        )
    repeated = simulate(c)
    for table in ("events", "oviposition", "offspring", "gamete_transmissions"):
        assert result[table] == repeated[table]
    with TemporaryDirectory() as temporary:
        path = Path(temporary) / "run"
        save_result(result, path)
        assert audit_run(path)["genetic_audit"]
        for table, column in (
            ("larvae", "residual_1"),
            ("larvae", "dose_2"),
            ("larvae", "W_insect"),
            ("plants", "defoliation"),
            ("plants", "recovery_rate"),
        ):
            reject_rehashed_tampering(
                path,
                table,
                lambda frame, field=column: frame.__setitem__(field, frame[field] + 0.1),
            )
        reject_rehashed_tampering(
            path, "oviposition", lambda frame: frame.__setitem__("contact_probabilities", "[]")
        )
        reject_rehashed_tampering(
            path,
            "oviposition",
            lambda frame: frame.__setitem__(
                "reason", frame.reason.replace("outside_option", "rejected")
            ),
        )
