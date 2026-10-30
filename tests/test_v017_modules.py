"""Expression/tuning, inherited encounters and native saved-data adversarial checks."""

import copy
import json
from dataclasses import asdict
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from evo_sim.core.config import merge_config
from evo_sim.models.plant_insect import Config, get_runner, simulate
from evo_sim.models.plant_insect.audit import audit_run
from evo_sim.models.plant_insect.batch import save_result
from evo_sim.models.plant_insect.helpers import ecology, modules, voc
from evo_sim.models.plant_insect.helpers.genetics import GeneticState
from evo_sim.models.plant_insect.helpers.inputs import load_inputs
from tests.support import raises
from tests.test_v016_audit import reject_rehashed_tampering

CONFIG = (
    Path(__file__).resolve().parents[1]
    / "evo_sim/models/plant_insect/configs/v017_single_ancestor.yaml"
)


def module_config(**changes):
    runner = get_runner(CONFIG)
    values = copy.deepcopy(runner.config["parameters"])
    values.update(generations=1, plant_sizes=[8], insect_sizes=[16])
    values.update(changes)
    c = Config(**values)
    c.spatial.update(
        initial_plant_cells=[[0, 0]] * sum(c.plant_sizes),
        initial_insect_cells=[[0, 0]] * sum(c.insect_sizes),
    )
    return c


def test_expression_tuning_capacity_saturation_and_no_flowering():
    c = module_config()
    f = c.functional_modules
    x = np.array([[0.6, 0.4, 0.0, 1.0, 0.5, 0.2, 0.0, 1.0]])
    I, _ = modules.effective(x, f)
    assert 0 < I.sum() < x[:, :2].sum()
    zero = x.copy()
    zero[:, :2] = 0
    assert not modules.effective(zero, f)[0].any()
    mismatch = dict(f, processing_width=1e-6)
    away = x.copy()
    away[:, 2:4] = 0.5
    assert not modules.effective(away, mismatch)[0].any()
    V = np.array([[0.1, 0.8], [0.8, 0.1]])
    score = modules.scores(V, x[0], c)
    amplified = x[0].copy()
    amplified[4:6] *= 2
    np.testing.assert_allclose(modules.scores(V, amplified, c) - c.b, 2 * (score - c.b))
    tuned = x[0].copy()
    tuned[6:8] = tuned[6:8][::-1]
    assert np.sign(np.diff(score))[0] == -np.sign(np.diff(modules.scores(V, tuned, c)))[0]
    fecundity, cp, cs = modules.maternal_fecundity(x, c)
    np.testing.assert_allclose(fecundity, 1 - cp - cs)
    g = GeneticState(c)
    phenotype, raw = g.phenotype(g.founders(c.plant_sizes, "plant"), "plant")
    assert raw.shape == (8, 3)
    assert set(phenotype) == {"D", "plant_chemical_traits"}
    assert len(g.loci["plant"]) == 6 and len(g.loci["insect"]) == 16
    assert not any("flower" in locus["functional_module"] for locus in c.genetics["loci"])
    invalid = copy.deepcopy(c)
    invalid.functional_modules["flowering_kappa"] = 0.0
    raises(invalid.validate)
    invalid = copy.deepcopy(c)
    invalid.genetics["background"]["plant"].append(0.0)
    raises(invalid.validate)
    for bad in (0, -1, np.inf):
        invalid = copy.deepcopy(f)
        invalid["processing_idle"] = bad
        raises(lambda: modules.validate(invalid))
    c.functional_modules = None
    raises(c.validate)


def test_v017_sensory_processing_and_zero_expression_cost_causal_controls():
    c = module_config()
    x = np.array([[0.6, 0.4, 0.0, 1.0, 0.5, 0.2, 0.0, 1.0]] * 3)
    x[1, 6:8] = x[1, 6:8][::-1]  # Sensory tuning only.
    x[2, :2] = x[2, :2][::-1]  # Processing allocation only; same total expression.
    V = np.array([[0.1, 0.8], [0.8, 0.1]])
    choices = [voc.conditional_probabilities(V, row, c, None) for row in x]
    assert np.argmax(choices[0]) != np.argmax(choices[1])
    np.testing.assert_array_equal(choices[0], choices[2])
    I, _ = modules.effective(x, c.functional_modules)
    # One fixed host and ample resources isolate processing from sensory host choice.
    feeding = ecology.settle_feeding(
        np.array([[1.0, 0.0]]),
        I,
        np.zeros(1, int),
        np.zeros(3, int),
        np.zeros((1, 2)),
        np.zeros(3, int),
        c,
        np.array([10.0]),
        np.ones(3),
        chemical=np.array([[1.0, 0.0, 0.0]]),
        module_traits=x,
    )
    np.testing.assert_array_equal(I[0], I[1])
    assert feeding["q"][0] == feeding["q"][1] == feeding["g"][0] == feeding["g"][1]
    assert feeding["q"][0] != feeding["q"][2] and feeding["g"][0] != feeding["g"][2]
    costly, cp, cs = modules.maternal_fecundity(x, c)
    assert (cp > 0).all() and (cs > 0).all() and (costly < 1).all()
    c.insect_cost = c.functional_modules["sensory_cost"] = 0.0
    free, cp, cs = modules.maternal_fecundity(x, c)
    np.testing.assert_array_equal(free, np.ones(3))
    assert not cp.any() and not cs.any()
    np.testing.assert_array_equal(I, modules.effective(x, c.functional_modules)[0])


def test_v017_monomorphic_no_mutation_fixes_both_genomes_and_module_expression():
    c = module_config(generations=3)
    for loci in c.genetics["founder_frequencies"].values():
        for locus, frequencies in loci.items():
            allele = min(item["allele_id"] for item in frequencies)
            loci[locus] = [dict(allele_id=allele, frequency=1.0)]
    c.genetics["mutation_enabled"] = dict(plant=False, insect=False)
    result = simulate(c)
    assert result["metadata"]["status"]["status"] == "completed"
    assert not result["mutations"] and len(result["snapshots"]) == 4
    for kind, phenotype in (
        ("plant", "plant_chemical_traits"),
        ("insect", "insect_module_traits"),
    ):
        initial = result["snapshots"][0]
        founder = initial[kind + "_genotypes"][0]
        assert np.all(founder == founder[0])
        for snapshot in result["snapshots"]:
            assert np.all(snapshot[kind + "_genotypes"] == founder)
            assert np.all(snapshot[phenotype] == initial[phenotype][0])


def test_v017_current_contract_rejects_removed_settlement_and_preserves_formal_design():
    for removed_value in (True, False):
        old = module_config()
        old.functional_modules["adult_settlement_enabled"] = removed_value
        old.validate(historical_audit=True)
        raises(old.validate)
        raises(lambda: simulate(old))
    loaded = get_runner(CONFIG.with_name("v017_200gen_400insects.yaml")).config
    assert loaded["simulation_metadata"]["num_steps_per_episode"] == 200
    assert loaded["batch"]["replicates"] == 3
    base = loaded["parameters"]
    assert base["plant_sizes"] == [200, 200] and base["insect_sizes"] == [400]
    assert base["genetics"]["schema_version"] == "slim-qtl-4"
    assert "adult_settlement_enabled" not in base["functional_modules"]
    assert "flowering_kappa" not in base["functional_modules"]
    cases = loaded["batch"]["cases"]
    assert [case["name"] for case in cases] == ["qtl_mutating", "mutation_off"]
    assert base["seed"] == 20261021
    for case, mutation in zip(cases, (True, False), strict=True):
        values = merge_config(base, case["overrides"])
        assert values["genetics"]["mutation_enabled"] == dict(plant=mutation, insect=mutation)


def test_v017_native_lifecycle_independent_audit_and_tamper_rejection():
    c = module_config(generations=2)
    result = simulate(c)
    assert result["metadata"]["status"]["status"] == "completed"
    assert result["metadata"]["version"] == "0.1.7"
    assert result["metadata"]["functional_schema"] == "expression-tuning-2"
    first = result["snapshots"][0]
    assert first["insect_module_traits"].shape == (16, 8)
    assert result["metadata"]["genetics_version"] == "slim-qtl-4"
    assert not any("flower" in key or "theta" in key for key in result["metadata"])
    assert not any("flower" in key or "theta" in key for key in first)
    assert "flowering_kappa" not in c.functional_modules
    assert not any("flower" in row.get("trait", "") for row in result["genetic_summary"])
    assert not any("wrapped_fraction" in row for row in result["genetic_summary"])
    assert first["reference_p"].shape == (16, 2)
    assert first["reference_plant_species"].tolist() == ["chem_A", "chem_B"]
    assert "adult_settlement" not in result
    assert "settlement_semantics" not in result["metadata"]
    assert result["gamete_transmissions"] and result["mutations"]
    for row in result["parents"]:
        ids = json.loads(row["candidate_mate_ids"])
        probability = json.loads(row["candidate_mate_probabilities"])
        assert len(ids) == len(probability)
        if probability:
            assert np.isclose(sum(probability), 1)
        if row["kind"] == "plant" and ids:
            np.testing.assert_allclose(probability, 1 / len(ids))
        state = result["snapshots"][row["generation"]]
        index = np.flatnonzero(state[row["kind"] + "_ids"] == row["individual_id"])[0]
        np.testing.assert_array_equal(
            [row["mating_row"], row["mating_col"]], state[row["kind"] + "_cells"][index]
        )
    with TemporaryDirectory() as temporary:
        path = Path(temporary) / "run"
        save_result(result, path)
        assert not (path / "adult_settlement.csv").exists()
        with np.load(path / "phenotypes.npz") as exported:
            assert not any("flower" in key or "theta" in key for key in exported.files)
        assert "wrapped_fraction" not in (path / "genetic_summary.csv").read_text()
        checked = audit_run(path)
        assert checked["genetic_audit"] and checked["contribution_audit"]
        for column in ("mating_row", "conditional_mu"):
            reject_rehashed_tampering(
                path, "parents", lambda frame, col=column: frame.__setitem__(col, frame[col] + 1)
            )
        reject_rehashed_tampering(
            path,
            "parents",
            lambda frame: frame.__setitem__("candidate_mate_probabilities", "[]"),
        )
    unseen_c = copy.deepcopy(c)
    unseen_c.record_observations = False
    unseen = simulate(unseen_c)
    for observed, unobserved in zip(result["snapshots"], unseen["snapshots"], strict=True):
        for field in ("plant_genotypes", "insect_genotypes", "plant_cells", "insect_cells"):
            np.testing.assert_array_equal(observed[field], unobserved[field])
    assert result["offspring"] == unseen["offspring"]


def test_v017_adults_mate_at_existing_cells_and_failure_is_atomic():
    c = module_config(choice_mode="forced")
    c.spatial.update(search_radius=10)
    c.sexual_reproduction["insect_mate_radius"] = 0
    stream = np.random.default_rng(np.random.SeedSequence([c.seed, 130, 6, 1, 0]))
    male = (
        stream.random(sum(c.insect_sizes))
        >= c.sexual_reproduction["insect_founder_female_probability"]
    )
    c.spatial["initial_insect_cells"] = [[2, 2] if sex else [1, 1] for sex in male]
    assert male.any() and (~male).any()
    separated = simulate(c)
    assert not any(r["kind"] == "insect" for r in separated["offspring"])
    assert all(
        not json.loads(r["candidate_mate_ids"])
        and [r["mating_row"], r["mating_col"]]
        == c.spatial["initial_insect_cells"][r["individual_id"]]
        for r in separated["parents"]
        if r["kind"] == "insect"
    )
    together_c = copy.deepcopy(c)
    together_c.spatial["initial_insect_cells"] = [[1, 1]] * sum(c.insect_sizes)
    together = simulate(together_c)
    assert any(r["kind"] == "insect" for r in together["offspring"])
    for result in (separated, together):
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / "run"
            save_result(result, path)
            assert audit_run(path)["genetic_audit"]
    failed_c = copy.deepcopy(together_c)
    failed_c.insect_offspring_lambda = {"H1": 1000.0}
    failed_c.max_individuals_safety = sum(c.plant_sizes) + sum(c.insect_sizes)
    failed = simulate(failed_c)
    assert failed["metadata"]["status"]["status"] == "computational_limit"
    assert "adult_settlement" not in failed and not failed["offspring"]
    np.testing.assert_array_equal(
        failed["snapshots"][0]["insect_cells"], together_c.spatial["initial_insect_cells"]
    )


def test_superseded_flowering_protocol_is_readable_but_not_executable():
    values = asdict(module_config())
    values["genetics"]["schema_version"] = "slim-qtl-3"
    values["genetics"]["background"]["plant"].append(0.0)
    values["functional_modules"]["flowering_kappa"] = 2.0
    values["functional_modules"]["adult_settlement_enabled"] = True
    legacy, _ = load_inputs(
        dict(
            loci="loci.csv",
            alleles="alleles.csv",
            founder_frequencies="single_founder_frequencies.csv",
        ),
        CONFIG.parents[1] / "inputs/v017",
        values,
    )
    old = Config(**legacy)
    old.validate(historical_audit=True)
    raises(old.validate)
    assert len([l for l in old.genetics["loci"] if l["organism"] == "plant"]) == 8
    raises(lambda: GeneticState(old))
    raises(lambda: simulate(old))
    legacy["genetics"]["schema_version"] = "slim-qtl-4"
    raises(Config(**legacy).validate)
