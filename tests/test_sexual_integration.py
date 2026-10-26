import copy
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect import get_runner, simulate
from evo_sim.models.plant_insect.audit import audit_run
from evo_sim.models.plant_insect.batch import save_result, write_json
from evo_sim.models.plant_insect.helpers.genetics import GeneticState
from evo_sim.models.plant_insect.helpers.slim_backend import SlimError
from tests.support import config, raises


def test_contract_rejections():
    c = config()
    c.validate()
    for key, value in (
        ("schema_version", "qtl-1"),
        ("inheritance", "clonal_diploid_no_recombination"),
        ("mode", "direct_trait_legacy"),
        ("extra", 1),
    ):
        bad = copy.deepcopy(c)
        bad.genetics[key] = value
        raises(bad.validate)
    for version in ("1.3", "1.4", "1.5", "1.6"):
        raises(config(model_version=version).validate)
    raises(config(genetics=None).validate)
    raises(lambda: config(evolve_plants=True), TypeError)
    for key, value in (("position_cM", float("nan")), ("chromosome_id", ""), ("targets", [9])):
        bad = copy.deepcopy(c)
        bad.genetics["loci"][0][key] = value
        raises(bad.validate)
    for key in ("pollen_radius", "insect_mate_radius", "insect_founder_female_probability"):
        bad = copy.deepcopy(c)
        bad.sexual_reproduction[key] = -1
        raises(bad.validate)


def test_terminal_and_atomic_failure():
    runner = get_runner()
    runner.config["simulation_metadata"]["num_steps_per_episode"] = 0
    runner.init()
    runtime = runner.state["environment"]["runtime"]
    before = copy.deepcopy(runtime.genetic.rng)
    runner.step()
    assert (
        not runner.result["offspring"]
        and not runtime.genetic.transmissions
        and not runtime.genetic.events
    )
    for key, streams in before.items():
        for a, b in zip(streams, runtime.genetic.rng[key]):
            assert a.bit_generator.state == b.bit_generator.state
    runner = get_runner()
    runner.init()
    runtime = runner.state["environment"]["runtime"]
    alleles, ids = copy.deepcopy(runtime.genetic.alleles), runtime.next_ids.copy()
    genes = copy.deepcopy(runtime.genetic.rng)
    agents = copy.deepcopy(runner.state["agents"])
    original = GeneticState.fertilize

    def fail_second(self, *args, **kwargs):
        if args[5] == "insect":
            raise MemoryError("injected failure after plant gametes")
        return original(self, *args, **kwargs)

    with patch.object(GeneticState, "fertilize", fail_second):
        runner.step(1)
    assert runner.result["metadata"]["status"]["status"] == "computational_limit"
    assert (
        not runner.result["offspring"]
        and not runtime.genetic.events
        and not runtime.genetic.transmissions
    )
    assert runtime.genetic.alleles == alleles and runtime.next_ids == ids
    for kind in agents:
        for key in agents[kind]:
            np.testing.assert_array_equal(agents[kind][key], runner.state["agents"][kind][key])
    for key, streams in genes.items():
        for a, b in zip(streams, runtime.genetic.rng[key]):
            assert a.bit_generator.state == b.bit_generator.state
    with TemporaryDirectory() as tmp:
        save_result(runner.result, Path(tmp) / "injected")
        assert audit_run(Path(tmp) / "injected")["status"] == "computational_limit"
    for field, value in (("max_alleles_safety", 48), ("max_genotype_values_safety", 736)):
        c = config()
        c.genetics[field] = value
        r = simulate(c)
        assert r["metadata"]["status"]["status"] == "computational_limit" and not r["offspring"]
    assert (
        simulate(config(max_individuals_safety=1))["metadata"]["status"]["status"]
        == "computational_limit"
    )


def test_linked_run_audit_and_tampering():
    c = config(generations=3)
    for l, locus in enumerate(c.genetics["loci"]):
        locus.update(chromosome_id=locus["organism"] + "_linked", position_cM=float(10 * l))
    result = simulate(c)
    assert result["metadata"]["status"]["status"] == "completed"
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "run"
        save_result(result, path)
        assert audit_run(path)["genetic_audit"]
        for table, column in (
            ("offspring", "paternal_id"),
            ("gamete_transmissions", "inherited_allele_id"),
            ("parents", "expected_maternal_offspring"),
        ):
            file = path / (table + ".csv")
            original = file.read_bytes()
            frame = pd.read_csv(file)
            frame.loc[0, column] += 100000
            frame.to_csv(file, index=False)
            meta = json.loads((path / "metadata.json").read_text())
            meta["output_sha256"][file.name] = hashlib.sha256(file.read_bytes()).hexdigest()
            write_json(path / "metadata.json", meta)
            raises(lambda: audit_run(path), (AssertionError, KeyError))
            file.write_bytes(original)
            meta["output_sha256"][file.name] = hashlib.sha256(original).hexdigest()
            write_json(path / "metadata.json", meta)
        assert audit_run(path)["genetic_audit"]


def test_insect_biparental_births_and_audit():
    c = config(
        generations=2,
        insect_sizes=[40],
        choice_mode="forced",
        sP0=1.0,
        resource_capacity={"P1": 100.0, "P2": 100.0},
    )
    c.spatial["initial_plant_cells"] = [[0, 0]] * 40
    c.spatial["initial_insect_cells"] = [[0, 0]] * 40
    result = simulate(c)
    assert any(row["kind"] == "insect" for row in result["offspring"])
    with TemporaryDirectory() as tmp:
        save_result(result, Path(tmp) / "run")
        audit_run(Path(tmp) / "run")


def test_slim_failure_does_not_commit_or_fall_back():
    runner = get_runner()
    runner.config["parameters"]["genetics"]["slim"]["executable"] = "missing_evosim_slim_binary"
    runner.init()
    runtime = runner.state["environment"]["runtime"]
    before = copy.deepcopy(runtime.genetic.alleles)
    runner.step(1)
    assert runner.result["metadata"]["status"]["status"] == "backend_error"
    assert not runtime.genetic.engine_calls and not runner.result["offspring"]
    assert runtime.genetic.alleles == before and not runtime.genetic.events
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "failed"
        save_result(runner.result, path)
        assert audit_run(path)["status"] == "backend_error"

    # Reject malformed returned material after a real SLiM process, before publication.
    from evo_sim.models.plant_insect.helpers import slim_backend

    original = slim_backend.subprocess.run

    def corrupt(command, **kwargs):
        result = original(command, **kwargs)
        output = next((a for a in command if a.startswith("OUTPUT=")), None)
        if output:
            Path(json.loads(output[len("OUTPUT=") :])).write_text("invalid output\n")
        return result

    c = config()
    g = GeneticState(c)
    parents = g.founders([2, 0], "plant")
    before = g.alleles.copy()
    with patch.object(slim_backend.subprocess, "run", corrupt):
        raises(
            lambda: g.fertilize(parents, np.zeros(2, int), np.arange(2), [0], [1], "plant", 2, 1),
            SlimError,
        )
    assert g.alleles == before and not g.events and not g.engine_calls
