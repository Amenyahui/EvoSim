"""Native-backend, stage-key and cumulative-storage regression checks."""

import copy
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd
import yaml

from evo_sim.models.plant_insect import Config, get_runner, simulate
from evo_sim.models.plant_insect.analysis.common import Batch
from evo_sim.models.plant_insect.analysis.summaries.lifecycle import (
    prepare_diagnostics,
    summarize_sensitivity,
)
from evo_sim.models.plant_insect.audit.run import audit, audit_run
from evo_sim.models.plant_insect.batch import COLUMNS, run, save_result
from evo_sim.models.plant_insect.config.schema import config_payload

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "evo_sim/models/plant_insect/configs/current_single_ancestor.yaml"


def small_config(**changes):
    values = copy.deepcopy(get_runner(CONFIG).config["parameters"])
    values.update(generations=1, plant_sizes=[8], insect_sizes=[16], choice_mode="forced")
    values.update(changes)
    c = Config(**values)
    c.plant_demography.update(microsites_per_cell=max(1, 4 * sum(c.plant_sizes)))
    c.spatial.update(
        initial_plant_cells=[[0, 0]] * sum(c.plant_sizes),
        initial_insect_cells=[[0, 0]] * sum(c.insect_sizes),
    )
    c.genetics["mutation_enabled"] = {"plant": False, "insect": False}
    return c


def stored_values(result):
    snapshots = sum(
        state[kind + "_genotypes"].size
        for state in result["snapshots"]
        for kind in ("plant", "insect")
    )
    return snapshots + len(result["gamete_transmissions"])


def test_native_stage_keys_and_empty_table_schema():
    c = small_config(resource_capacity={"P1": 100.0})
    c.genetics["mutation_enabled"] = {"plant": True, "insect": True}
    for locus in c.genetics["loci"]:
        locus["mutation_rate"] = 1.0
    result = simulate(c)
    assert result["metadata"]["status"]["status"] == "completed"
    assert {call["kind"] for call in result["metadata"]["slim_calls"]} == {"plant", "insect"}
    assert result["mutations"] and result["events"]
    with TemporaryDirectory() as tmp:
        full, empty = Path(tmp) / "full", Path(tmp) / "empty"
        save_result(result, full)
        zero = copy.deepcopy(c)
        zero.generations = 0
        save_result(simulate(zero), empty)
        checked = audit_run(full)
        assert checked["execution_success"] and checked["transitions"] == 1
        for name, columns in COLUMNS.items():
            assert pd.read_csv(full / f"{name}.csv").columns.tolist() == columns.split()
            assert pd.read_csv(empty / f"{name}.csv").columns.tolist() == columns.split()
        children = pd.read_csv(full / "offspring.csv")
        larvae = pd.read_csv(full / "larvae.csv").assign(kind="insect")
        joined = larvae.merge(
            children,
            left_on=["reproduction_round", "adult_cohort", "kind", "insect_id"],
            right_on=["reproduction_round", "adult_cohort", "kind", "offspring_id"],
            validate="one_to_one",
        )
        assert len(joined) == len(larvae)
        tx_file = full / "gamete_transmissions.csv"
        tx = pd.read_csv(tx_file)
        tx.loc[0, "adult_cohort"] += 1
        tx.to_csv(tx_file, index=False)
        meta_file = full / "metadata.json"
        meta = json.loads(meta_file.read_text(encoding="utf-8"))
        meta["output_sha256"][tx_file.name] = hashlib.sha256(tx_file.read_bytes()).hexdigest()
        meta_file.write_text(json.dumps(meta), encoding="utf-8")
        try:
            audit_run(full)
        except AssertionError:
            pass
        else:
            raise AssertionError("rehashed cohort mismatch passed independent audit")


def test_storage_budget_includes_initial_terminal_and_zero_birth_snapshots():
    c = small_config()
    baseline = simulate(c)
    assert baseline["metadata"]["status"]["status"] == "completed"
    total = stored_values(baseline)
    initial = sum(
        baseline["snapshots"][0][kind + "_genotypes"].size for kind in ("plant", "insect")
    )
    for budget, expected in ((total, "completed"), (total - 1, "computational_limit")):
        limited = copy.deepcopy(c)
        limited.genetics["max_genotype_values_safety"] = budget
        result = simulate(limited)
        assert result["metadata"]["status"]["status"] == expected
        assert result["metadata"]["stored_genotype_values"] == stored_values(result) <= budget
        if expected == "computational_limit":
            assert len(result["snapshots"]) == 1 and not result["offspring"]
            assert not result["events"] and not result["gamete_transmissions"]
        with TemporaryDirectory() as tmp:
            save_result(result, Path(tmp) / "run")
            audit_run(Path(tmp) / "run")
    for budget, expected in ((initial, "completed"), (initial - 1, "computational_limit")):
        zero = copy.deepcopy(c)
        zero.generations = 0
        zero.genetics["max_genotype_values_safety"] = budget
        result = simulate(zero)
        assert result["metadata"]["status"]["status"] == expected
        assert result["metadata"]["stored_genotype_values"] <= budget
    # Full plant occupancy and a lone insect yield no newborns.
    zero = small_config(plant_sizes=[1], insect_sizes=[1])
    zero.plant_demography.update(
        background_mortality_rate=0, defoliation_mortality_rate=0, microsites_per_cell=1
    )
    zero.spatial.update(
        grid_shape=[1, 1], search_radius=0, plant_dispersal_radius=0, insect_dispersal_radius=0
    )
    zero.sexual_reproduction["insect_mate_radius"] = 0
    initial = 2 * len(zero.genetics["loci"])
    retained_plant_values = 2 * sum(l["organism"] == "plant" for l in zero.genetics["loci"])
    zero.genetics["max_genotype_values_safety"] = initial + retained_plant_values
    result = simulate(zero)
    assert result["metadata"]["status"]["status"] == "completed"
    assert not result["offspring"] and not result["metadata"]["slim_calls"]
    assert len(result["snapshots"]) == 2
    assert (
        stored_values(result)
        == result["metadata"]["stored_genotype_values"]
        == initial + retained_plant_values
    )


def test_missing_backend_cli_and_strict_execution_acceptance():
    config = get_runner(CONFIG).config
    config["parameters"] = config_payload(small_config())
    config["parameters"].pop("generations")
    config["simulation_metadata"]["num_steps_per_episode"] = 1
    config["batch"] = dict(replicates=1, cases=[dict(name="missing", overrides={})])
    env = dict(os.environ, EVOSIM_SLIM="missing_evosim_slim_binary")
    with TemporaryDirectory() as tmp:
        path = Path(tmp)
        config_file = path / "missing.yaml"
        config_file.write_text(yaml.safe_dump(config), encoding="utf-8")
        output = path / "batch"
        command = [sys.executable, "-m", "evo_sim.models.plant_insect.batch"]
        failed = subprocess.run(
            command + ["--config", str(config_file), "--out", str(output)],
            env=env,
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        assert failed.returncode != 0
        status = json.loads((output / "data/status.json").read_text(encoding="utf-8"))
        assert not status["execution_success"] and status["status"] == "failed"
        report = audit(output, min_transitions=1)
        assert report["integrity_status"] == "passed" and report["transitions"] == 0
        assert not report["execution_success"] and not report["acceptance_passed"]
        strict = subprocess.run(
            [
                sys.executable,
                "-m",
                "evo_sim.models.plant_insect.audit",
                "--input",
                str(output),
                "--require-success",
                "--min-transitions",
                "1",
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        assert strict.returncode != 0 and "acceptance failed" in strict.stderr
        preflight = subprocess.run(
            [sys.executable, "-m", "tests"],
            env=env,
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        assert preflight.returncode != 0 and "SLiM preflight failed" in preflight.stderr
        assert "PASS " not in preflight.stdout


def test_single_pool_has_no_ghost_P2_resource_fraction():
    spec = dict(
        base=config_payload(small_config()),
        replicates=1,
        cases=[dict(name="single", overrides={})],
    )
    with TemporaryDirectory() as tmp:
        path = run(spec, Path(tmp) / "batch")
        checked = audit(path, min_transitions=1)
        assert checked["acceptance_passed"]
        batch = Batch.open(path)
        prepare_diagnostics(batch)
        fractions = pd.read_csv(batch.summary_dir / "host_resource_by_generation.csv")
        assert fractions.total_actual_resource.gt(0).all()
        assert fractions.P2_resource_fraction.isna().all()


def test_cue_shuffle_does_not_consume_choice_order_stream():
    states = []
    for shuffle in (False, True):
        runner = get_runner(CONFIG)
        c = small_config(shuffle_voc=shuffle)
        runner.config["parameters"] = config_payload(c)
        runner.config["parameters"].pop("generations")
        runner.config["simulation_metadata"]["num_steps_per_episode"] = 1
        runner.init()
        plants = runner.state["agents"]["plants"]
        # Equal cues make permutation causally inert; only extra stream consumption could differ.
        plants["chemical_traits"][:] = plants["chemical_traits"][0]
        plants["defense"][:] = plants["defense"][0]
        runtime = runner.state["environment"]["runtime"]
        before = copy.deepcopy(runtime.rng[8][0][0].bit_generator.state)
        runner.step()
        assert runner.result["metadata"]["status"]["status"] == "completed"
        assert (runtime.rng[8][0][0].bit_generator.state != before) == shuffle
        states.append((runner.result, runtime.rng[1][0][0].bit_generator.state))
    assert states[0][0]["oviposition"] == states[1][0]["oviposition"]
    assert states[0][0]["offspring"] == states[1][0]["offspring"]
    assert states[0][1] == states[1][1]


def test_zero_generation_is_valid_diagnosis_but_not_transition_acceptance():
    spec = dict(
        base=config_payload(small_config(generations=0)),
        replicates=1,
        cases=[dict(name="initial", overrides={})],
    )
    with TemporaryDirectory() as tmp:
        path = run(spec, Path(tmp) / "zero")
        checked = audit(path, min_transitions=1)
        assert checked["execution_success"] and checked["transitions"] == 0
        assert not checked["acceptance_passed"]
        batch = Batch.open(path)
        try:
            prepare_diagnostics(batch)
        except ValueError as error:
            assert "committed transition" in str(error)
        else:
            raise AssertionError("zero-transition batch entered lifecycle plotting")
        tables, _ = summarize_sensitivity(path)
        assert tables["runs"].audited.all() and tables["runs"].complete.all()
        assert tables["trajectories"].generation.eq(0).all()
        assert tables["diagnostics"].empty
        absent = small_config(plant_sizes=[0], insect_sizes=[0])
        absent.spatial.update(initial_plant_cells=None, initial_insect_cells=None)
        spec["base"] = config_payload(absent)
        empty_path = run(spec, Path(tmp) / "absent")
        checked = audit(empty_path, min_transitions=1)
        assert checked["execution_success"] and checked["transitions"] == 0
        assert not checked["acceptance_passed"]
        tables, _ = summarize_sensitivity(empty_path)
        assert tables["runs"].extinction_generation.eq(0).all()
        assert tables["trajectories"].insects.eq(0).all()


if __name__ == "__main__":
    for name, check in list(globals().items()):
        if name.startswith("test_") and callable(check):
            check()
            print(f"PASS {name}", flush=True)
