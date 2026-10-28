import copy
import io
import json
from contextlib import redirect_stdout
from dataclasses import MISSING, fields
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
import pandas as pd
import yaml

from evo_sim.core.config import merge_config
from evo_sim.models.plant_insect import Config, get_runner, simulate
from evo_sim.models.plant_insect.analysis.diagnose_initial_conditions import main as diagnose
from evo_sim.models.plant_insect.audit import audit_run
from evo_sim.models.plant_insect.batch import save_result
from evo_sim.models.plant_insect.helpers.genetics import GeneticState
from evo_sim.models.plant_insect.helpers.inputs import read_config
from tests.support import config, raises

CONFIGS = Path(__file__).resolve().parents[1] / "evo_sim/models/plant_insect/configs"


def test_yaml_is_complete_and_missing_parameters_fail():
    for path in CONFIGS.glob("*.yaml"):
        loaded = read_config(path)
        values = dict(
            loaded["parameters"], generations=loaded["simulation_metadata"]["num_steps_per_episode"]
        )
        assert set(values) == {field.name for field in fields(Config)}, path
        for case in loaded["batch"]["cases"]:
            Config(**merge_config(values, case["overrides"])).validate()
    assert all(
        field.default is MISSING and field.default_factory is MISSING
        for field in fields(Config)
        if field.name != "eta_shared"
    )
    values = dict(get_runner().config["parameters"], generations=1)
    for key in values:
        missing = {k: v for k, v in values.items() if k != key}
        if key == "eta_shared":
            raises(lambda: Config(**missing).validate())
        else:
            raises(lambda: Config(**missing), TypeError)
    for key in ("eta_shared", "weights"):
        raises(config(**{key: None}).validate)
    for bounds in ([0.8, 0.2], [-0.1, 1], [0, 1.2], [0, float("nan")], [False, True]):
        raises(config(voc_bounds=bounds).validate)
        c = config()
        c.genetics["phenotype_map"]["bounds"] = bounds
        raises(c.validate)
    for scale in (-1, True, float("inf")):
        c = config()
        c.genetics["phenotype_map"]["effect_scale"] = scale
        raises(c.validate)
    c = config()
    del c.genetics["phenotype_map"]
    raises(c.validate)
    c = config()
    c.genetics["loci"][0]["mutation_mean"][1] = 0.01  # Outside this locus's target.
    raises(c.validate)


def test_nested_yaml_and_case_overrides_preserve_other_parameters():
    raises(lambda: merge_config({}, None))
    original = read_config(CONFIGS / "model_parameters.yaml")
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "experiment.yaml"
        path.write_text(
            yaml.safe_dump(
                dict(
                    extends=str(CONFIGS / "model_parameters.yaml"),
                    parameters=dict(
                        genetics=dict(mutation_enabled=dict(plant=False)),
                        spatial=dict(grid_shape=[8, 12], initial_plant_cells=[[1, 2]]),
                    ),
                    batch=dict(replicates=2),
                )
            ),
            encoding="utf-8",
        )
        loaded = read_config(path)
        assert loaded["parameters"]["genetics"]["mutation_enabled"] == dict(
            plant=False, insect=True
        )
        assert (
            loaded["parameters"]["genetics"]["loci"] == original["parameters"]["genetics"]["loci"]
        )
        assert loaded["parameters"]["spatial"]["grid_shape"] == [8, 12]
        assert loaded["batch"]["cases"] == original["batch"]["cases"]
        child = Path(tmp) / "child.yaml"
        child.write_text(
            "extends: experiment.yaml\nparameters:\n  spatial:\n    initial_plant_cells: null\n",
            encoding="utf-8",
        )
        assert read_config(child)["parameters"]["spatial"]["initial_plant_cells"] is None
    base = copy.deepcopy(original["parameters"])
    base["genetics"]["loci"][0]["mutation_rate"] = 0.123
    case = merge_config(base, original["batch"]["cases"][1]["overrides"])
    assert case["genetics"]["loci"][0]["mutation_rate"] == 0.123
    assert case["genetics"]["mutation_enabled"] == dict(plant=False, insect=False)
    case["genetics"]["loci"][0]["mutation_rate"] = 0
    assert base["genetics"]["loci"][0]["mutation_rate"] == 0.123
    assert base["genetics"]["mutation_enabled"] == dict(plant=True, insect=True)


def test_nondefault_mapping_flows_through_simulation_and_independent_audit():
    for eta in (0.0, 1.0):
        c = config(generations=1, eta_shared=eta)
        c.genetics["phenotype_map"] = dict(effect_scale=0.25, bounds=[0.23, 0.51])
        c.validate()
        result = simulate(c)
        assert result["metadata"]["status"]["status"] == "completed"
        state = result["snapshots"][0]
        genetic = GeneticState(c)
        values, raw = genetic.phenotype(state["plant_genotypes"], "plant")
        default = config()
        _, baseline = GeneticState(default).phenotype(state["plant_genotypes"], "plant")
        np.testing.assert_array_equal(raw, baseline / 2)
        np.testing.assert_array_equal(values["D"], np.clip(raw, 0.23, 0.51)[:, : c.n])
        assert any(
            row["clipped_fraction"] > 0
            for row in result["genetic_summary"]
            if row["record_type"] == "trait"
        )
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "run"
            save_result(result, path)
            assert audit_run(path)["genetic_audit"]


def test_configured_mutation_mean_is_applied_and_audited():
    c = config(generations=1)
    for locus in c.genetics["loci"]:
        locus["mutation_rate"] = 1.0
        locus["mutation_sd"] = [0.0] * len(locus["mutation_sd"])
        for target in locus["targets"]:
            locus["mutation_mean"][target] = -0.01
    result = simulate(c)
    assert result["mutations"]
    assert all(-0.01 in json.loads(row["delta"]) for row in result["mutations"])
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "run"
        save_result(result, path)
        assert audit_run(path)["genetic_audit"]


def test_initial_diagnostics_reads_scenarios_and_replicates_from_yaml():
    with TemporaryDirectory() as tmp:
        path, out = Path(tmp) / "initial.yaml", Path(tmp) / "result"
        scenario = dict(
            extends=str(CONFIGS / "initial_conditions.yaml"),
            batch=dict(
                replicates=2,
                cases=[dict(name="custom", overrides=dict(plant_sizes=[2, 3], insect_sizes=[0]))],
            ),
        )
        path.write_text(yaml.safe_dump(scenario), encoding="utf-8")
        args = ["diagnose", "--config", str(path), "--out", str(out)]
        with patch("sys.argv", args), redirect_stdout(io.StringIO()):
            diagnose()
        frame = pd.read_csv(out / "data/summary/initial_conditions_raw.csv")
        assert frame.replicate.tolist() == [1, 2]
        assert (frame.case == "custom").all()
        assert (frame.plants == 5).all() and (frame.insects == 0).all()
