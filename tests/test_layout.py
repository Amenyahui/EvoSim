"""Run a small batch and check output separation and overwrite protection."""

import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from evo_sim.models.plant_insect.audit import audit
from evo_sim.models.plant_insect.batch import create_output_directory, run

ROOT = Path(__file__).resolve().parents[1]


def test_versioned_output_directory():
    with TemporaryDirectory() as tmp, patch("evo_sim.models.plant_insect.batch.ROOT", Path(tmp)):
        first = create_output_directory()
        marker = first / "preserve.txt"
        marker.write_text("existing results")
        second = create_output_directory()
        assert first.relative_to(tmp).as_posix() == "codex/validation/v0.1.7/run1"
        assert not (Path(tmp) / "results").exists()
        assert second.name == "run2" and marker.read_text() == "existing results"
        try:
            create_output_directory(first)
        except FileExistsError:
            pass
        else:
            raise AssertionError("explicit output overwrote an existing directory")


def test_layout():
    from evo_sim.models.plant_insect import get_runner

    config = get_runner().config
    spec = dict(config["batch"], base=dict(config["parameters"], generations=1))
    spec["base"]["generations"] = 1
    with TemporaryDirectory() as tmp:
        batch = Path(tmp) / "batch"
        run(spec, batch)
        assert {p.name for p in batch.iterdir()} == {"settings.json", "data"}
        assert (batch / "data/run_status.csv").is_file()
        assert not list(batch.rglob("*.py"))
        assert not list(batch.rglob("baseline*"))
        report = Path(tmp) / "reports"
        checked = audit(batch, report)
        assert checked["runs"] == 2 and checked["source_matches_current"]
        assert checked["genetically_audited_runs"] == 2
        assert (report / "integrity_audit.json").is_file()
        assert not (batch / "integrity_audit.json").exists()
        import json

        metadata = next((batch / "data").glob("*/rep01/metadata.json"))
        original = json.loads(metadata.read_text())
        changed = json.loads(metadata.read_text())
        changed["source_hash"]["models/plant_insect/analysis/plot_diagnostics.py"] = "changed"
        metadata.write_text(json.dumps(changed))
        try:
            audit(batch, report)
        except AssertionError:
            pass
        else:
            raise AssertionError("strict audit accepted source drift")
        checked = audit(batch, report, allow_plot_source_drift=True)
        assert len(checked["batch_run_source_differences"]) == 1
        changed["source_hash"]["models/plant_insect/simulator.py"] = "changed"
        metadata.write_text(json.dumps(changed))
        try:
            audit(batch, report, allow_plot_source_drift=True)
        except AssertionError:
            pass
        else:
            raise AssertionError("plot exception accepted simulation source drift")
        metadata.write_text(json.dumps(original))
        try:
            run(spec, batch)
        except FileExistsError:
            pass
        else:
            raise AssertionError("Existing results overwritten")
    optimized = subprocess.run(
        [sys.executable, "-O", "-c", "import evo_sim.models.plant_insect.audit"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert optimized.returncode != 0 and "assertions are required" in optimized.stderr


if __name__ == "__main__":
    test_layout()
    print("PASS output layout, independent audit and overwrite protection")


def test_long_run_inherits_cases():
    from evo_sim.models.plant_insect.helpers.inputs import read_config

    c = read_config(ROOT / "evo_sim/models/plant_insect/configs/sexual_qtl_200gen.yaml")
    assert c["batch"]["replicates"] == 3
    assert [x["name"] for x in c["batch"]["cases"]] == ["qtl_mutating", "mutation_off"]
