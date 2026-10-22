"""Run a small batch and check output separation and overwrite protection."""

import json
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

from evo_sim.models.plant_insect.audit import audit
from evo_sim.models.plant_insect.batch import run

ROOT = Path(__file__).resolve().parents[1]


def test_layout():
    spec = json.loads((ROOT / "tests/fixtures/qtl_development.json").read_text())
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
