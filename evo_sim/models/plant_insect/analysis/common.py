"""Shared helpers for the audited-batch analysis and plotting evo_sim.models.plant_insect.analysis.

Every batch script works on a finished run batch whose integrity audit report
(`tests/reports/<batch>/integrity_audit.json`, written by `python -m evo_sim.models.plant_insect.audit`)
matches `data/run_status.csv`. These helpers centralize that loading and its
standard assertions, the batch directory layout, the script/audit/settings hash
metadata, the headless plotting setup, and the command-line skeleton, so each
script only contains its own analysis logic.
"""

import argparse
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from evo_sim.models.plant_insect.audit import report_directory

BATCH_OUT_HELP = "Figure directory (default: INPUT/figures)"


def audit_file(batch):
    """Path of the integrity audit report for a batch."""
    return report_directory(batch) / "integrity_audit.json"


def read_audit(batch):
    """Parsed integrity audit report for a batch."""
    return json.loads(audit_file(batch).read_text(encoding="utf-8"))


def read_settings(batch):
    """Parsed settings.json (actual parameters and source hashes) for a batch."""
    return json.loads((batch / "settings.json").read_text(encoding="utf-8"))


def sha256(path):
    """SHA-256 hex digest of a file's bytes."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def batch_layout(batch, out=None):
    """Standard batch paths: figures directory `out` and `data/summary`, created if missing.

    Returns the resolved `(batch, out, summary_dir)` tuple.
    """
    batch = Path(batch)
    out = Path(out) if out is not None else batch / "figures"
    summary_dir = batch / "data" / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    return batch, out, summary_dir


def run_path(batch, case, replicate):
    """Directory of one run: `<batch>/data/<case>/rep<NN>`."""
    return Path(batch) / "data" / case / f"rep{int(replicate):02d}"


@dataclass(frozen=True)
class Batch:
    """A run batch whose audit report was cross-checked against run_status.csv."""

    path: Path
    out: Path
    summary_dir: Path
    audit: dict
    status: pd.DataFrame

    @classmethod
    def open(
        cls,
        batch,
        out=None,
        audit_key="runs",
        require_completed=False,
        require_terminal_generation=None,
        expected_runs=None,
    ):
        """Load a batch and verify its audit report against run_status.csv.

        `audit_key` selects which audit counter to compare (`runs` for standard
        audits, `genetically_audited_runs` for the QTL batch audit); the extra
        arguments mirror the per-script acceptance checks.
        """
        path, out, summary_dir = batch_layout(batch, out)
        audit = read_audit(path)
        status = pd.read_csv(path / "data" / "run_status.csv")
        if expected_runs is None:
            assert audit[audit_key] == len(status)
        else:
            assert audit[audit_key] == len(status) == expected_runs
        if require_completed:
            assert status.status.eq("completed").all()
        if require_terminal_generation is not None:
            assert status.generation.eq(require_terminal_generation).all()
        return cls(path, out, summary_dir, audit, status)


def write_json(path, payload):
    """Write JSON with the project's standard formatting (2-space indent, UTF-8)."""
    Path(path).write_text(json.dumps(payload, indent=2), encoding="utf-8")


def pyplot():
    """Headless matplotlib: select the Agg backend before pyplot is imported."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def batch_parser(description):
    """Standard `--input`/`--out` command line shared by the batch evo_sim.models.plant_insect.analysis."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out", type=Path, help=BATCH_OUT_HELP)
    return parser
