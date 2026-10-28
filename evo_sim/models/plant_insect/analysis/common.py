"""Shared paths, audit loading and metadata utilities for batch analysis."""

import argparse
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from evo_sim.models.plant_insect.audit import report_directory


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


def run_path(batch, case, replicate):
    """Directory of one run: `<batch>/data/<case>/rep<NN>`."""
    return Path(batch) / "data" / case / f"rep{int(replicate):02d}"


@dataclass(frozen=True)
class Batch:
    """A run batch whose audit report was cross-checked against run_status.csv."""

    path: Path
    summary_dir: Path
    audit: dict
    status: pd.DataFrame

    @classmethod
    def open(cls, batch):
        """Load and check shared batch metadata once; create only the summary directory."""
        path = Path(batch)
        if (path / "settings.json").exists() and read_settings(path).get("version") == "0.1.6":
            raise ValueError(
                "Legacy summaries assume self-feeding parents; v0.1.6 uses larvae. Read the stage-specific raw ledgers or summarize_lifecycle instead."
            )
        audit = read_audit(path)
        status = pd.read_csv(path / "data" / "run_status.csv")
        assert audit["runs"] == len(status)
        summary_dir = path / "data" / "summary"
        summary_dir.mkdir(parents=True, exist_ok=True)
        return cls(path, summary_dir, audit, status)


def write_json(path, payload):
    """Write JSON with the project's standard formatting (2-space indent, UTF-8)."""
    Path(path).write_text(json.dumps(payload, indent=2), encoding="utf-8")


def batch_parser(description):
    """Input argument for data-only summary commands."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--input", type=Path, required=True)
    return parser
