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
    ordinary = report_directory(batch) / "integrity_audit.json"
    corrected = Path(batch) / "data/summary/corrected_integrity_audit.json"
    return corrected if not ordinary.exists() and corrected.exists() else ordinary


def read_audit(batch):
    """Parsed integrity audit report for a batch."""
    path = audit_file(batch)
    report = json.loads(path.read_text(encoding="utf-8"))
    if path.name != "corrected_integrity_audit.json":
        return report
    # Independently audited corrected runs retain their original report, not a synthetic CLI audit.
    assert report["status"] == "passed"
    rows = report["runs"]
    status = pd.read_csv(Path(batch) / "data/run_status.csv")
    assert report["expected_runs"] == report["audited_runs"] == len(rows) == len(status)
    keys = [(r["case"], r["replicate"], r["seed"]) for r in rows]
    assert len(set(keys)) == len(keys)
    assert set(keys) == set(
        status[["case", "replicate", "seed"]].itertuples(index=False, name=None)
    )
    collected = []
    for name, digest in report["source_reports"].items():
        source = path.with_name(name)
        assert sha256(source) == digest, source
        collected.extend(json.loads(source.read_text(encoding="utf-8"))["runs"])
    assert sorted(collected, key=lambda r: (r["case"], r["replicate"])) == sorted(
        rows, key=lambda r: (r["case"], r["replicate"])
    )
    design_path = path.with_name("corrected_design_consistency.json")
    assert sha256(design_path) == report["design_report_sha256"]
    design = json.loads(design_path.read_text(encoding="utf-8"))
    assert design["status"] == "passed"
    assert sha256(Path(batch) / "settings.json") == design["settings_sha256"]
    assert sha256(Path(batch) / "data/run_status.csv") == design["run_status_sha256"]
    assert {(r["case"], r["replicate"]) for r in design["runs"]} == {k[:2] for k in keys}
    for row in design["runs"]:
        assert all(
            row[key]
            for key in (
                "config_matches_formal",
                "source_matches_frozen",
                "input_sources_match_settings",
                "run_status_matches_metadata",
            )
        )
        assert (
            sha256(run_path(batch, row["case"], row["replicate"]) / "metadata.json")
            == row["metadata_sha256"]
        )
    for row in rows:
        checked = row["audit"]
        terminal = status.loc[
            status.case.eq(row["case"]) & status.replicate.eq(row["replicate"]), "generation"
        ].iloc[0]
        assert row["original_files_unchanged"] and row["frozen_sources_match"]
        assert checked["status"] == "completed"
        assert checked["genetic_audit"] and checked["contribution_audit"]
        assert checked["states"] == terminal + 1
        assert checked["transitions"] == checked["contribution_transitions"] == terminal
    return dict(
        report,
        run_reports=rows,
        runs=len(rows),
        genetically_audited_runs=len(rows),
        audit_method="independent audit_run aggregate",
    )


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
    def open(cls, batch, *, allow_lifecycle=False):
        """Load and check shared batch metadata once; create only the summary directory."""
        path = Path(batch).resolve()
        if (
            (path / "settings.json").exists()
            and read_settings(path).get("version") in {"0.1.6", "0.1.7"}
            and not allow_lifecycle
        ):
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
