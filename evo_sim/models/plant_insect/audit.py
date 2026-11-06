"""Independent ledger integrity and explicit execution acceptance for the current model."""

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from evo_sim.core.config import merge_config
from evo_sim.core.helpers import source_hashes

from .batch import SUCCESS_STATUSES, write_json
from .helpers.config import Config, config_payload

if not __debug__:
    raise RuntimeError("Run scientific audits without python -O; assertions are required")


def report_directory(batch):
    return Path(batch).resolve() / "audits"


def audit_run(path):
    """Check saved records; a valid failure ledger does not imply successful execution."""
    path = Path(path)
    meta = json.loads((path / "metadata.json").read_text(encoding="utf-8"))
    if meta["version"] != "0.1.11":
        raise ValueError("Use the frozen historical source to audit a superseded model protocol")
    status = meta["status"]["status"]
    assert json.loads((path / "status.json").read_text(encoding="utf-8")) == meta["status"]
    for name, digest in meta["output_sha256"].items():
        assert hashlib.sha256((path / name).read_bytes()).hexdigest() == digest, name
    if status == "invalid_parameter":
        try:
            Config(**meta["config"]).validate()
        except (ValueError, TypeError, KeyError):
            return dict(
                status=status,
                execution_success=False,
                states=0,
                transitions=0,
                max_overrun=0.0,
                max_balance_error=0.0,
            )
        raise AssertionError("invalid_parameter record has valid configuration")
    from .audit_lifecycle import audit_run_lifecycle

    return dict(audit_run_lifecycle(path, meta), execution_success=status in SUCCESS_STATUSES)


def audit(out, report_dir=None, allow_plot_source_drift=False, *, min_transitions=0):
    if type(min_transitions) is not int or min_transitions < 0:
        raise ValueError("min_transitions must be a nonnegative integer")
    out = Path(out)
    settings = json.loads((out / "settings.json").read_text(encoding="utf-8"))
    report_dir = Path(report_dir) if report_dir is not None else report_directory(out)
    report_dir.mkdir(parents=True, exist_ok=True)
    status = pd.read_csv(out / "data" / "run_status.csv")
    rows = []
    source_differences = []
    allowed_drift = (
        {
            "models/plant_insect/analysis/plot_diagnostics.py",
            "models/plant_insect/analysis/summarize_qtl_run.py",
        }
        if allow_plot_source_drift
        else set()
    )
    assert len(status) == len(settings["spec"]["cases"]) * settings["spec"]["replicates"]
    for case in settings["spec"]["cases"]:
        for rep in range(1, settings["spec"]["replicates"] + 1):
            path = out / "data" / case["name"] / f"rep{rep:02d}"
            meta = json.loads((path / "metadata.json").read_text(encoding="utf-8"))
            differences = {
                key: {
                    "batch": settings["source_sha256"].get(key),
                    "run": meta["source_hash"].get(key),
                }
                for key in set(settings["source_sha256"]) | set(meta["source_hash"])
                if settings["source_sha256"].get(key) != meta["source_hash"].get(key)
            }
            assert set(differences) <= allowed_drift, differences
            if differences:
                source_differences.append(dict(case=case["name"], replicate=rep, files=differences))
            expected = merge_config(settings["spec"]["base"], case["overrides"])
            expected["seed"] += rep - 1
            if meta["status"]["status"] != "invalid_parameter":
                expected = config_payload(Config(**expected))
            assert meta["config"] == expected
            checked = audit_run(path)
            row = status[(status["case"] == case["name"]) & (status.replicate == rep)].iloc[0]
            assert (
                checked["status"] == row.status and meta["status"]["generation"] == row.generation
            )
            rows.append(dict(case=case["name"], replicate=rep, **checked))
            print("AUDITED", case["name"], rep, checked["status"], flush=True)
    pd.DataFrame(rows).to_csv(report_dir / "development_checks.csv", index=False)
    summary = dict(
        integrity_status="passed",
        execution_success=all(r["execution_success"] for r in rows),
        successful_runs=sum(r["execution_success"] for r in rows),
        failed_runs=sum(not r["execution_success"] for r in rows),
        zero_transition_runs=sum(r["transitions"] == 0 for r in rows),
        minimum_transitions_per_run=min_transitions,
        meets_minimum_transitions=all(r["transitions"] >= min_transitions for r in rows),
        acceptance_passed=all(
            r["execution_success"] and r["transitions"] >= min_transitions for r in rows
        ),
        runs=len(rows),
        allow_plot_source_drift=allow_plot_source_drift,
        batch_run_source_differences=source_differences,
        outcomes=pd.Series([r["status"] for r in rows]).value_counts().to_dict(),
        batch=str(out.resolve()),
        source_matches_current=settings["source_sha256"] == source_hashes(),
        source_verification="comparison with current files only; historical source snapshots are not verified",
        states=sum(r["states"] for r in rows),
        transitions=sum(r["transitions"] for r in rows),
        event_rows=sum(r.get("event_rows", 0) for r in rows),
        plant_rows=sum(r.get("plant_rows", 0) for r in rows),
        offspring_rows=sum(r.get("offspring_rows", 0) for r in rows),
        max_resource_overrun=max(r["max_overrun"] for r in rows),
        max_relative_balance_error=max(r["max_balance_error"] for r in rows),
    )
    if any(r.get("genetic_audit", False) for r in rows):
        summary.update(
            genetically_audited_runs=sum(r.get("genetic_audit", False) for r in rows),
            alleles_checked=sum(r.get("alleles_checked", 0) for r in rows),
            mutations_checked=sum(r.get("mutations_checked", 0) for r in rows),
        )
    write_json(report_dir / "integrity_audit.json", summary)
    print(json.dumps(summary, indent=2))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument(
        "--allow-plot-source-drift",
        action="store_true",
        help="Record and allow only the explicitly listed historical analysis-script changes",
    )
    parser.add_argument(
        "--require-success",
        action="store_true",
        help="Exit nonzero unless every run completed or all species became extinct",
    )
    parser.add_argument(
        "--min-transitions",
        type=int,
        default=0,
        help="Require at least this many committed transitions in every run for acceptance",
    )
    args = parser.parse_args()
    summary = audit(
        args.input,
        allow_plot_source_drift=args.allow_plot_source_drift,
        min_transitions=args.min_transitions,
    )
    if (args.require_success or args.min_transitions > 0) and not summary["acceptance_passed"]:
        raise SystemExit("Execution acceptance failed; integrity report retained for diagnosis")


if __name__ == "__main__":
    main()
