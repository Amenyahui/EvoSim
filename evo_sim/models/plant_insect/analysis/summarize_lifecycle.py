"""Stage-explicit v0.1.6 ledger summaries; no new ecological indices."""

import pandas as pd

from ..audit import audit_run
from .common import batch_parser, read_audit, read_settings, run_path


def summarize(batch):
    settings = read_settings(batch)
    if settings.get("version") != "0.1.6":
        raise ValueError("summarize_lifecycle requires v0.1.6 adult/larval ledgers")
    status = pd.read_csv(batch / "data/run_status.csv")
    assert read_audit(batch)["runs"] == len(status)
    sources = dict(
        lifecycle_counts="populations",
        host_background_flows="host_flows",
        host_background_groups="host_groups",
        lifecycle_networks="networks",
        parent_host_contributions="lagged_networks",
    )
    collected = {name: [] for name in sources}
    for run in status.itertuples(index=False):
        path = run_path(batch, run.case, run.replicate)
        audit_run(path)
        for target, source in sources.items():
            frame = pd.read_csv(path / f"{source}.csv")
            frame.insert(0, "replicate", run.replicate)
            frame.insert(0, "case", run.case)
            collected[target].append(frame)
    output = batch / "data/summary"
    output.mkdir(exist_ok=True)
    result = {}
    for name, frames in collected.items():
        frame = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        frame.to_csv(output / f"{name}.csv", index=False)
        result[name] = frame
    return result


def main():
    args = batch_parser(__doc__).parse_args()
    summarize(args.input)


if __name__ == "__main__":
    main()
