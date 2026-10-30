"""Stage-explicit v0.1.6 ledger summaries; no new ecological indices."""

import numpy as np
import pandas as pd

from ..audit import audit_run
from .common import batch_parser, read_audit, read_settings, run_path


def prepare_diagnostics(batch):
    """Reuse stage-neutral trait/genetic summaries; add adult endpoints and larval resource shares."""
    from .summarize_qtl_run import summarize as summarize_genetics

    assert read_settings(batch.path)["version"] in {"0.1.6", "0.1.7"}
    assert batch.audit["genetically_audited_runs"] == len(batch.status)
    summarize_genetics(batch, print_summary=False)
    populations, extinctions, resources = [], [], []
    for run in batch.status.itertuples(index=False):
        folder = run_path(batch.path, run.case, run.replicate)
        pop = pd.read_csv(folder / "populations.csv")
        terminal = int(run.generation)
        with np.load(folder / "snapshots" / f"g{terminal:06d}.npz") as state:
            for (kind, species), rows in pop.groupby(["kind", "species"]):
                count = int(np.count_nonzero(state[kind + "_species"] == int(species[1:]) - 1))
                populations.append(
                    dict(
                        case=run.case,
                        replicate=run.replicate,
                        generation=terminal,
                        kind=kind,
                        species=species,
                        N_start=count,
                        update_applied=False,
                    )
                )
                lost = rows.extinct_at_generation.dropna()
                observed = len(lost) > 0
                extinctions.append(
                    dict(
                        case=run.case,
                        replicate=run.replicate,
                        species=species,
                        extinction_observed=observed,
                        event_or_censor_generation=int(lost.min()) if observed else terminal,
                    )
                )
        populations.extend(pop.assign(case=run.case, replicate=run.replicate).to_dict("records"))
        net = pd.read_csv(folder / "networks.csv")
        for t, rows in net[net.insect_species == "H1"].groupby("generation"):
            total = rows.Q_actual.sum()
            resources.append(
                dict(
                    case=run.case,
                    replicate=run.replicate,
                    generation=t,
                    total_actual_resource=total,
                    P2_resource_fraction=rows.loc[rows.plant_species == "P2", "Q_actual"].sum()
                    / total
                    if total > 0
                    else np.nan,
                )
            )
    tables = dict(
        population_timeseries=pd.DataFrame(populations),
        species_extinction=pd.DataFrame(extinctions),
        host_resource_by_generation=pd.DataFrame(resources),
    )
    for name, frame in tables.items():
        frame.to_csv(batch.summary_dir / f"{name}.csv", index=False)
    # Keep the reused run summary's endpoint fields consistent with the terminal snapshot.
    runs = pd.read_csv(batch.summary_dir / "run_summary.csv")
    for endpoint in populations:
        if endpoint.get("update_applied") is False:
            mask = (runs.case == endpoint["case"]) & (runs.replicate == endpoint["replicate"])
            runs.loc[mask, endpoint["species"] + "_final"] = endpoint["N_start"]
            peak = endpoint["species"] + "_peak"
            runs.loc[mask, peak] = runs.loc[mask, peak].clip(lower=endpoint["N_start"])
    for event in extinctions:
        mask = (runs.case == event["case"]) & (runs.replicate == event["replicate"])
        runs.loc[mask, event["species"] + "_extinction_generation"] = (
            event["event_or_censor_generation"] if event["extinction_observed"] else np.nan
        )
    runs.to_csv(batch.summary_dir / "run_summary.csv", index=False)
    print(runs.to_string(index=False))


def summarize(batch):
    settings = read_settings(batch)
    if settings.get("version") not in {"0.1.6", "0.1.7"}:
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
