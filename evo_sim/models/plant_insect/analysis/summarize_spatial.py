"""Summarize spatial populations and measured host resource use."""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect.analysis.common import (
    Batch,
    batch_parser,
    sha256,
    write_json,
)


def summarize(batch: Batch):
    summary_dir, audit = batch.summary_dir, batch.audit
    batch = batch.path
    hostrows = []
    summary = []
    initial_links = []
    for folder in sorted((batch / "data").glob("*/rep*")):
        case = folder.parent.name
        replicate = int(folder.name[3:])
        pop = pd.read_csv(folder / "populations.csv")
        net = pd.read_csv(folder / "networks.csv")
        ev = pd.read_csv(folder / "events.csv")
        meta = json.loads((folder / "metadata.json").read_text())
        for t, group in net.groupby("generation"):
            total = group.Q_actual.sum()
            p2 = group[group.plant_species == "P2"].Q_actual.sum()
            hostrows.append(
                dict(
                    case=case,
                    replicate=replicate,
                    generation=t,
                    P2_resource_fraction=p2 / total if total else np.nan,
                    P2_actual_resource=p2,
                    total_actual_resource=total,
                )
            )
        first = ev[ev.generation == 0]
        initial_links.extend(
            dict(
                case=case,
                replicate=replicate,
                insect_species=row.insect_species,
                plant_species=row.plant_species,
                N_exposed=int(row.candidate_exposure),
                N_accept=int(row.N_accept),
                Q_actual=float(row.Q_actual),
            )
            for row in net[net.generation == 0].itertuples()
        )
        p2 = net[(net.plant_species == "P2") & (net.Q_actual > 0)]
        end = pop[pop.generation == pop.generation.max()]
        row = dict(
            case=case,
            replicate=replicate,
            status=meta["status"]["status"],
            terminal_generation=int(end.generation.iloc[0]),
            initial_P1_hosts=int(first.plant_species.eq("P1").sum()),
            initial_P2_hosts=int(first.plant_species.eq("P2").sum()),
            initial_unmatched=int(first.host_row.eq(-1).sum()),
            initial_insects=len(first),
            first_P2_resource_generation=int(p2.generation.min()) if len(p2) else None,
        )
        row.update({f"final_{r.species}": int(r.N_start) for r in end.itertuples()})
        summary.append(row)
    host = pd.DataFrame(hostrows)
    pd.DataFrame(summary).to_csv(summary_dir / "spatial_summary.csv", index=False)
    pd.DataFrame(initial_links).to_csv(summary_dir / "initial_links.csv", index=False)
    host.to_csv(summary_dir / "host_resource_by_generation.csv", index=False)
    write_json(
        summary_dir / "summarize_spatial_metadata.json",
        dict(
            batch=str(batch.resolve()),
            audited_runs=audit["runs"],
            source_sha256=sha256(Path(__file__)),
            settings_sha256=sha256(batch / "settings.json"),
            note="Interface demonstration; P2 feeding is not by itself evidence of evolved host switching",
        ),
    )
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    args = batch_parser(__doc__).parse_args()
    summarize(Batch.open(args.input))
