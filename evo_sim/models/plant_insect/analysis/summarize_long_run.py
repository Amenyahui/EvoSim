"""Descriptive 200-generation summaries from audited spatial runs."""

from pathlib import Path

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect.analysis.common import (
    Batch,
    audit_file,
    batch_parser,
    read_settings,
    run_path,
    sha256,
    write_json,
)


def summarize(batch, out=None):
    batch = Batch.open(
        batch,
        out,
        expected_runs=8,
        require_completed=True,
        require_terminal_generation=200,
    )
    settings = read_settings(batch.path)
    runs, generations, traits = [], [], []
    for row in batch.status.itertuples():
        path = run_path(batch.path, row.case, row.replicate)
        pop = pd.read_csv(path / "populations.csv")
        net = pd.read_csv(path / "networks.csv")
        events = pd.read_csv(path / "events.csv")
        parents = pd.read_csv(path / "parents.csv", low_memory=False)
        diagnostics = pd.read_csv(path / "diagnostics.csv")
        counts = pop.pivot(index="generation", columns="species", values="N_start")
        resource = net.pivot(index="generation", columns="plant_species", values="Q_actual")
        p2_share = resource.P2 / (resource.P1 + resource.P2).replace(0, np.nan)
        lost = counts.index[counts.H1.eq(0)]
        assert len(lost) and counts.index[-1] == 200
        death = int(lost[0])
        assert counts.H1.loc[death:].eq(0).all()
        assert counts[["P1", "P2"]].gt(0).all().all()
        pre = np.arange(max(0, death - 10), death)
        recent_events = events[events.generation.isin(pre)]
        recent_parents = parents[(parents.kind == "insect") & parents.generation.isin(pre)]
        recent_resource = resource.loc[pre].sum()
        reference = diagnostics[
            (diagnostics.generation == death - 1) & (diagnostics.reference_plant_species == "P2")
        ].iloc[0]
        runs.append(
            dict(
                case=row.case,
                replicate=row.replicate,
                seed=row.seed,
                terminal_generation=row.generation,
                H1_extinction_generation=death,
                H1_peak=int(counts.H1.max()),
                H1_peak_generation=int(counts.H1.idxmax()),
                H1_previous=int(counts.H1.loc[death - 1]),
                P1_at_H1_extinction=int(counts.P1.loc[death]),
                P2_at_H1_extinction=int(counts.P2.loc[death]),
                P1_final=int(counts.P1.loc[200]),
                P2_final=int(counts.P2.loc[200]),
                H1_final=int(counts.H1.loc[200]),
                last10_mean_g=recent_events.g.mean(),
                last10_acceptance_rate=recent_events.accepted.mean(),
                last10_mean_expected_offspring=recent_parents.expected_offspring.mean(),
                last10_P2_resource_fraction=recent_resource.P2 / recent_resource.sum(),
                P2_reference_preference_before_loss=reference.preference,
                P2_reference_q_before_loss=reference.q_unlimited,
            )
        )
        for t in counts.index:
            generations.append(
                dict(
                    case=row.case,
                    replicate=row.replicate,
                    seed=row.seed,
                    generation=int(t),
                    P1=int(counts.P1.loc[t]),
                    P2=int(counts.P2.loc[t]),
                    H1=int(counts.H1.loc[t]),
                    P1_resource=resource.P1.loc[t],
                    P2_resource=resource.P2.loc[t],
                    P2_resource_fraction=p2_share.loc[t],
                    resource_edges=int((resource.loc[t] > 0).sum()),
                )
            )
        for t in (0, death - 1, 200):
            with np.load(path / "snapshots" / f"g{t:06d}.npz") as state:
                for kind, label_key, names in (
                    ("plant", "plant_species", ("D", "plant_voc_independent_trait")),
                    ("insect", "insect_species", ("I", "A")),
                ):
                    for species_id in range(2 if kind == "plant" else 1):
                        selected = state[label_key] == species_id
                        for name in names:
                            for axis in range(state[name].shape[1]):
                                traits.append(
                                    dict(
                                        case=row.case,
                                        replicate=row.replicate,
                                        generation=t,
                                        kind=kind,
                                        species=f"{'P' if kind == 'plant' else 'H'}{species_id + 1}",
                                        trait=f"{name}_{axis + 1}",
                                        N=int(selected.sum()),
                                        mean=float(state[name][selected, axis].mean())
                                        if selected.any()
                                        else np.nan,
                                    )
                                )
    runframe = pd.DataFrame(runs)
    paired = runframe.pivot(index="seed", columns="case", values="H1_extinction_generation")
    paired["evolving_minus_frozen"] = paired.evolving - paired.frozen
    runframe.to_csv(batch.summary_dir / "run_summary.csv", index=False)
    pd.DataFrame(generations).to_csv(batch.summary_dir / "generation_summary.csv", index=False)
    pd.DataFrame(traits).to_csv(batch.summary_dir / "trait_means.csv", index=False)
    paired.to_csv(batch.summary_dir / "paired_extinction.csv")
    write_json(
        batch.summary_dir / "summarize_long_run_metadata.json",
        dict(
            batch=str(batch.path.resolve()),
            audited_runs=batch.audit["runs"],
            audit_sha256=sha256(audit_file(batch.path)),
            settings_sha256=sha256(batch.path / "settings.json"),
            source_sha256=sha256(Path(__file__)),
            purpose=settings["purpose"],
            interpretation="descriptive; paired seeds are not independent evidence of a treatment effect",
        ),
    )
    print(runframe.to_string(index=False))


if __name__ == "__main__":
    args = batch_parser(__doc__).parse_args()
    summarize(args.input, args.out)
