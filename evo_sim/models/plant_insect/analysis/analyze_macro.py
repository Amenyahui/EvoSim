"""Analyze every pilot run, including extinction and computational stops."""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect.analysis.common import (
    Batch,
    batch_parser,
    read_settings,
    run_path,
    sha256,
    write_json,
)


def effective_hosts(probabilities):
    values = np.asarray(probabilities, float)
    if not np.isfinite(values).all() or (values < 0).any() or not np.isclose(values.sum(), 1):
        return np.nan
    positive = values[values > 0]
    return float(np.exp(-np.sum(positive * np.log(positive))))


def edge_turnover(before, after):
    union = before | after
    return 1 - len(before & after) / len(union) if union else np.nan


def analyze(batch: Batch):
    settings = read_settings(batch.path)
    generation_rows, runs, extinctions, growth = [], [], [], []
    for record in batch.status.itertuples():
        folder = run_path(batch.path, record.case, record.replicate)
        meta = json.loads((folder / "metadata.json").read_text())
        c = meta["config"]
        context = dict(
            case=record.case,
            replicate=record.replicate,
            seed=record.seed,
            insect_lambda=c["insect_offspring_lambda"]["H1"],
            mode=record.case,
        )
        valid = record.status in ("completed", "all_species_extinct")
        if not valid:
            runs.append(
                dict(
                    context,
                    status=record.status,
                    terminal_generation=record.generation,
                    requested_generations=c["generations"],
                    biologically_resolved=False,
                )
            )
            continue
        pop = pd.read_csv(folder / "populations.csv")
        net = pd.read_csv(folder / "networks.csv")
        diag = pd.read_csv(folder / "diagnostics.csv")
        parents = pd.read_csv(folder / "parents.csv", low_memory=False)
        initial = parents[parents.generation == 0]
        for (kind, species), group in initial.groupby(["kind", "species"]):
            growth.append(
                dict(
                    context,
                    kind=kind,
                    species=species,
                    initial_expected_maternal_offspring=group.expected_maternal_offspring.mean(),
                    initial_maternal_role_count=int(
                        group.expected_maternal_offspring.notna().sum()
                    ),
                    initial_expected_maternal_offspring_per_individual=group.expected_maternal_offspring.fillna(
                        0
                    ).mean(),
                    initial_survival=group.S_i.mean(),
                    initial_fecundity=group.F_i.mean(),
                )
            )
        previous = None
        for t, p in pop.groupby("generation", sort=True):
            edges = net[net.generation == t]
            active = edges[edges.Q_actual > 0]
            current = set(zip(active.insect_species, active.plant_species))
            plants = p[(p.kind == "plant") & (p.N_start > 0)]
            insects = p[(p.kind == "insect") & (p.N_start > 0)]
            breadth, capacity = [], []
            for _, h in diag[diag.generation == t].groupby("insect_species"):
                if h.insect_count.iloc[0] > 0:
                    breadth.append(effective_hosts(h.preference.to_numpy()))
                    capacity.append(h.q_unlimited.mean())
            generation_rows.append(
                dict(
                    context,
                    generation=int(t),
                    plant_richness=len(plants),
                    insect_richness=len(insects),
                    plant_population=int(plants.N_start.sum()),
                    insect_population=int(insects.N_start.sum()),
                    resource_edges=len(current),
                    connectance_initial=len(current)
                    / (len(c["plant_sizes"]) * len(c["insect_sizes"])),
                    connectance_extant=len(current) / (len(plants) * len(insects))
                    if len(plants) * len(insects)
                    else np.nan,
                    edge_turnover=edge_turnover(previous, current)
                    if previous is not None
                    else np.nan,
                    reference_population_effective_hosts=np.mean(breadth) if breadth else np.nan,
                    reference_mean_q=np.mean(capacity) if capacity else np.nan,
                    resource_transfer=edges.Q_actual.sum(),
                    observed=True,
                )
            )
            previous = current
        local = [
            row
            for row in generation_rows
            if row["case"] == context["case"] and row["replicate"] == context["replicate"]
        ]
        end = local[-1]
        losses = []
        for species, group in pop.groupby("species"):
            zero = group[group.N_start == 0]
            event = not zero.empty
            time = int(zero.generation.min()) if event else c["generations"]
            if event:
                losses.append(time)
            extinctions.append(
                dict(
                    context,
                    species=species,
                    kind=group.kind.iloc[0],
                    extinction_observed=event,
                    event_or_censor_generation=time,
                    right_censored=not event,
                    termination_status=record.status,
                )
            )
        runs.append(
            dict(
                context,
                status=record.status,
                terminal_generation=int(end["generation"]),
                requested_generations=c["generations"],
                biologically_resolved=True,
                first_species_loss=min(losses) if losses else np.nan,
                plant_richness_end=end["plant_richness"],
                insect_richness_end=end["insect_richness"],
                plant_population_end=end["plant_population"],
                insect_population_end=end["insect_population"],
                total_extinction_generation=int(end["generation"])
                if record.status == "all_species_extinct"
                else np.nan,
                all_initial_species_retained=end["plant_richness"] == len(c["plant_sizes"])
                and end["insect_richness"] == len(c["insect_sizes"]),
                observed_population_updates=int(end["generation"]),
            )
        )
    frames = dict(
        generations=pd.DataFrame(generation_rows),
        runs=pd.DataFrame(runs),
        species_extinction=pd.DataFrame(extinctions),
        initial_growth=pd.DataFrame(growth),
    )
    for name, frame in frames.items():
        frame.to_csv(batch.summary_dir / f"{name}.csv", index=False)
    runframe = frames["runs"]
    manifest = dict(
        batch=str(batch.path.resolve()),
        source_hash=sha256(Path(__file__)),
        settings_sha256=sha256(batch.path / "settings.json"),
        runs=len(runframe),
        outcomes=runframe.status.value_counts().to_dict(),
        purpose=settings["purpose"],
        inference="descriptive exploratory only",
        breadth_definition="entropy of species-mean preference; not mean individual specialization",
        extinct_extension="Plots show absorbing zero states through requested horizon; raw tables retain observed times",
    )
    write_json(batch.summary_dir / "analyze_macro_metadata.json", manifest)

    print(json.dumps(manifest, indent=2))
    return frames


if __name__ == "__main__":
    args = batch_parser(__doc__).parse_args()
    analyze(Batch.open(args.input))
