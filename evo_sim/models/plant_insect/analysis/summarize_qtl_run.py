"""Descriptive population and trait summaries from an audited QTL batch."""

from pathlib import Path

import pandas as pd

from evo_sim.models.plant_insect.analysis.common import (
    Batch,
    audit_file,
    batch_parser,
    pyplot,
    run_path,
    sha256,
    write_json,
)

plt = pyplot()


def summarize(batch, out=None):
    batch = Batch.open(
        batch,
        out,
        audit_key="genetically_audited_runs",
        require_completed=True,
        require_terminal_generation=200,
    )
    summaries, populations, traits = [], [], []
    for row in batch.status.itertuples():
        path = run_path(batch.path, row.case, row.replicate)
        pop = pd.read_csv(path / "populations.csv")
        populations.append(pop.assign(case=row.case, replicate=row.replicate))
        catalog = pd.read_csv(path / "allele_catalog.csv")
        mutations = pd.read_csv(path / "mutations.csv")
        parts = []
        for chunk in pd.read_csv(path / "genetic_summary.csv", chunksize=100000):
            parts.append(chunk[chunk.record_type == "trait"])
        trait = pd.concat(parts, ignore_index=True)
        traits.append(trait.assign(case=row.case, replicate=row.replicate))
        record = dict(
            case=row.case,
            replicate=row.replicate,
            seed=row.seed,
            terminal_generation=row.generation,
            allele_catalog_entries=len(catalog),
            mutation_events=len(mutations),
            max_clipped_fraction=trait.clipped_fraction.max(),
        )
        for species, group in pop.groupby("species"):
            lost = group.loc[group.N_start.eq(0), "generation"]
            record[f"{species}_extinction_generation"] = int(lost.min()) if len(lost) else None
            record[f"{species}_final"] = int(group.loc[group.generation.eq(200), "N_start"].iloc[0])
            record[f"{species}_peak"] = int(group.N_start.max())
        summaries.append(record)
    summary = pd.DataFrame(summaries)
    pop = pd.concat(populations, ignore_index=True)
    trait = pd.concat(traits, ignore_index=True)
    summary.to_csv(batch.summary_dir / "run_summary.csv", index=False)
    pop.to_csv(batch.summary_dir / "population_timeseries.csv", index=False)
    trait.to_csv(batch.summary_dir / "trait_timeseries.csv", index=False)
    species = sorted(pop.species.unique())
    colors = dict(zip(species, plt.get_cmap("tab10").colors))
    fig, axes = plt.subplots(1, len(species), figsize=(4.2 * len(species), 4), squeeze=False)
    for ax, name in zip(axes.ravel(), species):
        for i, (case, group) in enumerate(pop[pop.species == name].groupby("case", sort=False)):
            for replicate, values in group.groupby("replicate"):
                ax.plot(
                    values.generation,
                    values.N_start,
                    color=colors[name],
                    linestyle="-" if i == 0 else "--",
                    label=f"{case} / rep{replicate}",
                )
        ax.set(title=name, xlabel="Generation", ylabel="Individuals")
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
    fig.suptitle("200-generation QTL trajectories; one paired seed, descriptive only")
    fig.tight_layout()
    fig.savefig(batch.out / "population_dynamics.png", dpi=160)
    fig.savefig(batch.out / "population_dynamics.pdf")
    plt.close(fig)
    names, cases = trait.trait.unique(), batch.status.case.unique()
    fig, axes = plt.subplots(
        len(cases),
        len(names),
        figsize=(3 * len(names), 3 * len(cases)),
        sharey="col",
        squeeze=False,
    )
    for r, case in enumerate(cases):
        for col, name in enumerate(names):
            ax = axes[r, col]
            for sp, group in trait[(trait.case == case) & (trait.trait == name)].groupby("species"):
                ax.plot(group.generation, group["mean"], color=colors[sp], label=sp)
            ax.set(title=name.replace("plant_voc_independent_trait", "Z"), xlabel="Generation")
            if col == 0:
                ax.set_ylabel(case + "\nMean trait")
            ax.legend(fontsize=8)
            ax.grid(alpha=0.2)
    fig.suptitle("Mean traits; extinct populations have missing values; one seed per treatment")
    fig.tight_layout()
    fig.savefig(batch.out / "trait_means.png", dpi=150)
    fig.savefig(batch.out / "trait_means.pdf")
    plt.close(fig)
    write_json(
        batch.summary_dir / "summarize_qtl_run_metadata.json",
        dict(
            batch=str(batch.path.resolve()),
            source_sha256=sha256(Path(__file__)),
            audit_sha256=sha256(audit_file(batch.path)),
            settings_sha256=sha256(batch.path / "settings.json"),
            interpretation="descriptive only; no independent replicate-based treatment-effect estimate",
        ),
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    args = batch_parser(__doc__).parse_args()
    summarize(args.input, args.out)
