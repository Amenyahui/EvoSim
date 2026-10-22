"""Plot an audited short QTL fixture; colors encode species, never seed."""

import math
from pathlib import Path

import numpy as np
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


def render(batch, out=None):
    batch = Batch.open(batch, out, require_completed=True)
    # ponytail: one-seed development figures; replicate-level uncertainty requires a separate experiment.
    if not batch.status.replicate.eq(1).all():
        raise ValueError("This short-fixture plot supports one seed per treatment")
    frames = []
    input_hashes = {}
    for row in batch.status.itertuples():
        path = run_path(batch.path, row.case, 1) / "genetic_summary.csv"
        input_hashes[str(path.relative_to(batch.path))] = sha256(path)
        data = pd.read_csv(path)
        data["case"] = row.case
        frames.append(data)
    data = pd.concat(frames, ignore_index=True)
    # ponytail: locus statistics repeat on every allele row; keep one row per locus.
    columns = [
        "case",
        "generation",
        "kind",
        "species",
        "N",
        "locus_id",
        "observed_heterozygosity",
        "gene_diversity",
        "mutation_count",
    ]
    data.loc[data.record_type.eq("allele"), columns].drop_duplicates().to_csv(
        batch.summary_dir / "genetic_timeseries.csv", index=False
    )
    species = sorted(data.species.unique())
    colors = dict(zip(species, plt.get_cmap("tab10").colors))
    traits = data[data.record_type == "trait"]
    names = list(traits.trait.unique())
    cases = batch.status.case.tolist()
    fig, axes = plt.subplots(
        len(cases),
        len(names),
        figsize=(3 * len(names), 3.2 * len(cases)),
        squeeze=False,
        sharey="col",
    )
    for r, case in enumerate(cases):
        for col, name in enumerate(names):
            ax = axes[r, col]
            subset = traits[(traits.case == case) & (traits.trait == name)]
            for sp, group in subset.groupby("species"):
                group = group.sort_values("generation")
                x, y, sd = (
                    group.generation.to_numpy(),
                    group["mean"].to_numpy(),
                    np.sqrt(group.variance.to_numpy()),
                )
                ax.plot(x, y, color=colors[sp], label=sp)
                ax.fill_between(x, y - sd, y + sd, color=colors[sp], alpha=0.15)
            ax.set(title=name.replace("plant_voc_independent_trait", "Z"), xlabel="Generation")
            ax.grid(alpha=0.2)
            if col == 0:
                ax.set_ylabel(case + "\nMean trait")
            ax.legend(fontsize=8)
    fig.suptitle(
        "One seed per treatment; shading = within-population SD, not replicate uncertainty"
    )
    fig.tight_layout()
    fig.savefig(batch.out / "trait_trajectories.png", dpi=150)
    fig.savefig(batch.out / "trait_trajectories.pdf")
    plt.close(fig)
    alleles = data[data.record_type == "allele"]
    for (case, kind), table in alleles.groupby(["case", "kind"]):
        loci = list(table.locus_id.unique())
        fig, axes = plt.subplots(
            math.ceil(len(loci) / 4), 4, figsize=(14, 3 * math.ceil(len(loci) / 4)), squeeze=False
        )
        for ax, locus in zip(axes.ravel(), loci):
            subset = table[table.locus_id == locus]
            for (sp, allele), group in subset.groupby(["species", "allele_id"]):
                group = group.sort_values("generation")
                if group.frequency.fillna(0).max() > 0:
                    ax.plot(
                        group.generation, group.frequency, color=colors[sp], alpha=0.6, linewidth=1
                    )
            for sp in subset.species.unique():
                ax.plot([], [], color=colors[sp], label=sp)
            ax.set(title=locus, xlabel="Generation", ylabel="Allele frequency", ylim=(-0.02, 1.02))
            ax.legend(fontsize=8)
        for ax in axes.ravel()[len(loci) :]:
            ax.set_visible(False)
        fig.suptitle(f"{case} / {kind}: one line per allele; color identifies species")
        fig.tight_layout()
        fig.savefig(batch.out / f"alleles_{case}_{kind}.png", dpi=150)
        fig.savefig(batch.out / f"alleles_{case}_{kind}.pdf")
        plt.close(fig)
    write_json(
        batch.summary_dir / "plot_genetics_metadata.json",
        dict(
            batch=str(batch.path.resolve()),
            inputs=input_hashes,
            audit_sha256=sha256(audit_file(batch.path)),
            source_sha256=sha256(Path(__file__)),
            interpretation="implementation demonstration only; no treatment-effect inference",
        ),
    )


if __name__ == "__main__":
    args = batch_parser(__doc__).parse_args()
    render(args.input, args.out)
