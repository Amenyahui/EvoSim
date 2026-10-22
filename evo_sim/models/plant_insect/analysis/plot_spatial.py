"""Grid views and descriptive summaries for the spatial interface demonstration."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm

from evo_sim.models.plant_insect.analysis.common import (
    batch_layout,
    batch_parser,
    pyplot,
    read_audit,
    read_settings,
    sha256,
    write_json,
)

plt = pyplot()


def counts(cells, shape):
    grid = np.zeros(shape, int)
    np.add.at(grid, (cells[:, 0], cells[:, 1]), 1)
    return grid


def render(batch, out=None):
    batch, out, summary_dir = batch_layout(batch, out)
    audit = read_audit(batch)
    settings = read_settings(batch)
    free_initial = settings["spec"]["base"]["model_version"] == "1.6"
    shape = tuple(settings["spec"]["base"]["spatial"]["grid_shape"])
    sample = batch / "data/evolving/rep01"
    with np.load(sample / "snapshots/g000000.npz") as snap:
        grids = [counts(snap["plant_cells"][snap["plant_species"] == s], shape) for s in (0, 1)]
        grids.append(counts(snap["insect_cells"], shape))
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.5))
    for ax, grid, title, cmap in zip(
        axes, grids, ["P1 plants", "P2 plants", "H1 insects"], ["Blues", "Greens", "Oranges"]
    ):
        maxcount = max(1, int(grid.max()))
        norm = BoundaryNorm(np.arange(-0.5, maxcount + 1.5), plt.get_cmap(cmap).N)
        im = ax.imshow(grid, origin="lower", cmap=cmap, norm=norm)
        ax.set_title(f"{title} (N={int(grid.sum())})")
        ax.set_xlabel("Patch column")
        ax.set_ylabel("Patch row")
        ax.set_xticks(range(shape[1]))
        ax.set_yticks(range(shape[0]))
        ax.set_xticks(np.arange(-0.5, shape[1]), minor=True)
        ax.set_yticks(np.arange(-0.5, shape[0]), minor=True)
        ax.grid(which="minor", color="#dddddd", linewidth=0.5)
        fig.colorbar(im, ax=ax, ticks=range(maxcount + 1), fraction=0.045, pad=0.03)
    rows, cols = np.nonzero(grids[2])
    axes[0].scatter(
        cols,
        rows,
        s=50 + 30 * grids[2][rows, cols],
        facecolors="none",
        edgecolors="#D55E00",
        linewidths=1.5,
    )
    fig.suptitle(
        "Initial grid: independent plant and insect locations"
        if free_initial
        else "Initial grid: H1 starts on P1 individuals",
        y=0.98,
    )
    fig.text(
        0.5,
        0.01,
        "Orange circles mark H1 patches over the P1 map; they do not denote feeding. Patches may hold multiple individuals."
        if free_initial
        else "Orange circles mark initial H1 patches on the P1 map. A patch can hold multiple individuals.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.12, 1, 0.94))
    fig.savefig(out / "initial_grid.png", dpi=180, bbox_inches="tight")
    fig.savefig(out / "initial_grid.pdf", bbox_inches="tight")
    plt.close(fig)

    population = []
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
        population.append(pop.assign(case=case, replicate=replicate))
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
    pop = pd.concat(population, ignore_index=True)
    host = pd.DataFrame(hostrows)
    pd.DataFrame(summary).to_csv(summary_dir / "spatial_summary.csv", index=False)
    pd.DataFrame(initial_links).to_csv(summary_dir / "initial_links.csv", index=False)
    host.to_csv(summary_dir / "host_resource_by_generation.csv", index=False)
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.6))
    for case, color, style in [("evolving", "#0072B2", "-"), ("frozen", "#D55E00", "--")]:
        for ax, species in zip(axes[:3], ["P1", "P2", "H1"]):
            frame = pop[(pop.case == case) & (pop.species == species)]
            for _, group in frame.groupby("replicate"):
                ax.plot(
                    group.generation, group.N_start, color=color, alpha=0.2, lw=1, linestyle=style
                )
            average = frame.groupby("generation").N_start.mean()
            ax.plot(average.index, average, color=color, label=case, lw=2, linestyle=style)
            ax.set_title(species)
            ax.set_xlabel("Generation")
            ax.set_ylabel("Individuals")
            ax.grid(alpha=0.15)
        frame = host[host.case == case]
        for _, group in frame.groupby("replicate"):
            axes[3].plot(
                group.generation,
                group.P2_resource_fraction,
                color=color,
                alpha=0.2,
                lw=1,
                linestyle=style,
            )
        average = frame.groupby("generation").P2_resource_fraction.mean()
        axes[3].plot(average.index, average, color=color, lw=2, linestyle=style)
    axes[3].set_title("P2 share of actual feeding")
    axes[3].set_ylim(0, 1)
    axes[3].set_xlabel("Generation")
    axes[3].grid(alpha=0.15)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False)
    fig.text(
        0.5,
        0.01,
        "Thin lines: four seeds per treatment. Bold lines: means. Generation 0 links emerge from local choice."
        if free_initial
        else "Thin lines: four seeds per treatment. Bold lines: means. Generation 0 hosts are prescribed.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.9))
    fig.savefig(out / "spatial_dynamics.png", dpi=180, bbox_inches="tight")
    fig.savefig(out / "spatial_dynamics.pdf", bbox_inches="tight")
    plt.close(fig)
    write_json(
        summary_dir / "plot_spatial_metadata.json",
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
    render(args.input, args.out)
