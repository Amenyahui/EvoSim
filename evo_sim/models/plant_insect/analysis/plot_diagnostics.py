"""Generate A–E diagnostics and their source tables from a completed audited batch.

Run with python -m evo_sim.models.plant_insect.analysis.plot_diagnostics --input BATCH.
Trajectories describe recorded model behavior; they do not establish causal effects.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon

C = {
    "P1": "#326955",
    "P2": "#92B56D",
    "H1": "#BA6747",
    "off": "#607A8A",
    "on": "#BA6747",
    "gray": "#68757B",
}
CASES = ["mutation_off", "qtl_mutating"]
CASE_NAMES = {"mutation_off": "Mutation off", "qtl_mutating": "QTL mutation on"}
REPS = (1, 2, 3)
HORIZON = 200
GRID = (10, 10)
SPATIAL_REPS = (1, 2)


def theme():
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial", "Liberation Sans"],
            "font.size": 9.1,
            "axes.titlesize": 10.2,
            "axes.labelsize": 9.3,
            "xtick.labelsize": 8.1,
            "ytick.labelsize": 8.1,
            "legend.fontsize": 8.0,
            "axes.linewidth": 0.75,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def axstyle(ax, *, grid=True):
    if grid:
        ax.grid(axis="y", color="#E7ECEC", lw=0.55, zorder=0)
        ax.set_axisbelow(True)
    ax.tick_params(pad=3)


def label(fig, letter, title):
    fig.text(0.025, 0.987, letter, ha="left", va="top", fontsize=17, weight="bold")
    fig.text(0.077, 0.980, title, ha="left", va="top", fontsize=13.3, weight="semibold")


def export(fig, output, stem, dpi):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(output / f"{stem}.{ext}", dpi=dpi, bbox_inches="tight", pad_inches=0.14)
    plt.close(fig)
    print("Saved", stem)


def load_data(path):
    names = [
        "population_timeseries",
        "species_extinction",
        "trait_timeseries",
        "host_resource_by_generation",
        "genetic_timeseries",
    ]
    return {key: pd.read_csv(path / f"{key}.csv") for key in names}


def figure_A(d, out, dpi):
    """Population counts: extinct replicates contribute zero to the unconditional mean."""
    pop, ext = d["population_timeseries"], d["species_extinction"]
    fig, axes = plt.subplots(3, 2, figsize=(9.0, 8.0), sharex=True, sharey="row")
    fig.subplots_adjust(left=0.12, right=0.98, top=0.89, bottom=0.16, hspace=0.25, wspace=0.14)
    label(fig, "A", "Population dynamics")
    for col, case in enumerate(CASES):
        axes[0, col].set_title(CASE_NAMES[case], weight="semibold", pad=10)
        for row, sp in enumerate(("P1", "P2", "H1")):
            ax = axes[row, col]
            s = pop[(pop.case == case) & (pop.species == sp)]
            for j, rep in enumerate(REPS):
                r = s[s.replicate == rep].sort_values("generation")
                ax.plot(
                    r.generation, r.N_start, color=C[sp], alpha=max(0.15, 0.5 - 0.1 * j), lw=1.1
                )
                if sp == "H1":
                    ev = ext[
                        (ext.case == case)
                        & (ext.replicate == rep)
                        & (ext.species == "H1")
                        & (ext.extinction_observed.eq(True))
                    ]
                    if not ev.empty:
                        ax.plot(
                            [ev.event_or_censor_generation.iloc[0]],
                            [0],
                            marker="x",
                            color=C[sp],
                            markersize=5,
                            ls="none",
                            zorder=4,
                        )
            mean = s.groupby("generation", as_index=False).N_start.mean()
            ax.plot(mean.generation, mean.N_start, color=C[sp], lw=2.0, zorder=3)
            ax.set_xlim(0, HORIZON)
            ax.set_xticks(np.linspace(0, HORIZON, 5))
            axstyle(ax)
            if col == 0:
                ax.set_ylabel(f"{sp} · individuals")
    for ax in axes.flat:
        ax.set_ylim(bottom=0)
    for a in axes[2, :]:
        a.set_xlabel("Generation")
    fig.legend(
        handles=[
            Line2D([0], [0], color=C["gray"], lw=2, label=f"Mean of {len(REPS)} runs"),
            Line2D([0], [0], color=C["gray"], lw=1, alpha=0.25, label="Individual run"),
            Line2D([0], [0], color=C["H1"], ls="none", marker="x", label="H1 extinction"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.54, 0.032),
        ncol=3,
        frameon=False,
    )
    fig.text(0.12, 0.085, "H1 means include extinct runs as N = 0.", size=8, color=C["gray"])
    export(fig, out, "A_population_dynamics", dpi)


def centered_xy(t, sp, case, rep, traits):
    """Per-run change since g=0, never fill missing phenotypes after extinction."""
    r = t[(t.species == sp) & (t.case == case) & (t.replicate == rep) & t.trait.isin(traits)]
    p = r.pivot_table(index="generation", columns="trait", values="mean", aggfunc="first")
    if any(k not in p for k in traits):
        return None
    p = p.loc[:, list(traits)].dropna().sort_index()
    if p.empty or 0 not in p.index:
        return None
    z = (p - p.loc[0]).to_numpy()
    return p.index.to_numpy(), z[:, 0], z[:, 1]


def evolution_path(ax, arr, color):
    """Pale to dark tracks simulation time; open circle is gen 0; diamond = last living record."""
    if arr is None:
        return
    g, x, y = arr
    if len(x) > 1:
        coords = np.column_stack([x, y])
        segs = np.stack([coords[:-1], coords[1:]], axis=1)
        from matplotlib.colors import to_rgba

        rgba = np.tile(np.array(to_rgba(color)), (len(segs), 1))
        rgba[:, 3] = 0.20 + 0.73 * (g[1:] / HORIZON)
        lc = LineCollection(segs, colors=rgba, linewidths=1.6, capstyle="round", zorder=2)
        ax.add_collection(lc)
    ax.plot(
        [0], [0], marker="o", mfc="white", mec=color, mew=1.0, markersize=4.0, ls="none", zorder=5
    )
    ax.plot(
        [x[-1]],
        [y[-1]],
        marker="D",
        mfc=color,
        mec="white",
        mew=0.55,
        markersize=4.2,
        ls="none",
        zorder=6,
    )
    # Sparse chronological checkpoints give direction without falsely implying a matching function.
    for k in (50, 100, 150):
        if k in g:
            ix = np.flatnonzero(g == k)[0]
            ax.scatter(x[ix], y[ix], s=7, color=color, alpha=0.4 + 0.5 * k / HORIZON, zorder=4)


def figure_B(d, out, dpi):
    """4 functional 2-D spaces × two treatments; baseline-centering reveals small changes.
    No unvalidated 'functional match' proxy: the model's interaction formulas are not in these CSVs.
    """
    t = d["trait_timeseries"].query("record_type=='trait'")
    chemical = "plant_chemical_traits_1" in set(t.trait)
    plant_pair = (
        ("plant_chemical_traits_1", "plant_chemical_traits_2") if chemical else ("D_1", "D_2")
    )
    odor_pair = (
        ("plant_chemical_traits_2", "plant_chemical_traits_3")
        if chemical
        else ("plant_voc_independent_trait_1", "plant_voc_independent_trait_2")
    )
    panels = [
        (
            "Shared / VOC-only chemical modules" if chemical else "Independent VOC trait",
            odor_pair,
            ("P1", "P2"),
            "Δ x_S" if chemical else "Δ VOC 1",
            "Δ x_V" if chemical else "Δ VOC 2",
        ),
        ("Insect detoxification", ("I_1", "I_2"), ("H1",), "Δ I₁", "Δ I₂"),
        ("Plant defense", plant_pair, ("P1", "P2"), "Δ D₁", "Δ D₂"),
        ("Insect recognition", ("A_1", "A_2"), ("H1",), "Δ A₁", "Δ A₂"),
    ]
    fig, axes = plt.subplots(4, 2, figsize=(9.5, 11.2))
    fig.subplots_adjust(left=0.125, right=0.965, top=0.908, bottom=0.115, hspace=0.43, wspace=0.27)
    label(fig, "B", "Functional trait evolution · changes relative to generation 0")
    for j, case in enumerate(CASES):
        axes[0, j].set_title(CASE_NAMES[case], weight="semibold", pad=10)
    for i, (name, pair, species, xlab, ylab) in enumerate(panels):
        allcoords = []
        for case in CASES:
            for sp in species:
                for rep in REPS:
                    dat = centered_xy(t, sp, case, rep, pair)
                    if dat is not None:
                        allcoords.append(dat)
        if not allcoords:
            continue
        xx = np.concatenate([a[1] for a in allcoords])
        yy = np.concatenate([a[2] for a in allcoords])
        # Identical axes within each treatment pair: critical for comparing evolutionary magnitude.
        xmax = max(abs(xx.min()), abs(xx.max()), 0.002) * 1.13
        ymax = max(abs(yy.min()), abs(yy.max()), 0.002) * 1.13
        for j, case in enumerate(CASES):
            ax = axes[i, j]
            for sp in species:
                for rep in REPS:
                    evolution_path(ax, centered_xy(t, sp, case, rep, pair), C[sp])
            ax.axhline(0, color="#C8D0D0", lw=0.65, zorder=0)
            ax.axvline(0, color="#C8D0D0", lw=0.65, zorder=0)
            ax.set_xlim(-xmax, xmax)
            ax.set_ylim(-ymax, ymax)
            ax.set_xlabel(xlab + " (model units)")
            if j == 0:
                ax.set_ylabel(ylab + " (model units)")
            ax.ticklabel_format(axis="both", style="plain", useOffset=False)
            axstyle(ax, grid=False)
            if j == 0:
                ax.text(
                    -0.26,
                    1.05,
                    name,
                    transform=ax.transAxes,
                    rotation=90,
                    va="top",
                    weight="semibold",
                    fontsize=9,
                    color="#3E4B4D",
                )
    handles = [Line2D([0], [0], color=C[s], lw=2, label=s) for s in ("P1", "P2", "H1")]
    handles += [
        Line2D(
            [0],
            [0],
            color=C["gray"],
            ls="none",
            marker="o",
            mfc="white",
            label="Start: generation 0",
        ),
        Line2D([0], [0], color=C["gray"], ls="none", marker="D", label="Last observed"),
    ]
    fig.legend(
        handles=handles, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.54, 0.047)
    )
    fig.text(
        0.125,
        0.018,
        "Tracks fade from early (pale) to late (dark). Individual runs per species/treatment;"
        " H1 stops at extinction. Axes are baseline-centered, not absolute trait values.",
        fontsize=7.9,
        color=C["gray"],
    )
    export(fig, out, "B_functional_trait_evolution", dpi)


def figure_C(d, out, dpi):
    """Host share observed only while total feeding >0; align real events at extinction."""
    fr = d["host_resource_by_generation"]
    ex = d["species_extinction"]
    pop = d["population_timeseries"]
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 7.9))
    fig.subplots_adjust(left=0.105, right=0.98, top=0.885, bottom=0.158, hspace=0.47, wspace=0.3)
    label(fig, "C", "Host use and extinction bottleneck")
    for j, case in enumerate(CASES):
        ax = axes[0, j]
        r = fr[fr.case == case]
        for rep in REPS:
            s = r[r.replicate == rep].sort_values("generation").copy()
            s.loc[s.total_actual_resource <= 0, "P2_resource_fraction"] = np.nan
            ax.plot(
                s.generation,
                s.P2_resource_fraction,
                color=C["off" if j == 0 else "on"],
                lw=1.45,
                alpha=0.65,
                label=f"Run {rep:02}",
            )
        ax.set_title(CASE_NAMES[case] + " · P2 resource fraction", weight="semibold", pad=9)
        ax.set_xlim(0, HORIZON)
        ax.set_ylim(0, 1)
        ax.set_xlabel("Generation")
        ax.set_ylabel("P2 / all actual plant resource" if j == 0 else "")
        ax.legend(frameon=False, fontsize=7.8, loc="upper left", ncol=3)
        axstyle(ax)

    # Extinction-aligned dynamics for 5 extinct H1 runs; surviving run is not an extinction event.
    ax = axes[1, 0]
    evs = ex[(ex.species == "H1") & (ex.extinction_observed.eq(True))].copy()
    runs = []
    for _, e in evs.iterrows():
        case, rep, at = e["case"], int(e["replicate"]), int(e["event_or_censor_generation"])
        sub = pop[(pop.case == case) & (pop.replicate == rep) & (pop.species == "H1")].copy()
        sub["before"] = sub.generation - at
        sub = sub[(sub.before >= -20) & (sub.before <= 0)].sort_values("before")
        runs.append((case, rep, at, sub))
        ax.plot(
            sub.before,
            sub.N_start,
            color=C["off" if case == "mutation_off" else "on"],
            lw=1.6,
            alpha=0.65,
            label=f"{'Off' if case == 'mutation_off' else 'On'} {rep:02}",
        )
    ax.set_xlim(-20, 0)
    ax.set_ylim(bottom=0)
    ax.axvline(0, lw=0.8, ls=":", color=C["gray"])
    ax.set_xlabel("Generations relative to H1 extinction")
    ax.set_ylabel("H1 population size")
    ax.set_title(f"Population collapse · {len(evs)} extinct runs", weight="semibold", pad=9)
    ax.legend(loc="upper right", frameon=False, fontsize=7.3, ncol=2)
    axstyle(ax)

    # Final generation's numbers; categories differ biologically (survivors vs mated mothers).
    ax = axes[1, 1]
    stages = ["N_start", "N_survived_to_reproduce", "N_mothers_with_mates", "N_offspring"]
    x = np.arange(4)
    for case, rep, at, _ in runs:
        last = pop[
            (pop.case == case)
            & (pop.replicate == rep)
            & (pop.species == "H1")
            & (pop.generation == at - 1)
        ]
        if last.empty:
            continue
        vals = last.iloc[0][stages].astype(float).to_numpy()
        col = C["off" if case == "mutation_off" else "on"]
        ax.plot(x, vals, color=col, alpha=0.65, lw=1.25, marker="o", ms=3.3)
    ax.set_xticks(x, ["Starting", "Surviving", "Mated\nmothers", "Offspring"])
    ax.set_ylabel("Individuals / mothers / offspring")
    ax.set_title("Last reproductive update before extinction", weight="semibold", pad=9)
    ax.set_ylim(bottom=0)
    axstyle(ax)
    fig.text(
        0.105,
        0.043,
        "Each line is a run; no across-run mean for host share (avoids survivor bias)."
        " Mated mothers are a subset, not an additional mortality stage.",
        fontsize=8,
        color=C["gray"],
    )
    export(fig, out, "C_host_use_extinction", dpi)


def cell_abundance(snapshot):
    counts = {}
    for sp in (0, 1):
        cells = snapshot["plant_cells"][snapshot["plant_species"] == sp]
        a = np.zeros(GRID, dtype=int)
        if len(cells):
            np.add.at(a, (cells[:, 0], cells[:, 1]), 1)
        counts[sp] = a
    ins = np.zeros(GRID, dtype=int)
    cells = snapshot["insect_cells"]
    if len(cells):
        np.add.at(ins, (cells[:, 0], cells[:, 1]), 1)
    return counts[0], counts[1], ins


def figure_D(snapshot_dir, out, dpi):
    """True patch-level snapshots; two runs share same mutation-off scenario."""
    gs = tuple(np.linspace(0, HORIZON, 5, dtype=int))
    dat = {}
    for rep in [f"rep{r:02d}" for r in SPATIAL_REPS]:
        for g in gs:
            f = snapshot_dir / "data" / "mutation_off" / rep / "snapshots" / f"g{g:06d}.npz"
            with np.load(f, allow_pickle=False) as a:
                dat[rep, g] = cell_abundance(a)
    maxima = [max(dat[k][i].max() for k in dat) for i in range(3)]
    maxima = [max(1, x) for x in maxima]
    print("D spatial per-cell maxima", maxima)
    cmap1 = LinearSegmentedColormap.from_list("P1", ["#F5F8F6", "#24654E"])
    cmap2 = LinearSegmentedColormap.from_list("P2", ["#F5F8F2", "#91B56A"])
    fig, axes = plt.subplots(2, 5, figsize=(12.5, 5.5))
    fig.subplots_adjust(left=0.11, right=0.985, top=0.83, bottom=0.18, hspace=0.30, wspace=0.13)
    label(fig, "D", "Recorded spatial dynamics · identical mutation treatment")
    for i, rep in enumerate([f"rep{r:02d}" for r in SPATIAL_REPS]):
        row_label = f"Mutation off · {rep}"
        for j, g in enumerate(gs):
            ax = axes[i, j]
            p1, p2, ins = dat[rep, g]
            for y in range(GRID[0]):
                for x in range(GRID[1]):
                    # Opposite halves of the same cell: no mixing P1/P2 abundances.
                    v1 = p1[y, x] / maxima[0]
                    v2 = p2[y, x] / maxima[1]
                    upper = [(x - 0.5, y + 0.5), (x - 0.5, y - 0.5), (x + 0.5, y + 0.5)]
                    lower = [(x + 0.5, y - 0.5), (x - 0.5, y - 0.5), (x + 0.5, y + 0.5)]
                    ax.add_patch(Polygon(upper, facecolor=cmap1(v1), ec="none"))
                    ax.add_patch(Polygon(lower, facecolor=cmap2(v2), ec="none"))
            ys, xs = np.where(ins > 0)
            if len(xs):
                v = ins[ys, xs]
                sizes = 16 * v
                ax.scatter(
                    xs, ys, s=sizes, c=C["H1"], ec="white", linewidths=0.38, alpha=0.93, zorder=3
                )
            ax.set_xlim(-0.5, GRID[1] - 0.5)
            ax.set_ylim(-0.5, GRID[0] - 0.5)
            ax.set_aspect("equal")
            ax.set_xticks([0, GRID[1] - 1])
            ax.set_yticks([0, GRID[0] - 1] if j == 0 else [])
            ax.tick_params(labelsize=7, length=2)
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_color("#CDD5D2")
                spine.set_linewidth(0.45)
            ax.set_title(f"Gen {g}", fontsize=9.6, weight="semibold", pad=4)
            ax.text(
                0.98,
                0.02,
                f"P1 {p1.sum():d} · P2 {p2.sum():d} · H1 {ins.sum():d}",
                transform=ax.transAxes,
                fontsize=6.5,
                ha="right",
                va="bottom",
                bbox=dict(facecolor="white", alpha=0.80, ec="none", pad=1.4),
            )
            if j == 0:
                ax.set_ylabel(row_label, fontsize=8.8, weight="semibold", labelpad=8)
    h = [
        Patch(facecolor=cmap1(0.85), label=f"P1 · upper-left triangle (max {maxima[0]}/cell)"),
        Patch(facecolor=cmap2(0.85), label=f"P2 · lower-right triangle (max {maxima[1]}/cell)"),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            mfc=C["H1"],
            mec="white",
            markersize=8,
            label="H1 · larger circles = more individuals",
        ),
    ]
    fig.legend(
        handles=h,
        loc="lower center",
        bbox_to_anchor=(0.55, 0.073),
        ncol=3,
        frameon=False,
        fontsize=8,
    )
    fig.text(
        0.12,
        0.018,
        "Patch colors represent count intensity on common species-specific scales across snapshots. "
        "These two trajectories are illustrative, not proof of a causal extinction effect.",
        fontsize=7.9,
        color=C["gray"],
    )
    export(fig, out, "D_spatial_dynamics", dpi)


def figure_E(d, out, dpi):
    g = d["genetic_timeseries"].copy()
    parsed = g.case.str.extract(r"^(mutation_off|qtl_mutating)_rep(\d+)$")
    g["condition"] = parsed[0]
    g["replicate"] = pd.to_numeric(parsed[1], errors="coerce")
    g = g.dropna(subset=["condition", "replicate", "gene_diversity"])
    g = g.groupby(
        ["condition", "replicate", "generation", "species"], as_index=False
    ).gene_diversity.mean()
    fig, axes = plt.subplots(3, 2, figsize=(9, 8), sharex=True, sharey="row")
    fig.subplots_adjust(left=0.12, right=0.98, top=0.89, bottom=0.16, hspace=0.25, wspace=0.14)
    label(fig, "E", "Genetic diversity through time")
    for col, case in enumerate(CASES):
        axes[0, col].set_title(CASE_NAMES[case], weight="semibold", pad=10)
        for row, sp in enumerate(("P1", "P2", "H1")):
            ax = axes[row, col]
            s = g[(g.condition == case) & (g.species == sp)]
            for j, rep in enumerate(REPS):
                r = s[s.replicate == rep].sort_values("generation")
                ax.plot(
                    r.generation,
                    r.gene_diversity,
                    color=C[sp],
                    alpha=max(0.2, 0.6 - 0.1 * j),
                    lw=1.0,
                )
            if sp != "H1" and s.groupby("generation").replicate.nunique().eq(len(REPS)).all():
                m = s.groupby("generation", as_index=False).gene_diversity.mean()
                ax.plot(m.generation, m.gene_diversity, color=C[sp], lw=2.1)
            ax.set_xlim(0, HORIZON)
            ax.set_ylim(0, 1.01)
            ax.set_xticks(np.linspace(0, HORIZON, 5))
            axstyle(ax)
            if col == 0:
                ax.set_ylabel(f"{sp} · gene diversity ($H_e$)")
    for ax in axes[-1, :]:
        ax.set_xlabel("Generation")
    fig.legend(
        handles=[
            Line2D([0], [0], color=C["gray"], lw=1, label="Individual run"),
            Line2D([0], [0], color=C["gray"], lw=2.1, label=f"Plant mean (n = {len(REPS)})"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.54, 0.029),
        ncol=2,
        frameon=False,
    )
    fig.text(
        0.12,
        0.080,
        "Within-run mean across recorded QTL loci; H1 is shown by run only (no surviving-only mean).",
        size=8,
        color=C["gray"],
    )
    export(fig, out, "E_genetic_diversity", dpi)


def main():
    import json

    from .common import Batch, sha256, write_json

    global REPS, HORIZON, GRID, SPATIAL_REPS
    ap = argparse.ArgumentParser(
        description="A–E descriptive diagnostics from audited batch summaries"
    )
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--dpi", type=int, default=350)
    ap.add_argument(
        "--spatial-reps",
        type=int,
        nargs=2,
        default=[1, 2],
        help="Explicit illustrative mutation-off runs; not a random sample",
    )
    a = ap.parse_args()
    batch = Batch.open(a.input)
    out = a.out or batch.path / "figures"
    if not batch.status.status.eq("completed").all():
        raise ValueError("A–E comparison requires complete runs; do not average censored runs")
    REPS = tuple(sorted(batch.status.replicate.unique()))
    if set(batch.status.case) != set(CASES):
        raise ValueError("A–E requires mutation_off and qtl_mutating")
    for case in CASES:
        assert set(batch.status.loc[batch.status.case.eq(case), "replicate"]) == set(REPS)
    HORIZON = int(batch.status.generation.iloc[0])
    assert HORIZON >= 4 and batch.status.generation.eq(HORIZON).all()
    SPATIAL_REPS = tuple(a.spatial_reps)
    assert set(SPATIAL_REPS) <= set(REPS)
    meta = json.loads(
        (a.input / "data/mutation_off" / f"rep{REPS[0]:02d}" / "metadata.json").read_text(
            encoding="utf-8"
        )
    )
    GRID = tuple(meta["config"]["spatial"]["grid_shape"])
    from .analyze_macro import analyze
    from .summarize_qtl_run import summarize
    from .summarize_spatial import summarize as spatial_summary

    summarize(batch)
    spatial_summary(batch)
    analyze(batch)
    d = load_data(batch.summary_dir)
    for key in ["population_timeseries", "trait_timeseries"]:
        frame = d[key]
        columns = ["case", "replicate", "generation", "species"]
        if key == "trait_timeseries":
            columns += ["trait"]
        assert not frame.duplicated(columns).any(), key
    out.mkdir(parents=True, exist_ok=True)
    theme()
    figure_A(d, out, a.dpi)
    figure_B(d, out, a.dpi)
    figure_C(d, out, a.dpi)
    figure_D(batch.path, out, a.dpi)
    figure_E(d, out, a.dpi)
    write_json(
        batch.summary_dir / "plot_diagnostics_metadata.json",
        dict(
            source_sha256=sha256(Path(__file__)),
            inputs={
                str(p.relative_to(batch.path)): sha256(p) for p in batch.summary_dir.glob("*.csv")
            },
            spatial_replicates=list(SPATIAL_REPS),
            interpretation="Exploratory diagnostics; selected spatial examples; no causal inference",
        ),
    )


if __name__ == "__main__":
    main()
