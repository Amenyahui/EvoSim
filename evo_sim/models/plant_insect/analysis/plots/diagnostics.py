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
from matplotlib.ticker import MaxNLocator

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
    handles = [
        Line2D([0], [0], color=C["gray"], lw=2, label=f"Mean of {len(REPS)} runs"),
        Line2D([0], [0], color=C["gray"], lw=1, alpha=0.25, label="Individual run"),
    ]
    if ((ext.species == "H1") & ext.extinction_observed.eq(True)).any():
        handles.append(
            Line2D([0], [0], color=C["H1"], ls="none", marker="x", label="H1 extinction")
        )
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.54, 0.032),
        ncol=len(handles),
        frameon=False,
    )
    fig.text(0.12, 0.085, "Means include extinct runs as N = 0.", size=8, color=C["gray"])
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
    panels = [
        (
            "Independent volatile signals",
            ("plant_chemical_traits_4", "plant_chemical_traits_5"),
            ("P1", "P2"),
            "Δ VOC1",
            "Δ VOC2",
        ),
        (
            "Leaf PSM concentrations",
            ("plant_chemical_traits_1", "plant_chemical_traits_2"),
            ("P1", "P2"),
            "Δ PSM1",
            "Δ PSM2",
        ),
    ] + [
        (name, tuple(f"insect_module_traits_{i}" for i in pair), ("H1",), xlab, ylab)
        for name, pair, xlab, ylab in (
            ("Processing expression", (1, 2), "Δ met E₁", "Δ met E₂"),
            ("Processing tuning", (3, 4), "Δ met μ₁", "Δ met μ₂"),
            ("Sensory expression", (5, 6), "Δ sense E₁", "Δ sense E₂"),
            ("Sensory tuning", (7, 8), "Δ sense μ₁", "Δ sense μ₂"),
            ("Contact response coordinates", (9, 10), "Δ G1 raw", "Δ G2 raw"),
        )
    ]
    fig, axes = plt.subplots(len(panels), 2, figsize=(9.5, 3.0 + 2.05 * len(panels)))
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
            short = (
                "Plant VOC",
                "Plant defense",
                "Metabolic expression",
                "Metabolic tuning",
                "Sensory expression",
                "Sensory tuning",
                "Contact response",
            )[i]
            ax.set_title(
                f"{CASE_NAMES[case]}\n{short}" if i == 0 else short,
                weight="semibold",
                pad=7,
                fontsize=9.2,
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
    evs = ex[(ex.species == "H1") & (ex.extinction_observed.eq(True))].copy()
    extinct = not evs.empty
    linestyles = {rep: ("-", "--", ":")[i % 3] for i, rep in enumerate(REPS)}
    label(fig, "C", "Host use and " + ("extinction bottleneck" if extinct else "adult recruitment"))
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
                ls=linestyles[rep],
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

    # Censored surviving runs must never be treated as extinction events.
    ax = axes[1, 0]
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
            ls=linestyles[rep],
            alpha=0.65,
            label=f"{'Off' if case == 'mutation_off' else 'On'} {rep:02}",
        )
    if not extinct:
        for case in CASES:
            for rep in REPS:
                sub = pop[
                    (pop.case == case)
                    & (pop.replicate == rep)
                    & (pop.species == "H1")
                    & (pop.generation >= HORIZON - 20)
                ].sort_values("generation")
                runs.append((case, rep, HORIZON, sub))
                ax.plot(
                    sub.generation,
                    sub.N_start,
                    color=C["off" if case == "mutation_off" else "on"],
                    lw=1.6,
                    ls=linestyles[rep],
                    alpha=0.65,
                    label=f"{'Off' if case == 'mutation_off' else 'On'} {rep:02}",
                )
    ax.set_xlim((-20, 0) if extinct else (HORIZON - 20, HORIZON))
    ax.set_ylim(bottom=0)
    if extinct:
        ax.axvline(0, lw=0.8, ls=":", color=C["gray"])
    ax.set_xlabel("Generations relative to H1 extinction" if extinct else "Generation")
    ax.set_ylabel("H1 population size")
    ax.set_title(
        f"Population collapse · {len(evs)} extinct runs"
        if extinct
        else "No H1 extinction · final 20 generations",
        weight="semibold",
        pad=9,
    )
    ax.legend(loc="upper right", frameon=False, fontsize=7.3, ncol=2)
    axstyle(ax)

    # Final generation's numbers; categories differ biologically (survivors vs mated mothers).
    ax = axes[1, 1]
    if d.get("lifecycle_schema") != "dose-tolerance-lifecycle-2":
        raise ValueError("Current recruitment diagnostics require stage-explicit larval ledgers")
    stages = ["N_start", "N_mothers_with_mates", "N_potential_eggs", "N_offspring", "N_next"]
    x = np.arange(len(stages))
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
        ax.plot(x, vals, color=col, alpha=0.65, lw=1.25, ls=linestyles[rep], marker="o", ms=3.3)
    ax.set_xticks(
        x,
        ["Adults", "Mated\nfemales", "Potential\neggs", "Larval\nbirths", "Adult\nrecruits"],
    )
    ax.set_ylabel("Individuals / eggs (stage-specific)")
    ax.set_title(
        "Last reproductive update before extinction"
        if extinct
        else f"Final reproductive update · {HORIZON - 1} → {HORIZON}",
        weight="semibold",
        pad=9,
    )
    ax.set_ylim(bottom=0)
    axstyle(ax)
    fig.text(
        0.105,
        0.043,
        "Each line is a run; no across-run mean for host share (avoids survivor bias)."
        " Females/eggs/larvae/recruits are distinct stages.",
        fontsize=8,
        color=C["gray"],
    )
    export(fig, out, "C_host_use_extinction" if extinct else "C_host_use_recruitment", dpi)


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


def lineage_plot_tables(tables):
    """Pool saved runs within a case using shared weights and complete copy denominators."""
    o = tables["matched_opportunities"]
    o = o[(o.target_group == "__all_cross__") & (o.weight_type == "conditional_mu")].copy()
    o["local_mass"] = o.p_local * o.matched_weight
    o["global_mass"] = o.p_global * o.matched_weight
    keys = ["case", "reproduction_round", "group_system"]
    opportunity = o.groupby(keys, as_index=False)[
        ["weight", "matched_weight", "no_mate_weight", "local_mass", "global_mass"]
    ].sum(min_count=1)
    denominator = opportunity.matched_weight.where(opportunity.matched_weight > 0)
    opportunity["p_local"] = opportunity.local_mass / denominator
    opportunity["p_global"] = opportunity.global_mass / denominator
    opportunity["coverage"] = opportunity.matched_weight / opportunity.weight.where(
        opportunity.weight > 0
    )
    opportunity["ri_relative"] = 1 - opportunity.p_local / opportunity.p_global.where(
        opportunity.p_global > 0
    )

    a = tables["qtl_ancestry"]
    a = a[a.stage.isin(["reference", "recruitment"]) & a.status.eq("available")].copy()
    if a.replicate.nunique() > 1:
        selected = a.source_system.eq("candidate_cluster") & ~a.source_group.isin(
            ["unknown", "unassigned"]
        )
        a.loc[selected, "source_group"] = (
            "r" + a.loc[selected, "replicate"].astype(str) + ":" + a.loc[selected, "source_group"]
        )
    keys = ["case", "reference_generation", "generation", "source_system"]
    # Each target's denominator repeats for every source; count it exactly once.
    totals = a.drop_duplicates(keys + ["replicate", "target_system", "target_group"])
    totals = totals.groupby(keys, as_index=False).total_qtl_copies.sum()
    a["source_group"] = a.source_group.replace({"unassigned": "unknown"})
    ancestry = a.groupby(keys + ["source_group"], as_index=False).qtl_copies.sum()
    ancestry = ancestry.merge(totals, on=keys, validate="many_to_one")
    ancestry["fraction"] = ancestry.qtl_copies / ancestry.total_qtl_copies.where(
        ancestry.total_qtl_copies > 0
    )

    n = tables["networks"]
    n = n[n.edge_semantics.eq("maternal_source")].copy()
    members = tables["memberships"]
    for kind, column in (("insect", "insect_cluster"), ("plant", "plant_cluster")):
        identity = members[members.kind.eq(kind)].drop_duplicates(
            ["case", "replicate", "generation", "cluster_id"]
        )
        identity = identity[["case", "replicate", "generation", "cluster_id", "lineage_id"]]
        identity = identity.rename(
            columns={"generation": "reproduction_round", "cluster_id": column}
        )
        n = n.merge(
            identity,
            on=["case", "replicate", "reproduction_round", column],
            how="left",
            validate="many_to_one",
        )
        if n.lineage_id.isna().any():
            raise ValueError("network identity missing from saved memberships")
        n[kind + "_lineage"] = n.pop("lineage_id")
        if n.replicate.nunique() > 1:
            n[kind + "_lineage"] = "r" + n.replicate.astype(str) + ":" + n[kind + "_lineage"]
    network = n.groupby(
        ["case", "reproduction_round", "insect_lineage", "plant_lineage"], as_index=False
    ).actual_consumption.sum()
    demographic = (
        tables["network_balances"]
        .groupby(["case", "reproduction_round"], as_index=False)[
            ["birth_count", "recruitment_count", "actual_consumption"]
        ]
        .sum()
    )
    expected = demographic.set_index(["case", "reproduction_round"]).actual_consumption
    actual = network.groupby(["case", "reproduction_round"]).actual_consumption.sum()
    if not np.allclose(actual.reindex(expected.index, fill_value=0), expected):
        raise ValueError("plotted network does not conserve consumption")
    return dict(
        opportunity=opportunity, ancestry=ancestry, network=network, demographic=demographic
    )


def plot_lineages(batch, out, dpi):
    """Render saved lineage observations, preserving unestimable probabilities as missing."""
    import json

    from evo_sim.core.io import atomic_csv
    from evo_sim.models.plant_insect.analysis.common import sha256, write_json

    folder = batch.summary_dir / "lineages"
    metadata = json.loads((folder / "metadata.json").read_text(encoding="utf-8"))
    for name, expected in metadata["input_sha256"].items():
        if sha256(batch.path / name) != expected:
            raise ValueError(f"lineage source changed: {name}")
    names = ["matched_opportunities", "qtl_ancestry", "networks", "memberships", "network_balances"]
    tables = {name: pd.read_csv(folder / f"{name}.csv") for name in names}
    data = lineage_plot_tables(tables)
    out.mkdir(parents=True, exist_ok=True)
    for name, frame in data.items():
        atomic_csv(frame, folder / f"plot_{name}.csv")
    theme()
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    cases = sorted(batch.status.case.unique())
    case_names = {"mutation_off": "关闭突变", "qtl_mutating": "开启 QTL 突变"}
    colors = ["#326955", "#BA6747", "#607A8A", "#9A79AA", "#C69D45"]

    def canvas(title, note):
        fig, axes = plt.subplots(len(cases), 2, figsize=(12, 3.5 * len(cases) + 1), squeeze=False)
        fig.subplots_adjust(left=0.10, right=0.97, top=0.88, bottom=0.21, hspace=0.70, wspace=0.38)
        fig.suptitle(title, fontsize=18, weight="bold", y=0.98)
        runs = batch.status.groupby("case").size()
        fig.text(
            0.10,
            0.92,
            "验证数据 · " + " / ".join(f"{case_names.get(c, c)}：{runs[c]} 次运行" for c in cases),
            color=C["gray"],
        )
        fig.text(0.10, 0.04, note, fontsize=10, color=C["gray"])
        return fig, axes

    fig, axes = canvas(
        "交配机会：局部与全局对照",
        "权重＝潜在卵数期望；全局限同一兼容池。覆盖率＝完整分类且有配偶的母本权重 / 全部母本权重。\n"
        "NA 不是零隔离：随机跨群概率为零或没有匹配母本时，RI 无定义。不同运行按权重合并。",
    )
    systems = [("birth_host", "出生寄主"), ("candidate_cluster", "候选遗传群")]
    for row, case in enumerate(cases):
        part = data["opportunity"].query("case == @case")
        if part.empty:
            for ax in axes[row]:
                ax.text(0.5, 0.5, "没有可绘制的交配机会记录", ha="center", transform=ax.transAxes)
            continue
        for j, (system, name) in enumerate(systems):
            s = part[part.group_system.eq(system)].sort_values("reproduction_round")
            for col, field in enumerate(["p_local", "coverage"]):
                ax = axes[row, col]
                ax.plot(s.reproduction_round, s[field], "o-", color=colors[j], label=name)
                if col == 0:
                    ax.plot(
                        s.reproduction_round,
                        s.p_global,
                        "x--",
                        color=colors[j],
                        label=name + "·全局",
                    )
                    missing = s[s.p_local.isna()]
                    for t in missing.reproduction_round:
                        ax.text(
                            t,
                            0.18 + j * 0.18,
                            name + " NA",
                            ha="center",
                            color=colors[j],
                            fontsize=8,
                        )
                ax.set_ylim(-0.05, 1.05)
                ax.set_xticks(sorted(part.reproduction_round.unique()))
                ax.set_xlim(
                    part.reproduction_round.min() - 0.4, part.reproduction_round.max() + 0.4
                )
                ax.set_xlabel("繁殖轮次 t")
                axstyle(ax)
        for col, title in enumerate(["跨群父本概率（实线＝局部）", "可估计风险集覆盖率"]):
            axes[row, col].set_title(case_names.get(case, case) + " · " + title)
            axes[row, col].legend(
                fontsize=8, loc="upper right" if col == 0 else "center right", frameon=False
            )
        if not part.ri_relative.notna().any():
            axes[row, 0].text(0.02, 0.90, "本组 RI 均不可估计", transform=axes[row, 0].transAxes)
    export(fig, out, "F_mating_opportunities", dpi)

    fig, axes = canvas(
        "跨代祖源：参考代来源如何留存在成年群体中",
        "每列为该来源分组在全体成年个体的 QTL 拷贝份额，分母包含未知来源。两种分组分别归一化，不能相加。\n"
        "仅追踪已观测 QTL；来源固定在参考代。未分类来源不等于新物种，来源份额不等于每代新增基因流。",
    )
    for row, case in enumerate(cases):
        for col, (system, name) in enumerate(systems):
            ax = axes[row, col]
            a = data["ancestry"]
            a = a[a.case.eq(case) & a.source_system.eq(system)]
            if a.empty:
                ax.text(0.5, 0.5, "未计算祖源或没有成年个体", ha="center", transform=ax.transAxes)
                continue
            p = a.pivot(index="generation", columns="source_group", values="fraction").fillna(0)
            bottom = np.zeros(len(p))
            for j, source in enumerate(p.columns):
                if not p[source].any():
                    continue
                name_source = "未分类来源" if source == "unknown" else source.replace("insect:", "")
                color = "#B8BEC3" if source == "unknown" else colors[j % len(colors)]
                ax.bar(
                    p.index, p[source], bottom=bottom, color=color, width=0.55, label=name_source
                )
                bottom += p[source].to_numpy()
            ax.set_title(
                f"{case_names.get(case, case)} · {name}（参考代 {int(a.reference_generation.iloc[0])}）"
            )
            ax.set_xticks(p.index)
            ax.set_ylim(0, 1.06)
            ax.set_xlabel("成年世代")
            ax.set_ylabel("QTL 来源份额")
            ax.legend(
                frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.19), ncol=2
            )
    export(fig, out, "G_qtl_ancestry", dpi)

    fig, axes = canvas(
        "取食网络与招募：资源流向和后代去向",
        "左图按已保存的操作谱系身份合并母本来源 → 植物来源的取食量，包含随后死亡的幼虫及未分类来源。\n"
        "右图为每轮出生和成年招募人数。不同运行计数相加；图中的谱系节点不是已确认物种。",
    )
    vmax = max(1, data["network"].actual_consumption.max())
    for row, case in enumerate(cases):
        n = data["network"].query("case == @case").copy()

        def short(value):
            return (
                value.replace("insect:", "")
                .replace("plant:", "")
                .replace("unassigned", "未分类")
                .replace("unknown", "未分类")
            )

        n["edge"] = n.insect_lineage.map(short) + " → " + n.plant_lineage.map(short)
        p = n.pivot(index="edge", columns="reproduction_round", values="actual_consumption").fillna(
            0
        )
        ax = axes[row, 0]
        empty_network = p.empty
        if empty_network:
            p = pd.DataFrame([[np.nan]], index=["无取食记录"], columns=["NA"])
        im = ax.imshow(p, cmap="YlGnBu", vmin=0, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(p.columns)), p.columns)
        ax.set_yticks(range(len(p)), p.index, fontsize=8)
        for i in range(len(p)):
            for j in range(len(p.columns)):
                value = p.iloc[i, j]
                ax.text(
                    j,
                    i,
                    "无取食记录" if empty_network else f"{value:.1f}",
                    ha="center",
                    va="center",
                    color="white" if value > vmax * 0.55 else "#16324A",
                )
        fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03, label="实际资源量")
        ax.set_xlabel("繁殖轮次 t")
        ax.set_title(case_names.get(case, case) + " · 母本来源网络")
        d = data["demographic"].query("case == @case")
        ax = axes[row, 1]
        for offset, field, name, color in [
            (-0.18, "birth_count", "出生幼虫", colors[2]),
            (0.18, "recruitment_count", "招募成虫", colors[0]),
        ]:
            bars = ax.bar(
                d.reproduction_round + offset, d[field], width=0.34, label=name, color=color
            )
            if len(d) <= 20:
                ax.bar_label(bars, padding=3, fontsize=9)
        ax.set_ylim(0, max(1, d.birth_count.max()) * 1.22)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
        ax.set_xlabel("繁殖轮次 t（招募进入 t+1 代）")
        ax.set_ylabel("个体数")
        ax.set_title(case_names.get(case, case) + " · 出生与招募")
        ax.legend(frameon=False)
        axstyle(ax)
    export(fig, out, "H_lineage_network_recruitment", dpi)
    write_json(
        folder / "plot_metadata.json",
        dict(
            plot_source_sha256=sha256(Path(__file__)),
            input_sha256={f"{name}.csv": sha256(folder / f"{name}.csv") for name in names},
            lineage_metadata_sha256=sha256(folder / "metadata.json"),
            aggregation="Within-case weight/copy/count pooling; descriptive, not causal or species evidence",
        ),
    )


def plot_sensitivity(path, out, dpi):
    """Plot audited per-run trajectories and fixed-horizon uncertainty."""
    from evo_sim.models.plant_insect.analysis.common import sha256, write_json
    from evo_sim.models.plant_insect.analysis.summaries.lifecycle import summarize_sensitivity

    tables, settings = summarize_sensitivity(path)
    cases = [c["name"] for c in settings["spec"]["cases"]]
    horizon = int(settings["spec"]["base"]["generations"])
    out.mkdir(parents=True, exist_ok=True)
    theme()
    plt.rcParams["font.sans-serif"] = [
        "Microsoft YaHei",
        "SimHei",
        "Noto Sans CJK SC",
        "DejaVu Sans",
    ]
    plt.rcParams["axes.unicode_minus"] = False
    chinese = {
        "baseline": "基线",
        "resource2": "资源 R=2",
        "resource4": "资源 R=4",
        "resource10": "资源 R=10",
        "resource30": "资源 R=30",
        "resource100": "资源 R=100",
        "fecundity8": "潜在卵 λ=8",
        "search3": "搜索半径 3",
        "mate3": "交配半径 3",
        "dispersal2": "昆虫扩散半径 2",
        "resource2_search3": "资源 R=2＋搜索半径 3",
    }
    labels = {case: chinese.get(case, case) for case in cases}
    fig, axes = plt.subplots(
        (len(cases) + 1) // 2, 2, figsize=(13, 2.8 * ((len(cases) + 1) // 2)), squeeze=False
    )
    for ax, case in zip(axes.flat, cases):
        frame = tables["trajectories"]
        if not frame.empty:
            for rep, rows in frame[frame.case.eq(case)].groupby("replicate"):
                ax.plot(rows.generation, rows.insects, lw=1, label=f"重复 {rep}")
        ax.set_title(labels[case], loc="left", fontsize=8)
        ax.set_ylabel("昆虫数")
        ax.set_xlabel("世代")
        ax.set_xlim(0, horizon)
        ax.set_ylim(bottom=0)
        axstyle(ax)
    for ax in list(axes.flat)[len(cases) :]:
        ax.set_visible(False)
    fig.suptitle("各情景全部重复：昆虫数量", y=1.0)
    fig.tight_layout()
    export(fig, out, "sensitivity_trajectories", dpi)

    fig, axes = plt.subplots(1, 2, figsize=(12, max(4, len(cases) * 0.55)))
    end = tables["survival"].query("generation == @horizon").set_index("case").loc[cases]
    for y, (_, row) in enumerate(end.iterrows()):
        axes[0].plot([row.lower, row.upper], [y, y], lw=6, color=C["off"])
        if not np.isnan(row.wilson_low):
            axes[0].plot([row.wilson_low, row.wilson_high], [y, y], color="black", lw=1.3)
        axes[0].scatter([row.lower, row.upper], [y, y], color=C["off"], s=12)
        axes[0].text(
            1.02,
            y,
            f"存续 {int(row.alive)} / 灭绝 {int(row.extinct)} / 未知 {int(row.unknown)}",
            va="center",
            fontsize=8,
            transform=axes[0].get_yaxis_transform(),
        )
    axes[0].set_yticks(range(len(cases)), [labels[case] for case in cases])
    axes[0].set_xlim(0, 1)
    axes[0].set_xlabel("存续比例：粗线为上下界；黑线为 Wilson 95%")
    axes[0].invert_yaxis()
    runs = tables["runs"]
    for y, case in enumerate(cases):
        sub = runs[(runs.case == case) & runs.complete & runs.final_insects.notna()]
        axes[1].scatter(
            sub.final_insects,
            y + (sub.replicate - (settings["spec"]["replicates"] + 1) / 2) * 0.06,
            alpha=0.65,
        )
    axes[1].set_yticks(range(len(cases)), [labels[case] for case in cases])
    axes[1].invert_yaxis()
    axes[1].set_xlabel("终点昆虫数（完成运行）")
    fig.suptitle(f"第 {horizon} 代存续与终点分布")
    fig.tight_layout()
    export(fig, out, "sensitivity_endpoints", dpi)

    fig, axes = plt.subplots(
        len(cases), 3, figsize=(13, max(4, len(cases) * 1.8)), squeeze=False, sharex=True
    )
    diag = tables["diagnostics"]
    metrics = ("mate_fraction", "acceptance", "resource_gain")
    for row, case in enumerate(cases):
        if not diag.empty:
            for _, records in diag[diag.case.eq(case)].groupby("replicate"):
                for ax, metric in zip(axes[row], metrics, strict=True):
                    ax.plot(records.generation, records[metric], alpha=0.55, lw=0.8)
        axes[row, 0].set_ylabel(labels[case])
        for ax in axes[row]:
            ax.set_ylim(0, 1.05)
            axstyle(ax)
    for ax, title in zip(
        axes[0], ("有配偶雌虫 / 雌虫", "出生幼虫 / 潜在卵", "实际摄入 / 幼虫最大需求"), strict=True
    ):
        ax.set_title(title)
    for ax in axes[-1]:
        ax.set_xlabel("繁殖轮次")
    fig.suptitle("阶段诊断：每条线为一个重复，分母为零时缺失")
    fig.tight_layout()
    export(fig, out, "sensitivity_stages", dpi)
    fig, axes = plt.subplots(
        len(cases), 3, figsize=(13, max(4, len(cases) * 1.8)), squeeze=False, sharex=True
    )
    for row, case in enumerate(cases):
        if not diag.empty:
            for _, records in diag[diag.case.eq(case)].groupby("replicate"):
                for ax, metric in zip(
                    axes[row],
                    (
                        "plant_retention_fraction",
                        "plant_establishment_fraction",
                        "plant_resource_fraction",
                    ),
                    strict=True,
                ):
                    ax.plot(records.generation, records[metric], alpha=0.55, lw=0.8)
        axes[row, 0].set_ylabel(labels[case])
        for ax in axes[row]:
            ax.set_ylim(0, 1.05)
            axstyle(ax)
    for ax, title in zip(
        axes[0],
        ("老株保留 / 轮初植物", "新建立株 / 繁殖提案", "轮初资源存量 / 最大容量"),
        strict=True,
    ):
        ax.set_title(title)
    for ax in axes[-1]:
        ax.set_xlabel("昆虫轮次")
    fig.suptitle("植物：保留、建立与资源存量；无提案时建立比例缺失")
    fig.tight_layout()
    export(fig, out, "sensitivity_resource_plants", dpi)
    write_json(
        path / "data/summary/sensitivity_plot_metadata.json",
        dict(
            source_sha256=sha256(Path(__file__)),
            settings_sha256=sha256(path / "settings.json"),
            audit_sha256=sha256(path / "data/summary/sensitivity_audit.json"),
            summary_sha256={
                p.name: sha256(p) for p in (path / "data/summary").glob("sensitivity_*.csv")
            },
        ),
    )


def main():
    import json

    from evo_sim.models.plant_insect.analysis.common import (
        Batch,
        audit_file,
        read_settings,
        run_path,
        sha256,
        write_json,
    )

    global REPS, HORIZON, GRID, SPATIAL_REPS
    ap = argparse.ArgumentParser(
        description="A–E descriptive diagnostics from audited batch summaries"
    )
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--dpi", type=int, default=350)
    ap.add_argument("--lineages", action="store_true", help="Plot saved lineage summaries as F–H")
    ap.add_argument(
        "--sensitivity", action="store_true", help="Plot multi-case sensitivity diagnostics"
    )
    ap.add_argument(
        "--spatial-reps",
        type=int,
        nargs=2,
        default=[1, 2],
        help="Explicit illustrative mutation-off runs; not a random sample",
    )
    a = ap.parse_args()
    if a.sensitivity:
        if a.lineages:
            ap.error("--sensitivity and --lineages are mutually exclusive")
        plot_sensitivity(a.input.resolve(), a.out or a.input / "figures/sensitivity", a.dpi)
        return
    if a.lineages:
        batch = Batch.open(a.input)
        plot_lineages(batch, a.out or batch.path / "figures/lineages", a.dpi)
        return
    version = read_settings(a.input).get("version")
    if version != "0.1.12":
        raise ValueError(
            "Current diagnostics require model protocol 0.1.12; use frozen history for old results"
        )
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
    raw_inputs = {}
    if version == "0.1.12":
        for run in batch.status.itertuples(index=False):
            folder = run_path(batch.path, run.case, run.replicate)
            metadata = folder / "metadata.json"
            recorded = json.loads(metadata.read_text(encoding="utf-8"))
            assert recorded["functional_schema"] == "expression-tuning-contact-3"
            assert recorded["genetics_version"] == "slim-qtl-6"
            assert "adult_settlement_enabled" not in recorded["config"]["functional_modules"]
            assert "flowering_kappa" not in recorded["config"]["functional_modules"]
            names = [
                "populations.csv",
                "genetic_summary.csv",
                "allele_catalog.csv",
                "mutations.csv",
                "networks.csv",
            ]
            generations = {HORIZON}
            if run.case == "mutation_off" and run.replicate in SPATIAL_REPS:
                generations.update(np.linspace(0, HORIZON, 5, dtype=int))
            names += [f"snapshots/g{g:06d}.npz" for g in sorted(generations)]
            for name in names:
                source = folder / name
                digest = sha256(source)
                assert digest == recorded["output_sha256"][name], source
                raw_inputs[str(source.relative_to(batch.path))] = digest
            raw_inputs[str(metadata.relative_to(batch.path))] = sha256(metadata)
    from evo_sim.models.plant_insect.analysis.summaries.lifecycle import prepare_diagnostics

    prepare_diagnostics(batch)
    d = load_data(batch.summary_dir)
    d["lifecycle_schema"] = "dose-tolerance-lifecycle-2"
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
            analysis_source_sha256={
                "plot_diagnostics.py": sha256(Path(__file__)),
                "common.py": sha256(Path(__file__).parents[1] / "common.py"),
                "summarize_qtl_run.py": sha256(Path(__file__).parents[1] / "summaries/qtl.py"),
                "summarize_lifecycle.py": sha256(
                    Path(__file__).parents[1] / "summaries/lifecycle.py"
                ),
            },
            inputs={
                str(p.relative_to(batch.path)): sha256(p) for p in batch.summary_dir.glob("*.csv")
            },
            raw_inputs=raw_inputs,
            audit_path=str(audit_file(batch.path).relative_to(batch.path)),
            audit_sha256=sha256(audit_file(batch.path)),
            audit_method=batch.audit.get("audit_method", "batch integrity audit"),
            functional_schema=meta.get("functional_schema"),
            genetics_schema=meta.get("genetics_version"),
            spatial_replicates=list(SPATIAL_REPS),
            lifecycle_schema=d["lifecycle_schema"],
            adapter_sha256=sha256(Path(__file__).parents[1] / "summaries/lifecycle.py"),
            interpretation="Exploratory diagnostics; selected spatial examples; no causal inference",
        ),
    )


if __name__ == "__main__":
    main()
