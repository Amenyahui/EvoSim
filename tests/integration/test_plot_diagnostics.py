import numpy as np
import pandas as pd

from evo_sim.models.plant_insect.analysis.plots.diagnostics import cell_abundance, centered_xy


def test_lineage_plot_weighting_missingness_and_copy_denominators():
    from evo_sim.models.plant_insect.analysis.plots.diagnostics import lineage_plot_tables

    base = dict(case="x", replicate=1)
    tables = {
        "matched_opportunities": pd.DataFrame(
            [
                dict(
                    **base,
                    reproduction_round=0,
                    group_system="birth_host",
                    target_group="__all_cross__",
                    weight_type="conditional_mu",
                    weight=4,
                    matched_weight=w,
                    no_mate_weight=0,
                    p_local=p,
                    p_global=0.5,
                )
                for w, p in [(1, 0.2), (3, 0.8)]
            ]
        ),
        "qtl_ancestry": pd.DataFrame(
            [
                dict(
                    **base,
                    reference_generation=1,
                    generation=1,
                    source_system="birth_host",
                    target_system="candidate_cluster",
                    stage="reference",
                    status="available",
                    target_group=target,
                    source_group=source,
                    qtl_copies=count,
                    total_qtl_copies=total,
                )
                for target, source, count, total in [
                    ("a", "P1", 1, 2),
                    ("a", "unknown", 1, 2),
                    ("b", "P1", 4, 4),
                    ("b", "unknown", 0, 4),
                ]
            ]
        ),
        "networks": pd.DataFrame(
            [
                dict(
                    **base,
                    reproduction_round=0,
                    edge_semantics="maternal_source",
                    insect_cluster=label,
                    plant_cluster="p",
                    actual_consumption=x,
                )
                for label, x in [("a", 2), ("unassigned", 1)]
            ]
        ),
        "memberships": pd.DataFrame(
            [
                dict(**base, generation=0, kind=kind, cluster_id=label, lineage_id=label)
                for kind, label in [("insect", "a"), ("insect", "unassigned"), ("plant", "p")]
            ]
        ),
        "network_balances": pd.DataFrame(
            [
                dict(
                    **base,
                    reproduction_round=0,
                    birth_count=3,
                    recruitment_count=1,
                    actual_consumption=3,
                )
            ]
        ),
    }
    data = lineage_plot_tables(tables)
    o = data["opportunity"].iloc[0]
    assert np.isclose(o.p_local, 0.65) and o.coverage == 0.5
    a = data["ancestry"].set_index("source_group")
    assert np.isclose(a.loc["P1", "fraction"], 5 / 6)
    assert np.isclose(a.fraction.sum(), 1)
    assert data["network"].actual_consumption.sum() == 3
    tables["matched_opportunities"]["matched_weight"] = 0
    tables["matched_opportunities"][["p_local", "p_global"]] = np.nan
    o = lineage_plot_tables(tables)["opportunity"].iloc[0]
    assert np.isnan(o.p_local) and np.isnan(o.ri_relative) and o.coverage == 0


def test_plot_coordinates_and_extinction():
    p1, p2, h = cell_abundance(
        dict(
            plant_cells=np.array([[2, 7]]),
            plant_species=np.array([0]),
            insect_cells=np.array([[3, 8]]),
        )
    )
    assert p1[2, 7] == 1 and h[3, 8] == 1 and p2.sum() == 0
    data = pd.DataFrame(
        [
            dict(species="H1", case="x", replicate=1, generation=g, trait=t, mean=v)
            for g, v in [(0, 0.2), (1, 0.3), (2, np.nan)]
            for t in ["I_1", "I_2"]
        ]
    )
    g, x, y = centered_xy(data, "H1", "x", 1, ("I_1", "I_2"))
    assert list(g) == [0, 1]
    np.testing.assert_allclose(x, [0, 0.1])


def test_summary_batch_has_no_figure_side_effect():
    import json
    from pathlib import Path
    from tempfile import TemporaryDirectory

    from evo_sim.models.plant_insect.analysis.common import Batch, batch_parser

    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "data").mkdir()
        (root / "settings.json").write_text(json.dumps({"version": "0.1.12"}))
        (root / "audits").mkdir()
        (root / "audits/integrity_audit.json").write_text(json.dumps({"runs": 1}))
        (root / "data/run_status.csv").write_text("case,replicate,status\nx,1,completed\n")
        batch = Batch.open(root)
        assert batch.summary_dir.is_dir()
        assert not (root / "figures").exists()
        assert not hasattr(batch_parser("summary").parse_args(["--input", tmp]), "out")


def test_v016_plot_tables_include_terminal_adults_and_larval_resource():
    from pathlib import Path
    from tempfile import TemporaryDirectory
    from unittest.mock import patch

    from evo_sim.models.plant_insect.analysis.common import Batch
    from evo_sim.models.plant_insect.analysis.plots.diagnostics import figure_C
    from evo_sim.models.plant_insect.analysis.summaries.lifecycle import prepare_diagnostics
    from evo_sim.models.plant_insect.audit.run import audit
    from evo_sim.models.plant_insect.batch import run
    from evo_sim.models.plant_insect.config.schema import config_payload
    from tests.support import config

    c = config(generations=4)
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "batch"
        run(
            dict(
                base=config_payload(c),
                cases=[
                    dict(name="mutation_off", overrides={}),
                    dict(name="qtl_mutating", overrides={}),
                ],
                replicates=1,
            ),
            root,
        )
        checked = audit(root)
        summary = root / "data/summary"
        summary.mkdir()
        batch = Batch(root, summary, checked, pd.read_csv(root / "data/run_status.csv"))
        prepare_diagnostics(batch)
        run_summary = pd.read_csv(summary / "run_summary.csv")
        for species in ("P1", "P2", "H1"):
            assert (run_summary[species + "_peak"] >= run_summary[species + "_final"]).all()
        pop = pd.read_csv(summary / "population_timeseries.csv")
        assert (
            pop.generation.max() == 4
            and not pop.duplicated(["case", "replicate", "generation", "species"]).any()
        )
        for case in ("mutation_off", "qtl_mutating"):
            with np.load(root / "data" / case / "rep01/snapshots/g000004.npz") as state:
                last = pop[(pop.case == case) & (pop.generation == 4)]
                assert last.loc[last.species == "H1", "N_start"].iloc[0] == len(state["insect_ids"])
        data = {
            name: pd.read_csv(summary / f"{name}.csv")
            for name in (
                "population_timeseries",
                "species_extinction",
                "host_resource_by_generation",
            )
        }
        data["lifecycle_schema"] = "dose-tolerance-lifecycle-2"
        assert data["host_resource_by_generation"].generation.max() == 3
        captured = []
        with patch(
            "evo_sim.models.plant_insect.analysis.plots.diagnostics.export",
            side_effect=lambda fig, *args: captured.append(fig),
        ):
            figure_C(data, Path(tmp), 80)
        ticks = [t.get_text() for t in captured[0].axes[-1].get_xticklabels()]
        assert ticks == [
            "Adults",
            "Mated\nfemales",
            "Potential\neggs",
            "Larval\nbirths",
            "Adult\nrecruits",
        ]
        from tests.support import raises

        legacy = dict(data)
        legacy.pop("lifecycle_schema")
        raises(lambda: figure_C(legacy, Path(tmp), 80))
        import matplotlib.pyplot as plt

        for fig in captured:
            plt.close(fig)


def test_v017_module_panels_and_surviving_recruitment():
    from pathlib import Path
    from tempfile import TemporaryDirectory
    from unittest.mock import patch

    import matplotlib.pyplot as plt

    from evo_sim.models.plant_insect.analysis.plots import diagnostics as plot

    traits = (
        [f"plant_chemical_traits_{i}" for i in range(1, 6)]
        + [f"insect_module_traits_{i}" for i in range(1, 11)]
        + ["I_1", "I_2", "A_1", "A_2"]
    )
    frame = pd.DataFrame(
        dict(
            record_type="trait",
            case=case,
            replicate=rep,
            species=sp,
            generation=g,
            trait=t,
            mean=i * g / 1000,
        )
        for case in plot.CASES
        for rep in plot.REPS
        for sp in ("P1", "P2", "H1")
        for g in (0, 200)
        for i, t in enumerate(traits, 1)
    )
    pop = pd.DataFrame(
        dict(
            case=case,
            replicate=rep,
            species="H1",
            generation=g,
            N_start=400,
            N_mothers_with_mates=100,
            N_potential_eggs=500,
            N_offspring=450,
            N_next=400,
        )
        for case in plot.CASES
        for rep in plot.REPS
        for g in range(201)
    )
    ext = pop.loc[pop.generation.eq(200), ["case", "replicate", "species"]].assign(
        extinction_observed=False, event_or_censor_generation=200
    )
    resource = pop.loc[pop.generation.lt(200), ["case", "replicate", "generation"]].assign(
        total_actual_resource=200, P2_resource_fraction=0.5
    )
    captured = []
    with (
        TemporaryDirectory() as tmp,
        patch.object(
            plot, "export", side_effect=lambda fig, out, stem, dpi: captured.append((fig, stem))
        ),
    ):
        plot.figure_B(dict(trait_timeseries=frame), Path(tmp), 60)
        plot.figure_C(
            dict(
                population_timeseries=pop,
                species_extinction=ext,
                host_resource_by_generation=resource,
                lifecycle_schema="dose-tolerance-lifecycle-2",
            ),
            Path(tmp),
            60,
        )
        plot.figure_A(dict(population_timeseries=pop, species_extinction=ext), Path(tmp), 60)
    fig = captured[0][0]
    assert len(fig.axes) == 14
    assert all(not ax.texts for ax in fig.axes)
    assert [ax.get_xlabel() for ax in fig.axes[::2]][2:] == [
        "Δ met E₁ (model units)",
        "Δ met μ₁ (model units)",
        "Δ sense E₁ (model units)",
        "Δ sense μ₁ (model units)",
        "Δ G1 raw (model units)",
    ]
    # Endpoint diamonds must follow the eight inherited module axes, never derived I/A.
    for ax, index in zip(fig.axes[4::2], (6, 8, 10, 12, 14), strict=True):
        diamonds = [line for line in ax.lines if line.get_marker() == "D"]
        assert diamonds
        np.testing.assert_allclose(diamonds[0].get_xdata(), [index * 0.2])
        np.testing.assert_allclose(diamonds[0].get_ydata(), [(index + 1) * 0.2])
    recruitment, stem = captured[1]
    assert stem == "C_host_use_recruitment"
    assert "No H1 extinction" in recruitment.axes[2].get_title()
    assert recruitment.axes[2].get_xlim() == (180, 200)
    assert len(recruitment.axes[2].lines) == 6
    assert len(recruitment.axes[3].lines) == 6
    for ax in recruitment.axes:
        assert [line.get_linestyle() for line in ax.lines[:3]] == ["-", "--", ":"]
    np.testing.assert_array_equal(
        recruitment.axes[3].lines[0].get_ydata(), [400, 100, 500, 450, 400]
    )
    assert "H1 extinction" not in [t.get_text() for t in captured[2][0].legends[0].get_texts()]
    for fig, _ in captured:
        plt.close(fig)


def test_current_batch_rejects_historical_protocols_and_report_fallback():
    import json
    from pathlib import Path
    from tempfile import TemporaryDirectory

    from evo_sim.models.plant_insect.analysis.common import Batch, audit_file, read_audit
    from tests.support import raises

    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        summary = root / "data/summary"
        summary.mkdir(parents=True)
        (root / "settings.json").write_text(json.dumps({"version": "0.1.6"}))
        (summary / "corrected_integrity_audit.json").write_text(json.dumps({"runs": 1}))
        raises(lambda: Batch.open(root))
        assert audit_file(root) == root / "audits/integrity_audit.json"
        raises(lambda: read_audit(root), FileNotFoundError)


def test_sensitivity_extinction_and_incomplete_bounds():
    from evo_sim.models.plant_insect.analysis.summaries.lifecycle import sensitivity_survival

    runs = pd.DataFrame(
        [
            dict(case="x", audited=True, last_generation=2, extinction_generation=2),
            dict(case="x", audited=True, last_generation=5, extinction_generation=np.nan),
            dict(case="x", audited=True, last_generation=1, extinction_generation=np.nan),
            dict(case="x", audited=False, last_generation=np.nan, extinction_generation=0),
        ]
    )
    at_two = sensitivity_survival(runs, 5).query("generation == 2").iloc[0]
    assert (at_two.alive, at_two.extinct, at_two.unknown) == (1, 1, 2)
    assert at_two.lower == 0.25 and at_two.upper == 0.75
    assert np.isnan(at_two.wilson_low)
    all_known = sensitivity_survival(runs.iloc[:2], 5).query("generation == 5").iloc[0]
    assert (all_known.alive, all_known.extinct, all_known.unknown) == (1, 1, 0)
    assert all_known.wilson_low < 0.5 < all_known.wilson_high
