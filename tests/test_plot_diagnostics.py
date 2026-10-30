import numpy as np
import pandas as pd

from evo_sim.models.plant_insect.analysis.plot_diagnostics import cell_abundance, centered_xy


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
    from evo_sim.models.plant_insect.analysis.plot_diagnostics import figure_C
    from evo_sim.models.plant_insect.analysis.summarize_lifecycle import prepare_diagnostics
    from evo_sim.models.plant_insect.audit import audit
    from evo_sim.models.plant_insect.batch import run
    from evo_sim.models.plant_insect.helpers.config import config_payload
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
        data["lifecycle_schema"] = "adult-oviposition-larval-recruitment-1"
        assert data["host_resource_by_generation"].generation.max() == 3
        captured = []
        with patch(
            "evo_sim.models.plant_insect.analysis.plot_diagnostics.export",
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
        legacy = dict(data)
        legacy.pop("lifecycle_schema")
        with patch(
            "evo_sim.models.plant_insect.analysis.plot_diagnostics.export",
            side_effect=lambda fig, *args: captured.append(fig),
        ):
            figure_C(legacy, Path(tmp), 80)
        assert [t.get_text() for t in captured[1].axes[-1].get_xticklabels()] == [
            "Starting",
            "Surviving",
            "Mated\nmothers",
            "Offspring",
        ]
        import matplotlib.pyplot as plt

        for fig in captured:
            plt.close(fig)


def test_v017_module_panels_and_surviving_recruitment():
    from pathlib import Path
    from tempfile import TemporaryDirectory
    from unittest.mock import patch

    import matplotlib.pyplot as plt

    from evo_sim.models.plant_insect.analysis import plot_diagnostics as plot

    traits = (
        [f"plant_chemical_traits_{i}" for i in range(1, 4)]
        + [f"insect_module_traits_{i}" for i in range(1, 9)]
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
                lifecycle_schema="adult-oviposition-larval-recruitment-1",
            ),
            Path(tmp),
            60,
        )
        plot.figure_A(dict(population_timeseries=pop, species_extinction=ext), Path(tmp), 60)
    fig = captured[0][0]
    assert len(fig.axes) == 12
    assert all(not ax.texts for ax in fig.axes)
    assert [ax.get_xlabel() for ax in fig.axes[::2]][2:] == [
        "Δ met E₁ (model units)",
        "Δ met μ₁ (model units)",
        "Δ sense E₁ (model units)",
        "Δ sense μ₁ (model units)",
    ]
    # Endpoint diamonds must follow the eight inherited module axes, never derived I/A.
    for ax, index in zip(fig.axes[4::2], (4, 6, 8, 10), strict=True):
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


def test_corrected_audit_checks_evidence_and_run_identity():
    import json
    from pathlib import Path
    from tempfile import TemporaryDirectory

    from evo_sim.models.plant_insect.analysis.common import audit_file, read_audit, sha256
    from tests.support import raises

    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        summary = root / "data/summary"
        summary.mkdir(parents=True)
        (root / "settings.json").write_text("{}")
        status = root / "data/run_status.csv"
        status.write_text("case,replicate,seed,generation,status\nx,1,123,200,completed\n")
        folder = root / "data/x/rep01"
        folder.mkdir(parents=True)
        (folder / "metadata.json").write_text("{}")
        row = dict(
            case="x",
            replicate=1,
            seed=123,
            original_files_unchanged=True,
            frozen_sources_match=True,
            audit=dict(
                status="completed",
                genetic_audit=True,
                contribution_audit=True,
                states=201,
                transitions=200,
                contribution_transitions=200,
            ),
        )
        source = summary / "corrected_run_audits_x.json"
        source.write_text(json.dumps(dict(runs=[row])))
        design = summary / "corrected_design_consistency.json"
        design.write_text(
            json.dumps(
                dict(
                    status="passed",
                    settings_sha256=sha256(root / "settings.json"),
                    run_status_sha256=sha256(status),
                    runs=[
                        dict(
                            case="x",
                            replicate=1,
                            metadata_sha256=sha256(folder / "metadata.json"),
                            config_matches_formal=True,
                            source_matches_frozen=True,
                            input_sources_match_settings=True,
                            run_status_matches_metadata=True,
                        )
                    ],
                )
            )
        )
        report = dict(
            status="passed",
            expected_runs=1,
            audited_runs=1,
            runs=[row],
            source_reports={source.name: sha256(source)},
            design_report_sha256=sha256(design),
        )
        target = summary / "corrected_integrity_audit.json"
        target.write_text(json.dumps(report))
        assert audit_file(root) == target
        checked = read_audit(root)
        assert checked["runs"] == checked["genetically_audited_runs"] == 1
        assert checked["run_reports"] == [row]
        (folder / "metadata.json").write_text('{"changed": true}')
        raises(lambda: read_audit(root), AssertionError)
        (folder / "metadata.json").write_text("{}")
        source.write_text(json.dumps(dict(runs=[])))
        raises(lambda: read_audit(root), AssertionError)
        source.write_text(json.dumps(dict(runs=[row])))
        row["seed"] = 124
        target.write_text(json.dumps(report))
        raises(lambda: read_audit(root), AssertionError)
