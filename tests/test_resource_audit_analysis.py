"""Independent plant resource, inheritance, and analysis checks."""

from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect import get_runner
from evo_sim.models.plant_insect.analysis.common import Batch
from evo_sim.models.plant_insect.analysis.summarize_lifecycle import (
    summarize,
    summarize_sensitivity,
)
from evo_sim.models.plant_insect.analysis.summarize_lineages import summarize as summarize_lineages
from evo_sim.models.plant_insect.audit import audit, audit_run
from evo_sim.models.plant_insect.batch import run
from evo_sim.models.plant_insect.helpers.config import config_payload
from tests.support import config, raises
from tests.test_v016_audit import reject_rehashed_tampering


def test_sexual_plants_are_audited_and_summarized():
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "batch"
        run(
            dict(
                base=dict(get_runner().config["parameters"], generations=2),
                cases=[dict(name="current", overrides={})],
                replicates=1,
            ),
            root,
        )
        folder = root / "data/current/rep01"
        assert audit(root, min_transitions=2)["acceptance_passed"]
        assert audit_run(folder)["genetic_audit"]
        births = pd.read_csv(folder / "offspring.csv").query("kind == 'plant'")
        assert len(births) and births.paternal_id.notna().all()
        assert births.maternal_id.ne(births.paternal_id).all()
        tx = pd.read_csv(folder / "gamete_transmissions.csv").query("kind == 'plant'")
        assert set(tx.role) == {"maternal", "paternal"} and set(tx["copy"]) == {0, 1}
        for table, field in (("plants", "defoliation"), ("parents", "expected_propagules")):
            reject_rehashed_tampering(
                folder,
                table,
                lambda frame, column=field: frame.__setitem__(column, frame[column] + 1),
            )
        raw = summarize(root)
        assert raw["plant_resource_cycles"].resource_capacity.gt(0).all()
        tables, _ = summarize_sensitivity(root)
        assert tables["diagnostics"].plant_reproduction_active.all()
        lineages = summarize_lineages(
            Batch.open(root),
            dict(distance_threshold=1.0, minimum_cluster_size=2, persistence_generations=2),
        )
        links = lineages["transitions"].query("kind == 'plant'")
        assert links.recruited_genetic_contribution.sum() == len(births)
        assert lineages["mating"].kind.eq("plant").any()


def test_audit_replays_cross_label_plant_parents():
    c = config(plant_sizes=[1, 1], insect_sizes=[0], generations=2)
    c.plant_demography.update(
        background_mortality_rate=0, defoliation_mortality_rate=0, microsites_per_cell=6
    )
    c.plant_offspring_lambda = {"P1": 100.0, "P2": 100.0}
    c.spatial.update(
        grid_shape=[1, 1],
        initial_plant_cells=[[0, 0], [0, 0]],
        plant_dispersal_radius=0,
        search_radius=0,
        insect_dispersal_radius=0,
    )
    c.sexual_reproduction["insect_mate_radius"] = 0
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "batch"
        run(
            dict(base=config_payload(c), cases=[dict(name="cross", overrides={})], replicates=1),
            root,
        )
        folder = root / "data/cross/rep01"
        assert audit_run(folder)["contribution_audit"]
        births = pd.read_csv(folder / "offspring.csv")
        first = births[births.generation.eq(1)]
        assert len(first) and first.maternal_species.ne(first.paternal_species).all()
        reject_rehashed_tampering(
            folder, "offspring", lambda f: f.__setitem__("paternal_id", f.maternal_id)
        )


def test_r30_resource_residual_replays_shared_budget_roundoff_without_relaxing_checks():
    from evo_sim.models.plant_insect.audit_lifecycle import check_rows, resource_shares

    # Actual R=30/rep01, round 7, plant 52: 25 larvae previously exposed an omitted
    # common-factor roundoff repair in the independent audit, not a biological error.
    demand = np.array(
        [
            0.9960833980957101,
            0.9599287362493887,
            1.0,
            1.0,
            1.0,
            1.0,
            0.9607031891431552,
            1.0,
            1.0,
            0.9956916870823135,
            0.9596623730298408,
            0.9994711510511312,
            0.8819639142548346,
            1.0,
            0.921381845502621,
            1.0,
            1.0,
            1.0,
            0.9564566560185253,
            0.9410962942690759,
            1.0,
            0.9970442219752372,
            1.0,
            0.9511782936242251,
            1.0,
        ]
    )
    stock = np.array([17.242060111183612])
    _, gamma, shares, consumed = resource_shares(np.zeros(len(demand), int), demand, stock)
    residual = stock[0] - consumed[0]
    assert residual == 1.0658141036401503e-14
    assert consumed[0] <= stock[0]
    np.testing.assert_allclose(shares / demand, gamma[0], rtol=1e-15, atol=0)
    frame = pd.DataFrame(dict(plant_id=[52], resource_after_feeding=[residual]))
    expected = {(52,): dict(resource_after_feeding=residual)}
    check_rows(frame, ("plant_id",), expected, "plants")
    frame.loc[0, "resource_after_feeding"] += 1e-6
    raises(lambda: check_rows(frame, ("plant_id",), expected, "plants"), AssertionError)
    for values in (np.array([]), np.array([0.0, 0.0]), np.array([0.5, 1.0])):
        _, _, allocated, used = resource_shares(np.zeros(len(values), int), values, np.array([0.0]))
        assert not allocated.any() and used.tolist() == [0.0]
