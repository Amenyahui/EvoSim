"""Independent perennial identity/resource audit and read-only cohort summaries."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect import Config, get_runner
from evo_sim.models.plant_insect.analysis.common import Batch
from evo_sim.models.plant_insect.analysis.summarize_lifecycle import (
    summarize,
    summarize_sensitivity,
)
from evo_sim.models.plant_insect.analysis.summarize_lineages import summarize as summarize_lineages
from evo_sim.models.plant_insect.audit import audit, audit_run
from evo_sim.models.plant_insect.batch import run, write_json
from evo_sim.models.plant_insect.helpers.config import config_payload
from tests.support import raises
from tests.test_v016_audit import reject_rehashed_tampering


def perennial_config():
    c = Config(**dict(get_runner().config["parameters"], generations=6))
    c.plant_sizes, c.insect_sizes = [4], [0]
    c.plant_offspring_lambda = {"P1": 20.0}
    c.plant_demography.update(
        background_mortality_rate=0.0,
        defoliation_mortality_rate=0.0,
        reproduction_interval=5,
        microsites_per_cell=1,
    )
    c.spatial.update(
        grid_shape=[3, 3],
        search_radius=1,
        plant_dispersal_radius=1,
        insect_dispersal_radius=1,
        initial_plant_cells=[[0, 0], [0, 1], [1, 0], [1, 1]],
        initial_insect_cells=None,
    )
    c.sexual_reproduction.update(pollen_radius=2, insect_mate_radius=2)
    return c


def test_perennial_saved_audit_rejects_rehashed_identity_stock_and_lottery_tampering():
    from evo_sim.models.plant_insect.analysis.common import sha256

    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "batch"
        c = perennial_config()
        run(
            dict(
                base=config_payload(c), cases=[dict(name="vacancies", overrides={})], replicates=1
            ),
            root,
        )
        checked = audit(root, min_transitions=6)
        assert checked["acceptance_passed"]
        folder = root / "data/vacancies/rep01"
        assert audit_run(folder)["transitions"] == 6
        proposals = pd.read_csv(folder / "plant_proposals.csv")
        assert proposals.reproduction_round.eq(4).all()
        assert len(proposals) > 5 and proposals.established.sum() == 5
        births = pd.read_csv(folder / "offspring.csv")
        assert len(births) == 5 and births.kind.eq("plant").all()
        with np.load(folder / "snapshots/g000006.npz") as snapshot:
            assert snapshot["plant_ids"][:4].tolist() == [0, 1, 2, 3]
            assert len(snapshot["plant_ids"]) == len(np.unique(snapshot["plant_microsite"])) == 9
            np.testing.assert_allclose(snapshot["plant_resource"][4:], 6.25)
        parents = pd.read_csv(folder / "parents.csv")
        assert parents.expected_maternal_offspring.isna().all()
        assert parents.maternal_offspring_variance.isna().all()
        assert parents.loc[parents.generation.ne(4), "expected_propagules"].eq(0).all()
        assert parents.loc[parents.generation.eq(4), "expected_propagules"].sum() > len(births)
        for table, field in (
            ("plants", "resource_after_feeding"),
            ("plants", "resource_next"),
            ("plant_proposals", "microsite"),
            ("parents", "expected_propagules"),
        ):
            reject_rehashed_tampering(
                folder,
                table,
                lambda frame, column=field: frame.__setitem__(column, frame[column] + 1),
            )
        file = folder / "snapshots/g000001.npz"
        original, original_meta = file.read_bytes(), (folder / "metadata.json").read_bytes()
        try:
            for field in ("plant_ids", "plant_resource"):
                with np.load(file) as saved:
                    state = {key: saved[key] for key in saved.files}
                if field == "plant_ids":
                    state[field][0] += 100000
                else:
                    state[field][0] -= 0.25
                    state["R"][0] = state[field][0]
                np.savez_compressed(file, **state)
                metadata = json.loads(original_meta)
                metadata["output_sha256"]["snapshots/g000001.npz"] = sha256(file)
                write_json(folder / "metadata.json", metadata)
                raises(lambda: audit_run(folder), (AssertionError, ValueError, KeyError))
                file.write_bytes(original)
        finally:
            file.write_bytes(original)
            (folder / "metadata.json").write_bytes(original_meta)
        assert audit_run(folder)["genetic_audit"]
        raw = summarize(root)
        assert len(raw["plant_establishment_proposals"]) == len(proposals)
        assert raw["plant_resource_cycles"].resource_capacity.gt(0).all()
        tables, _ = summarize_sensitivity(root)
        diagnostics = tables["diagnostics"]
        before = diagnostics[diagnostics.generation.lt(4)]
        assert (
            before.plant_newborns.eq(0).all() and before.plant_establishment_fraction.isna().all()
        )
        pulse = diagnostics[diagnostics.generation.eq(4)].iloc[0]
        assert pulse.plant_retained == 4 and pulse.plant_newborns == 5
        assert pulse.plant_propagules > pulse.plant_newborns
        lineages = summarize_lineages(
            Batch.open(root),
            dict(distance_threshold=1.0, minimum_cluster_size=2, persistence_generations=2),
        )
        links = lineages["transitions"].query("kind == 'plant'")
        assert links.retained_individuals.sum() == 29
        assert links.recruited_genetic_contribution.sum() == 5
        assert not lineages["lineage_events"].query("kind == 'plant'").event.eq("termination").any()


def test_retained_plants_without_any_births_accept_empty_event_ledgers():
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "batch"
        c = perennial_config()
        c.generations = 2
        run(
            dict(base=config_payload(c), cases=[dict(name="retained", overrides={})], replicates=1),
            root,
        )
        assert audit(root, min_transitions=2)["acceptance_passed"]
        frames = summarize_lineages(
            Batch.open(root),
            dict(distance_threshold=1.0, minimum_cluster_size=2, persistence_generations=2),
        )
        links = frames["transitions"].query("kind == 'plant'")
        assert links.retained_individuals.sum() == 8
        assert links.recruited_genetic_contribution.eq(0).all()
        assert frames["mating"].offspring_weighted_pairs.eq(0).all()


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
