"""New analysis consumes audited records without altering simulation or historical ledgers."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect import simulate
from evo_sim.models.plant_insect.analysis.summaries.allele_change import (
    COMPONENTS,
    decompose_run,
)
from evo_sim.models.plant_insect.audit.run import audit_run
from evo_sim.models.plant_insect.batch import save_result
from tests.support import config, raises


def test_saved_frequency_accounting_matches_births_recruits_and_actual_next_state():
    c = config(plant_sizes=[4, 4], insect_sizes=[12], generations=2)
    c.spatial.update(initial_plant_cells=[[0, 0]] * 8, initial_insect_cells=[[0, 0]] * 12)
    for locus in c.genetics["loci"]:
        locus["mutation_rate"] = 1.0
    result = simulate(c)
    assert result["metadata"]["status"]["status"] == "completed"
    with TemporaryDirectory() as temporary:
        path = Path(temporary) / "run"
        save_result(result, path)
        assert audit_run(path)["genetic_audit"]
        changes, hashes = decompose_run(path)
        assert set(changes.scope) == {"all", "type"}
        active = changes[changes.N_next.gt(0) & changes.N_parent.gt(0)]
        np.testing.assert_allclose(active[COMPONENTS].sum(axis=1), active.delta_total, atol=1e-12)
        assert changes.delta_mutation.abs().sum() > 0
        legacy = pd.DataFrame(result["allele_contributions"])
        for row in changes.query("scope == 'type'").itertuples():
            selection = legacy[
                legacy.generation.eq(row.generation)
                & legacy.kind.eq(row.kind)
                & legacy.species.eq(row.species)
                & legacy.locus_id.eq(row.locus_id)
                & legacy.allele_id.eq(row.allele_id)
                & legacy.stage.eq("birth")
            ]
            if not selection.empty:
                old = selection.iloc[0]
                np.testing.assert_allclose(
                    [row.p_parent, row.p_parent_weighted, row.p_birth],
                    [old.p_parent, old.p_parent_weighted, old.p_offspring],
                    equal_nan=True,
                )
        for row in changes.query("scope == 'all'").itertuples():
            snapshot = result["snapshots"][row.offspring_generation]
            names = [l["locus_id"] for l in result["genetic_loci"] if l["organism"] == row.kind]
            g = snapshot[row.kind + "_genotypes"][:, :, names.index(row.locus_id)]
            np.testing.assert_allclose(
                row.p_next, (g == row.allele_id).mean() if g.size else np.nan, equal_nan=True
            )
        again, same_hashes = decompose_run(path)
        pd.testing.assert_frame_equal(changes, again)
        assert hashes == same_hashes
        file = path / "gamete_transmissions.csv"
        original = file.read_bytes()
        file.write_bytes(file.read_bytes() + b"\n")
        raises(lambda: decompose_run(path))
        file.write_bytes(original)
        snapshot = path / "snapshots/g000001.npz"
        original = snapshot.read_bytes()
        snapshot.unlink()
        raises(lambda: decompose_run(path))
        snapshot.write_bytes(original)


def test_zero_round_and_extinct_runs_do_not_fabricate_frequency_change():
    with TemporaryDirectory() as temporary:
        for name, changes in (
            ("zero", dict(generations=0)),
            ("extinction", dict(generations=2, plant_sizes=[0, 0], insect_sizes=[1])),
        ):
            result = simulate(config(**changes))
            path = Path(temporary) / name
            save_result(result, path)
            table, _ = decompose_run(path)
            if name == "zero":
                assert table.empty
            else:
                assert not table.empty and table.N_next.eq(0).all()
                assert table[COMPONENTS + ["delta_total"]].isna().all().all()
            meta = json.loads((path / "metadata.json").read_text())
            assert meta["version"] == "0.1.12"
