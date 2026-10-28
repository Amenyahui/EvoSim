"""Saved-data adversarial checks for the independent v0.1.6 auditor."""

import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect import simulate
from evo_sim.models.plant_insect.audit import audit_run
from evo_sim.models.plant_insect.batch import save_result, write_json
from evo_sim.models.plant_insect.helpers.genetics import GeneticState
from tests.support import config, raises


def fixture(**changes):
    values = dict(
        generations=2,
        plant_sizes=[4, 4],
        insect_sizes=[12],
        choice_mode="forced",
        resource_capacity={"P1": 0.5, "P2": 0.5},
    )
    values.update(changes)
    c = config(**values)
    c.spatial.update(
        initial_plant_cells=[[0, 0]] * sum(c.plant_sizes) if sum(c.plant_sizes) else None,
        initial_insect_cells=[[0, 0]] * sum(c.insect_sizes) if sum(c.insect_sizes) else None,
    )
    return c


def reject_rehashed_tampering(path, table, corrupt):
    """A hash-valid file must still obey the biological event contract."""
    file = path / f"{table}.csv"
    original, original_meta = file.read_bytes(), (path / "metadata.json").read_bytes()
    frame = pd.read_csv(file, float_precision="round_trip")
    corrupt(frame)
    frame.to_csv(file, index=False)
    metadata = json.loads(original_meta)
    metadata["output_sha256"][file.name] = hashlib.sha256(file.read_bytes()).hexdigest()
    write_json(path / "metadata.json", metadata)
    try:
        raises(lambda: audit_run(path), (AssertionError, KeyError, ValueError))
    finally:
        file.write_bytes(original)
        (path / "metadata.json").write_bytes(original_meta)


def test_v016_audit_dead_larvae_and_rehashed_tampering():
    result = simulate(fixture())
    assert result["metadata"]["status"]["status"] == "completed"
    assert any(not r["survived_to_adult"] for r in result["offspring"] if r["kind"] == "insect")
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "run"
        save_result(result, path)
        assert audit_run(path)["genetic_audit"]
        for table, column in (
            ("parents", "realized_recruitment_contribution"),
            ("events", "x"),
            ("larvae", "g"),
            ("host_flows", "genetic_contribution"),
            ("networks", "successful_recruits"),
            ("allele_contributions", "p_parent_weighted"),
            ("genetic_summary", "frequency"),
        ):
            reject_rehashed_tampering(
                path, table, lambda f, col=column: f.__setitem__(col, f[col] + 0.25)
            )

        def corrupt_dead_transmission(frame):
            dead = next(
                r
                for r in result["offspring"]
                if r["kind"] == "insect" and not r["survived_to_adult"]
            )
            selected = (
                (frame.kind == "insect")
                & (frame.generation == dead["generation"])
                & (frame.offspring_id == dead["offspring_id"])
            )
            index = frame[selected].index[0]
            frame.loc[index, "inherited_allele_id"] += 10000

        reject_rehashed_tampering(path, "gamete_transmissions", corrupt_dead_transmission)
        reject_rehashed_tampering(
            path, "offspring", lambda f: f.__setitem__("paternal_id", f.paternal_id + 10000)
        )
        reject_rehashed_tampering(
            path, "oviposition", lambda f: f.__setitem__("event_id", "duplicate")
        )
        assert audit_run(path)["contribution_audit"]


def test_v016_audit_empty_stages_terminal_and_failure():
    cases = (
        fixture(generations=0),
        fixture(choice_mode="outside", outside_utility=1e6),
        fixture(resource_capacity={"P1": 0.0, "P2": 0.0}),
        fixture(insect_sizes=[0]),
        fixture(plant_sizes=[0, 0]),
        fixture(plant_sizes=[0, 0], insect_sizes=[0]),
    )
    with TemporaryDirectory() as tmp:
        for i, c in enumerate(cases):
            result = simulate(c)
            assert result["metadata"]["status"]["status"] in ("completed", "all_species_extinct")
            if i == 1:
                assert result["oviposition"]
                assert all(r["kind"] == "plant" for r in result["offspring"])
                assert not any(r["accepted"] for r in result["oviposition"])
            if i == 2:
                assert result["larvae"] and not any(
                    r["survived_to_adult"] for r in result["larvae"]
                )
            path = Path(tmp) / str(i)
            save_result(result, path)
            checked = audit_run(path)
            assert checked["contribution_audit"]
            assert not any(
                r["generation"] == len(result["snapshots"]) - 1 for r in result["events"]
            )
        original = GeneticState.fertilize

        def fail_after_larvae(self, *args, **kwargs):
            if args[5] == "plant":
                raise MemoryError("injected after larval survival")
            return original(self, *args, **kwargs)

        with patch.object(GeneticState, "fertilize", fail_after_larvae):
            result = simulate(fixture())
        assert result["metadata"]["status"]["status"] == "computational_limit"
        assert (
            not result["offspring"] and not result["events"] and not result["gamete_transmissions"]
        )
        path = Path(tmp) / "failed"
        save_result(result, path)
        assert audit_run(path)["transitions"] == 0


def test_v016_audit_choice_modes_and_snapshot_exports():
    with TemporaryDirectory() as tmp:
        for i, changes in enumerate(
            (
                dict(shuffle_voc=True),
                dict(tau=0.0),
                dict(choice_mode="outside", tau=0.5),
                dict(olfaction_mode="neutral"),
            )
        ):
            result = simulate(fixture(**changes))
            path = Path(tmp) / str(i)
            save_result(result, path)
            assert audit_run(path)["genetic_audit"]
            with np.load(path / "snapshots/g000000.npz") as data:
                assert (data["insect_birth_host_species"] == -1).all()
