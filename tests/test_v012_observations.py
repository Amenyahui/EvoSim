"""Scientific ledger regressions: joint responses, failed hosts and committed ancestry."""

import copy
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect import simulate
from evo_sim.models.plant_insect.audit import audit_run
from evo_sim.models.plant_insect.batch import save_result, write_json
from evo_sim.models.plant_insect.helpers.observations import standardized
from tests.support import config, raises


def test_joint_response_preserves_individual_associations():
    c = config(insect_sizes=[2])
    B, weights = c.validate()
    insects = SimpleNamespace(
        preference=np.array([[0.6, 0.6], [0.2, 0.2]]),
        resistance=np.array([[0.6, 0.2], [0.2, 0.6]]),
        species=np.zeros(2, int),
    )
    # The production population exposes its size through __len__.
    from evo_sim.models.plant_insect.helpers.agents import InsectPopulation

    state = SimpleNamespace(
        insects=InsectPopulation(
            ids=np.arange(2),
            species=insects.species,
            resistance=insects.resistance,
            preference=insects.preference,
        )
    )
    rows, p, q = standardized(state, c, B, weights, 0)
    assert p.shape == q.shape == (2, 2)
    for j, row in enumerate(rows):
        expected = sum(p[i, j] * q[i, j] for i in range(2)) / 2
        assert np.isclose(row["joint_pq"], expected)
        assert not np.isclose(expected, row["preference"] * row["q_unlimited"])
        assert np.isclose(row["cov_pq"], np.cov(p[:, j], q[:, j], ddof=0)[0, 1])
    state.insects = InsectPopulation(
        ids=np.empty(0, int),
        species=np.empty(0, int),
        resistance=np.empty((0, 2)),
        preference=np.empty((0, 2)),
    )
    rows, p, q = standardized(state, c, B, weights, 1)
    assert p.shape == q.shape == (0, 2)
    assert all(np.isnan(row["joint_pq"]) and row["insect_count"] == 0 for row in rows)


def test_unmatched_female_is_zero_not_excluded_from_maternal_mean():
    c = config(plant_sizes=[0, 0], insect_sizes=[8], generations=1)
    c.sexual_reproduction["insect_founder_female_probability"] = 1.0
    result = simulate(c)
    females = [r for r in result["parents"] if r["generation"] == 0]
    assert len(females) == 8
    assert all(r["is_maternal_role"] and r["reason"] == "no_host" for r in females)
    assert all(r["expected_maternal_offspring"] == 0 for r in females)
    assert all(r["maternal_offspring_variance"] == 0 for r in females)
    assert all(r["realized_genetic_contribution"] == 0 for r in females)
    assert all(r["contribution_committed"] for r in females)
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "no_host"
        save_result(result, path)
        assert audit_run(path)["contribution_transitions"] == 1


def test_contribution_audit_boundaries_and_tampering():
    c = config(
        generations=1,
        insect_sizes=[12],
        choice_mode="forced",
        sP0=1.0,
        resource_capacity={"P1": 100.0, "P2": 100.0},
    )
    c.spatial["initial_plant_cells"] = [[0, 0]] * 40
    c.spatial["initial_insect_cells"] = [[0, 0]] * 12
    for locus in c.genetics["loci"]:
        locus["mutation_rate"] = 1.0
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "run"
        result = simulate(c)
        save_result(result, path)
        assert audit_run(path)["contribution_transitions"] == 1
        assert any(r["kind"] == "insect" for r in result["offspring"])
        assert any(
            r["p_parent"] == 0 and r["p_offspring"] > 0 for r in result["allele_contributions"]
        )
        for table, column in (
            ("parents", "realized_genetic_contribution"),
            ("networks", "genetic_contribution"),
            ("allele_contributions", "p_parent_weighted"),
            ("diagnostics", "joint_pq"),
        ):
            file = path / (table + ".csv")
            original = file.read_bytes()
            metadata = json.loads((path / "metadata.json").read_text())
            frame = pd.read_csv(file)
            frame.loc[0, column] += 0.25
            frame.to_csv(file, index=False)
            metadata["output_sha256"][file.name] = hashlib.sha256(file.read_bytes()).hexdigest()
            write_json(path / "metadata.json", metadata)
            raises(lambda: audit_run(path), AssertionError)
            file.write_bytes(original)
            metadata["output_sha256"][file.name] = hashlib.sha256(original).hexdigest()
            write_json(path / "metadata.json", metadata)
        assert audit_run(path)["contribution_audit"]
        for label, changes in (
            ("terminal", {"generations": 0}),
            ("no_observations", {"record_observations": False}),
            ("initial_limit", {"max_individuals_safety": 1}),
            ("extinction", {"plant_sizes": [1, 0], "insect_sizes": [0]}),
        ):
            case = config(**changes)
            r = simulate(case)
            folder = Path(tmp) / label
            save_result(r, folder)
            assert audit_run(folder)["contribution_audit"]
            for row in r["parents"]:
                if not row["contribution_committed"]:
                    assert row["realized_genetic_contribution"] is None
        off = copy.deepcopy(c)
        off.record_observations = False
        unseen = simulate(off)
        for a, b in zip(result["snapshots"], unseen["snapshots"]):
            for field in ("plant_genotypes", "insect_genotypes", "D", "I", "A"):
                np.testing.assert_array_equal(a[field], b[field])
        assert result["offspring"] == unseen["offspring"]
