"""Maternal-only cost boundaries, biparental transmission and independent reconstruction."""

import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import numpy as np

from evo_sim.models.plant_insect import get_runner, simulate
from evo_sim.models.plant_insect.audit import audit_run
from evo_sim.models.plant_insect.batch import run, save_result, write_json
from evo_sim.models.plant_insect.helpers.config import config_payload
from evo_sim.models.plant_insect.helpers.ecology import settle_feeding
from evo_sim.models.plant_insect.substeps.demography import SampleMaternallyProducedOffspring
from tests.support import config, raises


def test_v015_roles_zero_cost_multiplier_and_mate_boundaries():
    runner = get_runner()
    runner.config["parameters"].update(plant_sizes=[3, 1], insect_sizes=[4])
    runner.init()
    runtime = runner.state["environment"]["runtime"]
    runtime.update = True
    runtime.plant_phi = np.ones(4)
    runtime.c.sexual_reproduction.update(pollen_radius=0, insect_mate_radius=0)
    runner.state["agents"]["plants"]["cells"][:] = [[0, 0], [0, 0], [5, 5], [0, 0]]
    runner.state["agents"]["insects"]["sex"][:] = [0, 1, 0, 1]
    runtime.e = dict(
        S_plant=np.ones(4),
        fp=np.array([1.0, 0.0, 1.0, 1.0]),
        S_insect=np.array([1.0, 1.0, 0.0, 0.0]),
        fi=np.array([1.0, 0.0, 0.0, 0.0]),
        host=np.array([0, 0, -1, 0]),
    )
    # Fixed integer Poisson outcomes isolate role/mean rules from finite-sample noise.
    for kind in (0, 1):
        runtime.rng[3][kind][:] = [
            SimpleNamespace(poisson=lambda mu: mu.astype(int)) for _ in runtime.rng[3][kind]
        ]
        component = SampleMaternallyProducedOffspring(
            runner.config,
            {"runtime": "environment/runtime", "agents": "agents"},
            ["offspring"],
            {"fixed": {"taxon": kind}},
        )
        runtime.result["parents"] = []
        sampled = component.forward(runner.state, {})["offspring"]
        rows = runtime.result["parents"]
        assert sampled["candidates"][0].tolist() == [1]  # Surviving F_m=0 father.
        assert rows[1]["eligible_as_donor"] and rows[1]["F_i"] == 0
        assert sampled["offspring_count"][0] == (2 if kind == 0 else 4)
        assert not sampled["offspring_count"][1:].any()
        assert rows[0]["expected_maternal_offspring"] == (2 if kind == 0 else 4)
        if kind == 0:
            assert rows[0]["eligible_as_donor"]  # Same plant can serve both roles.
            assert rows[1]["eligible_as_mother"] and rows[1]["conditional_mu"] == 0
            assert not sampled["candidates"][2].size  # Outside mating radius.
            assert not sampled["candidates"][3].size  # Different species cannot substitute.
        else:
            assert not rows[0]["eligible_as_donor"] and not rows[1]["eligible_as_mother"]
            assert np.isnan(rows[1]["expected_maternal_offspring"])
            assert rows[2]["expected_maternal_offspring"] == 0  # No-host female, not NA.
            assert not rows[3]["eligible_as_donor"]  # Dead male cannot donate.
        runtime.e["fp" if kind == 0 else "fi"][:] = 0
        assert not component.forward(runner.state, {})["offspring"]["offspring_count"].any()
        runtime.e["fp" if kind == 0 else "fi"][0] = 1
        runtime.e["S_plant" if kind == 0 else "S_insect"][1] = 0
        no_father = component.forward(runner.state, {})["offspring"]
        assert not no_father["candidates"][0].size
        assert not no_father["offspring_count"].any()
        assert runtime.result["parents"][-4]["expected_maternal_offspring"] == 0


def test_v015_cost_functions_and_maternal_mean():
    c = config(plant_cost=0.5, voc_cost=0.25, insect_cost=0.5)
    D, I = np.array([[1.0, 1.0], [0.0, 0.0]]), np.array([[1.0, 1.0], [0.0, 0.0]])
    result = settle_feeding(
        D,
        I,
        np.array([0, 0]),
        np.array([0, 0]),
        np.ones((2, 2)),
        np.array([0, 1]),
        c,
        np.full(2, 10.0),
        np.ones(2),
        chemical=np.column_stack((D, np.ones(2))),
    )
    np.testing.assert_allclose(result["fp"], [0, 0.75])
    np.testing.assert_allclose(result["fi"], [0, 1])
    # Costs are maternal multipliers; survival remains a function of the resource ledger.
    np.testing.assert_allclose(result["g"], [1, 1])
    assert np.isclose(result["h"].sum(), 0.2)


def test_v015_native_zero_multiplier_fathers_and_contribution_audit():
    for kind in ("plant", "insect"):
        c = config(
            generations=1,
            plant_sizes=[20, 0] if kind == "plant" else [1, 0],
            insect_sizes=[0] if kind == "plant" else [40],
            sP0=1.0,
            choice_mode="forced",
            plant_cost=1.2 if kind == "plant" else 0.0,
            voc_cost=0.0,
            insect_cost=1.2,
            resource_capacity={"P1": 1000.0, "P2": 1000.0},
        )
        c.plant_offspring_lambda["P1"] = 100.0
        c.insect_offspring_lambda["H1"] = 100.0
        c.spatial.update(
            initial_plant_cells=[[0, 0]] * sum(c.plant_sizes),
            initial_insect_cells=[[0, 0]] * sum(c.insect_sizes) if sum(c.insect_sizes) else None,
        )
        c.genetics["mutation_enabled"] = dict(plant=False, insect=False)
        if kind == "insect":
            # Zero plant defense makes feeding survival exactly one for this role check.
            for allele in c.genetics["alleles"]:
                if any(
                    l["locus_id"] == allele["locus_id"] and l["organism"] == "plant"
                    for l in c.genetics["loci"]
                ):
                    allele["effect"][: c.n] = [0.0] * c.n
        result = simulate(c)
        assert result["metadata"]["status"]["status"] == "completed"
        parents = {
            r["individual_id"]: r
            for r in result["parents"]
            if r["generation"] == 0 and r["kind"] == kind
        }
        zero = [r for r in parents.values() if r["F_i"] == 0]
        assert zero and all(r["offspring_count_as_mother"] == 0 for r in zero)
        assert any(r["offspring_count_sired"] > 0 and r["eligible_as_donor"] for r in zero)
        offspring = [r for r in result["offspring"] if r["kind"] == kind]
        assert offspring and sum(r["offspring_count_as_mother"] for r in parents.values()) == len(
            offspring
        )
        assert sum(r["offspring_count_sired"] for r in parents.values()) == len(offspring)
        assert sum(r["realized_genetic_contribution"] for r in parents.values()) == len(offspring)
        for child in offspring:
            assert child["maternal_id"] != child["paternal_id"]
            if kind == "insect":
                assert parents[child["maternal_id"]]["sex_or_reproductive_role"] == "female"
                assert parents[child["paternal_id"]]["sex_or_reproductive_role"] == "male"
        start = result["snapshots"][0]
        start_lookup = {int(i): row for row, i in enumerate(start[f"{kind}_ids"])}
        for row, child in enumerate(offspring):
            for copy, role in enumerate(("maternal_id", "paternal_id")):
                source = start[f"{kind}_genotypes"][start_lookup[child[role]]]
                for locus, allele in enumerate(json.loads(child["genotype"])[copy]):
                    assert allele in source[:, locus]
        with TemporaryDirectory() as tmp:
            folder = Path(tmp) / kind
            save_result(result, folder)
            assert audit_run(folder)["genetic_audit"]
            # Metadata relabeling must not silently apply the new pool rule to old versions.
            path = folder / "metadata.json"
            original = path.read_bytes()
            meta = json.loads(original)
            meta["version"] = meta["config"]["model_version"] = "0.1.4"
            meta["config_sha256"] = hashlib.sha256(
                json.dumps(meta["config"], sort_keys=True, separators=(",", ":")).encode()
            ).hexdigest()
            write_json(path, meta)
            raises(lambda: audit_run(folder), (AssertionError, ValueError))
            path.write_bytes(original)
            assert audit_run(folder)["contribution_audit"]


def test_v015_historical_configuration_execution_boundary():
    old = config(model_version="0.1.4", generations=0)
    # A new chemical protocol cannot be relabeled as a pre-chemical model version.
    raises(old.validate)
    raises(lambda: simulate(old))
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "old_config"
        run(
            dict(base=config_payload(old), cases=[dict(name="check", overrides={})], replicates=1),
            path,
        )
        assert audit_run(path / "data/check/rep01")["status"] == "invalid_parameter"
