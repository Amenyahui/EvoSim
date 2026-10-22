"""T01–T16 acceptance checks; statistical checks validate code distributions only."""

import copy
import hashlib
import json
from dataclasses import asdict, replace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect.audit import audit, audit_run
from evo_sim.models.plant_insect.batch import run, save_result, write_json
from evo_sim.models.plant_insect.helpers import runtime as model
from evo_sim.models.plant_insect.helpers import traits as trait_functions
from evo_sim.models.plant_insect.helpers import voc
from evo_sim.models.plant_insect.helpers.config import Config
from evo_sim.models.plant_insect.helpers.demography import (
    alpha_from_K,
    density_phi,
    draw_offspring,
    inherit,
)
from evo_sim.models.plant_insect.helpers.runtime import ecology, streams
from evo_sim.models.plant_insect.simulator import simulate

ROOT = Path(__file__).resolve().parents[1]
BASE = json.loads(
    (ROOT / "tests/fixtures/development_fixture_v13.json").read_text(encoding="utf-8")
)["base"]


def tiny(**kwargs):
    data = copy.deepcopy(BASE)
    data.update(plant_sizes=[4, 4], insect_sizes=[4, 4], generations=3)
    data.update(kwargs)
    return Config(**data)


def state_for(c, D, I, A=None, trait=None):
    B, w = c.validate()
    if A is None:
        A = np.zeros((len(I), len(B)))
    if trait is None:
        trait = D @ B.T
    state = dict(
        D=D,
        I=I,
        A=A,
        plant_voc_independent_trait=trait,
        plant_species=np.repeat(np.arange(len(c.plant_sizes)), c.plant_sizes),
        insect_species=np.repeat(np.arange(len(c.insect_sizes)), c.insect_sizes),
        plant_ids=np.arange(len(D)),
        insect_ids=np.arange(len(I)),
    )
    return model.CommunityState.from_arrays(state), B, w


def test_T01_single_pair():
    c = tiny(
        plant_sizes=[1, 0],
        insect_sizes=[1, 0],
        choice_mode="forced",
        plant_cost=0.0,
        insect_cost=0.0,
        resource_capacity={"P1": 2.0, "P2": 1.0},
    )
    state, B, w = state_for(c, np.array([[0.4, 0.8]]), np.array([[0.2, 0.4]]))
    e = ecology(state, c, B, w, np.random.default_rng(1), streams(c)[1][1])
    assert np.allclose([e["q"][0], e["g"][0], e["S_insect"][0], e["h"][0]], [0.5, 0.5, 0.5, 0.25])
    assert np.allclose(e["x"], [0.5])
    assert np.allclose(e["gamma"], [1.0])
    assert np.allclose(e["S_plant"], [0.775])


def test_T02_zero_food():
    for c, reason in [
        (tiny(choice_mode="outside", outside_utility=10000.0), "rejected"),
        (tiny(plant_sizes=[0, 0]), "no_candidate"),
        (tiny(choice_mode="forced", insect_I_initial_max=[0.0, 0.0]), "accepted"),
    ]:
        r = simulate(c)
        ev = pd.DataFrame(r["events"]).query("generation==0")
        assert ev.reason.eq(reason).all() and (ev[["d", "x", "g", "W_insect"]] == 0).all().all()
        assert ev.q.isna().all() if reason != "accepted" else ev.q.eq(0).all()
        parent = pd.DataFrame(r["parents"]).query('generation==0 and kind=="insect"')
        assert (
            parent.S_i.eq(0).all()
            and not parent.survived_to_reproduce.any()
            and parent.offspring_count.eq(0).all()
        )


def test_T03_survival_endpoints():
    a = draw_offspring(
        np.array([0.0, 0.0, 1.0, 1.0]),
        np.ones(4),
        np.full(4, 2.0),
        np.ones(4),
        np.random.default_rng(8),
        np.random.default_rng(9),
    )
    assert a["survived_to_reproduce"].tolist() == [False, False, True, True]
    assert a["offspring_count"][:2].tolist() == [0, 0]
    assert np.array_equal(a["offspring_count"][2:], np.random.default_rng(9).poisson([2.0, 2.0]))


def test_T04_zero_poisson():
    a = draw_offspring(
        np.ones(100), np.zeros(100), 2.0, 1.0, np.random.default_rng(2), np.random.default_rng(3)
    )
    assert a["survived_to_reproduce"].all() and not a["offspring_count"].any()


def test_T05_two_stage_moments():
    a = draw_offspring(
        np.full(200000, 0.5),
        np.ones(200000),
        2.0,
        1.0,
        np.random.default_rng(20),
        np.random.default_rng(21),
    )
    assert np.all(a["conditional_offspring_mean_mu"] == 2)
    assert np.all(a["expected_offspring"] == 1) and np.all(a["offspring_variance"] == 2)
    b = a["offspring_count"]
    assert abs(b.mean() - 1) < 0.015 and abs(b.var() - 2) < 0.04
    print("T05 implementation moments:", float(b.mean()), float(b.var()))


def test_T06_BH_equilibrium_not_cap():
    alpha = alpha_from_K(2.0, 80.0)
    phi = density_phi([[alpha]], [80])[0]
    assert phi == 0.5 and 80 * 2 * phi == 80
    a = draw_offspring(
        np.ones(80 * 2000),
        np.ones(80 * 2000),
        2.0,
        phi,
        np.random.default_rng(5),
        np.random.default_rng(6),
    )
    totals = a["offspring_count"].reshape(-1, 80).sum(axis=1)
    assert (totals > 80).any() and (totals < 80).any() and abs(totals.mean() - 80) < 1
    for lam, K in [(1.0, 80.0), (0.5, 80.0), (2.0, 0.0)]:
        try:
            alpha_from_K(lam, K)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid ideal K parameterization accepted")


def test_T07_competition_rows_and_species_keys():
    c = tiny(
        plant_competition_alpha={"P2": {"P2": 0.03, "P1": 0.04}, "P1": {"P2": 0.02, "P1": 0.01}}
    )
    _, _, a = c.resolved_demography()
    assert np.allclose(density_phi(a, [80, 20]), [1 / (1 + 0.8 + 0.4), 1 / (1 + 3.2 + 0.6)])
    reordered = replace(
        c, plant_competition_alpha={"P1": {"P1": 0.01, "P2": 0.02}, "P2": {"P1": 0.04, "P2": 0.03}}
    )
    assert np.array_equal(reordered.resolved_demography()[2], a)
    assert np.array_equal(density_phi(np.zeros((3, 3)), [1, 0, 5]), np.ones(3))
    data = asdict(tiny())
    ps = ["P1", "P2", "P3"]
    hs = ["H1", "H2", "H3"]
    entry = data["parameter_provenance"]["plant_offspring_lambda"]["P1"]
    data.update(
        plant_sizes=[2, 1, 0],
        insect_sizes=[1, 2, 0],
        plant_centers=None,
        insect_I_initial_max=None,
        insect_A_initial_bounds=None,
        B=[[1.0, 0.0], [0.0, 1.0], [0.5, 0.5]],
        weights=[1 / 3] * 3,
        resource_capacity={k: 1.0 for k in ps},
        max_feeding_demand={k: 1.0 for k in hs},
        plant_offspring_lambda={k: 2.0 for k in ps},
        insect_offspring_lambda={k: 2.0 for k in hs},
        plant_competition_alpha={k: {j: 0.01 for j in ps} for k in ps},
    )
    data["parameter_provenance"].update(
        plant_offspring_lambda={k: entry for k in ps},
        insect_offspring_lambda={k: entry for k in hs},
        plant_competition_alpha={k: {j: entry for j in ps} for k in ps},
    )
    r = simulate(Config(**data))
    assert len(pd.DataFrame(r["networks"]).query("generation==0")) == 9
    assert r["snapshots"][0]["A"].shape == (3, 3)


def test_T08_no_extra_insect_density():
    c = tiny(plant_sizes=[1, 0], insect_sizes=[4, 4], choice_mode="forced")
    state, B, w = state_for(c, np.zeros((1, 2)), np.zeros((8, 2)))
    e = ecology(state, c, B, w, np.random.default_rng(9), streams(c)[1][1])
    assert np.isclose(e["x"].sum(), 1.0) and np.allclose(e["S_insect"], np.full(8, 1 / 8))
    r = simulate(
        c,
        initial=(
            state.plants.defense,
            state.insects.resistance,
            state.insects.preference,
            state.plants.voc_trait,
        ),
    )
    parents = pd.DataFrame(r["parents"]).query('generation==0 and kind=="insect"')
    assert parents.phi_i.eq(1).all()
    assert np.allclose(parents.conditional_offspring_mean_mu, parents.lambda_s * parents.F_i)


def test_T09_partial_extinction_absorbing():
    r = simulate(tiny(insect_I_initial_max=[0.4, 0.0], sP0=1.0, plant_cost=0.0))
    pop = pd.DataFrame(r["populations"])
    assert pop.query('species=="H2" and generation>=1').N_start.eq(0).all()
    assert r["metadata"]["status"]["generation"] == 3
    assert pop.query('kind=="plant" and generation==1').N_start.sum() > 0
    assert (
        pd.DataFrame(r["networks"])
        .query('insect_species=="H2" and generation>=1')
        .N_accept.eq(0)
        .all()
    )


def test_T10_no_plants():
    r = simulate(tiny(plant_sizes=[0, 0]))
    assert (
        r["metadata"]["status"]["status"] == "all_species_extinct"
        and r["metadata"]["status"]["generation"] == 1
    )
    ev = pd.DataFrame(r["events"])
    assert ev.reason.eq("no_candidate").all() and ev.q.isna().all()
    assert len(r["snapshots"][-1]["insect_ids"]) == 0


def test_T11_all_extinct_initial():
    r = simulate(tiny(plant_sizes=[0, 0], insect_sizes=[0, 0]))
    assert r["metadata"]["status"]["status"] == "all_species_extinct"
    assert pd.DataFrame(r["populations"]).N_start.eq(0).all()
    assert len(r["networks"]) == 4 and r["events"] == r["plants"] == []
    assert all(pd.isna(row["mean_g"]) and pd.isna(row["mean_h"]) for row in r["populations"])


def test_T12_inheritance_and_mutation():
    c = tiny(mutation_rate=0.0, mu_A=0.0, mu_Z=0.0)
    traits = {"I": np.array([[0.1, 0.2], [0.8, 0.9]]), "A": np.array([[0.2, 0.3], [0.9, 0.8]])}
    children, labels, ids, rows, nxt = inherit(
        traits,
        np.array([0, 1]),
        np.array([10, 20]),
        np.array([2, 3]),
        c,
        "insect",
        [np.random.default_rng(2), np.random.default_rng(3)],
        30,
        1,
    )
    assert labels.tolist() == [0, 0, 1, 1, 1] and ids.tolist() == list(range(30, 35))
    assert [row["parent_id"] for row in rows] == [10, 10, 20, 20, 20]
    for key in traits:
        assert np.array_equal(children[key], traits[key][[0, 0, 1, 1, 1]])
    c = replace(c, mutation_rate=1.0, mu_Z=1.0)
    traits = {
        "D": np.array([[0.1, 0.2], [0.8, 0.9]]),
        "plant_voc_independent_trait": np.array([[0.2, 0.3], [0.9, 0.8]]),
    }
    rngs = [np.random.default_rng(7), np.random.default_rng(8)]
    reference = copy.deepcopy(rngs)
    children, *_ = inherit(
        traits, np.array([0, 1]), np.array([10, 20]), np.array([1, 1]), c, "plant", rngs, 30, 1
    )
    for s in range(2):
        assert np.array_equal(
            children["D"][s : s + 1],
            trait_functions.mutate(
                traits["D"][s : s + 1], c.mutation_rate, c.mutation_sd, reference[s]
            ),
        )
        assert np.array_equal(
            children["plant_voc_independent_trait"][s : s + 1],
            voc.mutate_Z(traits["plant_voc_independent_trait"][s : s + 1], c, reference[s])[0],
        )


def test_T13_mass_units_and_streams():
    c = tiny(resource_capacity={"P1": 0.5, "P2": 2.0}, max_feeding_demand={"H1": 1.0, "H2": 2.0})
    a = simulate(c)
    b = simulate(c)
    scaled = simulate(
        replace(
            c,
            resource_capacity={"P1": 500.0, "P2": 2000.0},
            max_feeding_demand={"H1": 1000.0, "H2": 2000.0},
            resource_unit="mg_edible_dry_biomass",
            demand_unit="mg_edible_dry_biomass",
        )
    )
    for key in ("events", "plants", "parents", "offspring", "populations", "networks"):
        pd.testing.assert_frame_equal(pd.DataFrame(a[key]), pd.DataFrame(b[key]), check_exact=True)
    pa, pb = pd.DataFrame(a["parents"]), pd.DataFrame(scaled["parents"])
    for key in ("S_i", "F_i", "phi_i", "W_i", "offspring_count"):
        assert np.allclose(pa[key], pb[key], equal_nan=True)
    for x, y in zip(a["snapshots"], scaled["snapshots"]):
        for key in (
            "D",
            "I",
            "A",
            "plant_voc_independent_trait",
            "plant_species",
            "insect_species",
        ):
            assert np.array_equal(x[key], y[key])
    assert np.allclose(pd.DataFrame(a["events"]).x * 1000, pd.DataFrame(scaled["events"]).x)
    x, y = streams(c), streams(c)
    x[2][0][0].random(500)
    assert np.array_equal(x[1][0][0].random(10), y[1][0][0].random(10))


def test_T14_parameters_and_safety():
    invalid = [
        dict(plant_offspring_lambda={"P1": -1.0, "P2": 2.0}),
        dict(insect_offspring_lambda=None),
        dict(plant_offspring_lambda={"P1": float("inf"), "P2": 2.0}),
        dict(plant_competition_alpha={"P1": {"P1": 1.0}}),
        dict(parameter_provenance=None),
        dict(sI0=0.1),
        dict(demand_unit="wrong"),
        dict(illustrative_only=False),
        dict(plant_sizes=[-1, 4]),
    ]
    for change in invalid:
        try:
            tiny(**change).validate()
        except (ValueError, TypeError):
            pass
        else:
            raise AssertionError(change)
    c = tiny(plant_offspring_lambda={"P1": 1000.0, "P2": 1000.0}, max_individuals_safety=20)
    r = simulate(c)
    assert r["metadata"]["status"]["status"] == "computational_limit"
    pop = pd.DataFrame(r["populations"])
    assert pop.N_offspring.sum() > 20 and not pop.update_applied.any() and pop.N_next.isna().all()
    assert not r["offspring"]
    initial = simulate(tiny(max_individuals_safety=2))
    assert initial["metadata"]["status"]["status"] == "computational_limit"
    with TemporaryDirectory() as t:
        for name, result in [("sampled", r), ("initial", initial)]:
            save_result(result, Path(t) / name)
            assert (
                json.loads((Path(t) / name / "status.json").read_text())["status"]
                == "computational_limit"
            )
        spec = dict(
            replicates=1, base=asdict(tiny()), cases=[dict(name="bad", overrides=invalid[0])]
        )
        out = run(spec, Path(t) / "invalid")
        assert audit(out, Path(t) / "reports")["outcomes"] == {"invalid_parameter": 1}
    # An unrepresentable resource total is arithmetic_error, never extinction.
    arithmetic = simulate(
        tiny(
            plant_sizes=[1, 0],
            insect_sizes=[8, 0],
            resource_capacity={"P1": 8e307, "P2": 1.0},
            max_feeding_demand={"H1": 8e307, "H2": 1.0},
            choice_mode="forced",
            initial_spread=0.0,
            plant_centers=None,
            insect_I_initial_max=[0.0, 0.0],
        )
    )
    assert arithmetic["metadata"]["status"]["status"] == "arithmetic_error"
    with TemporaryDirectory() as t:
        save_result(arithmetic, Path(t) / "arithmetic")
        assert json.loads((Path(t) / "arithmetic/status.json").read_text())["message"]


def test_T15_resource_conservation_and_order():
    rng = np.random.default_rng(901)
    host = rng.integers(0, 7, 1000)
    q, capacity, demand = (
        rng.uniform(0, 1, 1000),
        rng.uniform(0.01, 1, 7),
        rng.uniform(0.1, 2, 1000),
    )
    a = model.shared.allocate_resources(host, q, capacity, demand)
    assert (a["x"] >= 0).all() and (a["x"] <= a["d"]).all()
    assert (a["X"] <= capacity).all()
    np.testing.assert_allclose(a["X"], np.minimum(a["total_demand"], capacity))
    order = rng.permutation(len(host))
    b = model.shared.allocate_resources(host[order], q[order], capacity, demand[order])
    np.testing.assert_allclose(a["x"][order], b["x"])


def test_T16_saved_audit_and_damage_detection():
    spec = dict(
        replicates=1,
        base=asdict(tiny()),
        cases=[
            dict(name="normal", overrides={}),
            dict(name="empty", overrides=dict(plant_sizes=[0, 0])),
        ],
    )
    with TemporaryDirectory() as t:
        out = run(spec, Path(t) / "run")
        summary = audit(out, Path(t) / "reports")
        assert summary["runs"] == 2 and summary["max_resource_overrun"] == 0
        path = out / "data/normal/rep01"
        file = path / "parents.csv"
        table = pd.read_csv(file)
        table.loc[0, "W_i"] += 0.1
        table.to_csv(file, index=False)
        meta = json.loads((path / "metadata.json").read_text())
        meta["output_sha256"]["parents.csv"] = hashlib.sha256(file.read_bytes()).hexdigest()
        write_json(path / "metadata.json", meta)
        try:
            audit_run(path)
        except AssertionError:
            pass
        else:
            raise AssertionError("independent audit accepted corrupted individual W")


if __name__ == "__main__":
    tests = [f for name, f in globals().copy().items() if name.startswith("test_T")]
    for test in tests:
        test()
        print("PASS", test.__name__, flush=True)
    print(f"{len(tests)} T01–T16 acceptance groups passed; no ecological hypothesis tests")
