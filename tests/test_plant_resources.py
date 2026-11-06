"""Resource stock, finite establishment and atomic-storage contracts."""

import copy
from types import SimpleNamespace

import numpy as np

from evo_sim.models.plant_insect import get_runner, simulate
from evo_sim.models.plant_insect.helpers.demography import plant_resource_update
from evo_sim.models.plant_insect.helpers.spatial import establish_propagules
from tests.support import config, raises


def plant_config(**changes):
    c = config(**dict(dict(plant_sizes=[2, 0], insect_sizes=[0], generations=10), **changes))
    c.plant_demography.update(
        background_mortality_rate=0.0,
        defoliation_mortality_rate=0.0,
        microsites_per_cell=4,
    )
    c.spatial.update(
        grid_shape=[1, 1],
        search_radius=0,
        plant_dispersal_radius=0,
        insect_dispersal_radius=0,
        initial_plant_cells=[[0, 0]] * sum(c.plant_sizes),
        initial_insect_cells=None,
    )
    c.sexual_reproduction.update(insect_mate_radius=0)
    c.genetics["mutation_enabled"] = {"plant": False, "insect": False}
    return c


def runner_for(c):
    from evo_sim.models.plant_insect.config.schema import config_payload

    r = get_runner()
    r.config["parameters"] = config_payload(c)
    r.config["parameters"].pop("generations")
    r.config["simulation_metadata"]["num_steps_per_episode"] = c.generations
    r.init()
    return r


def test_resource_equations_and_invalid_stock():
    c = plant_config()
    spec = c.plant_demography
    spec.update(background_mortality_rate=np.log(2), defoliation_mortality_rate=np.log(4))
    post, deficit, survival, recovered = plant_resource_update(
        [10.0, 4.0, 0.0],
        [0.0, 2.0, 0.0],
        [10.0, 8.0, 10.0],
        spec,
    )
    np.testing.assert_allclose(post, [10, 2, 0])
    np.testing.assert_allclose(deficit, [0, 0.75, 1])
    np.testing.assert_allclose(survival, 0.5 * 4**-deficit)
    np.testing.assert_allclose(recovered, [10, 5, 5])
    for stock, consumed, capacity in (
        ([-1], [0], [2]),
        ([3], [0], [2]),
        ([1], [2], [2]),
        ([0], [0], [0]),
    ):
        raises(lambda: plant_resource_update(stock, consumed, capacity, spec))
    for key in ("sP0", "plant_competition_alpha"):
        values = dict(get_runner().config["parameters"], generations=1, **{key: 1})
        from evo_sim.models.plant_insect import Config

        raises(lambda: Config(**values), TypeError)
    bad = copy.deepcopy(c)
    bad.resource_capacity["P1"] = 0
    raises(bad.validate)


def test_resource_deficit_depends_only_on_current_stock_and_recovers():
    spec = plant_config().plant_demography
    spec.update(background_mortality_rate=-np.log(0.99), defoliation_mortality_rate=0.3)
    post, deficit, survival, recovered = plant_resource_update([2.5], [0], [10], spec)
    np.testing.assert_allclose(deficit, 0.75)
    np.testing.assert_allclose(survival, 0.99 * np.exp(-0.3 * 0.75))
    np.testing.assert_allclose(recovered, 6.25)
    _, following, _, _ = plant_resource_update(recovered, [0], [10], spec)
    assert following[0] < deficit[0]


def test_initial_sites_are_unique_and_explicit_overflow_is_rejected():
    c = plant_config()
    c.plant_demography["microsites_per_cell"] = 1
    raises(c.validate)
    c.spatial.update(grid_shape=[2, 2], initial_plant_cells=None)
    c.plant_sizes = [2, 2]
    first, second = runner_for(c), runner_for(c)
    sites = first.state["agents"]["plants"]["microsite"]
    assert len(sites) == len(set(sites)) == 4
    np.testing.assert_array_equal(sites, second.state["agents"]["plants"]["microsite"])


def test_full_sites_retain_adults_without_manufacturing_newborns():
    c = plant_config()
    c.plant_demography["microsites_per_cell"] = 2
    result = simulate(c)
    assert result["metadata"]["status"]["status"] == "completed"
    assert result["plant_proposals"] and not result["offspring"]
    assert not result["gamete_transmissions"] and not result["metadata"]["slim_calls"]
    for state in result["snapshots"]:
        np.testing.assert_array_equal(state["plant_ids"], [0, 1])
        np.testing.assert_array_equal(state["plant_microsite"], [0, 1])
    plants = [p for p in result["populations"] if p["kind"] == "plant" and p["species"] == "P1"]
    assert all(p["N_retained"] == p["N_next"] == 2 and p["N_offspring"] == 0 for p in plants)


def test_each_round_establishment_and_newborn_resource():
    c = plant_config()
    c.plant_offspring_lambda = {"P1": 20.0, "P2": 20.0}
    result = simulate(c)
    assert result["metadata"]["status"]["status"] == "completed"
    births = [p for p in result["offspring"] if p["kind"] == "plant"]
    assert len(births) == 2 and {p["generation"] for p in births} == {1}
    new_ids = {p["offspring_id"] for p in births}
    for parent in result["parents"]:
        if parent["kind"] == "plant":
            assert parent["expected_maternal_offspring"] is None
            assert parent["reproduction_active"]
            if parent["individual_id"] in new_ids:
                assert parent["eligible_as_mother"]
                assert parent["conditional_mu"] > 0
    assert {row["reproduction_round"] for row in result["plant_proposals"]} == set(range(10))
    for newborn in births:
        state = result["snapshots"][1]
        i = np.flatnonzero(state["plant_ids"] == newborn["offspring_id"])[0]
        assert state["plant_resource"][i] == c.plant_demography["newborn_resource_fraction"]
        assert newborn["paternal_id"] != newborn["maternal_id"]


def test_shared_lottery_is_unbiased_and_never_fills_untargeted_sites():
    c = plant_config()
    c.plant_demography["microsites_per_cell"] = 1
    c.spatial["grid_shape"] = [2, 2]
    rng = np.random.default_rng(771)
    counts = np.zeros(4, int)
    for _ in range(2000):
        winners, sites = establish_propagules(np.zeros((4, 2), int), [], c, rng)
        assert len(winners) == 1 and sites[winners[0]] == 0
        counts[winners[0]] += 1
    assert np.all(np.abs(counts / counts.sum() - 0.25) < 0.04), counts
    winners, sites = establish_propagules(np.zeros((4, 2), int), [0], c, rng)
    assert not len(winners) and (sites == -1).all()


def test_survivor_snapshot_budget_and_failure_rollback():
    c = plant_config(generations=4)
    c.plant_demography["microsites_per_cell"] = 2
    initial = 2 * sum(c.plant_sizes) * 6
    c.genetics["max_genotype_values_safety"] = initial * 5
    passed = simulate(c)
    assert passed["metadata"]["status"]["status"] == "completed"
    assert passed["metadata"]["stored_genotype_values"] == initial * 5
    c.genetics["max_genotype_values_safety"] = 2 * initial - 1
    r = runner_for(c)
    before = copy.deepcopy(r.state["agents"])
    rt = r.state["environment"]["runtime"]
    rng_before = copy.deepcopy(rt.rng)
    r.step(1)
    assert r.result["metadata"]["status"]["status"] == "computational_limit"
    assert r.result["metadata"]["stored_genotype_values"] == initial
    assert not any(r.result[name] for name in ("plants", "parents", "offspring", "plant_proposals"))
    for kind, values in before.items():
        for key, value in values.items():
            np.testing.assert_array_equal(value, r.state["agents"][kind][key])
    for phase in rng_before:
        for old_taxon, new_taxon in zip(rng_before[phase], rt.rng[phase], strict=True):
            for old, new in zip(old_taxon, new_taxon, strict=True):
                assert old.bit_generator.state == new.bit_generator.state


def test_only_established_seedlings_get_native_genotypes():
    c = plant_config(generations=1)
    c.plant_demography.update(microsites_per_cell=3)
    c.plant_offspring_lambda = {"P1": 20.0, "P2": 20.0}
    result = simulate(c)
    assert result["metadata"]["status"]["status"] == "completed"
    assert len(result["offspring"]) == 1
    assert len(result["plant_proposals"]) > 1
    assert sum(row["established"] for row in result["plant_proposals"]) == 1
    assert len(result["gamete_transmissions"]) == 12
    assert [call["offspring_count"] for call in result["metadata"]["slim_calls"]] == [1]
    assert result["snapshots"][1]["plant_ids"].tolist() == [0, 1, 2]


def test_feeding_uses_carried_stock_then_survivors_regrow():
    c = plant_config(generations=2, insect_sizes=[2], choice_mode="forced")
    r = runner_for(c)
    r.state["agents"]["insects"]["sex"][:] = [0, 1]
    rt = r.state["environment"]["runtime"]
    rt.rng[3][1][0] = SimpleNamespace(poisson=lambda mu: np.full(len(mu), 12, int))
    r.step()
    assert r.result["metadata"]["status"]["status"] == "completed"
    first = [row for row in r.result["plants"] if row["generation"] == 0]
    second = [row for row in r.result["plants"] if row["generation"] == 1]
    assert sum(row["X"] for row in first) > 0
    parents = {
        row["individual_id"]: row
        for row in r.result["parents"]
        if row["kind"] == "plant" and row["generation"] == 0
    }
    for row in first:
        expected_f = max(0, 1 - row["defense_cost"] - row["voc_cost"])
        expected_f *= row["resource_after_feeding"] / row["resource_capacity"]
        assert np.isclose(parents[row["plant_id"]]["F_i"], expected_f)
        assert np.isclose(row["W_plant"], parents[row["plant_id"]]["W_i"])
        if row["resource_after_feeding"] == 0:
            assert parents[row["plant_id"]]["propagule_count"] == 0
    for before, after in zip(first, second, strict=True):
        expected = 1 - (1 - before["resource_after_feeding"]) * np.exp(
            -c.plant_demography["resource_recovery_rate"]
        )
        assert before["resource_next"] == after["R"] == expected
    assert any(row["R"] < row["resource_capacity"] for row in second)


if __name__ == "__main__":
    for name, check in list(globals().items()):
        if name.startswith("test_") and callable(check):
            check()
            print(f"PASS {name}", flush=True)
