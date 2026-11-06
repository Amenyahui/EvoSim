"""Resource-coupled plant reproduction and cross-label inheritance."""

from types import SimpleNamespace

import numpy as np

from evo_sim.models.plant_insect import simulate
from evo_sim.models.plant_insect.biology.plants.demography import plant_resource_update
from evo_sim.models.plant_insect.simulation.stages.adult_reproduction import sample_parents
from tests.support import config


def test_cross_label_parents_and_retained_genotypes():
    c = config(plant_sizes=[1, 1], insect_sizes=[0], generations=2)
    c.plant_demography.update(background_mortality_rate=0, defoliation_mortality_rate=0)
    c.plant_offspring_lambda = {"P1": 100.0, "P2": 100.0}
    c.spatial.update(
        grid_shape=[1, 1],
        initial_plant_cells=[[0, 0], [0, 0]],
        plant_dispersal_radius=0,
        search_radius=0,
        insect_dispersal_radius=0,
    )
    c.plant_demography["microsites_per_cell"] = 6
    c.sexual_reproduction["insect_mate_radius"] = 0
    result = simulate(c)
    assert result["metadata"]["status"]["status"] == "completed"
    births = [r for r in result["offspring"] if r["kind"] == "plant" and r["generation"] == 1]
    assert {r["maternal_species"] for r in births} == {"P1", "P2"}
    assert all(r["maternal_species"] != r["paternal_species"] for r in births)
    assert all(r["maternal_id"] != r["paternal_id"] for r in births)
    initial = result["snapshots"][0]
    for snapshot in result["snapshots"][1:]:
        for i, ident in enumerate(initial["plant_ids"]):
            j = np.flatnonzero(snapshot["plant_ids"] == ident)[0]
            np.testing.assert_array_equal(
                snapshot["plant_genotypes"][j], initial["plant_genotypes"][i]
            )
    assert {r["generation"] for r in result["parents"] if r["kind"] == "plant"} == {0, 1}
    transmissions = [r for r in result["gamete_transmissions"] if r["kind"] == "plant"]
    assert {r["copy"] for r in transmissions} == {0, 1}
    assert all(r["role"] == ("maternal" if r["copy"] == 0 else "paternal") for r in transmissions)
    for row in result["allele_contributions"]:
        if row["kind"] == "plant" and row["N_offspring"]:
            np.testing.assert_allclose(
                row["p_parent"] + row["delta_parent_contribution"] + row["delta_transmission"],
                row["p_offspring"],
                atol=1e-12,
            )


def test_resource_limits_plant_reproduction():
    c = config()
    post, deficit, survival, recovered = plant_resource_update(
        [10, 4, 0], [0, 2, 0], [10] * 3, c.plant_demography
    )
    np.testing.assert_allclose(post, [10, 2, 0])
    np.testing.assert_allclose(deficit, [0, 0.8, 1])
    assert np.all(recovered >= post) and np.all(survival <= 1)

    class Population(SimpleNamespace):
        def __len__(self):
            return len(self.ids)

    pop = Population(
        ids=np.array([0, 1, 2]), species=np.array([0, 0, 0]), cells=np.zeros((3, 2), int)
    )
    rng = [[np.random.default_rng(1)], [np.random.default_rng(2)]]
    runtime = SimpleNamespace(c=SimpleNamespace(), lp=np.array([100.0]), rng=[None, None, rng, rng])
    _, _, rows, survived = sample_parents(
        pop, "plant", np.ones(3), np.array([1.0, 0.2, 0.0]), runtime, 0
    )
    assert survived.all()
    assert [row["conditional_mu"] for row in rows] == [100, 20, 0]
    assert rows[2]["propagule_count"] == 0
    assert all(row["candidate_mate_count"] == 2 for row in rows)
