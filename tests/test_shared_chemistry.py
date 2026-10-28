"""Shared chemical axes: analytical paths, real QTL widths and checked persistence."""

from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from evo_sim.models.plant_insect import simulate
from evo_sim.models.plant_insect.audit import audit_run
from evo_sim.models.plant_insect.batch import save_result
from evo_sim.models.plant_insect.helpers import chemistry, voc
from evo_sim.models.plant_insect.helpers.genetics import GeneticState
from evo_sim.models.plant_insect.helpers.inputs import load_inputs
from evo_sim.models.plant_insect.helpers.traits import consumption
from tests.support import config, raises


def test_shared_module_paths_and_hand_calculation():
    c = config(eta_shared=0.5)
    x = np.array([[0.2, 0.6, 0.4]])
    D, V = chemistry.outputs(x, c.eta_shared)
    np.testing.assert_array_equal(D, [[0.2, 0.6]])
    np.testing.assert_allclose(V, [[0.3, 0.4]])
    fm, cd, cv = voc.fecundity(D, V, c, x)
    np.testing.assert_allclose([fm[0], cd[0], cv[0]], [0.88, 0.10, 0.02])
    for j in range(3):
        changed = x.copy()
        changed[:, j] += 0.1
        d, v = chemistry.outputs(changed, c.eta_shared)
        new_fm = voc.fecundity(d, v, c, changed)[0]
        np.testing.assert_allclose(fm - new_fm, 0.1 * (c.plant_cost if j < 2 else c.voc_cost))
        np.testing.assert_allclose(d - D, [[0.1 if j == 0 else 0, 0.1 if j == 1 else 0]])
        np.testing.assert_allclose(v - V, [[0.05 if j == 1 else 0, 0.1 if j == 2 else 0]])
    for eta in (0.0, 1.0):
        d, v = chemistry.outputs(x, eta)
        np.testing.assert_array_equal(d, D)
        np.testing.assert_allclose(v, [[eta * 0.6, 0.4]])
        np.testing.assert_array_equal(voc.fecundity(d, v, config(eta_shared=eta), x)[0], fm)


def test_shared_module_boundaries_and_configuration_contract():
    c = config()
    for value in (0.0, 1.0):
        x = np.full((2, 3), value)
        D, V = chemistry.outputs(x, c.eta_shared)
        assert np.isfinite(D).all() and np.isfinite(V).all()
        expected = max(0, 1 - 2 * c.plant_cost * value - c.voc_cost * value)
        np.testing.assert_allclose(voc.fecundity(D, V, c, x)[0], expected)
    x = np.ones((1, 3))
    D, V = chemistry.outputs(x, 1.0)
    assert voc.fecundity(D, V, config(plant_cost=1.0), x)[0][0] == 0
    np.testing.assert_array_equal(consumption(np.zeros((1, 2)), np.zeros((1, 2))), [1.0])
    for bad in (
        np.full((1, 4), 0.5),
        [[np.nan, 0, 0]],
        [[np.inf, 0, 0]],
        [[-0.1, 0, 0]],
        [[0, 1.1, 0]],
    ):
        raises(lambda: chemistry.outputs(bad, 1.0))
    for bad in (-0.1, 1.1, np.nan, np.inf, True, None):
        raises(config(eta_shared=bad).validate)
    raises(config(rho_voc=0.5).validate)
    raises(config(B=[[1.0, 0.0], [0.0, 1.0]]).validate)
    c.genetics["alleles"][0]["effect"].append(0.0)
    raises(c.validate)


def test_shared_qtl_shapes_native_inheritance_and_chemical_audit():
    c = config(generations=1)
    genetic = GeneticState(c)
    for kind, sizes, width, count in (
        ("plant", c.plant_sizes, 3, 6),
        ("insect", c.insect_sizes, 4, 8),
    ):
        g = genetic.founders(sizes, kind)
        values, raw = genetic.phenotype(g, kind)
        assert g.shape == (sum(sizes), 2, count) and raw.shape == (sum(sizes), width)
        if kind == "plant":
            assert values["plant_chemical_traits"].shape[1] == 3
    result = simulate(c)
    assert result["metadata"]["genetics_version"] == "slim-qtl-2"
    assert result["metadata"]["plant_trait_order"] == ["x_N", "x_S", "x_V"]
    for state in result["snapshots"]:
        assert "plant_voc_independent_trait" not in state
        D, V = chemistry.outputs(state["plant_chemical_traits"], c.eta_shared)
        np.testing.assert_array_equal(D, state["D"])
        np.testing.assert_array_equal(V, state["V_emit"])
        assert state["I"].shape[1] == D.shape[1] == state["A"].shape[1] == V.shape[1] == 2
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "run"
        save_result(result, path)
        checked = audit_run(path)
        assert checked["genetic_audit"] and checked["contribution_audit"]
        assert checked["max_overrun"] == 0


def test_chemical_csv_requires_empty_unused_plant_coordinate():
    c = config()
    parameters = dict(n=2, B=None, mapping="shared_chemical_modules", genetics=c.genetics)
    locus = c.genetics["loci"][0]["locus_id"]
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "a.csv"
        path.write_text(
            f"allele_id,locus_id,effect_0,effect_1,effect_2,effect_3\n0,{locus},0.1,0,0,\n"
        )
        loaded, _ = load_inputs(dict(alleles="a.csv"), path.parent, parameters)
        assert loaded["genetics"]["alleles"][0]["effect"] == [0.1, 0.0, 0.0]
        path.write_text(path.read_text().replace("0.1,0,0,\n", "0.1,0,0,0\n"))
        raises(lambda: load_inputs(dict(alleles="a.csv"), path.parent, parameters))
