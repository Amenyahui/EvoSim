"""Independent PSM/VOC axes, costs, QTL widths and CSV boundaries."""

from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from evo_sim.models.plant_insect import simulate
from evo_sim.models.plant_insect.audit.run import audit_run
from evo_sim.models.plant_insect.batch import save_result
from evo_sim.models.plant_insect.biology.genetics.qtl import GeneticState
from evo_sim.models.plant_insect.biology.plants import chemistry
from evo_sim.models.plant_insect.biology.plants import defense as _defense
from evo_sim.models.plant_insect.config.input_loader import load_inputs
from tests.support import config, raises


def test_independent_module_paths_and_costs():
    c = config()
    values = np.array([[0.2, 0.6, 0.5, 0.3, 0.4]])
    defense, odor = chemistry.outputs(values)
    np.testing.assert_array_equal(defense, [[0.2, 0.6]])
    np.testing.assert_array_equal(odor, [[0.3, 0.4]])
    fecundity = _defense.fecundity(defense, odor, c, values)[0]
    np.testing.assert_allclose(fecundity, [0.865])
    for axis in range(5):
        altered = values.copy()
        altered[:, axis] += 0.1
        new_defense, new_odor = chemistry.outputs(altered)
        cost = _defense.fecundity(new_defense, new_odor, c, altered)[0]
        np.testing.assert_allclose(
            fecundity - cost, 0.1 * (c.plant_cost if axis < 2 else 0 if axis == 2 else c.voc_cost)
        )
        if axis < 3:
            np.testing.assert_array_equal(odor, new_odor)
        else:
            np.testing.assert_array_equal(defense, new_defense)


def test_chemical_boundaries_and_configuration_contract():
    c = config()
    for value in (0.0, 1.0):
        values = np.full((2, 5), value)
        defense, odor = chemistry.outputs(values)
        np.testing.assert_allclose(
            _defense.fecundity(defense, odor, c, values)[0],
            max(0, 1 - 2 * (c.plant_cost + c.voc_cost) * value),
        )
    for bad in (
        np.zeros((1, 3)),
        np.full((1, 5), np.nan),
        np.full((1, 5), -0.1),
        np.full((1, 5), 1.1),
    ):
        raises(lambda: chemistry.outputs(bad))
    raises(config(rho_voc=0.5).validate)
    raises(config(B=[[1.0, 0.0], [0.0, 1.0]]).validate)
    c.genetics["alleles"][0]["effect"].append(0.0)
    raises(c.validate)


def test_qtl_shapes_native_inheritance_and_chemical_audit():
    c = config(generations=1)
    genetic = GeneticState(c)
    for kind, sizes, width, count in (
        ("plant", c.plant_sizes, 5, 10),
        ("insect", c.insect_sizes, 10, 20),
    ):
        genome = genetic.founders(sizes, kind)
        _, raw = genetic.phenotype(genome, kind)
        assert genome.shape == (sum(sizes), 2, count) and raw.shape == (sum(sizes), width)
    result = simulate(c)
    assert result["metadata"]["plant_trait_order"] == ["PSM1", "PSM2", "Z", "VOC1", "VOC2"]
    for state in result["snapshots"]:
        defense, odor = chemistry.outputs(state["plant_chemical_traits"])
        np.testing.assert_array_equal(defense, state["D"])
        np.testing.assert_array_equal(odor, state["V_emit"])
    with TemporaryDirectory() as temporary:
        path = Path(temporary) / "run"
        save_result(result, path)
        checked = audit_run(path)
        assert (
            checked["genetic_audit"]
            and checked["contribution_audit"]
            and checked["max_overrun"] == 0
        )


def test_chemical_csv_requires_empty_unused_plant_coordinates():
    c = config()
    parameters = dict(n=2, B=None, mapping="independent_psm_voc", genetics=c.genetics)
    locus = c.genetics["loci"][0]["locus_id"]
    with TemporaryDirectory() as temporary:
        path = Path(temporary) / "a.csv"
        header = "allele_id,locus_id," + ",".join(f"effect_{axis}" for axis in range(10)) + "\n"
        path.write_text(header + f"0,{locus},0.1,0,0,0,0,,,,,\n")
        loaded, _ = load_inputs(dict(alleles="a.csv"), path.parent, parameters)
        assert loaded["genetics"]["alleles"][0]["effect"] == [0.1, 0, 0, 0, 0]
        path.write_text(header + f"0,{locus},0.1,0,0,0,0,0,,,,\n")
        raises(lambda: load_inputs(dict(alleles="a.csv"), path.parent, parameters))
