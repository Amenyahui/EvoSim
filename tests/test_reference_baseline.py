"""The resident/mismatched profiles differ only in the intended ancestral A effects."""

from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from evo_sim.models.plant_insect import Config, simulate
from evo_sim.models.plant_insect.audit import audit_run
from evo_sim.models.plant_insect.batch import save_result
from evo_sim.models.plant_insect.helpers.genetics import GeneticState
from evo_sim.models.plant_insect.helpers.inputs import read_config
from evo_sim.models.plant_insect.helpers.voc import conditional_probabilities

CONFIGS = Path(__file__).resolve().parents[1] / "evo_sim/models/plant_insect/configs"


def test_matched_reference_is_encoded_in_inherited_alleles():
    configs = []
    for name in ("v016_recognition_mismatch", "v016_ancestral_matched"):
        loaded = read_config(CONFIGS / f"{name}.yaml")
        c = Config(**dict(loaded["parameters"], generations=1))
        c.validate()
        configs.append(c)
    original, matched = configs
    assert original.insect_offspring_lambda == matched.insect_offspring_lambda == {"H1": 4.0}
    assert original.genetics["background"] == matched.genetics["background"]
    assert original.genetics["founder_frequencies"] == matched.genetics["founder_frequencies"]
    changed = []
    for a, b in zip(original.genetics["alleles"], matched.genetics["alleles"]):
        delta = np.asarray(b["effect"]) - a["effect"]
        if np.any(delta):
            changed.append(a["allele_id"])
            np.testing.assert_allclose(delta, [0, 0, 0.2, 0])
    assert changed == [32, 33, 34, 35]
    states = [GeneticState(c) for c in configs]
    genotypes = [g.founders(c.insect_sizes, "insect") for g, c in zip(states, configs)]
    np.testing.assert_array_equal(*genotypes)
    first, _ = states[0].phenotype(genotypes[0], "insect")
    second, raw = states[1].phenotype(genotypes[1], "insect")
    np.testing.assert_array_equal(first["I"], second["I"])
    np.testing.assert_allclose(second["A"], first["A"] + [0.4, 0])
    assert ((raw >= 0) & (raw <= 1)).all()
    V = np.array([[0.62, 0.62], [0.22, 0.22]])
    p = conditional_probabilities(V, np.array([0.22, 0.62]), original, np.array([0.5, 0.5]))
    np.testing.assert_allclose(p, [0.5, 0.5])
    p = conditional_probabilities(V, V[0], matched, np.array([0.5, 0.5]))
    np.testing.assert_allclose(p[0], 1 / (1 + np.exp(-2)))
    # Native inherited offspring must preserve the intended mapping and auditable provenance.
    matched.plant_sizes, matched.insect_sizes = [10, 10], [8]
    matched.spatial.update(initial_plant_cells=[[0, 0]] * 20, initial_insect_cells=[[0, 0]] * 8)
    result = simulate(matched)
    assert result["metadata"]["status"]["status"] == "completed"
    assert result["larvae"]
    with TemporaryDirectory() as tmp:
        folder = Path(tmp) / "run"
        save_result(result, folder)
        assert audit_run(folder)["genetic_audit"]
