from pathlib import Path

import numpy as np

from evo_sim.core.config import merge_config
from evo_sim.models.plant_insect.helpers.config import Config
from evo_sim.models.plant_insect.helpers.inputs import read_config


def test_ten_plant_type_design():
    path = Path(__file__).parents[1] / "evo_sim/models/plant_insect/configs/ten_plant_types.yaml"
    design = read_config(path)
    params = dict(design["parameters"], generations=30)
    cases = design["batch"]["cases"]
    assert design["simulation_metadata"]["num_steps_per_episode"] == 30
    assert design["batch"]["replicates"] == 4
    assert [case["name"] for case in cases] == ["ten_types", "ten_types_no_insects"]
    for case in cases:
        Config(**merge_config(params, case["overrides"])).validate()
    assert params["plant_sizes"] == [20] * 10 and params["insect_sizes"] == [200]
    assert cases[1]["overrides"]["insect_sizes"] == [0]
    assert set(params["resource_capacity"].values()) == {10.0}
    assert params["max_feeding_demand"] == {"H1": 1.0}
    assert params["seed"] == 20261102
    assert params["genetics"]["mutation_enabled"] == {"plant": False, "insect": False}
    assert "plant_selfing" not in params["sexual_reproduction"]
    assert "pollen_radius" not in params["sexual_reproduction"]
    assert params["spatial"]["grid_shape"] == [10, 10]
    assert params["plant_demography"]["microsites_per_cell"] == 3

    genetics = params["genetics"]
    alleles = genetics["alleles"]
    centers = []
    for i in range(1, 11):
        loci = genetics["founder_frequencies"][f"P{i}"]
        mean = np.zeros(3)
        for entries in loci.values():
            assert len(entries) == 2 and all(e["frequency"] == 0.5 for e in entries)
            effects = [np.array(alleles[e["allele_id"]]["effect"]) for e in entries]
            assert not np.array_equal(*effects)
            mean += sum(effects) / 2
        centers.append(mean)
    assert len({tuple(np.round(center, 8)) for center in centers}) == 10
    panel = params["observation_reference_traits"]
    assert panel["reference_ids"] == [f"P{i}" for i in range(1, 11)]
    assert np.allclose(centers, panel["plant_chemical_traits"])
