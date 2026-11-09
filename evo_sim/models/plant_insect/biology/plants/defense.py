"""Defense mechanisms; equations preserved from v0.1.12."""

import numpy as np


def fecundity(D, V_emit, c, chemical=None):
    defense_cost = c.plant_cost * D.sum(axis=1)
    from evo_sim.models.plant_insect.biology.plants.chemistry import outputs

    chemical = np.asarray(chemical, float)
    defense, _ = outputs(chemical)
    if not np.array_equal(defense, D):
        raise ValueError("PSM costs require aligned genetic concentrations and defense")
    voc_cost = c.voc_cost * chemical[:, 3:5].sum(axis=1)
    # Preserve exact arithmetic order of v1.1 when the extra cost is disabled.
    fp = (
        np.maximum(0, 1 - defense_cost)
        if c.voc_cost == 0
        else np.maximum(0, 1 - defense_cost - voc_cost)
    )
    return fp, defense_cost, voc_cost
