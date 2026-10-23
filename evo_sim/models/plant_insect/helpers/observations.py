"""Deterministic observations on a fixed, equally available ancestral panel."""

import numpy as np

from .traits import consumption
from .voc import conditional_probabilities, map_voc


def standardized(state, c, B, w, generation):
    """Each reference species has one equal-weight option; no ecological RNG."""
    D = np.asarray(c.observation_reference_traits["D"], float)
    Z = np.asarray(c.observation_reference_traits["plant_voc_independent_trait"], float)
    V = map_voc(D, Z, B, c.rho_voc)
    probability = np.empty((len(state.insects), len(D)))
    for i, preference in enumerate(state.insects.preference):
        probability[i] = conditional_probabilities(V, preference, c, w)
    q = consumption(
        np.broadcast_to(D, (len(state.insects), *D.shape)),
        state.insects.resistance[:, None, :],
    )
    rows = []
    for h in range(len(c.insect_sizes)):
        mask = state.insects.species == h
        n = int(mask.sum())
        for p in range(len(D)):
            rows.append(
                dict(
                    generation=generation,
                    insect_species=f"H{h + 1}",
                    reference_plant_species=f"P{p + 1}",
                    insect_count=n,
                    preference=float(probability[mask, p].mean()) if n else np.nan,
                    q_unlimited=float(q[mask, p].mean()) if n else np.nan,
                    panel="equal_ancestral_species",
                )
            )
    return rows
