"""Fixed stage/taxon/species random streams, unchanged from v0.1.12."""

import numpy as np


def streams(c):
    """Fixed stage/taxon/species IDs, so demographic draws never consume ecological draws."""

    def rng(phase, kind, s):
        return np.random.default_rng(np.random.SeedSequence([c.seed, 130, phase, kind, s]))

    return {
        phase: (
            [rng(phase, 0, s) for s in range(len(c.plant_sizes))],
            [rng(phase, 1, s) for s in range(len(c.insect_sizes))],
        )
        for phase in range(10)
    }
