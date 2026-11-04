"""Fixed, normalized shared chemical modules; no emission or storage dynamics."""

import numpy as np


def outputs(x, eta):
    x = np.asarray(x, float)
    if (
        x.ndim != 2
        or x.shape[1] != 3
        or not np.isfinite(x).all()
        or ((x < 0) | (x > 1)).any()
        or type(eta) not in (int, float)
        or not np.isfinite(eta)
        or not 0 <= eta <= 1
    ):
        raise ValueError("chemical traits require N x 3 values in [0,1] and eta in [0,1]")
    return x[:, :2].copy(), np.column_stack((eta * x[:, 1], x[:, 2]))


def plant_voc(D, independent, c, B=None, chemical=None):
    defense, V = outputs(chemical, c.eta_shared)
    if not np.array_equal(defense, D):
        raise ValueError("plant D is inconsistent with its genetic chemical modules")
    return V


def reference(c, B=None):
    return outputs(c.observation_reference_traits["plant_chemical_traits"], c.eta_shared)
