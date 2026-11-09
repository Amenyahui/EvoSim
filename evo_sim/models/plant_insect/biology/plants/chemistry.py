"""Independent normalized leaf PSM and VOC coordinates, with a separate tolerance axis."""

import numpy as np


def outputs(x):
    x = np.asarray(x, float)
    if x.ndim != 2 or x.shape[1] != 5 or not np.isfinite(x).all() or ((x < 0) | (x > 1)).any():
        raise ValueError("plant traits require N x 5 values in [0,1]")
    return x[:, :2].copy(), x[:, 3:5].copy()


def plant_voc(D, independent, c, B=None, chemical=None):
    defense, V = outputs(chemical)
    if not np.array_equal(defense, D):
        raise ValueError("plant D is inconsistent with its genetic chemical modules")
    return V


def reference(c, B=None):
    return outputs(c.observation_reference_traits["plant_chemical_traits"])
