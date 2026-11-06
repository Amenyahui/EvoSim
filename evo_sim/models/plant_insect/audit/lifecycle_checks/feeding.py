"""Independent feeding checks; equations preserved from v0.1.11 audit."""

import numpy as np


def module_profiles(values, c):
    """Independent saved-data reconstruction, deliberately separate from runtime kernels."""
    f = c.functional_modules
    z = np.asarray(f["chemical_features"])
    processing = np.exp(-0.5 * ((values[:, 2:4, None] - z) / f["processing_width"]) ** 2)
    processing /= f["processing_idle"] + processing.sum(axis=2, keepdims=True)
    I = np.sum(values[:, :2, None] * processing, axis=1)
    sensing = np.exp(-0.5 * ((values[:, 6:8, None] - z) / f["sensory_width"]) ** 2)
    return I, np.sum(values[:, 4:6, None] * sensing, axis=1)


def module_scores(V, x, c):
    f = c.functional_modules
    if c.olfaction_mode == "neutral":
        return np.full(len(V), c.b)
    affinity = np.exp(
        -0.5 * ((x[6:8, None] - np.asarray(f["chemical_features"])) / f["sensory_width"]) ** 2
    )
    activation = V @ affinity.T
    return c.b + (activation / (f["sensory_half"] + activation)) @ (
        x[4:6] * np.asarray(f["sensory_valence"])
    )


def module_choice(V, x, c, forced=False):
    scores = module_scores(V, x, c)
    if forced or c.choice_mode == "forced":
        logits = (scores - scores.max()) * c.tau
        mass = np.exp(logits)
        return np.append(mass / mass.sum(), 0.0)
    logits = (
        np.append(scores, c.outside_utility / c.tau)
        if c.tau >= 1
        else np.append(scores * c.tau if c.tau else np.zeros(len(scores)), c.outside_utility)
    )
    logits = (logits - logits.max()) * (c.tau if c.tau >= 1 else 1)
    mass = np.exp(logits)
    probability = mass / mass.sum()
    return probability / probability.sum()
