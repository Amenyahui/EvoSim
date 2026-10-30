"""Inherited VOC mapping, preference probabilities and plant costs."""

import numpy as np

from . import traits


def map_voc(D, Z, B, rho, bounds):
    if rho == 1:
        return traits.map_voc(D, B, bounds)
    return np.clip(rho * (D @ B.T) + (1 - rho) * Z, *bounds)


def choice_probabilities(scores, tau, outside, mode):
    if mode == "outside":
        return traits.choice_probabilities(scores, tau, outside)
    scores = np.asarray(scores, float)
    if (
        mode != "forced"
        or scores.ndim != 1
        or not len(scores)
        or np.isnan(scores).any()
        or np.isposinf(scores).any()
        or not np.isfinite(tau)
        or tau < 0
    ):
        raise ValueError("invalid forced choice inputs")
    if tau == 0:
        return np.append(np.ones(len(scores)) / len(scores), 0.0)
    if not np.isfinite(scores).any():
        raise ValueError("forced choice needs a finite maximum score")
    with np.errstate(over="ignore", under="ignore"):
        values = np.exp((scores - scores.max()) * tau)
    return np.append(values / values.sum(), 0.0)


def conditional_probabilities(V, A, c, w):
    if c.genetics["schema_version"] == "slim-qtl-4":
        return choice_probabilities(
            traits.chemical_scores(V, A, c, w), c.tau, c.outside_utility, mode="forced"
        )[:-1]
    if c.olfaction_mode == "neutral" or c.tau == 0:
        relative = np.zeros(len(V))
    else:
        distance = (V - A) ** 2 @ w
        # Center before sigma division so the nearest candidate stays finite.
        with np.errstate(over="ignore"):
            relative = -((distance - distance.min()) / c.sigma) / c.sigma / 2
    return choice_probabilities(relative, c.tau, c.outside_utility, mode="forced")[:-1]


def fecundity(D, V_emit, c, chemical=None):
    defense_cost = c.plant_cost * D.sum(axis=1)
    if c.mapping == "shared_chemical_modules":
        from .chemistry import outputs

        chemical = np.asarray(chemical, float)
        defense, _ = outputs(chemical, c.eta_shared)
        if not np.array_equal(defense, D):
            raise ValueError("shared chemical costs require aligned genetic modules and defense")
        voc_cost = c.voc_cost * chemical[:, 2]
    else:
        voc_cost = c.voc_cost * V_emit.mean(axis=1)
    # Preserve exact arithmetic order of v1.1 when the extra cost is disabled.
    fp = (
        np.maximum(0, 1 - defense_cost)
        if c.voc_cost == 0
        else np.maximum(0, 1 - defense_cost - voc_cost)
    )
    return fp, defense_cost, voc_cost
