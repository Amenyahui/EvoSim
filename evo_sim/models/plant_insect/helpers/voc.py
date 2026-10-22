"""Inherited VOC mapping, preference probabilities and plant costs."""

import numpy as np

from . import traits


def map_voc(D, Z, B, rho):
    if rho == 1:
        return traits.map_voc(D, B)
    return np.clip(rho * (D @ B.T) + (1 - rho) * Z, 0, 1)


def choice_probabilities(scores, tau=1.0, outside=0.0, mode="outside"):
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
    if c.olfaction_mode == "neutral" or c.tau == 0:
        relative = np.zeros(len(V))
    else:
        distance = (V - A) ** 2 @ w
        # Center before sigma division so the nearest candidate stays finite.
        with np.errstate(over="ignore"):
            relative = -((distance - distance.min()) / c.sigma) / c.sigma / 2
    return choice_probabilities(relative, c.tau, mode="forced")[:-1]


def fecundity(D, V_emit, c):
    defense_cost = c.plant_cost * D.sum(axis=1)
    voc_cost = c.voc_cost * V_emit.mean(axis=1)
    # Preserve exact arithmetic order of v1.1 when the extra cost is disabled.
    fp = (
        np.maximum(0, 1 - defense_cost)
        if c.voc_cost == 0
        else np.maximum(0, 1 - defense_cost - voc_cost)
    )
    return fp, defense_cost, voc_cost


def mutate_Z(Z, c, rng):
    result = Z.copy()
    mask = rng.random(result.shape) < c.mu_Z
    changed = (result[mask] + rng.normal(0, c.sigma_Z, int(mask.sum()))) % 2
    result[mask] = np.where(changed <= 1, changed, 2 - changed)
    return result, int(mask.sum())
