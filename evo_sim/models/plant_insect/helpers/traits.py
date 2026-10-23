"""Inherited quantitative traits, encounter, preference and feeding functions."""

import numpy as np


def map_voc(D, B):
    return np.clip(D @ B.T, 0, 1)


def encounter_candidates(pool, size, rng):
    return rng.choice(pool, min(size, len(pool)), replace=False)


def chemical_scores(V, A, c, weights):
    if c.olfaction_mode == "neutral":
        return np.full(len(V), c.b)
    distance = (V - A) ** 2 @ weights
    # Divide sequentially: sigma**2 can underflow even when a zero distance is exact.
    with np.errstate(over="ignore"):
        return c.b - (distance / c.sigma) / c.sigma / 2


def choice_probabilities(scores, tau=1.0, outside=0.0):
    scores = np.asarray(scores, float)
    if (
        np.isnan(scores).any()
        or np.isposinf(scores).any()
        or not np.isfinite(tau)
        or tau < 0
        or not np.isfinite(outside)
    ):
        raise ValueError("invalid choice parameters")
    # Center BEFORE large tau multiplication; negative infinities mean zero mass.
    with np.errstate(over="ignore", under="ignore"):
        if tau >= 1:
            logits = np.append(scores, outside / tau)
            logits = (logits - logits.max()) * tau
        else:
            logits = np.append(scores * tau if tau > 0 else np.zeros(len(scores)), outside)
            logits = logits - logits.max()
        values = np.exp(logits)
    probabilities = values / values.sum()
    return probabilities / probabilities.sum()


def draw_host(candidates, probabilities, rng):
    selected = rng.choice(len(probabilities), p=probabilities)
    return -1 if selected == len(candidates) else int(candidates[selected])


def consumption(D, I):
    ratios = np.divide(I, D, out=np.full_like(D, np.inf, dtype=float), where=D > 0)
    return np.minimum(1, ratios.min(axis=-1))
