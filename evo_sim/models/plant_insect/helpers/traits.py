"""Inherited quantitative traits, encounter, preference and feeding functions."""

import numpy as np


def map_voc(D, B, bounds):
    return np.clip(D @ B.T, *bounds)


def encounter_candidates(pool, size, rng):
    return rng.choice(pool, min(size, len(pool)), replace=False)


def chemical_scores(V, A, c, weights):
    if c.genetics["schema_version"] == "slim-qtl-4":
        from .modules import scores

        return scores(V, A, c)
    if c.olfaction_mode == "neutral":
        return np.full(len(V), c.b)
    distance = (V - A) ** 2 @ weights
    # Divide sequentially: sigma**2 can underflow even when a zero distance is exact.
    with np.errstate(over="ignore"):
        return c.b - (distance / c.sigma) / c.sigma / 2


def choice_probabilities(scores, tau, outside):
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
    """Most limiting positive-defense component; absent defenses impose no limit."""
    D, I = np.broadcast_arrays(np.asarray(D, float), np.asarray(I, float))
    if D.ndim == 0 or D.shape[-1] == 0:
        raise ValueError("defense and processing require a nonempty component axis")
    if not np.isfinite(D).all() or not np.isfinite(I).all() or (D < 0).any() or (I < 0).any():
        raise ValueError("defense and processing must be finite and nonnegative")
    # A finite capacity divided by a tiny defense may overflow: q still saturates at 1.
    with np.errstate(over="ignore"):
        ratios = np.divide(I, D, out=np.full_like(D, np.inf, dtype=float), where=D > 0)
    return np.minimum(1, ratios.min(axis=-1))
