"""Sensory mechanisms; equations preserved from v0.1.12."""

import numpy as np


def chemical_scores(V, A, c, weights):
    from evo_sim.models.plant_insect.biology.traits import scores

    return scores(V, A, c)


def outside_choice_probabilities(scores, tau, outside):
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


def choice_probabilities(scores, tau, outside, mode):
    if mode == "outside":
        return outside_choice_probabilities(scores, tau, outside)
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
    return choice_probabilities(
        chemical_scores(V, A, c, w), c.tau, c.outside_utility, mode="forced"
    )[:-1]
