"""Host selection mechanisms; equations preserved from v0.1.11."""

import numpy as np

from evo_sim.models.plant_insect.biology.insects import sensory


def encounter_candidates(pool, size, rng):
    return rng.choice(pool, min(size, len(pool)), replace=False)


def draw_host(candidates, probabilities, rng):
    selected = rng.choice(len(probabilities), p=probabilities)
    return -1 if selected == len(candidates) else int(candidates[selected])


def choose_hosts(A, p_labels, h_labels, V_cue, c, w, rng_p, rng_H, candidate_pools=None):
    """The pool is fixed during choice: acceptance never removes a plant."""
    host = np.full(len(A), -1, int)
    reason = np.full(len(A), "rejected", dtype="<U12")
    rank = np.empty(len(A), int)
    candidate_count = np.zeros(len(A), int)
    exposures = np.zeros((len(c.insect_sizes), len(p_labels)), int)
    candidate_share = np.full((len(A), len(c.plant_sizes)), np.nan)
    conditional_share = np.full_like(candidate_share, np.nan)
    pool = np.arange(len(p_labels))
    recorded_candidates = [np.empty(0, int) for _ in A] if candidate_pools is not None else None
    for order, i in enumerate(rng_p.permutation(len(A))):
        rank[i] = order
        rng = rng_H[h_labels[i]]
        available = pool if candidate_pools is None else candidate_pools[i]
        candidates = encounter_candidates(available, c.candidate_size, rng)
        candidate_count[i] = len(candidates)
        if recorded_candidates is not None:
            recorded_candidates[i] = candidates.copy()
        if not len(candidates):
            reason[i] = "no_candidate"
            continue
        exposures[h_labels[i], candidates] += 1
        conditional = sensory.conditional_probabilities(V_cue[candidates], A[i], c, w)
        if c.choice_mode == "forced":
            probabilities = np.append(conditional, 0.0)
        else:
            scores = sensory.chemical_scores(V_cue[candidates], A[i], c, w)
            probabilities = sensory.choice_probabilities(
                scores, c.tau, c.outside_utility, c.choice_mode
            )
        for s in range(len(c.plant_sizes)):
            mask = p_labels[candidates] == s
            candidate_share[i, s] = mask.mean()
            conditional_share[i, s] = conditional[mask].sum()
        host[i] = draw_host(candidates, probabilities, rng)
        if host[i] >= 0:
            reason[i] = "accepted"
    result = dict(
        host=host,
        host_found=host >= 0,
        reason=reason,
        choice_order=rank,
        candidate_count=candidate_count,
        candidate_exposure=exposures,
        candidate_share=candidate_share,
        conditional_share=conditional_share,
    )
    if recorded_candidates is not None:
        result["candidate_rows"] = recorded_candidates
    return result
