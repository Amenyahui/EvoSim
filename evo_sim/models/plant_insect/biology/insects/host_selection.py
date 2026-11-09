"""Host selection mechanisms; equations preserved from v0.1.12."""

import numpy as np

from evo_sim.models.plant_insect.biology.insects import sensory


def encounter_candidates(pool, size, rng):
    return rng.choice(pool, min(size, len(pool)), replace=False)


def draw_host(candidates, probabilities, rng):
    selected = rng.choice(len(probabilities), p=probabilities)
    return -1 if selected == len(candidates) else int(candidates[selected])


def contact_acceptance(PSM, traits, spec):
    """Adult female contact-response trait during oviposition, not a GR gene model.

    Each plant has a conditional Bernoulli acceptance probability after contact.
    These probabilities are not normalized across plants. G does not control
    larval host rejection or feeding in the current model.
    """
    PSM, traits = np.asarray(PSM, float), np.asarray(traits, float)
    if (
        PSM.ndim != 2
        or PSM.shape[1] != 2
        or traits.shape != (10,)
        or not np.isfinite(PSM).all()
        or not np.isfinite(traits).all()
        or ((PSM < 0) | (PSM > 1)).any()
        or ((traits < 0) | (traits > 1)).any()
    ):
        raise ValueError("normalized PSM and ten insect coordinates required")
    weights = 2 * traits[8:10] - 1
    utility = spec["intercept"] + spec["slope"] * (PSM @ weights)
    return np.exp(-np.logaddexp(0, -utility))


def choose_hosts(A, p_labels, h_labels, V_cue, c, w, rng_p, rng_H, candidate_pools=None, PSM=None):
    """The pool is fixed during choice: acceptance never removes a plant."""
    host = np.full(len(A), -1, int)
    reason = np.full(len(A), "rejected", dtype="<U16")
    rank = np.empty(len(A), int)
    candidate_count = np.zeros(len(A), int)
    exposures = np.zeros((len(c.insect_sizes), len(p_labels)), int)
    candidate_share = np.full((len(A), len(c.plant_sizes)), np.nan)
    conditional_share = np.full_like(candidate_share, np.nan)
    pool = np.arange(len(p_labels))
    recorded_candidates = [np.empty(0, int) for _ in A] if candidate_pools is not None else None
    if PSM is None or np.asarray(PSM).shape != (len(p_labels), 2):
        raise ValueError("plant PSM concentrations required for contact checks")
    contact_rows = [[] for _ in A]
    contact_probabilities = [[] for _ in A]
    contact_accepted = [[] for _ in A]
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
        for s in range(len(c.plant_sizes)):
            mask = p_labels[candidates] == s
            candidate_share[i, s] = mask.mean()
            conditional_share[i, s] = conditional[mask].sum()
        remaining = candidates.copy()
        while len(remaining):
            scores = sensory.chemical_scores(V_cue[remaining], A[i], c, w)
            probabilities = sensory.choice_probabilities(
                scores, c.tau, c.outside_utility, c.choice_mode
            )
            selected = rng.choice(len(probabilities), p=probabilities)
            if selected == len(remaining):
                reason[i] = "outside_option"
                break
            plant = int(remaining[selected])
            probability = float(contact_acceptance(PSM[[plant]], A[i], c.contact_response)[0])
            accepted = bool(rng.random() < probability)
            contact_rows[i].append(plant)
            contact_probabilities[i].append(probability)
            contact_accepted[i].append(accepted)
            if accepted:
                host[i] = plant
                break
            remaining = np.delete(remaining, selected)
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
        contact_rows=contact_rows,
        contact_probabilities=contact_probabilities,
        contact_accepted=contact_accepted,
    )
    if recorded_candidates is not None:
        result["candidate_rows"] = recorded_candidates
    return result
