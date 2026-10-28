"""Simultaneous host choices and a conserved shared resource ledger."""

import numpy as np

from . import traits, voc


def allocate_resources(host, q, capacity, dmax):
    """Deterministic proportional sharing. Rejected q is NaN, its demand/x/g are zero."""
    host = np.asarray(host)
    q, capacity, dmax = (np.asarray(x, dtype=float) for x in (q, capacity, dmax))
    if (
        host.ndim != 1
        or host.dtype.kind not in "iu"
        or q.shape != host.shape
        or dmax.shape != host.shape
        or capacity.ndim != 1
    ):
        raise ValueError("host must be integer vector; q/dmax match it; capacity is a plant vector")
    if (
        not np.isfinite(capacity).all()
        or (capacity <= 0).any()
        or not np.isfinite(dmax).all()
        or (dmax <= 0).any()
    ):
        raise ValueError("resource and demand units must be positive and finite")
    if ((host < -1) | (host >= len(capacity))).any():
        raise ValueError("host indices outside [-1,n_plants)")
    found = host >= 0
    if (
        not np.isfinite(q[found]).all()
        or ((q[found] < 0) | (q[found] > 1)).any()
        or not np.isnan(q[~found]).all()
    ):
        raise ValueError("accepted q must be in [0,1]; unmatched q must be NaN")
    demand = np.zeros(len(host))
    demand[found] = dmax[found] * q[found]
    with np.errstate(over="ignore"):
        total = np.bincount(host[found], weights=demand[found], minlength=len(capacity))
    if not np.isfinite(total).all():
        raise ValueError("total feeding demand exceeds floating-point range")
    gamma = np.ones(len(capacity))
    limited = total > capacity
    gamma[limited] = capacity[limited] / total[limited]
    x = np.zeros(len(host))
    x[found] = gamma[host[found]] * demand[found]
    consumed = np.bincount(host[found], weights=x[found], minlength=len(capacity))
    # Roundoff correction is shared by all consumers on the same plant, never first-come.
    for _ in range(3):
        excess = consumed > capacity
        if not excess.any():
            break
        correction = np.ones(len(capacity))
        correction[excess] = np.nextafter(capacity[excess] / consumed[excess], 0.0)
        gamma *= correction
        x[found] *= correction[host[found]]
        consumed = np.bincount(host[found], weights=x[found], minlength=len(capacity))
    if (consumed > capacity).any() or not np.allclose(
        consumed, np.minimum(total, capacity), rtol=1e-12, atol=0
    ):
        raise ValueError("resource allocation cannot be represented within the budget")
    load = np.bincount(host[found], minlength=len(capacity))
    return dict(
        d=demand,
        x=x,
        g=x / dmax,
        total_demand=total,
        gamma=gamma,
        X=consumed,
        h=consumed / capacity,
        herbivore_load=load,
        resource_limited=limited,
        resource_exhausted=(consumed >= capacity * (1 - 1e-12)),
    )


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
        candidates = traits.encounter_candidates(available, c.candidate_size, rng)
        candidate_count[i] = len(candidates)
        if recorded_candidates is not None:
            recorded_candidates[i] = candidates.copy()
        if not len(candidates):
            reason[i] = "no_candidate"
            continue
        exposures[h_labels[i], candidates] += 1
        conditional = voc.conditional_probabilities(V_cue[candidates], A[i], c, w)
        if c.choice_mode == "forced":
            probabilities = np.append(conditional, 0.0)
        else:
            scores = traits.chemical_scores(V_cue[candidates], A[i], c, w)
            probabilities = voc.choice_probabilities(
                scores, c.tau, c.outside_utility, c.choice_mode
            )
        for s in range(len(c.plant_sizes)):
            mask = p_labels[candidates] == s
            candidate_share[i, s] = mask.mean()
            conditional_share[i, s] = conditional[mask].sum()
        host[i] = traits.draw_host(candidates, probabilities, rng)
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


def settle_feeding(D, I, p_labels, h_labels, V_emit, host, c, capacity, dmax, chemical=None):
    """Fitness uses the same x ledger on both sides; species only annotate competition."""
    found = host >= 0
    q = np.full(len(I), np.nan)
    q[found] = traits.consumption(D[host[found]], I[found])
    allocation = allocate_resources(host, q, capacity, dmax)
    fp, defense_cost, voc_cost = voc.fecundity(D, V_emit, c, chemical)
    fi = np.maximum(0, 1 - c.insect_cost * I.sum(axis=1))
    wp = (c.sP0 + (1 - c.sP0) * (1 - allocation["h"])) * fp
    wi = np.zeros(len(I))
    wi[found] = (c.sI0 + (1 - c.sI0) * allocation["g"][found]) * fi[found]
    loads = np.zeros((len(c.insect_sizes), len(D)), int)
    demands = np.zeros_like(loads, dtype=float)
    np.add.at(loads, (h_labels[found], host[found]), 1)
    np.add.at(demands, (h_labels[found], host[found]), allocation["d"][found])
    same_load = np.zeros(len(I), int)
    other_load = np.zeros(len(I), int)
    same_demand = np.zeros(len(I))
    other_demand = np.zeros(len(I))
    same_load[found] = loads[h_labels[found], host[found]] - 1
    other_load[found] = (
        allocation["herbivore_load"][host[found]] - loads[h_labels[found], host[found]]
    )
    same_demand[found] = np.maximum(
        0, demands[h_labels[found], host[found]] - allocation["d"][found]
    )
    other_demand[found] = np.maximum(
        0, allocation["total_demand"][host[found]] - demands[h_labels[found], host[found]]
    )
    return dict(
        allocation,
        q=q,
        wp=wp,
        wi=wi,
        fp=fp,
        fi=fi,
        defense_cost=defense_cost,
        voc_cost=voc_cost,
        same_species_others=same_load,
        other_species_others=other_load,
        same_species_codemand=same_demand,
        other_species_codemand=other_demand,
    )
