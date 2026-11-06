"""Feeding mechanisms; equations preserved from v0.1.11."""

import numpy as np

from evo_sim.models.plant_insect.biology.plants import defense


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


def larval_survival_probability(g):
    """Baseline survival is x/dmax, with no maternal expression cost applied again."""
    probability = np.asarray(g, dtype=float)
    if not np.isfinite(probability).all() or ((probability < 0) | (probability > 1)).any():
        raise ValueError("larval resource fraction must be finite and in [0,1]")
    return probability


def founder_feeding_diagnostics(D, I, capacity, dmax):
    """Probe validated internal founder arrays, without encounters, competition or RNG."""
    D, I, capacity, dmax = (np.asarray(value, float) for value in (D, I, capacity, dmax))
    count = len(D) * len(I)
    result = dict(scope="all_initial_pairs_without_encounter_or_choice", pair_count=count)
    if not count:
        return dict(result, status="no_initial_pairs")
    delta, tolerance = 1e-4, 1e-12
    perturbations = []
    for axis in range(D.shape[1]):
        lower, upper = D.copy(), D.copy()
        lower[:, axis] = np.maximum(0, lower[:, axis] - delta)
        upper[:, axis] = np.minimum(1, upper[:, axis] + delta)
        perturbations.append((lower, upper))
    q_min, q_max, q_sum, unsaturated, limited = 1.0, 0.0, 0.0, 0, 0
    q_sensitive = np.zeros(D.shape[1], dtype=np.int64)
    g_sensitive = np.zeros_like(q_sensitive)
    resource_masked = np.zeros_like(q_sensitive)
    # ponytail: one insect against all plants bounds memory by the initial plant count.
    for response, demand in zip(I, dmax, strict=True):
        q = consumption(D, response)
        q_min, q_max = min(q_min, float(q.min())), max(q_max, float(q.max()))
        q_sum += float(q.sum())
        unsaturated += int(np.count_nonzero(q < 1 - tolerance))
        with np.errstate(over="ignore"):
            resource_ceiling = capacity / demand
        limited += int(np.count_nonzero(resource_ceiling < q))
        for axis, (lower, upper) in enumerate(perturbations):
            q_lower = consumption(lower, response)
            q_upper = consumption(upper, response)
            q_changed = q_lower - q_upper > tolerance
            g_changed = (
                np.minimum(q_lower, resource_ceiling) - np.minimum(q_upper, resource_ceiling)
                > tolerance
            )
            q_sensitive[axis] += np.count_nonzero(q_changed)
            g_sensitive[axis] += np.count_nonzero(g_changed)
            resource_masked[axis] += np.count_nonzero(q_changed & ~g_changed)
    return dict(
        result,
        status="counterfactual_initial_pairs",
        q_min=q_min,
        q_max=q_max,
        q_mean=q_sum / count,
        q_unsaturated_fraction=unsaturated / count,
        defense_delta=delta,
        response_tolerance=tolerance,
        q_sensitive_fraction=(q_sensitive / count).tolist(),
        single_larva_g_sensitive_fraction=(g_sensitive / count).tolist(),
        single_larva_resource_masked_q_sensitive_fraction=(resource_masked / count).tolist(),
        single_larva_resource_limited_fraction=limited / count,
    )


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
        or (capacity < 0).any()
        or not np.isfinite(dmax).all()
        or (dmax <= 0).any()
    ):
        raise ValueError("resource must be nonnegative; demand must be positive; both finite")
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
    excess = consumed > capacity
    if excess.any():
        # ponytail: bounded common-factor search only for hosts still over budget.
        low = np.zeros(len(capacity))
        high = np.ones(len(capacity))
        for _ in range(64):
            middle = (low + high) / 2
            trial = np.bincount(
                host[found], weights=x[found] * middle[host[found]], minlength=len(capacity)
            )
            feasible = trial <= capacity
            low[excess & feasible] = middle[excess & feasible]
            high[excess & ~feasible] = middle[excess & ~feasible]
        correction = np.ones(len(capacity))
        correction[excess] = low[excess]
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
        h=np.divide(
            consumed, capacity, out=np.zeros_like(consumed, dtype=float), where=capacity > 0
        ),
        herbivore_load=load,
        resource_limited=limited,
        resource_exhausted=(consumed >= capacity * (1 - 1e-12)),
    )


def settle_feeding(
    D, I, p_labels, h_labels, V_emit, host, c, capacity, dmax, chemical=None, module_traits=None
):
    """Return larval survival plus separate maternal output and plant fitness diagnostics."""
    found = host >= 0
    q = np.full(len(I), np.nan)
    q[found] = consumption(D[host[found]], I[found])
    allocation = allocate_resources(host, q, capacity, dmax)
    fp, defense_cost, voc_cost = defense.fecundity(D, V_emit, c, chemical)
    from evo_sim.models.plant_insect.biology.traits import maternal_fecundity

    if module_traits is None:
        raise ValueError("module feeding diagnostics require original expression/tuning traits")
    fi = maternal_fecundity(module_traits, c)[0]
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
        survival_probability=larval_survival_probability(allocation["g"]),
        q=q,
        fp=fp,
        fi=fi,
        defense_cost=defense_cost,
        voc_cost=voc_cost,
        same_species_others=same_load,
        other_species_others=other_load,
        same_species_codemand=same_demand,
        other_species_codemand=other_demand,
    )
