"""Actual food sharing followed by explicit post-ingestion chemical processing."""

import numpy as np

from evo_sim.models.plant_insect.biology.plants import defense


def chemical_survival(PSM, processing, intake, maximum, spec):
    """Post-ingestion loads in abstract mass units and dimensionless survival components.

    Processing is a per-window capacity, independent of demand. Antinutrition
    acts on the ingested dose independently of chemical processing. All weights
    and the nutritional exponent are explicit theoretical configuration values.
    """
    PSM, processing = np.broadcast_arrays(np.asarray(PSM, float), np.asarray(processing, float))
    intake, maximum = np.broadcast_arrays(np.asarray(intake, float), np.asarray(maximum, float))
    if (
        PSM.shape != intake.shape + (2,)
        or not np.isfinite(PSM).all()
        or not np.isfinite(processing).all()
        or (PSM < 0).any()
        or (processing < 0).any()
    ):
        raise ValueError("two aligned finite nonnegative chemical axes required")
    if (
        not np.isfinite(intake).all()
        or not np.isfinite(maximum).all()
        or (maximum <= 0).any()
        or (intake < 0).any()
        or (intake > maximum).any()
    ):
        raise ValueError("intake must satisfy 0 <= intake <= positive maximum")
    dose = intake[..., None] * PSM * np.asarray(spec["dose_scale"])
    capacity = processing * np.asarray(spec["processing_scale"])
    if not np.isfinite(dose).all() or not np.isfinite(capacity).all():
        raise ValueError("chemical load exceeds floating-point range")
    processed = np.minimum(dose, capacity)
    residual = dose - processed
    antinutrient = dose.copy()
    g = intake / maximum
    nutrition = np.zeros_like(g)
    np.power(g, spec["nutrition_exponent"], out=nutrition, where=g > 0)
    with np.errstate(over="ignore"):
        toxicity = np.exp(-(residual * np.asarray(spec["toxicity_weights"])).sum(axis=-1))
        antinutrition = np.exp(
            -(antinutrient * np.asarray(spec["antinutrition_weights"])).sum(axis=-1)
        )
    return dict(
        dose=dose,
        processing_capacity=capacity,
        processed=processed,
        residual=residual,
        antinutrient=antinutrient,
        nutrition_survival=nutrition,
        toxicity_survival=toxicity,
        antinutrition_survival=antinutrition,
        survival_probability=nutrition * toxicity * antinutrition,
    )


def founder_feeding_diagnostics(D, I, capacity, dmax, spec):
    """Counterfactual single-larva exposure, without encounters, competition or RNG."""
    probabilities = []
    for processing, maximum in zip(I, dmax, strict=True):
        intake = np.minimum(capacity, maximum)
        probabilities.extend(
            chemical_survival(D, processing, intake, np.full(len(D), maximum), spec)[
                "survival_probability"
            ]
        )
    return dict(
        scope="initial_single_larva_pairs_without_choice",
        pair_count=len(probabilities),
        status="counterfactual_initial_pairs" if probabilities else "no_initial_pairs",
        survival_mean=float(np.mean(probabilities)) if probabilities else None,
    )


def allocate_resources(host, capacity, dmax):
    """Proportional sharing of actual food demand; chemical traits do not affect demand."""
    host = np.asarray(host)
    capacity, dmax = (np.asarray(x, dtype=float) for x in (capacity, dmax))
    if (
        host.ndim != 1
        or host.dtype.kind not in "iu"
        or dmax.shape != host.shape
        or capacity.ndim != 1
    ):
        raise ValueError("host must be integer vector; dmax matches it; capacity is a plant vector")
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
    demand = np.zeros(len(host))
    demand[found] = dmax[found]
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
    allocation = allocate_resources(host, capacity, dmax)
    concentrations = np.zeros((len(I), 2))
    concentrations[found] = D[host[found]]
    exposure = chemical_survival(concentrations, I, allocation["x"], dmax, c.chemical_processing)
    unlimited = chemical_survival(
        D[host[found]], I[found], dmax[found], dmax[found], c.chemical_processing
    )["survival_probability"]
    q = np.full(len(I), np.nan)
    q[found] = unlimited
    chemical = np.asarray(chemical)
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
        **exposure,
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
