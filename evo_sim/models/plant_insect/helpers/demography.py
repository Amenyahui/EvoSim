"""F07–F14: explicit demographic assumptions, not fixed-size resampling."""

import copy

import numpy as np

from . import traits as core
from . import voc


class ComputationalLimit(RuntimeError):
    pass


def density_phi(alpha, counts):
    """F07 [R3,R4]: rows are recipients; counts are pre-feeding numbers."""
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        return 1 / (1 + np.asarray(alpha, float) @ np.asarray(counts, float))


def alpha_from_K(lam, K):
    """F08 [R3,R4], ideal single species only; K is never a hard ceiling."""
    if not np.isfinite([lam, K]).all() or lam <= 1 or K <= 0:
        raise ValueError("positive ideal K and lambda>1 required")
    return (lam - 1) / K


def draw_offspring(S, F, lam, phi, survival_rng, offspring_rng):
    """F09–F12 [R1,R4]: Bernoulli survival THEN Poisson for survivors."""
    S, F, lam, phi = np.broadcast_arrays(*[np.asarray(x, float) for x in (S, F, lam, phi)])
    if S.ndim != 1 or not all(np.isfinite(x).all() for x in (S, F, lam, phi)):
        raise ValueError("finite one-dimensional demographic inputs required")
    if ((S < 0) | (S > 1) | (F < 0) | (F > 1) | (lam <= 0) | (phi < 0) | (phi > 1)).any():
        raise ValueError("S/F/phi in [0,1], lambda positive required")
    with np.errstate(over="raise", invalid="raise"):
        mu = lam * F * phi
        expected = S * mu
        variance = S * mu + S * (1 - S) * mu**2
    survived = survival_rng.random(S.shape) < S
    births = np.zeros(S.shape, dtype=np.int64)
    # No Poisson draw for dead parents; there is no hidden second factor g.
    try:
        births[survived] = offspring_rng.poisson(mu[survived])
    except ValueError as error:
        raise ComputationalLimit(
            "conditional Poisson mean exceeds integer sampler range"
        ) from error
    return dict(
        survived_to_reproduce=survived,
        offspring_count=births,
        conditional_offspring_mean_mu=mu,
        expected_offspring=expected,
        offspring_variance=variance,
    )


def inherit(traits, labels, ids, births, c, kind, rngs, next_id, generation):
    """F13–F14 [R1,R2,R4]: exactly B children, same parent and species."""
    # Caller checks the complete sampled total before allocating ANY descendants.
    parent = np.repeat(np.arange(len(ids)), births)
    child_labels = labels[parent]
    child_ids = np.arange(next_id, next_id + len(parent), dtype=np.int64)
    children = {name: values[parent].copy() for name, values in traits.items()}
    masks = {name: np.zeros(value.shape, bool) for name, value in children.items()}
    for s in range(len(rngs)):
        # Frozen taxa start homogeneous within species: no segregating variation
        # and no mutation, so selection cannot change their trait distribution.
        if not (c.evolve_plants if kind == "plant" else c.evolve_insects):
            continue
        selected = np.flatnonzero(child_labels == s)
        rng = rngs[s]
        for name in children:
            rate, sd = (
                (c.mu_Z, c.sigma_Z)
                if name == "plant_voc_independent_trait"
                else ((c.mu_A, c.sigma_A) if name == "A" else (c.mutation_rate, c.mutation_sd))
            )
            # Replay only the selection mask on a clone; actual mutations use untouched old functions.
            probe = copy.deepcopy(rng)
            masks[name][selected] = probe.random(children[name][selected].shape) < rate
            if name == "plant_voc_independent_trait":
                children[name][selected] = voc.mutate_Z(children[name][selected], c, rng)[0]
            else:
                children[name][selected] = core.mutate(children[name][selected], rate, sd, rng)
    records = []
    prefix = "P" if kind == "plant" else "H"
    for j, index in enumerate(parent):
        row = dict(
            generation=generation,
            kind=kind,
            species=f"{prefix}{child_labels[j] + 1}",
            offspring_id=int(child_ids[j]),
            parent_id=int(ids[index]),
            parent_genotype_id=f"{prefix}:{ids[index]}",
            mutation_stream=f"phase4/{kind}/species{child_labels[j]}",
        )
        for name, value in children.items():
            for k in range(value.shape[1]):
                row[f"{name}_{k + 1}"] = value[j, k]
                row[f"{name}_mutated_{k + 1}"] = bool(masks[name][j, k])
        records.append(row)
    return children, child_labels, child_ids, records, next_id + len(parent)
