"""F07–F14: explicit demographic assumptions, not fixed-size resampling."""

import numpy as np

from .spatial import distance


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


def validate_sexual(c):
    spec = c.sexual_reproduction
    fixed = dict(
        plant_mating_system="hermaphroditic_obligate_outcrossing",
        plant_selfing=False,
        insect_mating_system="dioecious",
        father_sampling="uniform_with_replacement",
    )
    probabilities = ("insect_founder_female_probability", "insect_offspring_female_probability")
    radii = ("pollen_radius", "insect_mate_radius")
    if not isinstance(spec, dict) or set(spec) != set(fixed) | set(probabilities) | set(radii) | {
        "provenance"
    }:
        raise ValueError("complete sexual_reproduction contract required; unknown keys rejected")
    for key, value in fixed.items():
        if type(spec[key]) is not type(value) or spec[key] != value:
            raise ValueError(f"unsupported {key}")
    for key in probabilities:
        x = spec[key]
        if type(x) not in (int, float) or not np.isfinite(x) or not 0 <= x <= 1:
            raise ValueError(f"{key} must be a probability")
    for key in radii:
        if type(spec[key]) is not int or not 0 <= spec[key] <= max(c.spatial["grid_shape"]):
            raise ValueError(f"{key} must be a bounded integer radius")
    provenance = spec["provenance"]
    if (
        not isinstance(provenance, dict)
        or provenance.get("source_kind") != "theoretical_placeholder"
        or not isinstance(provenance.get("description"), str)
        or not provenance["description"].strip()
        or not c.illustrative_only
    ):
        raise ValueError("sexual reproduction assumptions require theoretical provenance")


def mating_candidates(labels, cells, donor, mother, radius, shape):
    # ponytail: quadratic scan for bounded prototypes; patch indexing if profiling warrants it.
    candidates = []
    for i in range(len(labels)):
        mask = donor & (labels == labels[i]) & (distance(cells[i], cells, shape) <= radius)
        mask[i] = False
        candidates.append(np.flatnonzero(mask) if mother[i] else np.empty(0, dtype=int))
    return candidates


def maternal_moments(S, mu, candidates, mother):
    q = np.array([1 - np.prod(1 - S[pool]) for pool in candidates])
    p = S * q
    return (
        q,
        np.where(mother, p * mu, np.nan),
        np.where(mother, p * mu + p * (1 - p) * mu**2, np.nan),
    )
