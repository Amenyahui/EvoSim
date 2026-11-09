"""Actual feeding damage and genetically variable compensatory leaf recovery."""

import numpy as np


def leaf_recovery_rate(tolerance, defoliation, spec):
    """Base recovery plus damage-dependent compensation, per insect generation.

    Z and actual feeding damage H are aligned dimensionless coordinates in [0,1].
    The theoretical rate r0 + rZ*Z*H responds only to current-round damage;
    it does not represent constitutive growth or persistent induced recovery.
    """
    z, h = np.asarray(tolerance, float), np.asarray(defoliation, float)
    if (
        z.shape != h.shape
        or not np.isfinite(z).all()
        or not np.isfinite(h).all()
        or ((z < 0) | (z > 1)).any()
        or ((h < 0) | (h > 1)).any()
    ):
        raise ValueError("aligned tolerance and defoliation in [0,1] required")
    return spec["resource_recovery_rate"] + spec["tolerance_recovery_rate"] * z * h


def plant_resource_update(stock, consumed, capacity, spec, tolerance):
    """Return post-feeding stock, actual damage fraction, survival and recovered stock.

    Stock, consumption and capacity share the configured edible biomass unit.
    Z is a normalized inherited tolerance coordinate; configured rates are per
    insect generation. Low initial stock alone is not herbivore damage.
    Recovery is conditional on survival; callers apply it only to retained plants.
    """
    stock, consumed, capacity = (np.asarray(a, float) for a in (stock, consumed, capacity))
    if (
        stock.shape != consumed.shape
        or stock.shape != capacity.shape
        or not all(np.isfinite(a).all() for a in (stock, consumed, capacity))
        or (capacity <= 0).any()
        or (stock < 0).any()
        or (consumed < 0).any()
        or (consumed > stock).any()
        or (stock > capacity).any()
    ):
        raise ValueError("plant stock must cover consumption and satisfy 0 <= stock <= K")
    post = stock - consumed
    defoliation = np.divide(consumed, stock, out=np.zeros_like(stock), where=stock > 0)
    survival = np.exp(
        -spec["background_mortality_rate"] - spec["defoliation_mortality_rate"] * defoliation
    )
    rate = leaf_recovery_rate(tolerance, defoliation, spec)
    recovered = capacity - (capacity - post) * np.exp(-rate)
    return post, defoliation, survival, recovered
