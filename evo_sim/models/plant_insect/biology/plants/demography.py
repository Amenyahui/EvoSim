"""Actual feeding damage and genetically variable compensatory leaf recovery."""

import numpy as np


def leaf_recovery_rate(tolerance, spec):
    """Map inherited normalized Z to a leaf recovery rate per insect generation."""
    return spec["resource_recovery_rate"] + spec["tolerance_recovery_rate"] * np.asarray(
        tolerance, float
    )


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
    tolerance = np.asarray(tolerance, float)
    if (
        tolerance.shape != stock.shape
        or not np.isfinite(tolerance).all()
        or ((tolerance < 0) | (tolerance > 1)).any()
    ):
        raise ValueError("aligned tolerance Z in [0,1] required")
    post = stock - consumed
    defoliation = np.divide(consumed, stock, out=np.zeros_like(stock), where=stock > 0)
    survival = np.exp(
        -spec["background_mortality_rate"] - spec["defoliation_mortality_rate"] * defoliation
    )
    rate = leaf_recovery_rate(tolerance, spec)
    recovered = capacity - (capacity - post) * np.exp(-rate)
    return post, defoliation, survival, recovered
