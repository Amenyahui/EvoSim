"""Demography mechanisms; equations preserved from v0.1.11."""

import numpy as np


def plant_resource_update(stock, consumed, capacity, spec):
    """Update resource and survival from the post-feeding capacity deficit."""
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
    deficit = 1 - post / capacity
    survival = np.exp(
        -spec["background_mortality_rate"] - spec["defoliation_mortality_rate"] * deficit
    )
    recovered = capacity - (capacity - post) * np.exp(-spec["resource_recovery_rate"])
    return post, deficit, survival, recovered
