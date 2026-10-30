"""Abstract expression/tuning modules; every numerical response parameter is a model assumption."""

import numpy as np


def enabled(c):
    return c.genetics["schema_version"] == "slim-qtl-4"


def validate(spec, *, historical_flowering=False, historical_audit=False):
    keys = {
        "chemical_features",
        "processing_width",
        "sensory_width",
        "processing_idle",
        "sensory_half",
        "sensory_valence",
        "sensory_cost",
    }
    if historical_flowering:
        keys.add("flowering_kappa")
    if historical_audit:
        keys.add("adult_settlement_enabled")
    if not isinstance(spec, dict) or set(spec) != keys:
        raise ValueError(
            "functional_modules requires the current module contract; settlement removed"
        )
    for key in ("chemical_features", "sensory_valence"):
        values = np.asarray(spec[key])
        if values.shape != (2,) or values.dtype.kind not in "iuf" or not np.isfinite(values).all():
            raise ValueError(f"{key} requires two finite numbers")
    if (
        (np.asarray(spec["chemical_features"]) < 0) | (np.asarray(spec["chemical_features"]) > 1)
    ).any():
        raise ValueError("chemical_features must be in [0,1]")
    for key in (
        "processing_width",
        "sensory_width",
        "processing_idle",
        "sensory_half",
        "sensory_cost",
    ):
        value = spec[key]
        if type(value) not in (int, float) or not np.isfinite(value) or value < 0:
            raise ValueError(f"{key} must be finite and nonnegative")
        if (
            key in ("processing_width", "sensory_width", "processing_idle", "sensory_half")
            and value == 0
        ):
            raise ValueError(f"{key} must be positive")
    if historical_flowering:
        kappa = spec["flowering_kappa"]
        if type(kappa) not in (int, float) or not np.isfinite(kappa) or not 0 <= kappa <= 300:
            raise ValueError("historical flowering_kappa must be finite and in [0,300]")
    if historical_audit and type(spec["adult_settlement_enabled"]) is not bool:
        raise ValueError("adult_settlement_enabled must be boolean")


def profile(tuning, features, width):
    with np.errstate(over="ignore", under="ignore"):
        return np.exp(-0.5 * ((tuning[..., None] - np.asarray(features)) / width) ** 2)


def effective(values, spec):
    """Conserved processing capacity and sensory expression profile on the fixed catalog."""
    values = np.asarray(values, float)
    if values.ndim != 2 or values.shape[1] != 8 or not np.isfinite(values).all():
        raise ValueError("insect module traits require N x 8 finite coordinates")
    # Assumption: two fixed-width Gaussian modules allocate a fixed expression capacity,
    # including an idle channel; this allocation rule is not an empirical law.
    processing = profile(values[:, 2:4], spec["chemical_features"], spec["processing_width"])
    processing /= spec["processing_idle"] + processing.sum(axis=2, keepdims=True)
    capacity = (values[:, :2, None] * processing).sum(axis=1)
    sensory = profile(values[:, 6:8], spec["chemical_features"], spec["sensory_width"])
    return capacity, (values[:, 4:6, None] * sensory).sum(axis=1)


def scores(V, values, c):
    """Receptor activation saturates before expression gain and fixed neural valence."""
    if c.olfaction_mode == "neutral":
        return np.full(len(V), c.b)
    spec = c.functional_modules
    tuning = profile(np.asarray(values)[6:8], spec["chemical_features"], spec["sensory_width"])
    activation = np.asarray(V) @ tuning.T
    response = activation / (spec["sensory_half"] + activation)
    # Assumption: expression is gain after saturation; neural valence and tuning width are fixed.
    return c.b + response @ (np.asarray(spec["sensory_valence"]) * np.asarray(values)[4:6])


def maternal_fecundity(values, c):
    # Assumption: maternal expression costs are additive and linear, not an empirical law.
    processing_cost = c.insect_cost * values[:, :2].sum(axis=1)
    sensory_cost = c.functional_modules["sensory_cost"] * values[:, 4:6].sum(axis=1)
    return np.maximum(0, 1 - processing_cost - sensory_cost), processing_cost, sensory_cost
