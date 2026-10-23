"""Ecological resource conservation retained through the genetic migration."""

import numpy as np

from evo_sim.models.plant_insect.helpers import runtime as model


def test_T15_resource_conservation_and_order():
    rng = np.random.default_rng(901)
    host = rng.integers(0, 7, 1000)
    q, capacity, demand = (
        rng.uniform(0, 1, 1000),
        rng.uniform(0.01, 1, 7),
        rng.uniform(0.1, 2, 1000),
    )
    a = model.shared.allocate_resources(host, q, capacity, demand)
    assert (a["x"] >= 0).all() and (a["x"] <= a["d"]).all()
    assert (a["X"] <= capacity).all()
    np.testing.assert_allclose(a["X"], np.minimum(a["total_demand"], capacity))
    order = rng.permutation(len(host))
    b = model.shared.allocate_resources(host[order], q[order], capacity, demand[order])
    np.testing.assert_allclose(a["x"][order], b["x"])
