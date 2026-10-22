"""Patch grid, local encounters, local competition and bounded dispersal."""

import numpy as np


def validate_spatial(c):
    s = c.spatial
    required = {
        "grid_shape",
        "boundary",
        "search_radius",
        "competition_radius",
        "plant_dispersal_radius",
        "insect_dispersal_radius",
        "initial_plant_cells",
    }
    required.add(
        "initial_host_by_insect_species" if c.model_version == "1.5" else "initial_insect_cells"
    )
    if not isinstance(s, dict) or set(s) != required:
        raise ValueError(f"spatial requires exactly {sorted(required)}")
    shape = s["grid_shape"]
    if (
        not isinstance(shape, (tuple, list))
        or len(shape) != 2
        or any(type(x) is not int or x < 1 for x in shape)
    ):
        raise ValueError("grid_shape requires two positive integers")
    if s["boundary"] != "periodic":
        raise ValueError("Only explicit periodic boundaries are implemented")
    for name in (
        "search_radius",
        "competition_radius",
        "plant_dispersal_radius",
        "insect_dispersal_radius",
    ):
        if type(s[name]) is not int or not 0 <= s[name] <= max(shape):
            raise ValueError(f"{name} must be an integer between 0 and max(grid_shape)")
    if c.model_version == "1.5":
        hosts = s["initial_host_by_insect_species"]
        if not isinstance(hosts, dict) or set(hosts) != {
            f"H{i + 1}" for i in range(len(c.insect_sizes))
        }:
            raise ValueError("Explicit initial host required for every insect species")
        for i, n in enumerate(c.insect_sizes):
            value = hosts[f"H{i + 1}"]
            keys = [f"P{p + 1}" for p in range(len(c.plant_sizes))]
            if value not in keys or (n > 0 and c.plant_sizes[keys.index(value)] == 0):
                raise ValueError("Initial insect host must name an initially present plant species")
    cell_fields = [("initial_plant_cells", sum(c.plant_sizes))]
    if c.model_version == "1.6":
        cell_fields.append(("initial_insect_cells", sum(c.insect_sizes)))
    for name, count in cell_fields:
        cells = s[name]
        if cells is not None:
            cells = np.asarray(cells)
            if (
                cells.shape != (count, 2)
                or cells.dtype.kind not in "iu"
                or ((cells < 0) | (cells >= shape)).any()
            ):
                raise ValueError(f"{name} must provide a valid integer cell per individual")


def distance(origin, targets, shape):
    delta = np.abs(np.asarray(targets) - np.asarray(origin))
    return np.minimum(delta, np.asarray(shape) - delta).max(axis=-1)


def initial_positions(c, plant_species, insect_species):
    shape = np.asarray(c.spatial["grid_shape"])
    rng = np.random.default_rng(np.random.SeedSequence([c.seed, 150, 0]))
    cells = c.spatial["initial_plant_cells"]
    plants = (
        rng.integers(0, shape, size=(len(plant_species), 2))
        if cells is None
        else np.asarray(cells, int).reshape(-1, 2).copy()
    )
    if c.model_version == "1.6":
        insect_cells = c.spatial["initial_insect_cells"]
        insects = (
            np.random.default_rng(np.random.SeedSequence([c.seed, 150, 2])).integers(
                0, shape, size=(len(insect_species), 2)
            )
            if insect_cells is None
            else np.asarray(insect_cells, int).reshape(-1, 2).copy()
        )
        return plants, insects, None
    hosts = np.empty(len(insect_species), int)
    for h in range(len(c.insect_sizes)):
        p = int(c.spatial["initial_host_by_insect_species"][f"H{h + 1}"][1:]) - 1
        available = np.flatnonzero(plant_species == p)
        index = np.flatnonzero(insect_species == h)
        if len(index):
            hosts[index] = rng.choice(available, len(index), replace=len(index) > len(available))
    return plants, plants[hosts].copy(), hosts


def candidate_pools(plant_cells, insect_cells, c):
    # ponytail: O(Nplant*Ninsect) scan for bounded prototypes; index patch occupants for large runs.
    return [
        np.flatnonzero(
            distance(cell, plant_cells, c.spatial["grid_shape"]) <= c.spatial["search_radius"]
        )
        for cell in insect_cells
    ]


def local_competition(plant_cells, labels, alpha, c):
    """Count each patch once under periodic boundaries; include the focal plant."""
    shape = tuple(c.spatial["grid_shape"])
    occupancy = np.zeros((*shape, len(c.plant_sizes)), dtype=np.int64)
    np.add.at(occupancy, (plant_cells[:, 0], plant_cells[:, 1], labels), 1)
    radius = c.spatial["competition_radius"]
    offsets = [np.unique(np.arange(-radius, radius + 1) % length) for length in shape]
    neighbors = np.zeros_like(occupancy)
    for row in offsets[0]:
        for col in offsets[1]:
            neighbors += np.roll(occupancy, (int(row), int(col)), axis=(0, 1))
    counts = neighbors[plant_cells[:, 0], plant_cells[:, 1]]
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        return 1 / (1 + np.sum(alpha[labels] * counts, axis=1))


def disperse(origins, radius, shape, rng):
    steps = rng.integers(-radius, radius + 1, size=origins.shape)
    return (origins + steps) % np.asarray(shape)


def offspring_cells(state, host, births, kind, c, rng):
    parents = np.repeat(np.arange(len(births)), births)
    if kind == "plant":
        origins = state.plants.cells[parents]
    else:
        selected = host[parents]
        if (selected < 0).any():
            raise ValueError("Obligate-host offspring cannot originate from an unmatched parent")
        origins = state.plants.cells[selected]
    return disperse(origins, c.spatial[f"{kind}_dispersal_radius"], c.spatial["grid_shape"], rng)
