"""Patch grid, local encounters, finite plant microsites and bounded dispersal."""

import numpy as np


def validate_spatial(c):
    s = c.spatial
    required = {
        "grid_shape",
        "boundary",
        "search_radius",
        "plant_dispersal_radius",
        "insect_dispersal_radius",
        "initial_plant_cells",
    }
    required.add("initial_insect_cells")
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
        "plant_dispersal_radius",
        "insect_dispersal_radius",
    ):
        if type(s[name]) is not int or not 0 <= s[name] <= max(shape):
            raise ValueError(f"{name} must be an integer between 0 and max(grid_shape)")
    cell_fields = [
        ("initial_plant_cells", sum(c.plant_sizes)),
        ("initial_insect_cells", sum(c.insect_sizes)),
    ]
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
    capacity = int(shape[0]) * int(shape[1]) * c.plant_demography["microsites_per_cell"]
    if capacity >= 2**63 or sum(c.plant_sizes) > capacity:
        raise ValueError("initial plants exceed finite int64 microsite capacity")
    cells = s["initial_plant_cells"]
    if cells is not None and len(cells):
        _, counts = np.unique(np.asarray(cells), axis=0, return_counts=True)
        if (counts > c.plant_demography["microsites_per_cell"]).any():
            raise ValueError("initial plant cells exceed microsites_per_cell; no plants truncated")


def distance(origin, targets, shape):
    delta = np.abs(np.asarray(targets) - np.asarray(origin))
    return np.minimum(delta, np.asarray(shape) - delta).max(axis=-1)


def initial_positions(c, plant_species, insect_species):
    shape = np.asarray(c.spatial["grid_shape"])
    rng = np.random.default_rng(np.random.SeedSequence([c.seed, 150, 0]))
    cells = c.spatial["initial_plant_cells"]
    m = c.plant_demography["microsites_per_cell"]
    if cells is None:
        sites = rng.choice(int(np.prod(shape)) * m, len(plant_species), replace=False)
        plants = np.column_stack((sites // m // shape[1], sites // m % shape[1]))
    else:
        plants = np.asarray(cells, int).reshape(-1, 2).copy()
        assigned = {}
        sites = np.empty(len(plants), dtype=np.int64)
        for i, (row, col) in enumerate(plants):
            cell = int(row * shape[1] + col)
            slot = assigned.get(cell, 0)
            sites[i] = cell * m + slot
            assigned[cell] = slot + 1
    insect_cells = c.spatial["initial_insect_cells"]
    insects = (
        np.random.default_rng(np.random.SeedSequence([c.seed, 150, 2])).integers(
            0, shape, size=(len(insect_species), 2)
        )
        if insect_cells is None
        else np.asarray(insect_cells, int).reshape(-1, 2).copy()
    )
    return plants, insects, np.asarray(sites, dtype=np.int64)


def candidate_pools(plant_cells, insect_cells, c):
    # ponytail: O(Nplant*Ninsect) scan for bounded prototypes; index patch occupants for large runs.
    return [
        np.flatnonzero(
            distance(cell, plant_cells, c.spatial["grid_shape"]) <= c.spatial["search_radius"]
        )
        for cell in insect_cells
    ]


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


def establish_propagules(cells, occupied, c, rng):
    """Shared random proposal order claims each target cell's first vacant microsite."""
    m = c.plant_demography["microsites_per_cell"]
    occupied = set(map(int, occupied))
    sites = np.full(len(cells), -1, dtype=np.int64)
    for i in rng.permutation(len(cells)):
        base = (int(cells[i, 0]) * c.spatial["grid_shape"][1] + int(cells[i, 1])) * m
        for site in range(base, base + m):
            if site not in occupied:
                occupied.add(site)
                sites[i] = site
                break
    return np.flatnonzero(sites >= 0), sites
