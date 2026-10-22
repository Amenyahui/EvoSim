"""Grid invariants, prescribed initial host and unrestricted later local choice."""

import hashlib
import json
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect import Config, simulate
from evo_sim.models.plant_insect.audit import audit_run
from evo_sim.models.plant_insect.batch import save_result, write_json
from evo_sim.models.plant_insect.helpers.agents import CommunityState
from evo_sim.models.plant_insect.helpers.config import config_payload
from evo_sim.models.plant_insect.helpers.runtime import ecology, streams
from evo_sim.models.plant_insect.helpers.spatial import candidate_pools, distance, local_competition

ROOT = Path(__file__).resolve().parents[1]


def config(**changes):
    values = json.loads((ROOT / "tests/fixtures/spatial_demo_v15.json").read_text())["base"]
    values.update(plant_sizes=[4, 4], insect_sizes=[2], generations=3)
    values.update(changes)
    return Config(**values)


def test_grid_and_trajectory():
    c = config()
    assert distance([0, 0], [[9, 9], [0, 1], [5, 5]], [10, 10]).tolist() == [1, 1, 5]
    cells = np.array([[0, 0], [0, 0], [1, 1], [5, 5]])
    labels = np.array([0, 1, 0, 1])
    alpha = np.array([[0.2, 0.1], [0.3, 0.4]])
    phi = local_competition(cells, labels, alpha, c)
    np.testing.assert_allclose(phi, [1 / 1.5, 1 / 2.0, 1 / 1.5, 1 / 1.4])
    tiny = replace(c, spatial=dict(c.spatial, grid_shape=[1, 1], competition_radius=1))
    np.testing.assert_allclose(
        local_competition(np.zeros((4, 2), int), labels, alpha, tiny),
        [1 / 1.6, 1 / 2.4, 1 / 1.6, 1 / 2.4],
    )
    a = simulate(c)
    b = simulate(c)
    assert a["metadata"]["status"]["status"] == "completed"
    for name in (
        "events",
        "plants",
        "parents",
        "offspring",
        "populations",
        "networks",
        "diagnostics",
    ):
        pd.testing.assert_frame_equal(
            pd.DataFrame(a[name]), pd.DataFrame(b[name]), check_exact=True
        )
    first = pd.DataFrame(a["events"]).query("generation==0")
    assert first.plant_species.eq("P1").all() and first.reason.eq("initial_host").all()
    assert first.host_id.nunique() == len(first)
    silent = simulate(replace(c, record_observations=False))
    for x, y in zip(a["snapshots"], silent["snapshots"]):
        for key in x:
            np.testing.assert_array_equal(x[key], y[key])
    with TemporaryDirectory() as tmp:
        save_result(a, Path(tmp) / "run")
        audit_run(Path(tmp) / "run")
    # At a later state only P2 lies within this insect's radius: ancestral P1
    # designation must not restrict current host candidates or force P1.
    state = CommunityState.from_arrays(a["snapshots"][0])
    state.plants.cells[:4] = [5, 5]
    state.plants.cells[4:] = [0, 0]
    state.insects.cells[:] = [0, 0]
    local = replace(c, spatial=dict(c.spatial, search_radius=0), choice_mode="forced")
    pools = candidate_pools(state.plants.cells, state.insects.cells, local)
    assert all(np.array_equal(p, np.arange(4, 8)) for p in pools)
    B, w = local.validate()
    e = ecology(state, local, B, w, np.random.default_rng(8), streams(local)[1][1])
    assert (state.plants.species[e["host"]] == 1).all()
    print(
        "PASS periodic geometry, no duplicate patches, reproducibility, P1 initial hosts, local P2 access, observation isolation"
    )


def test_boundaries_and_tampering():
    for params in (
        dict(boundary="reflect"),
        dict(search_radius=-1),
        dict(grid_shape=[0, 10]),
        dict(initial_host_by_insect_species={"H1": "P3"}),
        dict(initial_plant_cells=[[0, 0]]),
    ):
        c = config()
        c.spatial.update(params)
        try:
            c.validate()
        except ValueError:
            pass
        else:
            raise AssertionError(params)
    for changes in (
        dict(insect_sizes=[0]),
        dict(plant_sizes=[0, 0], insect_sizes=[0]),
        dict(plant_sizes=[0, 4], insect_sizes=[0]),
    ):
        r = simulate(config(**changes))
        with TemporaryDirectory() as tmp:
            save_result(r, Path(tmp) / "run")
            audit_run(Path(tmp) / "run")
    zero = config(
        spatial=dict(config().spatial, plant_dispersal_radius=0, insect_dispersal_radius=0)
    )
    r = simulate(zero)
    with TemporaryDirectory() as tmp:
        save_result(r, Path(tmp) / "run")
        audit_run(Path(tmp) / "run")
    for name, edit in [("plants.csv", "density_phi"), ("events.csv", "host_cell_row")]:
        with TemporaryDirectory() as tmp:
            folder = Path(tmp) / "run"
            save_result(r, folder)
            frame = pd.read_csv(folder / name)
            frame[edit] = frame[edit].astype(float)
            frame.loc[0, edit] += 0.2
            frame.to_csv(folder / name, index=False)
            meta = json.loads((folder / "metadata.json").read_text())
            meta["output_sha256"][name] = hashlib.sha256((folder / name).read_bytes()).hexdigest()
            write_json(folder / "metadata.json", meta)
            try:
                audit_run(folder)
            except AssertionError:
                pass
            else:
                raise AssertionError("Spatial corruption undetected")
    assert config_payload(config())["spatial"] is not None
    print(
        "PASS spatial validation, empty taxa, zero dispersal and independent corruption detection"
    )


if __name__ == "__main__":
    test_grid_and_trajectory()
    test_boundaries_and_tampering()
