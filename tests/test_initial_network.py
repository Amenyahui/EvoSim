"""Initial links emerge from the same local encounter and choice rules as later links."""

import json
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect import Config, simulate
from evo_sim.models.plant_insect.audit import audit_run
from evo_sim.models.plant_insect.batch import save_result

ROOT = Path(__file__).resolve().parents[1]


def config(**changes):
    values = json.loads((ROOT / "tests/fixtures/network_emergence.json").read_text())["base"]
    values.update(plant_sizes=[4, 4], insect_sizes=[2], generations=0)
    values.update(changes)
    return Config(**values)


def test_initial_search():
    c = config(
        choice_mode="forced",
        spatial=dict(
            config().spatial,
            search_radius=0,
            initial_plant_cells=[[0, 0]] * 4 + [[5, 5]] * 4,
            initial_insect_cells=[[0, 0], [5, 5]],
        ),
    )
    result = simulate(c)
    events = pd.DataFrame(result["events"])
    assert events.plant_species.tolist() == ["P1", "P2"]
    assert events.reason.eq("accepted").all()
    assert events.host_assignment.eq("local_choice").all()
    assert events.candidate_count.eq(4).all()
    assert events.q.iloc[0] == 1 and np.isclose(events.q.iloc[1], 1 / 3)
    assert not np.array_equal(
        result["snapshots"][0]["insect_cells"], result["snapshots"][0]["plant_cells"][:2]
    )
    with TemporaryDirectory() as tmp:
        folder = Path(tmp) / "run"
        save_result(result, folder)
        audit_run(folder)

    empty = replace(
        c,
        generations=1,
        spatial=dict(c.spatial, initial_insect_cells=[[3, 3], [3, 3]]),
    )
    result = simulate(empty)
    first = pd.DataFrame(result["events"])
    first = first[first.generation == 0]
    assert first.reason.eq("no_candidate").all() and first.host_row.eq(-1).all()
    assert (
        pd.DataFrame(result["populations"])
        .query("generation == 1 and species == 'H1'")
        .N_start.iloc[0]
        == 0
    )
    with TemporaryDirectory() as tmp:
        folder = Path(tmp) / "run"
        save_result(result, folder)
        audit_run(folder)


def test_validation_and_reproducibility():
    for change in (
        dict(initial_insect_cells=[[0, 0]]),
        dict(initial_insect_cells=[[0, 10], [0, 0]]),
        dict(initial_host_by_insect_species={"H1": "P1"}),
    ):
        c = config()
        c.spatial.update(change)
        try:
            c.validate()
        except ValueError:
            pass
        else:
            raise AssertionError(change)
    a, b = simulate(config()), simulate(config())
    assert a["metadata"]["status"]["status"] == "completed"
    for left, right in zip(a["snapshots"], b["snapshots"]):
        for key in left:
            np.testing.assert_array_equal(left[key], right[key])
    pd.testing.assert_frame_equal(
        pd.DataFrame(a["events"]), pd.DataFrame(b["events"]), check_exact=True
    )
    c = config(record_observations=False)
    silent = simulate(c)
    for left, right in zip(a["snapshots"], silent["snapshots"]):
        for key in left:
            np.testing.assert_array_equal(left[key], right[key])


if __name__ == "__main__":
    test_initial_search()
    test_validation_and_reproducibility()
    print(
        "PASS free initial search, P1/P2 choice, empty neighborhoods, validation and reproducibility"
    )
