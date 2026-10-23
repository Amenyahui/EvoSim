"""Configured lifecycle must preserve trajectories, including incremental stepping."""

from copy import deepcopy
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from evo_sim.core.config import read_config
from evo_sim.models.plant_insect import get_runner, simulate
from evo_sim.models.plant_insect.helpers.config import Config
from tests.support import config as tiny


def same_results(actual, expected):
    assert actual.keys() == expected.keys()
    for key in actual:
        if key == "metadata":
            assert actual[key]["status"] == expected[key]["status"]
        elif key == "snapshots":
            assert len(actual[key]) == len(expected[key])
            for a, b in zip(actual[key], expected[key]):
                assert a.keys() == b.keys()
                for name in a:
                    np.testing.assert_array_equal(a[name], b[name], err_msg=name)
        else:
            pd.testing.assert_frame_equal(
                pd.DataFrame(actual[key]), pd.DataFrame(expected[key]), check_exact=True
            )


def test_lifecycle():
    runner = get_runner()
    runner.config["simulation_metadata"]["num_steps_per_episode"] = 3
    runner.init()
    runner.step(0)
    assert not runner.result["snapshots"]
    runner.step(1)
    assert runner.state["current_step"] == 1 and not runner.state["terminated"]
    runner.step(2)
    assert runner.state["terminated"] and len(runner.state_trajectory) == 4
    expected = deepcopy(runner.result)
    runner.step(2)
    same_results(runner.result, expected)
    runner.reset()
    runner.step()
    same_results(runner.result, expected)
    config = Config(**dict(runner.config["parameters"], generations=3))
    same_results(simulate(config), expected)
    assert simulate(replace(tiny(), generations=0))["metadata"]["status"]["generation"] == 0


def test_config_validation():
    config = get_runner().config
    broken = deepcopy(config)
    broken["substeps"]["0"]["policy"]["insects"]["choose"]["generator"] = "missing"
    try:
        get_runner(config=broken)
    except ValueError:
        pass
    else:
        raise AssertionError("Unregistered policy accepted")
    with TemporaryDirectory() as tmp:
        p = Path(tmp) / "cycle.yaml"
        p.write_text("extends: cycle.yaml")
        try:
            read_config(p)
        except ValueError:
            pass
        else:
            raise AssertionError("Cyclic configuration accepted")
