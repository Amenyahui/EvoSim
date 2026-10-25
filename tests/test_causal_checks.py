"""Small reproducible checks; ecological outcomes are not success conditions."""

import copy
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from evo_sim.core.config import merge_config
from evo_sim.models.plant_insect import Config
from evo_sim.models.plant_insect.analysis.check_causal_chain import (
    MODEL,
    check,
    clean,
    mirror_check,
    mutation_scale,
)
from evo_sim.models.plant_insect.helpers.inputs import read_config
from tests.support import raises


def scenario(name):
    loaded = read_config(MODEL / "configs" / name)
    return Config(
        **dict(
            loaded["parameters"], generations=loaded["simulation_metadata"]["num_steps_per_episode"]
        )
    )


def test_v012_inputs_and_mutation_scale():
    for path in (MODEL / "configs").glob("v012*.yaml"):
        loaded = read_config(path)
        base = dict(
            loaded["parameters"], generations=loaded["simulation_metadata"]["num_steps_per_episode"]
        )
        for case in loaded["batch"]["cases"]:
            Config(**merge_config(base, case["overrides"])).validate()
    base, linked = scenario("v012_base.yaml"), scenario("v012_linked.yaml")
    assert mutation_scale(base) == mutation_scale(linked)
    np.testing.assert_allclose([x["raw_increment_variance"] for x in mutation_scale(base)], 4e-6)
    off = copy.deepcopy(base)
    off.genetics["mutation_enabled"] = dict(plant=False, insect=False)
    assert all(x["raw_increment_variance"] == 0 for x in mutation_scale(off))
    mirror = scenario("v012_mirror.yaml")
    assert mirror_check(base, mirror)["input_permutation_verified"]
    mirror.resource_capacity["P1"] += 1
    raises(lambda: mirror_check(base, mirror), AssertionError)
    formation = scenario("v012_formation.yaml")
    effects = {}
    for allele in formation.genetics["alleles"]:
        if allele["locus_id"].startswith("insect"):
            effects.setdefault(allele["locus_id"], []).append(sum(allele["effect"]))
    assert all(value == [0.19, 0.23] for value in effects.values())
    assert scenario("v012_multi.yaml").insect_sizes == [12, 12]
    assert clean(dict(missing=np.nan, empty=np.inf)) == dict(missing=None, empty=None)


def test_v012_causal_suite_reproducible_and_self_contained():
    with TemporaryDirectory() as tmp:
        folder = Path(tmp)
        a = check(MODEL / "configs/v012_causal.yaml", folder / "first")
        b = check(MODEL / "configs/v012_causal.yaml", folder / "second")
        assert a == b
        assert len(a["runs"]) == 20 and len(a["paired_checks"]) == 8
        saved = json.loads((folder / "first/causal_checks.json").read_text(encoding="utf-8"))
        assert saved == a and a["mirror"]["input_permutation_verified"]
        assert all(len(run["initial"]["founder_IA_associations"]) > 0 for run in a["runs"])
        off = [run for run in a["runs"] if run["case"] == "mutation_off"]
        assert all(run["mutation_events"] == 0 for run in off)
        assert all(run["contribution_checks"] for run in a["runs"])
        assert all(
            len(run["initial"]["panels"]) == 4 for run in a["runs"] if run["case"] == "multi"
        )
        raises(lambda: check(MODEL / "configs/v012_causal.yaml", folder / "first"), FileExistsError)
