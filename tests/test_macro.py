"""v1.4 controls, non-invasive observations and independent record audit."""

import copy
import hashlib
import json
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect import Config, simulate
from evo_sim.models.plant_insect.analysis.analyze_macro import edge_turnover, effective_hosts
from evo_sim.models.plant_insect.audit import audit_run
from evo_sim.models.plant_insect.batch import save_result, write_json
from evo_sim.models.plant_insect.helpers.config import config_payload

ROOT = Path(__file__).resolve().parents[1]


def config(**changes):
    args = json.loads((ROOT / "tests/fixtures/macro_pilot_v14.json").read_text())["base"]
    args.update(
        plant_sizes=[4, 4, 4],
        insect_sizes=[4, 4, 4],
        generations=3,
        mutation_rate=1.0,
        mu_Z=1.0,
        mu_A=1.0,
    )
    args.update(changes)
    return Config(**args)


def test_controls_and_observation():
    checked = 0
    for plant in (False, True):
        for insect in (False, True):
            c = config(evolve_plants=plant, evolve_insects=insect)
            result = simulate(c)
            assert result["metadata"]["status"]["status"] == "completed"
            with TemporaryDirectory() as tmp:
                save_result(result, Path(tmp) / "run")
                audit_run(Path(tmp) / "run")
            for kind, evolving, keys in [
                ("plant", plant, ("D", "plant_voc_independent_trait")),
                ("insect", insect, ("I", "A")),
            ]:
                changed = False
                for snapshot in result["snapshots"]:
                    for key in keys:
                        reference = np.asarray(c.founder_traits[key])[snapshot[kind + "_species"]]
                        if not evolving:
                            np.testing.assert_array_equal(snapshot[key], reference)
                        changed |= not np.array_equal(snapshot[key], reference)
                if evolving:
                    assert changed
            checked += 1
    c = config()
    observed = simulate(c)
    silent = simulate(replace(c, record_observations=False))
    for x, y in zip(observed["snapshots"], silent["snapshots"]):
        for key in x:
            np.testing.assert_array_equal(x[key], y[key])
    for name in ("events", "plants", "parents", "offspring", "populations", "networks"):
        base = pd.DataFrame(silent[name])
        pd.testing.assert_frame_equal(
            pd.DataFrame(observed[name])[base.columns], base, check_exact=True
        )
    print(f"PASS {checked} freeze modes; identical trajectories with observations on/off")


def test_boundaries_and_tampering():
    for changes in (
        dict(evolve_plants=False, founder_traits=None),
        dict(evolve_insects=0),
        dict(founder_traits={"D": [[1, 2]]}),
        dict(model_version="1.3"),
    ):
        try:
            config(**changes).validate()
        except ValueError:
            pass
        else:
            raise AssertionError(changes)
    try:
        simulate(config(), initial=([], [], [], []))
    except ValueError:
        pass
    else:
        raise AssertionError("Founder initial overridden")
    for changes in (
        dict(plant_sizes=[0, 0, 0]),
        dict(insect_sizes=[0, 0, 0]),
        dict(plant_sizes=[0, 0, 0], insect_sizes=[0, 0, 0]),
        dict(olfaction_mode="neutral"),
        dict(generations=0),
    ):
        result = simulate(config(**changes))
        with TemporaryDirectory() as tmp:
            save_result(result, Path(tmp) / "run")
            audit_run(Path(tmp) / "run")
        if changes.get("olfaction_mode") == "neutral":
            diag = pd.DataFrame(result["diagnostics"])
            assert np.allclose(diag.loc[diag.insect_count > 0, "preference"], 1 / 3)
    result = simulate(config())
    for target, column in [
        ("diagnostics.csv", "q_unlimited"),
        ("events.csv", "P1_candidate_share"),
    ]:
        with TemporaryDirectory() as tmp:
            folder = Path(tmp) / "run"
            save_result(copy.deepcopy(result), folder)
            frame = pd.read_csv(folder / target)
            frame.loc[0, column] += 0.2
            frame.to_csv(folder / target, index=False)
            meta = json.loads((folder / "metadata.json").read_text())
            meta["output_sha256"][target] = hashlib.sha256(
                (folder / target).read_bytes()
            ).hexdigest()
            write_json(folder / "metadata.json", meta)
            try:
                audit_run(folder)
            except AssertionError:
                pass
            else:
                raise AssertionError("Tampered diagnostic passed independent audit")
    assert config_payload(config())["evolve_plants"] is True
    print("PASS invalid controls, empty taxa, neutral reference, final-only state and tampering")


if __name__ == "__main__":
    assert np.isclose(effective_hosts([1, 0, 0]), 1)
    assert np.isclose(effective_hosts([1 / 3] * 3), 3)
    assert np.isnan(effective_hosts([np.nan] * 3))
    assert np.isnan(edge_turnover(set(), set()))
    assert edge_turnover({("H1", "P1")}, {("H1", "P1")}) == 0
    assert edge_turnover({("H1", "P1")}, {("H1", "P2")}) == 1
    test_controls_and_observation()
    test_boundaries_and_tampering()
