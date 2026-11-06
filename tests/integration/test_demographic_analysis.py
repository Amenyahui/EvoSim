"""Prespecified demographic cases and read-only finite-horizon interpretation."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
import pandas as pd

from evo_sim.core.config import merge_config
from evo_sim.models.plant_insect.analysis.summaries.lifecycle import (
    coexistence_outcome,
    summarize_sensitivity,
)
from evo_sim.models.plant_insect.audit.run import audit_run
from evo_sim.models.plant_insect.config.loader import read_config
from evo_sim.models.plant_insect.config.schema import Config, config_payload

CONFIGS = Path(__file__).resolve().parents[2] / "evo_sim/models/plant_insect/configs"


def test_demographic_design_preserves_baseline_and_single_factor_cases():
    baseline = read_config(CONFIGS / "current.yaml")
    design = read_config(CONFIGS / "demographic_sensitivity.yaml")
    assert design["parameters"] == baseline["parameters"]
    assert design["simulation_metadata"]["num_steps_per_episode"] == 200
    assert design["batch"]["replicates"] == 4
    cases = {case["name"]: case["overrides"] for case in design["batch"]["cases"]}
    assert len(cases) == 7 and cases["baseline"] == {}
    assert cases["damage_off"] == {"plant_demography": {"defoliation_mortality_rate": 0.0}}
    demography = baseline["parameters"]["plant_demography"]
    for case, field, factor in [
        ("damage_half", "defoliation_mortality_rate", 0.5),
        ("recovery_half", "resource_recovery_rate", 0.5),
        ("recovery_double", "resource_recovery_rate", 2.0),
    ]:
        assert cases[case] == {"plant_demography": {field: demography[field] * factor}}
    assert cases["newborn_half"] == {"plant_demography": {"newborn_resource_fraction": 0.5}}
    assert cases["plants_only"] == {"insect_sizes": [0]}
    for overrides in cases.values():
        parameters = merge_config(design["parameters"], overrides)
        assert parameters["resource_capacity"] == baseline["parameters"]["resource_capacity"]
        assert (
            parameters["insect_offspring_lambda"]
            == baseline["parameters"]["insect_offspring_lambda"]
        )
        assert parameters["genetics"] == baseline["parameters"]["genetics"]
        for replicate in range(4):
            config = Config(**dict(parameters, generations=200, seed=20261102 + replicate))
            config.validate()


def test_finite_horizon_coexistence_keeps_loss_and_unknown_separate():
    record = dict(
        coexistence_applicable=True,
        audited=True,
        last_generation=200,
        extinction_generation=np.nan,
        plant_extinction_generation=np.nan,
    )
    assert coexistence_outcome(record, 200) == "coexistence_at_horizon"
    assert coexistence_outcome(dict(record, last_generation=30), 200) == "unknown"
    assert coexistence_outcome(dict(record, plant_extinction_generation=20), 200) == "known_loss"
    assert coexistence_outcome(dict(record, extinction_generation=40), 200) == "known_loss"
    assert coexistence_outcome(dict(record, extinction_generation=201), 200) == (
        "coexistence_at_horizon"
    )
    assert coexistence_outcome(dict(record, audited=False), 200) == "unknown"
    assert coexistence_outcome(dict(record, coexistence_applicable=False), 200) == "not_applicable"


def test_demographic_summary_separates_reference_deficit_from_size_and_counts_cycle_losses():
    parameters = read_config(CONFIGS / "current.yaml")["parameters"]
    parameters["generations"] = 2
    parameters = config_payload(Config(**parameters))
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        folder = root / "data/baseline/rep01"
        (folder / "snapshots").mkdir(parents=True)
        (root / "settings.json").write_text(
            json.dumps(
                dict(
                    version="0.1.11",
                    spec=dict(
                        base=parameters,
                        cases=[dict(name="baseline", overrides={})],
                        replicates=1,
                    ),
                    source_sha256={},
                )
            ),
            encoding="utf-8",
        )
        (folder / "metadata.json").write_text(
            json.dumps(
                dict(
                    config=parameters,
                    status=dict(status="completed", generation=2),
                    source_hash={},
                )
            ),
            encoding="utf-8",
        )
        pd.DataFrame(
            [
                dict(
                    case="baseline",
                    replicate=1,
                    seed=parameters["seed"],
                    status="completed",
                    generation=2,
                )
            ]
        ).to_csv(root / "data/run_status.csv", index=False)
        np.savez(folder / "snapshots/g000002.npz", insect_ids=[1], plant_ids=[1, 2])
        populations = []
        for generation in range(2):
            populations.extend(
                [
                    dict(
                        generation=generation,
                        kind="plant",
                        N_start=2 if generation == 0 else 1,
                        N_retained=1,
                        N_propagules=2 * generation,
                        N_offspring=generation,
                    ),
                    dict(
                        generation=generation,
                        kind="insect",
                        N_start=1,
                        N_female=1,
                        N_mothers_with_mates=1,
                        N_potential_eggs=2,
                        N_offspring=2,
                    ),
                ]
            )
        pd.DataFrame(populations).to_csv(folder / "populations.csv", index=False)
        pd.DataFrame(
            dict(generation=[0, 0, 1], q=[0.5, 1.0, 1.0], dmax=[1.0] * 3, x=[0.5] * 3)
        ).to_csv(folder / "larvae.csv", index=False)
        pd.DataFrame(
            dict(
                generation=[0, 1],
                resource_capacity=[10.0, 10.0],
                R=[2.5, 2.5],
                T=[0.0, 5.0],
                resource_after_feeding=[2.5, 0.0],
                resource_recovery=[3.75, 5.0],
                resource_next=[6.25, 5.0],
                defoliation=[0.75, 1.0],
                gamma=[1.0, 0.5],
            )
        ).to_csv(folder / "plants.csv", index=False)
        pd.DataFrame(
            dict(
                generation=[0, 1],
                kind=["plant"] * 2,
                expected_propagules=[0.0, 2.0],
                eligible_as_mother=[False, True],
                has_mate=[False, False],
            )
        ).to_csv(folder / "parents.csv", index=False)
        # Synthetic ledgers exercise summary semantics; integrity itself is tested elsewhere.
        with patch(
            "evo_sim.models.plant_insect.analysis.summaries.lifecycle.audit_run",
            return_value=dict(status="completed"),
        ):
            tables, _ = summarize_sensitivity(root)
        assert tables["runs"].iloc[0].coexistence_outcome == "coexistence_at_horizon"
        before, pulse = tables["diagnostics"].iloc[0], tables["diagnostics"].iloc[1]
        assert before.plant_resource_fraction == 0.25 and before.plant_defoliation_mean == 0.75
        assert before.plant_losses == before.plant_cycle_losses == 1
        assert before.plant_resource_limited_fraction == 0.0
        assert before.larval_q_saturated_fraction == 0.5 and before.larval_q_mean == 0.75
        assert pulse.plant_cycle_losses == 0 and np.isnan(pulse.plant_cycle_replacement_fraction)
        assert pulse.plant_reproductive_mothers == 1
        assert pulse.plant_resource_limited_fraction == 1.0


def test_readers_reject_old_invalid_run_and_empty_batch():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        status = dict(status="invalid_parameter")
        parameters = read_config(CONFIGS / "current.yaml")["parameters"]
        parameters["model_version"] = "0.1.9"
        parameters["generations"] = 2
        (root / "metadata.json").write_text(
            json.dumps(dict(version="0.1.9", status=status, output_sha256={}, config=parameters)),
            encoding="utf-8",
        )
        (root / "status.json").write_text(json.dumps(status), encoding="utf-8")
        (root / "settings.json").write_text(
            json.dumps(
                dict(
                    version="0.1.9",
                    spec=dict(
                        base=parameters, cases=[dict(name="baseline", overrides={})], replicates=1
                    ),
                )
            ),
            encoding="utf-8",
        )
        for reader in (audit_run, summarize_sensitivity):
            try:
                reader(root)
            except ValueError:
                pass
            else:
                raise AssertionError("old protocol must be rejected")
