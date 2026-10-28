import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect.audit import audit
from evo_sim.models.plant_insect.batch import run
from evo_sim.models.plant_insect.helpers.inputs import read_config


def main():
    parser = argparse.ArgumentParser(description="Zero-generation resource and density diagnosis")
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "configs/initial_conditions.yaml",
    )
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    config = read_config(args.config)
    generations = config["simulation_metadata"]["num_steps_per_episode"]
    if generations != 0 or any(
        case["overrides"].get("generations", generations) != 0 for case in config["batch"]["cases"]
    ):
        raise ValueError("initial-condition diagnosis requires zero generations")
    spec = dict(config["batch"], base=dict(config["parameters"], generations=generations))
    out = args.out
    run(spec, out, runtime_config=config)
    audit(out)
    summary_dir = out / "data" / "summary"
    summary_dir.mkdir()
    rows = []
    for case in spec["cases"]:
        for rep in range(1, spec["replicates"] + 1):
            folder = out / "data" / case["name"] / f"rep{rep:02d}"
            with np.load(folder / "snapshots/g000000.npz") as snapshot:
                plant_count = len(snapshot["plant_ids"])
                insect_count = len(snapshot["insect_ids"])
                total_resource = float(snapshot["R"].sum())
            rows.append(
                dict(
                    case=case["name"],
                    replicate=rep,
                    plants=plant_count,
                    insects=insect_count,
                    total_resource=total_resource,
                    observation_stage="adult_plant_snapshot_no_larvae",
                    mean_plant_damage=np.nan,
                    mean_phi=np.nan,
                    mean_plant_expected_maternal_offspring=np.nan,
                    mean_plant_potential_maternal_mu_before_survival_and_mate_limitation=np.nan,
                )
            )
    frame = pd.DataFrame(rows)
    frame.to_csv(summary_dir / "initial_conditions_raw.csv", index=False)
    summary = (
        frame.groupby(["case", "plants", "insects"], sort=False)
        .mean(numeric_only=True)
        .drop(columns="replicate")
        .reset_index()
    )
    summary.to_csv(summary_dir / "initial_conditions_summary.csv", index=False)
    print(
        summary[
            [
                "case",
                "mean_plant_damage",
                "mean_phi",
                "mean_plant_expected_maternal_offspring",
                "mean_plant_potential_maternal_mu_before_survival_and_mate_limitation",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
