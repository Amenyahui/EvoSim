import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from evo_sim.core.config import read_config
from evo_sim.models.plant_insect.audit import audit
from evo_sim.models.plant_insect.batch import run


def main():
    config = read_config(
        Path(__file__).resolve().parents[1] / "configs/sexual_qtl_development.yaml"
    )
    pilot = {"base": config["parameters"]}
    base = dict(pilot["base"], generations=0)
    cases = []
    for plants, insects in [(60, 60), (60, 6), (180, 60), (600, 60), (60, 0), (180, 0), (600, 0)]:
        cases.append(
            dict(
                name=f"P{plants}_H{insects}",
                overrides=dict(plant_sizes=[plants // 2] * 2, insect_sizes=[insects]),
            )
        )
    spec = dict(
        purpose="Post-pilot initial-state resource/density diagnosis; zero population updates",
        replicates=4,
        base=base,
        cases=cases,
    )
    parser = argparse.ArgumentParser(description="Zero-generation resource and density diagnosis")
    parser.add_argument("--out", type=Path, required=True)
    out = parser.parse_args().out
    run(spec, out)
    audit(out)
    summary_dir = out / "data" / "summary"
    summary_dir.mkdir()
    rows = []
    for case in cases:
        for rep in range(1, 5):
            folder = out / "data" / case["name"] / f"rep{rep:02d}"
            p = pd.read_csv(folder / "parents.csv")
            b = pd.read_csv(folder / "plants.csv")
            e = pd.read_csv(folder / "events.csv")
            plant = p[p.kind == "plant"]
            insect = p[p.kind == "insect"]
            rows.append(
                dict(
                    case=case["name"],
                    replicate=rep,
                    plants=len(plant),
                    insects=len(insect),
                    total_resource=b.R.sum(),
                    maximum_demand=e.dmax.sum(),
                    potential_demand=e.d.sum(),
                    actual_resource=e.x.sum(),
                    mean_plant_damage=b.h.mean(),
                    mean_phi=plant.phi_i.mean(),
                    mean_plant_expected_maternal_offspring=plant.expected_maternal_offspring.mean(),
                    mean_plant_expected_without_damage=(
                        plant.lambda_s * plant.F_i * plant.phi_i
                    ).mean(),
                    mean_insect_expected_maternal_offspring=insect.expected_maternal_offspring.mean()
                    if len(insect)
                    else np.nan,
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
                "mean_plant_expected_without_damage",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
