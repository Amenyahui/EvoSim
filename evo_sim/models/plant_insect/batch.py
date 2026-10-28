"""Controlled simulation batches; never overwrite an output directory."""

import argparse
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

from evo_sim import __version__
from evo_sim.core.config import merge_config
from evo_sim.core.helpers import source_hashes
from evo_sim.core.io import atomic_csv, atomic_text

from .helpers.config import Config
from .helpers.contributions import ALLELE_COLUMNS, NETWORK_CONTRIBUTION_COLUMNS
from .helpers.inputs import read_config
from .simulator import simulate

ROOT = Path(__file__).resolve().parents[3]
COLUMNS = {
    "events": "generation seed insect_id insect_row species species_id host_row host_id plant_species accepted reason candidate_count choice_order q d x g dmax R gamma host_h I_sum W_insect",
    "plants": "generation plant_id plant_row species species_id R T gamma X h herbivore_load D_sum V_emit_mean defense_cost voc_cost W_plant",
    "parents": "generation kind species species_id individual_id is_maternal_role S_i F_i W_i phi_i lambda_s conditional_mu expected_maternal_offspring maternal_offspring_variance survived_to_mating_stage offspring_count_as_mother offspring_count_sired mate_probability sampling_applied survival_stream offspring_stream realized_genetic_contribution contribution_committed",
    "offspring": "generation kind species offspring_id maternal_id paternal_id mutations_maternal mutations_paternal offspring_sex",
    "populations": "generation kind species N_start N_survived_to_reproduce N_offspring N_next update_applied extinct_at_generation zero_fitness_count mean_g mean_h total_resource_transfer density_phi population_status",
    "networks": "generation insect_species plant_species N_accept Q_actual L_per_insect resource_share "
    + NETWORK_CONTRIBUTION_COLUMNS,
    "allele_contributions": ALLELE_COLUMNS,
}
DIAGNOSTIC_COLUMNS = (
    "generation insect_species reference_plant_species insect_count preference q_unlimited "
    "joint_pq cov_pq variance_p variance_q variance_pq panel"
)


def safe_json(value):
    if isinstance(value, float) and not math.isfinite(value):
        return f"NONFINITE:{value}"
    if isinstance(value, dict):
        return {k: safe_json(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [safe_json(v) for v in value]
    return value


def write_json(path, value):
    atomic_text(path, json.dumps(safe_json(value), indent=2, ensure_ascii=False, allow_nan=False))


def save_result(result, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    write_json(out / "status.json", result["metadata"]["status"])
    columns_by_table = dict(COLUMNS)
    if "diagnostics" in result:
        columns_by_table["diagnostics"] = DIAGNOSTIC_COLUMNS
    if "allele_catalog" in result:
        columns_by_table.update(
            haplotype_summary="generation kind species locus_left locus_right allele_left allele_right frequency",
            genetic_loci="locus_id organism functional_module chromosome_id position_cM targets mutation_rate mutation_sd",
            gamete_transmissions="generation kind offspring_id role parent_id chromosome_id locus_id position_cM inherited_allele_id child_allele_id transmitted_homolog_before_mutation crossover_count switch_before_locus mutation_event_id",
            allele_catalog="allele_id locus_id parent_allele_id birth_generation effect",
            mutations="generation kind species offspring_id locus_id copy parent_allele_id allele_id delta mutation_stream",
            genetic_summary="generation kind species N record_type locus_id allele_id frequency observed_heterozygosity gene_diversity mutation_count trait mean variance clipped_fraction",
        )
    if result.get("metadata", {}).get("version") == "0.1.6":
        columns_by_table["allele_contributions"] += " stage"
        columns_by_table["offspring"] += (
            " reproduction_round adult_cohort survived_to_adult genotype"
        )
        columns_by_table["parents"] += (
            " potential_egg_count recruited_count_as_mother recruited_count_sired realized_recruitment_contribution"
        )
        columns_by_table["networks"] = (
            "generation reproduction_round adult_cohort insect_species plant_species N_accept Q_actual L_per_insect resource_share birth_count successful_recruits edge_semantics contribution_committed"
        )
        columns_by_table.update(
            oviposition="event_id offspring_id reproduction_round adult_cohort maternal_id egg_index host_id plant_species accepted reason candidate_count candidate_plant_ids choice_order",
            larvae=COLUMNS["events"]
            + " reproduction_round adult_cohort maternal_id paternal_id birth_host_id birth_host_species survived_to_adult insect_life_stage",
            host_groups="generation insect_species birth_host_species N status",
            host_flows="reproduction_round adult_cohort insect_species stage source_host_species destination_host_species genetic_contribution",
            lagged_networks="reproduction_round adult_cohort insect_species plant_species maternal_offspring paternal_offspring genetic_contribution recruitment_contribution edge_semantics",
        )
    for name, columns in columns_by_table.items():
        table = pd.DataFrame(result.get(name, []))
        if table.empty:
            table = pd.DataFrame(columns=columns.split())
        atomic_csv(table, out / f"{name}.csv")
    (out / "snapshots").mkdir()
    for state in result.get("snapshots", []):
        np.savez_compressed(out / "snapshots" / f"g{int(state['generation']):06d}.npz", **state)
    if "allele_catalog" in result:
        genotypes, phenotypes = {}, {}
        for state in result["snapshots"]:
            prefix = f"g{int(state['generation']):06d}_"
            for key, value in state.items():
                if key in {
                    kind + suffix
                    for kind in ("plant", "insect")
                    for suffix in ("_genotypes", "_ids", "_species")
                }:
                    genotypes[prefix + key] = value
                if key in ("D", "I", "A", "plant_voc_independent_trait", "plant_chemical_traits"):
                    phenotypes[prefix + key] = value
        np.savez_compressed(out / "genotypes.npz", **genotypes)
        np.savez_compressed(out / "phenotypes.npz", **phenotypes)
    hashes = {
        f.relative_to(out).as_posix(): hashlib.sha256(f.read_bytes()).hexdigest()
        for f in out.rglob("*")
        if f.is_file()
    }
    write_json(out / "metadata.json", dict(result["metadata"], output_sha256=hashes))


def create_output_directory(out=None):
    if out is not None:
        path = Path(out).resolve()
        path.mkdir(parents=True, exist_ok=False)
        return path
    root = ROOT / "results" / f"results_v{__version__}"
    root.mkdir(parents=True, exist_ok=True)
    numbers = [int(p.name[3:]) for p in root.iterdir() if re.fullmatch(r"run[0-9]+", p.name)]
    number = max(numbers, default=0) + 1
    while True:
        path = root / f"run{number}"
        try:
            path.mkdir()
            return path
        except FileExistsError:
            number += 1


def run(spec, out=None, *, runtime_config=None):
    out = create_output_directory(out)
    write_json(
        out / "settings.json",
        dict(
            version=spec.get("base", {}).get("model_version", "unspecified"),
            purpose=spec.get("purpose", "implementation, boundary and ledger checks only"),
            spec=spec,
            source_sha256=source_hashes(),
            pipeline=runtime_config,
        ),
    )
    (out / "data").mkdir()
    statuses = []
    try:
        if type(spec["replicates"]) is not int or spec["replicates"] < 1:
            raise ValueError("replicates must be a positive integer")
        cases = spec["cases"]
        case_names = [x["name"] for x in cases]
        if (
            not cases
            or len(set(case_names)) != len(case_names)
            or any(
                not isinstance(n, str) or not re.fullmatch(r"[A-Za-z0-9_-]+", n) for n in case_names
            )
        ):
            raise ValueError("case names must be unique simple names")
        for case in cases:
            for rep in range(spec["replicates"]):
                args = merge_config(spec["base"], case["overrides"])
                try:
                    if type(args.get("seed")) is not int:
                        raise ValueError("explicit integer seed required")
                    args["seed"] += rep
                    c = Config(**args)
                    c.validate()
                    result = simulate(c, runtime_config=runtime_config)
                except (ValueError, TypeError, KeyError) as error:
                    result = dict(
                        metadata=dict(
                            version=args.get("model_version", "unspecified"),
                            config=args,
                            source_hash=source_hashes(),
                            status=dict(
                                status="invalid_parameter", generation=0, message=str(error)
                            ),
                        )
                    )
                folder = out / "data" / case["name"] / f"rep{rep + 1:02d}"
                save_result(result, folder)
                status = dict(
                    case=case["name"],
                    replicate=rep + 1,
                    seed=args.get("seed"),
                    **result["metadata"]["status"],
                )
                statuses.append(status)
                atomic_csv(pd.DataFrame(statuses), out / "data" / "run_status.csv")
                print(case["name"], rep + 1, status["status"], status["generation"], flush=True)
    except (ValueError, TypeError, KeyError) as error:
        write_json(
            out / "data" / "status.json", dict(status="invalid_parameter", message=str(error))
        )
        return out
    write_json(
        out / "data" / "status.json",
        dict(
            status="finished",
            runs=len(statuses),
            outcome_counts=pd.Series([s["status"] for s in statuses]).value_counts().to_dict(),
        ),
    )
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--out", type=Path, help="Default: results/results_v<version>/run<N>")
    args = parser.parse_args()
    config = read_config(args.config)
    spec = dict(
        config["batch"],
        base=dict(
            config["parameters"], generations=config["simulation_metadata"]["num_steps_per_episode"]
        ),
    )
    output = run(spec, args.out, runtime_config=config)
    print(f"Output: {output}")


if __name__ == "__main__":
    main()
