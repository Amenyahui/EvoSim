"""Descriptive population and trait summaries from an audited QTL batch."""

from pathlib import Path

import pandas as pd

from evo_sim.models.plant_insect.analysis.common import (
    Batch,
    audit_file,
    batch_parser,
    run_path,
    sha256,
    write_json,
)


def summarize(batch: Batch):
    assert batch.audit["genetically_audited_runs"] == len(batch.status)
    genetic = []
    input_hashes = {}
    columns = "case generation kind species N locus_id observed_heterozygosity gene_diversity mutation_count".split()
    summaries, populations, traits = [], [], []
    for row in batch.status.itertuples():
        path = run_path(batch.path, row.case, row.replicate)
        pop = pd.read_csv(path / "populations.csv")
        populations.append(pop.assign(case=row.case, replicate=row.replicate))
        catalog = pd.read_csv(path / "allele_catalog.csv")
        mutations = pd.read_csv(path / "mutations.csv")
        input_hashes[str((path / "genetic_summary.csv").relative_to(batch.path))] = sha256(
            path / "genetic_summary.csv"
        )
        parts = []
        for chunk in pd.read_csv(path / "genetic_summary.csv", chunksize=100000, low_memory=False):
            parts.append(chunk[chunk.record_type == "trait"])
            loci = chunk[chunk.record_type.eq("allele")].assign(
                case=f"{row.case}_rep{row.replicate:02d}"
            )
            genetic.append(loci[columns].drop_duplicates())
        trait = pd.concat(parts, ignore_index=True)
        traits.append(trait.assign(case=row.case, replicate=row.replicate))
        record = dict(
            case=row.case,
            replicate=row.replicate,
            seed=row.seed,
            terminal_generation=row.generation,
            allele_catalog_entries=len(catalog),
            mutation_events=len(mutations),
            max_clipped_fraction=trait.clipped_fraction.max(),
        )
        for species, group in pop.groupby("species"):
            lost = group.loc[group.N_start.eq(0), "generation"]
            record[f"{species}_extinction_generation"] = int(lost.min()) if len(lost) else None
            record[f"{species}_final"] = int(
                group.loc[group.generation.eq(group.generation.max()), "N_start"].iloc[0]
            )
            record[f"{species}_peak"] = int(group.N_start.max())
        summaries.append(record)
    pd.concat(genetic, ignore_index=True).drop_duplicates().to_csv(
        batch.summary_dir / "genetic_timeseries.csv", index=False
    )
    summary = pd.DataFrame(summaries)
    pop = pd.concat(populations, ignore_index=True)
    trait = pd.concat(traits, ignore_index=True)
    summary.to_csv(batch.summary_dir / "run_summary.csv", index=False)
    pop.to_csv(batch.summary_dir / "population_timeseries.csv", index=False)
    trait.to_csv(batch.summary_dir / "trait_timeseries.csv", index=False)
    write_json(
        batch.summary_dir / "summarize_qtl_run_metadata.json",
        dict(
            batch=str(batch.path.resolve()),
            inputs=input_hashes,
            source_sha256=sha256(Path(__file__)),
            audit_sha256=sha256(audit_file(batch.path)),
            settings_sha256=sha256(batch.path / "settings.json"),
            interpretation="descriptive only; no independent replicate-based treatment-effect estimate",
        ),
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    args = batch_parser(__doc__).parse_args()
    summarize(Batch.open(args.input))
