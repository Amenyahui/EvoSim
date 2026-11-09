"""Descriptive population and trait summaries from an audited QTL batch."""

from pathlib import Path

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect.analysis.common import (
    Batch,
    audit_file,
    batch_parser,
    run_path,
    sha256,
    write_json,
)


def partition_trait_variance(traits: pd.DataFrame) -> pd.DataFrame:
    """Partition population variance by recorded type, weighted by individual counts.

    Input variances use ddof=0. Outputs are descriptive coordinate variances,
    in squared trait units, not estimates of additive genetic variance.
    """
    keys = ["generation", "kind", "trait"]
    columns = keys + [
        "N",
        "type_count",
        "mean",
        "within_type_variance",
        "between_type_variance",
        "total_variance",
    ]
    if traits.duplicated(keys + ["species"]).any():
        raise ValueError("trait variance requires one row per generation, kind, trait and type")
    if not np.isfinite(traits.N).all() or (traits.N < 0).any() or (traits.N % 1 != 0).any():
        raise ValueError("trait counts must be finite nonnegative integers")
    rows = []
    for identity, group in traits.groupby(keys, sort=False):
        active = group[group.N > 0]
        if not np.isfinite(active[["mean", "variance"]]).all().all() or (active.variance < 0).any():
            raise ValueError("nonempty types require finite means and nonnegative variances")
        n = int(active.N.sum())
        mean = float((active.N * active["mean"]).sum() / n) if n else np.nan
        within = float((active.N * active.variance).sum() / n) if n else np.nan
        between = float((active.N * (active["mean"] - mean) ** 2).sum() / n) if n else np.nan
        rows.append(
            dict(
                zip(keys, identity),
                N=n,
                type_count=len(active),
                mean=mean,
                within_type_variance=within,
                between_type_variance=between,
                total_variance=within + between,
            )
        )
    return pd.DataFrame(rows, columns=columns)


def functional_trait_moments(path: Path) -> pd.DataFrame:
    """Reduce recorded D/I snapshots by type; do not recalculate biological outputs."""
    rows = []
    for file in sorted((path / "snapshots").glob("g*.npz")):
        with np.load(file) as state:
            generation = int(state["generation"])
            for kind, key, prefix in (("plant", "D", "P"), ("insect", "I", "H")):
                labels, values = state[kind + "_species"], state[key]
                for label in np.unique(labels):
                    selected = values[labels == label]
                    for axis in range(values.shape[1]):
                        rows.append(
                            dict(
                                generation=generation,
                                kind=kind,
                                species=f"{prefix}{int(label) + 1}",
                                N=len(selected),
                                trait=f"{key}_{axis + 1}",
                                mean=selected[:, axis].mean(),
                                variance=selected[:, axis].var(),
                            )
                        )
    return pd.DataFrame(
        rows, columns=["generation", "kind", "species", "N", "trait", "mean", "variance"]
    )


def summarize(batch: Batch, *, print_summary=True):
    assert batch.audit["genetically_audited_runs"] == len(batch.status)
    genetic = []
    input_hashes = {}
    columns = "case generation kind species N locus_id observed_heterozygosity gene_diversity mutation_count".split()
    summaries, populations, traits, partitions = [], [], [], []
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
        moments = pd.concat([trait, functional_trait_moments(path)], ignore_index=True)
        partitions.append(
            partition_trait_variance(moments).assign(case=row.case, replicate=row.replicate)
        )
        for file in sorted((path / "snapshots").glob("g*.npz")):
            input_hashes[str(file.relative_to(batch.path))] = sha256(file)
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
    pd.concat(partitions, ignore_index=True).to_csv(
        batch.summary_dir / "trait_variance_partition.csv", index=False
    )
    write_json(
        batch.summary_dir / "summarize_qtl_run_metadata.json",
        dict(
            batch=str(batch.path.resolve()),
            inputs=input_hashes,
            source_sha256=sha256(Path(__file__)),
            audit_sha256=sha256(audit_file(batch.path)),
            settings_sha256=sha256(batch.path / "settings.json"),
            interpretation="descriptive only; no independent replicate-based treatment-effect estimate",
            variance_partition=dict(
                weighting="individual counts within each run, generation, kind and trait",
                variance="ddof=0; total = within_type + between_type; squared trait units",
                type_semantics="plant maternal-source labels; insect preset compatibility pools",
                traits="inherited coordinates plus recorded defense D and processing capacity I",
                empty_types="N=0 has zero weight; absent functional types have no snapshot rows",
                limitation="not additive genetic variance, evolvability or a cross-trait speed ranking",
            ),
        ),
    )
    if print_summary:
        print(summary.to_string(index=False))


def main():
    """Summarize an audited batch through the current command-line interface."""
    args = batch_parser(__doc__).parse_args()
    summarize(Batch.open(args.input))


if __name__ == "__main__":
    main()
