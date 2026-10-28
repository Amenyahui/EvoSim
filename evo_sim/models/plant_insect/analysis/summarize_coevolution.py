"""Read-only functional time crosses, contribution differentials and host-edge histories."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from evo_sim.models.plant_insect.analysis.analyze_macro import effective_hosts
from evo_sim.models.plant_insect.analysis.common import (
    Batch,
    audit_file,
    batch_parser,
    run_path,
    sha256,
    write_json,
)
from evo_sim.models.plant_insect.helpers import chemistry, traits, voc
from evo_sim.models.plant_insect.helpers.config import Config
from evo_sim.models.plant_insect.helpers.ecology import allocate_resources


def read_design(batch):
    """Load a complete, prespecified assay design, without hidden scientific defaults."""
    default = Path(__file__).with_name("v014_analysis.yaml")
    path = batch / "analysis.yaml"
    if not path.exists():
        path = default
    design = yaml.safe_load(path.read_text(encoding="utf-8"))
    keys = set(yaml.safe_load(default.read_text(encoding="utf-8")))
    if not isinstance(design, dict) or set(design) != keys:
        raise ValueError(f"analysis YAML requires exactly {sorted(keys)}")
    for key in ("early_generation", "persistence_generations"):
        if type(design[key]) is not int or design[key] < (key == "persistence_generations"):
            raise ValueError(f"invalid {key}")
    late = design["late_generation"]
    if late != "last_recorded" and (type(late) is not int or late < design["early_generation"]):
        raise ValueError("late_generation must be last_recorded or an integer >= early_generation")
    for key in keys - {"early_generation", "late_generation", "persistence_generations"}:
        value = design[key]
        if type(value) not in (int, float) or not np.isfinite(value) or value < 0:
            raise ValueError(f"invalid {key}")
        if key.startswith("assay_") and value == 0:
            raise ValueError(f"{key} must be positive")
    return design, path


def cross_assay(plants, insects, c, design):
    """Exact equal-individual means of independent single-host trials; no mating or RNG.

    Each trial offers one plant and the configured rejection option. It is NOT a
    simultaneous multi-host choice assay. SF is an ecological proxy, not offspring.
    """
    w = np.asarray(c.weights, float)
    D = plants["D"]
    V = chemistry.plant_voc(
        D,
        plants.get("plant_voc_independent_trait"),
        c,
        chemical=plants.get("plant_chemical_traits"),
    )
    fp = voc.fecundity(D, V, c, plants.get("plant_chemical_traits"))[0]
    metrics = (
        "acceptance",
        "q",
        "resource",
        "insect_survival",
        "plant_damage",
        "plant_survival",
        "insect_SF",
        "plant_SF",
    )
    rows = []
    for p in range(len(c.plant_sizes)):
        pm = plants["plant_species"] == p
        for h in range(len(c.insect_sizes)):
            indices = np.flatnonzero(insects["insect_species"] == h)
            count = int(pm.sum())
            total = np.zeros(len(metrics))
            # ponytail: exact O(Np*Nh) trials, O(Np) memory; stratified sampling only if needed.
            for i in indices if count else []:
                score = traits.chemical_scores(V[pm], insects["A"][i], c, w)
                if c.choice_mode == "forced":
                    acceptance = np.ones(count)
                else:
                    with np.errstate(over="ignore", under="ignore"):
                        logits = score * c.tau if c.tau else np.zeros(count)
                        acceptance = np.exp(-np.logaddexp(0, c.outside_utility - logits))
                q = traits.consumption(D[pm], insects["I"][i])
                allocation = allocate_resources(
                    np.arange(count),
                    q,
                    np.full(count, design["assay_resource_capacity"]),
                    np.full(count, design["assay_max_feeding_demand"]),
                )
                gain, damage = acceptance * allocation["g"], acceptance * allocation["h"]
                survival = 1 - (1 - c.sP0) * damage
                fi = max(0, 1 - c.insect_cost * insects["I"][i].sum())
                values = (
                    acceptance,
                    q,
                    acceptance * allocation["x"],
                    gain,
                    damage,
                    survival,
                    gain * fi,
                    survival * fp[pm],
                )
                total += [value.sum() for value in values]
            pairs = count * len(indices)
            rows.append(
                dict(
                    plant_species=f"P{p + 1}",
                    insect_species=f"H{h + 1}",
                    plant_count=count,
                    insect_count=len(indices),
                    pair_count=pairs,
                    available=bool(pairs),
                    **dict(zip(metrics, total / pairs if pairs else np.full(len(metrics), np.nan))),
                )
            )
    return rows


def contribution_differentials(snapshot, following, parents):
    """Realized biparental weighting of phenotype, with a separately reported residual.

    cov(z,w)/mean(w) = weighted parental mean - unweighted parental mean.
    The residual includes segregation, mutation and clipping; no causal selection claim.
    """
    rows = []
    plant_fields = (
        ("plant_chemical_traits",)
        if "plant_chemical_traits" in snapshot
        else ("D", "plant_voc_independent_trait")
    )
    for kind, fields in (("plant", plant_fields), ("insect", ("I", "A"))):
        ledger = parents[parents.kind == kind].set_index("individual_id")
        ids, species = snapshot[f"{kind}_ids"], snapshot[f"{kind}_species"]
        if ledger.index.has_duplicates or set(ledger.index) != set(ids):
            raise ValueError("parent/snapshot identity mismatch")
        ledger = ledger.reindex(ids)
        for s in np.unique(species):
            mask = species == s
            committed = bool(ledger.loc[mask, "contribution_committed"].all())
            weights = ledger.loc[mask, "realized_genetic_contribution"].to_numpy(float)
            if committed and (not np.isfinite(weights).all() or (weights < 0).any()):
                raise ValueError("invalid committed contribution")
            total = weights.sum() if committed else np.nan
            for field in fields:
                for coordinate in range(snapshot[field].shape[1]):
                    z = snapshot[field][mask, coordinate]
                    mean = float(z.mean())
                    weighted = float(z @ weights / total) if total > 0 else np.nan
                    nxt = np.empty(0)
                    if following is not None and committed:
                        nxt = following[field][following[f"{kind}_species"] == s, coordinate]
                    next_mean = float(nxt.mean()) if len(nxt) else np.nan
                    rows.append(
                        dict(
                            kind=kind,
                            species=f"{'P' if kind == 'plant' else 'H'}{s + 1}",
                            trait=field,
                            coordinate=coordinate,
                            N=len(z),
                            mean=mean,
                            variance=float(z.var()),
                            contribution_committed=committed,
                            mean_contribution=total / len(z),
                            contribution_weighted_mean=weighted,
                            realized_differential=weighted - mean,
                            next_mean=next_mean,
                            transmission_residual=next_mean - weighted,
                        )
                    )
    return rows


def edge_histories(networks, design):
    """Per-edge consecutive persistence; terminal uncommitted contributions are unknown."""
    rows = []
    for (h, p), group in networks.groupby(["insect_species", "plant_species"], sort=True):
        for layer, column, threshold in (
            ("resource", "Q_actual", design["resource_edge_threshold"]),
            ("contribution", "genetic_contribution", design["contribution_edge_threshold"]),
        ):
            streak, previous, first = 0, None, None
            for record in group.sort_values("generation").itertuples():
                known = layer == "resource" or record.contribution_committed
                value = float(getattr(record, column)) if known else np.nan
                if known and (not np.isfinite(value) or value < 0):
                    raise ValueError("invalid network edge weight")
                present = value > threshold if known else None
                streak = (streak + 1 if previous == record.generation - 1 else 1) if present else 0
                if present and first is None:
                    first = record.generation
                rows.append(
                    dict(
                        generation=record.generation,
                        insect_species=h,
                        plant_species=p,
                        layer=layer,
                        weight=value,
                        observed=bool(known),
                        present=present,
                        consecutive_generations=streak if known else np.nan,
                        persistent=streak >= design["persistence_generations"] if known else None,
                        first_positive_generation=first,
                    )
                )
                previous = record.generation
    return rows


def summarize(batch: Batch):
    """Summarize existing saved runs; validate every consumed raw file against its manifest."""
    design, design_path = read_design(batch.path)
    outputs = {
        name: []
        for name in (
            "cross_assay",
            "contribution_differentials",
            "host_relations",
            "host_edges",
            "host_breadth",
        )
    }
    sources, runs = {}, []
    for record in batch.status.itertuples():
        folder = run_path(batch.path, record.case, record.replicate)
        meta = json.loads((folder / "metadata.json").read_text(encoding="utf-8"))
        c = Config(**meta["config"])
        c.validate()
        context = dict(
            case=record.case, replicate=record.replicate, seed=record.seed, run_status=record.status
        )
        sources[(folder / "metadata.json").relative_to(batch.path).as_posix()] = sha256(
            folder / "metadata.json"
        )

        def checked(relative):
            path = folder / relative
            digest = sha256(path)
            if meta["output_sha256"].get(relative) != digest:
                raise ValueError(f"raw input hash mismatch: {path}")
            sources[path.relative_to(batch.path).as_posix()] = digest
            return path

        parents = pd.read_csv(
            checked("parents.csv"),
            usecols=[
                "generation",
                "kind",
                "individual_id",
                "contribution_committed",
                "realized_genetic_contribution",
            ],
        )
        networks = pd.read_csv(checked("networks.csv"))
        diagnostics = pd.read_csv(checked("diagnostics.csv"))
        events = pd.read_csv(
            checked("events.csv"), usecols=["generation", "species", "accepted", "reason"]
        )
        populations = pd.read_csv(checked("populations.csv"))
        generations = sorted(populations.generation.unique())

        def snapshot(g):
            with np.load(checked(f"snapshots/g{g:06d}.npz"), allow_pickle=False) as saved:
                return dict(saved)

        early = design["early_generation"]
        late = (
            generations[-1]
            if design["late_generation"] == "last_recorded"
            else design["late_generation"]
        )
        if early not in generations or late not in generations or early > late:
            raise ValueError(f"requested assay generations unavailable: {folder}")
        runs.append(
            dict(
                context,
                early_generation=int(early),
                late_generation=int(late),
                resource_unit=c.resource_unit,
                feeding_window=c.feeding_window,
            )
        )
        endpoints = {g: snapshot(g) for g in sorted({early, late})}
        for pg in endpoints:
            for hg in endpoints:
                outputs["cross_assay"].extend(
                    dict(context, plant_generation=int(pg), insect_generation=int(hg), **row)
                    for row in cross_assay(endpoints[pg], endpoints[hg], c, design)
                )
        current = snapshot(generations[0])
        for j, g in enumerate(generations):
            following = snapshot(generations[j + 1]) if j + 1 < len(generations) else None
            if following is not None and generations[j + 1] != g + 1:
                raise ValueError("missing generation in saved run")
            outputs["contribution_differentials"].extend(
                dict(context, generation=int(g), **row)
                for row in contribution_differentials(
                    current, following, parents[parents.generation == g]
                )
            )
            if "reference_p" in current:
                for h in range(len(c.insect_sizes)):
                    mask = current["insect_species"] == h
                    for key in ("reference_p", "reference_q"):
                        values = current[key][mask]
                        totals = values.sum(axis=1)
                        breadth = [effective_hosts(v / t) for v, t in zip(values, totals) if t > 0]
                        outputs["host_breadth"].append(
                            dict(
                                context,
                                generation=int(g),
                                insect_species=f"H{h + 1}",
                                metric=key,
                                N=len(values),
                                positive_response_individuals=len(breadth),
                                mean_individual_effective_hosts=np.mean(breadth)
                                if breadth
                                else np.nan,
                                variance_individual_effective_hosts=np.var(breadth)
                                if breadth
                                else np.nan,
                            )
                        )
            current = following
        outputs["host_edges"].extend(
            dict(context, **row) for row in edge_histories(networks, design)
        )
        relation = networks.merge(
            diagnostics.rename(columns={"reference_plant_species": "plant_species"}),
            on=["generation", "insect_species", "plant_species"],
            how="left",
            validate="one_to_one",
        )
        for kind, key in (("plant", "plant_species"), ("insect", "insect_species")):
            counts = populations.loc[populations.kind == kind, ["generation", "species", "N_start"]]
            relation = relation.merge(
                counts.rename(columns={"species": key, "N_start": f"{kind}_abundance"}),
                on=["generation", key],
                how="left",
                validate="many_to_one",
            )
        for (g, h), group in events.groupby(["generation", "species"]):
            mask = (relation.generation == g) & (relation.insect_species == h)
            relation.loc[mask, "insect_count_actual"] = len(group)
            relation.loc[mask, "acceptance_fraction_actual"] = group.accepted.mean()
            relation.loc[mask, "no_candidate_fraction_actual"] = (
                group.reason == "no_candidate"
            ).mean()
        outputs["host_relations"].extend(
            dict(context, **row) for row in relation.to_dict("records")
        )
    for name, rows in outputs.items():
        pd.DataFrame(rows).to_csv(batch.summary_dir / f"{name}.csv", index=False)
    write_json(
        batch.summary_dir / "coevolution_metadata.json",
        dict(
            schema="functional-observations-1",
            design=design,
            design_sha256=sha256(design_path),
            audit_sha256=sha256(audit_file(batch.path)),
            status_sha256=sha256(batch.path / "data" / "run_status.csv"),
            raw_sha256=sources,
            runs=runs,
            analysis_sha256={
                path.name: sha256(path)
                for path in (
                    Path(__file__),
                    Path(traits.__file__),
                    Path(voc.__file__),
                    Path(chemistry.__file__),
                    Path(__file__).parents[1] / "helpers" / "ecology.py",
                    Path(__file__).parents[1] / "helpers" / "config.py",
                    Path(__file__).with_name("analyze_macro.py"),
                    Path(__file__).with_name("common.py"),
                )
            },
            interpretation="Independent single-host assays; SF is not offspring. Realized differentials include drift. No coevolution or speciation classification.",
        ),
    )
    return outputs


def main():
    args = batch_parser(__doc__).parse_args()
    summarize(Batch.open(args.input))


if __name__ == "__main__":
    main()
