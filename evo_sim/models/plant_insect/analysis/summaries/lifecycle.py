"""Stage-explicit current ledger summaries; no new ecological indices."""

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect.analysis.common import (
    batch_parser,
    read_audit,
    read_settings,
    run_path,
)
from evo_sim.models.plant_insect.audit.run import audit_run


def prepare_diagnostics(batch):
    """Reuse stage-neutral trait/genetic summaries; add adult endpoints and larval resource shares."""
    from evo_sim.models.plant_insect.analysis.summaries.qtl import summarize as summarize_genetics

    assert read_settings(batch.path)["version"] == "0.1.12"
    assert batch.audit["genetically_audited_runs"] == len(batch.status)
    if batch.status.generation.lt(1).any():
        raise ValueError(
            "Lifecycle diagnostic plots require a committed transition in every run; "
            "use diagnose_initial_conditions for zero-generation snapshots"
        )
    summarize_genetics(batch, print_summary=False)
    populations, extinctions, resources = [], [], []
    for run in batch.status.itertuples(index=False):
        folder = run_path(batch.path, run.case, run.replicate)
        pop = pd.read_csv(folder / "populations.csv")
        terminal = int(run.generation)
        with np.load(folder / "snapshots" / f"g{terminal:06d}.npz") as state:
            for (kind, species), rows in pop.groupby(["kind", "species"]):
                count = int(np.count_nonzero(state[kind + "_species"] == int(species[1:]) - 1))
                populations.append(
                    dict(
                        case=run.case,
                        replicate=run.replicate,
                        generation=terminal,
                        kind=kind,
                        species=species,
                        N_start=count,
                        update_applied=False,
                    )
                )
                lost = rows.extinct_at_generation.dropna()
                observed = len(lost) > 0
                extinctions.append(
                    dict(
                        case=run.case,
                        replicate=run.replicate,
                        species=species,
                        extinction_observed=observed,
                        event_or_censor_generation=int(lost.min()) if observed else terminal,
                    )
                )
        populations.extend(pop.assign(case=run.case, replicate=run.replicate).to_dict("records"))
        net = pd.read_csv(folder / "networks.csv")
        for t, rows in net[net.insect_species == "H1"].groupby("generation"):
            total = rows.Q_actual.sum()
            resources.append(
                dict(
                    case=run.case,
                    replicate=run.replicate,
                    generation=t,
                    total_actual_resource=total,
                    P2_resource_fraction=rows.loc[rows.plant_species == "P2", "Q_actual"].sum()
                    / total
                    if total > 0 and (rows.plant_species == "P2").any()
                    else np.nan,
                )
            )
    tables = dict(
        population_timeseries=pd.DataFrame(populations),
        species_extinction=pd.DataFrame(extinctions),
        host_resource_by_generation=pd.DataFrame(resources),
    )
    for name, frame in tables.items():
        frame.to_csv(batch.summary_dir / f"{name}.csv", index=False)
    # Keep the reused run summary's endpoint fields consistent with the terminal snapshot.
    runs = pd.read_csv(batch.summary_dir / "run_summary.csv")
    for endpoint in populations:
        if endpoint.get("update_applied") is False:
            mask = (runs.case == endpoint["case"]) & (runs.replicate == endpoint["replicate"])
            runs.loc[mask, endpoint["species"] + "_final"] = endpoint["N_start"]
            peak = endpoint["species"] + "_peak"
            runs.loc[mask, peak] = runs.loc[mask, peak].clip(lower=endpoint["N_start"])
    for event in extinctions:
        mask = (runs.case == event["case"]) & (runs.replicate == event["replicate"])
        runs.loc[mask, event["species"] + "_extinction_generation"] = (
            event["event_or_censor_generation"] if event["extinction_observed"] else np.nan
        )
    runs.to_csv(batch.summary_dir / "run_summary.csv", index=False)
    print(runs.to_string(index=False))


def summarize(batch):
    settings = read_settings(batch)
    if settings.get("version") != "0.1.12":
        raise ValueError("summarize_lifecycle requires the current adult/larval ledgers")
    status = pd.read_csv(batch / "data/run_status.csv")
    assert read_audit(batch)["runs"] == len(status)
    sources = dict(
        lifecycle_counts="populations",
        plant_resource_cycles="plants",
        plant_establishment_proposals="plant_proposals",
        host_background_flows="host_flows",
        host_background_groups="host_groups",
        lifecycle_networks="networks",
        parent_host_contributions="lagged_networks",
    )
    collected = {name: [] for name in sources}
    for run in status.itertuples(index=False):
        path = run_path(batch, run.case, run.replicate)
        audit_run(path)
        for target, source in sources.items():
            frame = pd.read_csv(path / f"{source}.csv")
            frame.insert(0, "replicate", run.replicate)
            frame.insert(0, "case", run.case)
            collected[target].append(frame)
    output = batch / "data/summary"
    output.mkdir(exist_ok=True)
    result = {}
    for name, frames in collected.items():
        frame = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        frame.to_csv(output / f"{name}.csv", index=False)
        result[name] = frame
    return result


def main():
    parser = batch_parser(__doc__)
    parser.add_argument(
        "--sensitivity",
        action="store_true",
        help="Audit finite-horizon outcomes and stage diagnostics",
    )
    args = parser.parse_args()
    if args.sensitivity:
        summarize_sensitivity(args.input)
    else:
        summarize(args.input)


def sensitivity_survival(runs, horizon):
    """Known survival bounds; incomplete living runs are unknown after their last state."""
    from math import sqrt

    rows = []
    for case, group in runs.groupby("case", sort=False):
        for generation in range(horizon + 1):
            extinct = (
                group.audited
                & group.extinction_generation.notna()
                & group.extinction_generation.le(generation)
            )
            alive = ~extinct & group.last_generation.ge(generation) & group.audited
            n = len(group)
            a, e = int(alive.sum()), int(extinct.sum())
            unknown = n - a - e
            lo, hi = a / n, (a + unknown) / n
            ci_lo = ci_hi = np.nan
            if not unknown:
                z = 1.959963984540054
                center = (lo + z * z / (2 * n)) / (1 + z * z / n)
                half = z * sqrt(lo * (1 - lo) / n + z * z / (4 * n * n)) / (1 + z * z / n)
                ci_lo, ci_hi = center - half, center + half
            rows.append(
                dict(
                    case=case,
                    generation=generation,
                    n=n,
                    alive=a,
                    extinct=e,
                    unknown=unknown,
                    lower=lo,
                    upper=hi,
                    wilson_low=ci_lo,
                    wilson_high=ci_hi,
                )
            )
    return pd.DataFrame(rows)


def coexistence_outcome(record, horizon):
    """Classify observed finite-horizon coexistence, separately from execution success."""
    if not record["coexistence_applicable"]:
        return "not_applicable"
    if not record["audited"]:
        return "unknown"
    if any(
        pd.notna(record[key]) and record[key] <= horizon
        for key in ("extinction_generation", "plant_extinction_generation")
    ):
        return "known_loss"
    return "coexistence_at_horizon" if record["last_generation"] >= horizon else "unknown"


def summarize_sensitivity(path):
    """Audit recorded runs and summarize abundance and stage denominators without imputing failure."""
    import json
    from pathlib import Path

    from evo_sim.core.config import merge_config
    from evo_sim.models.plant_insect.analysis.common import sha256
    from evo_sim.models.plant_insect.config.schema import Config, config_payload

    path = Path(path)
    settings = read_settings(path)
    if settings.get("version") != "0.1.12":
        raise ValueError("summarize_sensitivity requires the current model protocol")
    spec = settings["spec"]
    horizon = int(spec["base"]["generations"])
    status_file = path / "data/run_status.csv"
    status = pd.read_csv(status_file) if status_file.exists() else pd.DataFrame()
    if not status.empty and status.duplicated(["case", "replicate"]).any():
        raise ValueError("duplicate case/replicate")
    runs, trajectories, diagnostics, evidence = [], [], [], {}
    for case in spec["cases"]:
        for rep in range(1, spec["replicates"] + 1):
            name = case["name"]
            expected = merge_config(spec["base"], case["overrides"])
            expected["seed"] += rep - 1
            match = (
                status[(status.case == name) & (status.replicate == rep)]
                if not status.empty
                else status
            )
            record = dict(
                case=name,
                replicate=rep,
                status="not_started",
                audited=False,
                complete=False,
                safety_stop=False,
                last_generation=np.nan,
                extinction_generation=np.nan,
                plant_extinction_generation=np.nan,
                coexistence_applicable=sum(expected["plant_sizes"]) > 0
                and sum(expected["insect_sizes"]) > 0,
                final_insects=np.nan,
                final_plants=np.nan,
                audit_error="",
            )
            if match.empty:
                runs.append(record)
                continue
            entry = match.iloc[0]
            record.update(
                status=entry.status,
                seed=entry.seed,
                safety_stop=entry.status == "computational_limit",
            )
            folder = run_path(path, name, rep)
            trajectory_start, diagnostic_start = len(trajectories), len(diagnostics)
            try:
                checked = audit_run(folder)
                if checked["status"] != entry.status:
                    raise ValueError("status mismatch")
                meta = json.loads((folder / "metadata.json").read_text(encoding="utf-8"))
                if int(meta["status"]["generation"]) != int(entry.generation):
                    raise ValueError("generation mismatch")
                if entry.status != "invalid_parameter":
                    expected = config_payload(Config(**expected))
                if meta["config"] != expected or int(entry.seed) != expected["seed"]:
                    raise ValueError("case configuration or seed mismatch")
                if meta["source_hash"] != settings["source_sha256"]:
                    raise ValueError("batch and run source fingerprints differ")
                if entry.status == "invalid_parameter":
                    raise ValueError("invalid parameter: no ecological state")
                terminal = int(entry.generation)
                if entry.status == "completed" and terminal != horizon:
                    raise ValueError("completed run did not reach horizon")
                with np.load(folder / "snapshots" / f"g{terminal:06d}.npz") as state:
                    final = len(state["insect_ids"])
                    final_plants = len(state["plant_ids"])
                pop = pd.read_csv(folder / "populations.csv")
                ins = pop[pop.kind == "insect"]
                counts = ins.groupby("generation").N_start.sum().to_dict()
                counts[terminal] = final
                plant_counts = (
                    pop[pop.kind == "plant"].groupby("generation").N_start.sum().to_dict()
                )
                plant_counts[terminal] = final_plants
                # All insect pools must be absent; a per-species extinction alone is insufficient.
                zero = [g for g, n in counts.items() if n == 0]
                extinction = min(zero) if zero else np.nan
                plant_zero = [g for g, n in plant_counts.items() if n == 0]
                if final == 0 and not pd.isna(extinction):
                    for g in range(int(extinction), horizon + 1):
                        counts[g] = 0
                record.update(
                    audited=True,
                    complete=entry.status in ("completed", "all_species_extinct"),
                    last_generation=terminal,
                    extinction_generation=extinction,
                    plant_extinction_generation=min(plant_zero) if plant_zero else np.nan,
                    final_insects=final,
                    final_plants=final_plants,
                )
                trajectories.extend(
                    dict(
                        case=name,
                        replicate=rep,
                        generation=g,
                        insects=n,
                        plants=plant_counts.get(
                            g, 0 if entry.status == "all_species_extinct" else np.nan
                        ),
                    )
                    for g, n in sorted(counts.items())
                )
                larvae = pd.read_csv(folder / "larvae.csv")
                plant_ledger = pd.read_csv(folder / "plants.csv")
                parent_ledger = pd.read_csv(folder / "parents.csv")
                oviposition = pd.read_csv(folder / "oviposition.csv")
                plant_populations = pop[pop.kind == "plant"]
                for g, frame in ins.groupby("generation"):
                    female = frame.N_female.sum()
                    mates = frame.N_mothers_with_mates.sum()
                    eggs = frame.N_potential_eggs.sum()
                    born = frame.N_offspring.sum()
                    diag = dict(
                        case=name,
                        replicate=rep,
                        generation=int(g),
                        females=female,
                        mated=mates,
                        eggs=eggs,
                        born=born,
                        mate_fraction=mates / female if female else np.nan,
                        acceptance=born / eggs if eggs else np.nan,
                    )
                    choices = oviposition[oviposition.reproduction_round == g]
                    diag.update(
                        females_without_mates=female - mates,
                        eggs_rejected=int(choices.reason.eq("rejected").sum()),
                        eggs_outside_option=int(choices.reason.eq("outside_option").sum()),
                        eggs_without_candidates=int(choices.reason.eq("no_candidate").sum()),
                    )
                    selected = larvae[larvae.generation == g]
                    demand = selected.dmax.sum()
                    consumed = selected.x.sum()
                    diag.update(
                        demand=demand,
                        consumed=consumed,
                        resource_gain=consumed / demand if demand else np.nan,
                        larval_count=len(selected),
                        larval_expected_recruits=selected.W_insect.sum(),
                        larval_realized_recruits=int(selected.survived_to_adult.eq(True).sum()),
                        larval_survival_fraction=selected.survived_to_adult.eq(True).mean()
                        if len(selected)
                        else np.nan,
                        larval_q_mean=selected.q.mean(),
                        larval_q_saturated_fraction=selected.q.ge(1.0).mean()
                        if len(selected)
                        else np.nan,
                    )
                    plants = plant_ledger[plant_ledger.generation == g]
                    plant_population = plant_populations[plant_populations.generation == g]
                    plant_parents = parent_ledger[
                        (parent_ledger.generation == g) & (parent_ledger.kind == "plant")
                    ]
                    k_total = plants.resource_capacity.sum()
                    b_total = plants.R.sum()
                    propagules = plant_population.N_propagules.sum()
                    newborns = plant_population.N_offspring.sum()
                    n_plants = plant_population.N_start.sum()
                    retained = plant_population.N_retained.sum()
                    losses = n_plants - retained
                    eligible = plant_parents.loc[plant_parents.eligible_as_mother.eq(True)]
                    reproductive = eligible[eligible.has_mate.eq(True)]
                    positive_stock = plants.R > 0
                    attacked = positive_stock & plants.herbivore_load.gt(0)
                    exhausted = plants.resource_after_feeding.le(plants.R * 1e-12)
                    diag.update(
                        plant_start=n_plants,
                        plant_retained=retained,
                        plant_propagules=propagules,
                        plant_newborns=newborns,
                        plant_losses=losses,
                        plant_net_change=retained + newborns - n_plants,
                        plant_cycle_losses=losses,
                        plant_cycle_replacement_fraction=newborns / losses if losses else np.nan,
                        plant_expected_propagules=plant_parents.expected_propagules.sum(),
                        plant_reproduction_active=True,
                        plant_reproductive_mothers=len(eligible),
                        plant_mothers_with_mates=len(reproductive),
                        plant_survivors_without_mates=int(eligible.has_mate.eq(False).sum()),
                        plant_zero_mean_mothers=int(reproductive.conditional_mu.eq(0).sum()),
                        plant_stochastic_zero_propagules=int(
                            (
                                reproductive.conditional_mu.gt(0)
                                & reproductive.propagule_count.eq(0)
                            ).sum()
                        ),
                        plant_retention_fraction=retained / n_plants if n_plants else np.nan,
                        plant_establishment_fraction=newborns / propagules
                        if propagules
                        else np.nan,
                        plant_resource_before=b_total,
                        plant_resource_capacity=k_total,
                        plant_resource_after_feeding=plants.resource_after_feeding.sum(),
                        plant_resource_recovery=plants.resource_recovery.sum(),
                        plant_resource_next_retained=plants.resource_next.sum(),
                        plant_resource_fraction=b_total / k_total if k_total else np.nan,
                        plant_defoliation_mean=plants.defoliation.mean(),
                        plant_post_resource_deficit_mean=(
                            1 - plants.resource_after_feeding / plants.resource_capacity
                        ).mean(),
                        plant_start_deficit_mean=(1 - plants.R / plants.resource_capacity).mean(),
                        plant_consumption_capacity_fraction_mean=(
                            (plants.R - plants.resource_after_feeding) / plants.resource_capacity
                        ).mean(),
                        plant_positive_stock_count=int(positive_stock.sum()),
                        plant_attacked_positive_stock_count=int(attacked.sum()),
                        plant_zero_stock_fraction=plants.R.eq(0).mean() if len(plants) else np.nan,
                        plant_resource_exhausted_fraction=exhausted[positive_stock].mean(),
                        plant_attacked_resource_exhausted_fraction=exhausted[attacked].mean(),
                        plant_resource_limited_fraction=plants["T"].gt(plants.R).mean()
                        if len(plants)
                        else np.nan,
                    )
                    diagnostics.append(diag)
                evidence[f"{name}/rep{rep:02d}"] = dict(
                    audit=checked, metadata_sha256=sha256(folder / "metadata.json")
                )
            except (AssertionError, OSError, ValueError, KeyError, TypeError) as error:
                del trajectories[trajectory_start:]
                del diagnostics[diagnostic_start:]
                record.update(
                    audited=False,
                    complete=False,
                    last_generation=np.nan,
                    extinction_generation=np.nan,
                    plant_extinction_generation=np.nan,
                    final_insects=np.nan,
                    final_plants=np.nan,
                    audit_error=str(error) or type(error).__name__,
                )
            runs.append(record)
    output = path / "data/summary"
    output.mkdir(parents=True, exist_ok=True)
    tables = dict(
        runs=pd.DataFrame(runs),
        trajectories=pd.DataFrame(trajectories),
        diagnostics=pd.DataFrame(diagnostics),
    )
    tables["runs"]["coexistence_outcome"] = tables["runs"].apply(
        coexistence_outcome, axis=1, horizon=horizon
    )
    tables["survival"] = sensitivity_survival(tables["runs"], horizon)
    for name, frame in tables.items():
        frame.to_csv(output / f"sensitivity_{name}.csv", index=False)
    from evo_sim.models.plant_insect.analysis.common import write_json

    write_json(
        output / "sensitivity_audit.json",
        dict(
            settings_sha256=sha256(path / "settings.json"),
            status_sha256=sha256(status_file) if status_file.exists() else None,
            source_sha256=sha256(Path(__file__)),
            runs=evidence,
            diagnostic_definitions=dict(
                generations="committed transitions only; no extrapolation of stage diagnostics",
                exhaustion="post-feeding resource <= 1e-12 * starting stock; denominator excludes zero starting stock",
                attacked_exhaustion="same threshold; denominator is positive-stock plants with herbivore_load > 0",
                resource_deficit="start deficit + consumption/capacity = post resource deficit; defoliation is instead consumption/starting resource",
                plant_reproduction="zero mean and positive-mean Poisson zeros among surviving mothers with mates",
                larval_recruits="sum of recorded W_insect versus count of survived_to_adult; zero larvae gives NA fraction",
                oviposition_reasons="rejected means all contacts rejected; outside_option is terminal search exit, possibly after contact rejection; no_candidate means no spatial candidates",
            ),
        ),
    )
    return tables, settings


if __name__ == "__main__":
    main()
