"""Small paired input/intervention checks, not a test for desirable ecological outcomes."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from evo_sim.core.config import merge_config
from evo_sim.models.plant_insect import Config, simulate
from evo_sim.models.plant_insect.audit import audit_run
from evo_sim.models.plant_insect.batch import save_result
from evo_sim.models.plant_insect.helpers.chemistry import reference
from evo_sim.models.plant_insect.helpers.inputs import read_config
from evo_sim.models.plant_insect.helpers.traits import consumption
from evo_sim.models.plant_insect.helpers.voc import conditional_probabilities

MODEL = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def clean(value):
    """Strict JSON: absent/undefined statistics are null, never zero or NaN."""
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def distribution(values):
    a = np.asarray(values, float)
    if not len(a):
        return dict(n=0, mean=None, variance=None, min=None, median=None, max=None)
    return dict(
        n=len(a), mean=a.mean(), variance=a.var(), min=a.min(), median=np.median(a), max=a.max()
    )


def reference_arrays(snapshot, c):
    """Rebuild independently of the persisted diagnostics aggregations."""
    B, w = c.validate()
    D, V = reference(c, B)
    p = np.array([conditional_probabilities(V, a, c, w) for a in snapshot["A"]])
    p = p.reshape(len(snapshot["A"]), len(D))
    q = consumption(np.broadcast_to(D, (len(p), *D.shape)), snapshot["I"][:, None, :])
    return p, q


def initial_report(result, c):
    snapshot = result["snapshots"][0]
    p, q = reference_arrays(snapshot, c)
    for field, expected in [("reference_p", p), ("reference_q", q)]:
        np.testing.assert_allclose(snapshot[field], expected, rtol=1e-12, atol=1e-12)
    panels = []
    for h in range(len(c.insect_sizes)):
        mask = snapshot["insect_species"] == h
        for plant in range(len(c.plant_sizes)):
            a, b = p[mask, plant], q[mask, plant]
            joint = (a * b).mean() if len(a) else np.nan
            cov = joint - a.mean() * b.mean() if len(a) else np.nan
            panels.append(
                dict(
                    insect_species=f"H{h + 1}",
                    plant_species=f"P{plant + 1}",
                    p=distribution(a),
                    q=distribution(b),
                    pq=distribution(a * b),
                    joint_pq=joint,
                    cov_pq=cov,
                )
            )
            row = next(
                x
                for x in result["diagnostics"]
                if x["generation"] == 0
                and x["insect_species"] == f"H{h + 1}"
                and x["reference_plant_species"] == f"P{plant + 1}"
            )
            for field, expected in [
                ("preference", a.mean() if len(a) else np.nan),
                ("q_unlimited", b.mean() if len(b) else np.nan),
                ("joint_pq", joint),
                ("cov_pq", cov),
            ]:
                np.testing.assert_allclose(row[field], expected, equal_nan=True, atol=1e-12)
    loci = [x for x in c.genetics["loci"] if x["organism"] == "insect"]
    indices = {x["locus_id"]: i for i, x in enumerate(loci)}
    # These same pairs are measured even in the independent map (not just linked groups).
    associations = []
    for coordinate in (1, 2):
        for copy in (1, 2):
            left, right = f"insect_{coordinate}_{copy}", f"insect_{coordinate + 2}_{copy}"
            for h in range(len(c.insect_sizes)):
                g = snapshot["insect_genotypes"][snapshot["insect_species"] == h]
                x, y = g[:, :, indices[left]].ravel(), g[:, :, indices[right]].ravel()
                if not len(x):
                    associations.append(
                        dict(
                            species=f"H{h + 1}",
                            left=left,
                            right=right,
                            allele_left=None,
                            allele_right=None,
                            D=None,
                            r2=None,
                        )
                    )
                    continue
                for a in np.unique(x):
                    for b in np.unique(y):
                        fx, fy = np.mean(x == a), np.mean(y == b)
                        d = np.mean((x == a) & (y == b)) - fx * fy
                        denominator = fx * (1 - fx) * fy * (1 - fy)
                        associations.append(
                            dict(
                                species=f"H{h + 1}",
                                left=left,
                                right=right,
                                allele_left=int(a),
                                allele_right=int(b),
                                D=d,
                                r2=d * d / denominator if denominator else None,
                            )
                        )
    return dict(panels=panels, founder_IA_associations=associations)


def mutation_scale(c):
    scale = c.genetics["phenotype_map"]["effect_scale"]
    rows = []
    for kind in ("plant", "insect"):
        for coordinate in range(len(c.genetics["background"][kind])):
            mean, variance = 0.0, 0.0
            for locus in c.genetics["loci"]:
                if locus["organism"] != kind:
                    continue
                mu = locus["mutation_rate"] if c.genetics["mutation_enabled"][kind] else 0.0
                m, sd = locus["mutation_mean"][coordinate], locus["mutation_sd"][coordinate]
                mean += 2 * scale * mu * m
                variance += 2 * scale**2 * (mu * sd**2 + mu * (1 - mu) * m**2)
            rows.append(
                dict(
                    kind=kind,
                    coordinate=coordinate,
                    raw_increment_mean=mean,
                    raw_increment_variance=variance,
                )
            )
    return rows


def contribution_checks(result):
    checks = []
    for pop in result["populations"]:
        if not pop["update_applied"]:
            continue
        t, kind, species = pop["generation"], pop["kind"], pop["species"]
        parents = [
            x
            for x in result["parents"]
            if x["generation"] == t and x["kind"] == kind and x["species"] == species
        ]
        offspring = [
            x
            for x in result["offspring"]
            if x["generation"] == t + 1 and x["kind"] == kind and x["species"] == species
        ]
        total = 0.0
        for parent in parents:
            maternal = sum(x["maternal_id"] == parent["individual_id"] for x in offspring)
            paternal = sum(x["paternal_id"] == parent["individual_id"] for x in offspring)
            require(parent["offspring_count_as_mother"] == maternal, "maternal ledger mismatch")
            require(parent["offspring_count_sired"] == paternal, "paternal ledger mismatch")
            contribution = (maternal + paternal) / 2
            require(parent["contribution_committed"], "committed population lacks contribution")
            require(
                parent["realized_genetic_contribution"] == contribution,
                "realized contribution mismatch",
            )
            total += contribution
        require(total == len(offspring) == pop["N_offspring"], "birth contribution sum mismatch")
        recruits = [x for x in offspring if x["survived_to_adult"]]
        require(
            sum(x["realized_recruitment_contribution"] for x in parents)
            == len(recruits)
            == pop["N_next"],
            "recruitment contribution sum mismatch",
        )
        checks.append(
            dict(
                generation=t,
                kind=kind,
                species=species,
                contribution_sum=total,
                offspring=len(offspring),
            )
        )
    return checks


def focal_outcomes(result, kind):
    """Unconditional mean and within-species relative realized contributions at t=0."""
    rows = []
    for species in sorted({x["species"] for x in result["parents"] if x["kind"] == kind}):
        parents = [
            x
            for x in result["parents"]
            if x["generation"] == 0 and x["kind"] == kind and x["species"] == species
        ]
        committed = next(
            x["update_applied"]
            for x in result["populations"]
            if x["generation"] == 0 and x["kind"] == kind and x["species"] == species
        )
        contributions = (
            np.array(
                [(x["offspring_count_as_mother"] + x["offspring_count_sired"]) / 2 for x in parents]
            )
            if committed
            else np.array([])
        )
        mean = contributions.mean() if len(contributions) else np.nan
        relative = contributions / mean if mean > 0 else np.full(len(contributions), np.nan)
        ecological = result["plants"] if kind == "plant" else result["events"]
        id_field = "plant_id" if kind == "plant" else "insect_id"
        by_id = {x[id_field]: x for x in ecological if x["generation"] == 0}
        ecology_fields = ("h", "X") if kind == "plant" else ("x", "g", "q", "candidate_count")
        if kind == "insect":
            # Maternal outcomes summarize her offspring; adults do not feed again.
            ecology = {}
            for field in ecology_fields:
                ecology[field] = []
                for parent in parents:
                    larvae = [
                        x
                        for x in ecological
                        if x["generation"] == 0 and x["maternal_id"] == parent["individual_id"]
                    ]
                    values = [x[field] for x in larvae]
                    ecology[field].append(float(np.mean(values)) if values else np.nan)
        else:
            ecology = {
                field: [by_id[x["individual_id"]][field] for x in parents]
                for field in ecology_fields
            }
        rows.append(
            dict(
                species=species,
                ids=[x["individual_id"] for x in parents],
                ecology=ecology,
                survival=[x["S_i"] for x in parents],
                contribution=contributions.tolist(),
                relative_contribution=relative.tolist(),
                mean_contribution=mean,
            )
        )
    return rows


def same_snapshot(a, b, fields):
    for field in fields:
        np.testing.assert_array_equal(a[field], b[field], err_msg=f"paired initial {field}")


def paired_checks(results, configs):
    rows = []
    seeds = sorted(seed for name, seed in results if name == "independent")
    for seed in seeds:
        base, linked = results["independent", seed], results["linked", seed]
        same_snapshot(
            base["snapshots"][0],
            linked["snapshots"][0],
            [
                "plant_genotypes",
                "insect_genotypes",
                "D",
                "I",
                "A",
                "plant_chemical_traits",
                "plant_cells",
                "insect_cells",
            ],
        )
        a, b = configs["independent", seed].genetics, configs["linked", seed].genetics

        def strip(locus):
            return {k: v for k, v in locus.items() if k not in ("chromosome_id", "position_cM")}

        require(
            [strip(x) for x in a["loci"]] == [strip(x) for x in b["loci"]],
            "linked comparison altered non-map locus inputs",
        )
        off = results["mutation_off", seed]
        same_snapshot(
            base["snapshots"][0],
            off["snapshots"][0],
            ["plant_genotypes", "insect_genotypes", "plant_cells", "insect_cells"],
        )
        require(not off["mutations"], "mutation-off emitted mutations")
        require(len(off["allele_catalog"]) == len(a["alleles"]), "mutation-off created alleles")
        rows.append(dict(seed=seed, comparison="architecture_and_mutation", controlled=True))
        ref = results["intervention_baseline", seed]
        for name, kind in [
            ("plant_intervention", "insect"),
            ("insect_intervention", "plant"),
            ("search_radius", "both"),
        ]:
            changed = results[name, seed]
            fields = ["plant_genotypes", "insect_genotypes", "plant_cells", "insect_cells"]
            fields += ["I", "A", "insect_sex"] if kind == "insect" else []
            fields += ["D", "plant_chemical_traits"] if kind == "plant" else []
            if kind == "both":
                fields += ["D", "I", "A", "plant_chemical_traits", "insect_sex"]
            same_snapshot(ref["snapshots"][0], changed["snapshots"][0], fields)
            ca, cb = configs["intervention_baseline", seed], configs[name, seed]
            for field in (
                "plant_sizes",
                "insect_sizes",
                "resource_capacity",
                "max_feeding_demand",
                "plant_competition_alpha",
            ):
                require(getattr(ca, field) == getattr(cb, field), f"intervention changed {field}")
            if kind == "both":
                da, db = dict(vars(ca)), dict(vars(cb))
                sa, sb = dict(da.pop("spatial")), dict(db.pop("spatial"))
                sa.pop("search_radius")
                sb.pop("search_radius")
                require(da == db and sa == sb, "search intervention changed another scale")
            outcomes = []
            for focal in ("plant", "insect") if kind == "both" else (kind,):
                before, after = focal_outcomes(ref, focal), focal_outcomes(changed, focal)
                for x, y in zip(before, after):
                    dx = np.asarray(y["relative_contribution"]) - x["relative_contribution"]
                    outcomes.append(
                        dict(
                            kind=focal,
                            species=x["species"],
                            baseline=x,
                            treatment=y,
                            mean_survival_difference=np.mean(y["survival"])
                            - np.mean(x["survival"]),
                            mean_absolute_relative_contribution_change=np.mean(np.abs(dx))
                            if np.isfinite(dx).all()
                            else None,
                        )
                    )
            rows.append(
                dict(seed=seed, comparison=name, focal=kind, controlled=True, outcomes=outcomes)
            )
    return rows


def mirror_check(a, b):
    """Check complete input relabeling, not identical stochastic trajectories."""
    swap = {"P1": "P2", "P2": "P1", "H1": "H1"}
    for field in ("resource_capacity", "plant_offspring_lambda"):
        require(
            {swap[k]: v for k, v in getattr(a, field).items()} == getattr(b, field),
            f"mirror {field}",
        )
    require(
        {
            swap[k]: {swap[j]: v for j, v in row.items()}
            for k, row in a.plant_competition_alpha.items()
        }
        == b.plant_competition_alpha,
        "mirror competition axes",
    )
    require(
        {swap[k]: v for k, v in a.genetics["founder_frequencies"].items()}
        == b.genetics["founder_frequencies"],
        "mirror founder frequencies",
    )
    for key in a.observation_reference_traits:
        require(
            a.observation_reference_traits[key][::-1] == b.observation_reference_traits[key],
            "mirror reference panel",
        )
    require(a.plant_sizes[::-1] == b.plant_sizes, "mirror counts")
    n = a.plant_sizes[0]
    cells = a.spatial["initial_plant_cells"]
    require(cells[n:] + cells[:n] == b.spatial["initial_plant_cells"], "mirror positions")
    require(a.genetics["alleles"] == b.genetics["alleles"], "mirror changed global alleles")
    for field in ("max_feeding_demand", "insect_offspring_lambda", "insect_sizes"):
        require(getattr(a, field) == getattr(b, field), f"mirror altered insects: {field}")
    for field in ("plant_offspring_lambda", "insect_offspring_lambda"):
        require(
            {swap[k]: v for k, v in a.parameter_provenance[field].items()}
            == b.parameter_provenance[field],
            f"mirror provenance {field}",
        )
    require(
        {
            swap[k]: {swap[j]: v for j, v in row.items()}
            for k, row in a.parameter_provenance["plant_competition_alpha"].items()
        }
        == b.parameter_provenance["plant_competition_alpha"],
        "mirror alpha provenance",
    )
    require(
        a.spatial["initial_insect_cells"] == b.spatial["initial_insect_cells"],
        "mirror insect positions",
    )
    return dict(
        input_permutation_verified=True,
        trajectory_equality_required=False,
        reason="Random streams are indexed by species slots; two seeds cannot test distributional equivalence.",
    )


def markdown(report):
    lines = [
        "# Small causal-chain checks",
        "",
        "Numerical consistency and paired-intervention design checks, not proof of coevolution.",
        "",
        "| Case | Seed | Status | Final generation | Mutations |",
        "|---|---:|---|---:|---:|",
    ]
    for run in report["runs"]:
        lines.append(
            f"| {run['case']} | {run['seed']} | {run['status']['status']} | "
            f"{run['status']['generation']} | {run['mutation_events']} |"
        )
    lines += [
        "",
        "## Interpretation limits",
        *["- " + x for x in report["limitations"]],
        "",
        "The JSON contains initial p/q/pq distributions, identical I/A locus-pair associations,",
        "raw mutation-input moments, clipping fractions and individual paired contribution outcomes.",
        "No specialization, coexistence, effect direction or significance threshold is required.",
    ]
    return "\n".join(lines) + "\n"


def check(config, out):
    config, out = Path(config).resolve(), Path(out).resolve()
    suite = read_config(config)
    out.mkdir(parents=True, exist_ok=False)
    results, configs, runs = {}, {}, []
    for entry in suite["causal_checks"]["scenarios"]:
        path = (config.parent / entry).resolve()
        loaded = read_config(path)
        base = dict(
            loaded["parameters"], generations=loaded["simulation_metadata"]["num_steps_per_episode"]
        )
        for case in loaded["batch"]["cases"]:
            for rep in range(loaded["batch"]["replicates"]):
                parameters = merge_config(base, case["overrides"])
                parameters["seed"] += rep
                c = Config(**parameters)
                c.validate()
                key = case["name"], c.seed
                require(key not in results, "duplicate scenario/seed")
                result = simulate(c, runtime_config=loaded)
                require(
                    result["metadata"]["status"]["status"] in ("completed", "all_species_extinct"),
                    f"run not interpretable: {key}: {result['metadata']['status']}",
                )
                folder = out / "data" / case["name"] / f"rep{rep + 1:02d}"
                save_result(result, folder)
                audit_result = audit_run(folder)
                results[key], configs[key] = result, c
                runs.append(
                    dict(
                        case=case["name"],
                        seed=c.seed,
                        status=result["metadata"]["status"],
                        config=str(path),
                        yaml_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                        input_sources=loaded.get("input_sources", {}),
                        initial=initial_report(result, c),
                        mutation_input=mutation_scale(c),
                        mutation_events=len(result["mutations"]),
                        contribution_checks=contribution_checks(result),
                        audit=audit_result,
                        boundaries=[
                            x for x in result["genetic_summary"] if x["record_type"] == "trait"
                        ],
                        populations=result["populations"],
                    )
                )
                print(case["name"], c.seed, result["metadata"]["status"]["status"], flush=True)
    seed = min(s for name, s in configs if name == "independent")
    report = dict(
        schema="v012-causal-checks-1",
        runs=runs,
        paired_checks=paired_checks(results, configs),
        mirror=mirror_check(configs["independent", seed], configs["mirror", seed]),
        limitations=[
            "Two fixed seeds and short runs are numerical/design checks, not inferential evidence.",
            "Mirror checks full input permutation, not pathwise equality or distributional equivalence.",
            "Formation input is symmetric in distribution; finite founder samples need not be symmetric.",
            "Formation standing variance differs from the adapted baseline and is not a pure mean intervention.",
            "Linked pairs use 10 cM; linkage does not impose initial LD. Both maps use the same measured pairs.",
            "Opponent interventions change total phenotype composition, not one isolated molecular mechanism.",
            "Plant D intervention also changes VOC under rho=0.5; insect intervention changes both I and A.",
            "Search-radius effects may propagate through host-associated mating; only the input scale is isolated.",
            "2P2H preserves total insect abundance but splits conspecific mate pools; it is not a neutral control.",
            "Raw mutation variance is pre-clipping; observed boundary fractions are reported separately.",
            "No effect direction, specialization, persistence or coexistence is a passing criterion.",
        ],
    )
    report = clean(report)
    (out / "causal_checks.json").write_text(
        json.dumps(report, indent=2, allow_nan=False), encoding="utf-8"
    )
    (out / "causal_checks.md").write_text(markdown(report), encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=MODEL / "configs/v012_causal.yaml")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = check(args.config, args.out)
    print(f"Checked {len(report['runs'])} runs; report: {args.out / 'causal_checks.json'}")


if __name__ == "__main__":
    main()
