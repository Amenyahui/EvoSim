"""Explicit QTL states; SLiM executes meiosis and the configured QTL mutation kernel."""

import copy
import json
from collections import Counter
from dataclasses import asdict, dataclass
from itertools import combinations

import numpy as np

from .config import normalized_bounds, numeric
from .demography import ComputationalLimit
from .slim_backend import SlimError, reproduce

KINDS = ("plant", "insect")
TRAITS = {"plant": ("D", "plant_voc_independent_trait"), "insect": ("I", "A")}


def validate_genome_map(loci):
    groups = {}
    for i, locus in enumerate(loci):
        group, position = locus.get("chromosome_id"), locus.get("position_cM")
        if (
            not isinstance(group, str)
            or not group.strip()
            or type(position) not in (int, float)
            or not np.isfinite(position)
            or position < 0
        ):
            raise ValueError("nonempty chromosome_id and finite nonnegative position_cM required")
        groups.setdefault(group, []).append(i)
    for indices in groups.values():
        indices.sort(key=lambda i: loci[i]["position_cM"])
        positions = [loci[i]["position_cM"] for i in indices]
        if len(set(positions)) != len(positions):
            raise ValueError("duplicate position within linkage group")
    return groups


def vector(value, width, name):
    a = np.asarray(value)
    if a.dtype.kind not in "iuf" or a.shape != (width,) or not np.isfinite(a).all():
        raise ValueError(f"{name}: finite numeric vector of length {width} required")
    return a.astype(float)


def validate(c, K):
    g = c.genetics
    required = {
        "mode",
        "schema_version",
        "ploidy",
        "inheritance",
        "recombination_model",
        "background",
        "phenotype_map",
        "loci",
        "alleles",
        "founder_frequencies",
        "mutation_enabled",
        "provenance",
        "max_alleles_safety",
        "max_genotype_values_safety",
    }
    slim = isinstance(g, dict) and g.get("schema_version") == "slim-qtl-1"
    if slim:
        required.add("slim")
    if not isinstance(g, dict) or set(g) != required:
        raise ValueError("genetics requires the complete inline qtl schema; unknown keys rejected")
    if (
        g["mode"] != "sexual_qtl"
        or g["schema_version"] not in {"sexual-qtl-1", "slim-qtl-1"}
        or type(g["ploidy"]) is not int
        or g["ploidy"] != 2
        or g["inheritance"] != "biparental_meiosis"
        or g["recombination_model"]
        != ("slim_haldane_marker_map" if slim else "haldane_poisson_no_interference")
    ):
        raise ValueError("explicit supported biparental QTL protocol required")
    if slim:
        settings = g["slim"]
        if not isinstance(settings, dict) or set(settings) != {"executable", "timeout_seconds"}:
            raise ValueError("slim requires executable and timeout_seconds")
        if not isinstance(settings["executable"], str) or not settings["executable"].strip():
            raise ValueError("nonempty SLiM executable required")
        numeric(settings["timeout_seconds"], "slim.timeout_seconds")
        if settings["timeout_seconds"] <= 0:
            raise ValueError("SLiM timeout must be positive")
    width = c.n + K
    mapping = g["phenotype_map"]
    if not isinstance(mapping, dict) or set(mapping) != {"effect_scale", "bounds"}:
        raise ValueError("phenotype_map requires explicit effect_scale and bounds")
    numeric(mapping["effect_scale"], "phenotype_map.effect_scale")
    normalized_bounds(mapping["bounds"], "phenotype_map.bounds")
    if set(g["background"]) != set(KINDS) or set(g["mutation_enabled"]) != set(KINDS):
        raise ValueError("background/mutation_enabled require plant and insect")
    for kind in KINDS:
        vector(g["background"][kind], width, "background")
        if type(g["mutation_enabled"][kind]) is not bool:
            raise ValueError("mutation_enabled must be boolean")
    for name in ("max_alleles_safety", "max_genotype_values_safety"):
        if type(g[name]) is not int or not 1 <= g[name] < 2**63:
            raise ValueError(f"{name}: positive int64 count required")
    if not isinstance(g["provenance"], dict) or set(g["provenance"]) != {
        "mapping",
        "genome_map",
        "mutation",
        "founders",
        "limits",
    }:
        raise ValueError("genetic mapping, mutation, founders and limits need provenance")
    for entry in g["provenance"].values():
        if (
            not isinstance(entry, dict)
            or entry.get("source_kind")
            not in {"theoretical_placeholder", "implementation_choice", "empirical"}
            or not isinstance(entry.get("description"), str)
            or not entry["description"].strip()
        ):
            raise ValueError("genetic parameter provenance requires source_kind and description")
        if entry["source_kind"] == "theoretical_placeholder" and not c.illustrative_only:
            raise ValueError("uncalibrated genetics requires illustrative_only")
    loci = {}
    for locus in g["loci"]:
        if set(locus) != {
            "locus_id",
            "organism",
            "functional_module",
            "chromosome_id",
            "position_cM",
            "targets",
            "mutation_rate",
            "mutation_mean",
            "mutation_sd",
        }:
            raise ValueError("complete locus record required")
        key, kind = locus["locus_id"], locus["organism"]
        if not isinstance(key, str) or not key or key in loci or kind not in KINDS:
            raise ValueError("unique locus_id and valid organism required")
        if not isinstance(locus["functional_module"], str) or not locus["functional_module"]:
            raise ValueError("functional_module required")
        targets = locus["targets"]
        if (
            not isinstance(targets, list)
            or not targets
            or len(set(targets)) != len(targets)
            or any(type(x) is not int or not 0 <= x < width for x in targets)
        ):
            raise ValueError("unique target indices within D/Z or I/A coordinates required")
        mu = locus["mutation_rate"]
        if type(mu) not in (int, float) or not np.isfinite(mu) or not 0 <= mu <= 1:
            raise ValueError("mutation_rate per copy/locus/generation must be in [0,1]")
        mean = vector(locus["mutation_mean"], width, "mutation_mean")
        sd = vector(locus["mutation_sd"], width, "mutation_sd")
        if (sd < 0).any() or any(
            mean[j] != 0 or sd[j] != 0 for j in set(range(width)) - set(targets)
        ):
            raise ValueError(
                "nonnegative SD; mutation mean/SD zero outside declared targets required"
            )
        loci[key] = locus
    for kind in KINDS:
        targets = {j for x in loci.values() if x["organism"] == kind for j in x["targets"]}
        if targets != set(range(width)):
            raise ValueError("all trait dimensions require declared QTL coverage")
    for kind in KINDS:
        validate_genome_map([x for x in loci.values() if x["organism"] == kind])
    alleles = g["alleles"]
    if not alleles or len(alleles) > g["max_alleles_safety"]:
        raise ValueError("nonempty founder allele catalog within safety limit required")
    for i, a in enumerate(alleles):
        if (
            set(a) != {"allele_id", "locus_id", "effect"}
            or type(a["allele_id"]) is not int
            or a["allele_id"] != i
        ):
            raise ValueError("founder allele IDs must be consecutive integers starting at zero")
        if a["locus_id"] not in loci:
            raise ValueError("allele references unknown locus")
        effect = vector(a["effect"], width, "allele effect")
        if any(effect[j] != 0 for j in set(range(width)) - set(loci[a["locus_id"]]["targets"])):
            raise ValueError("allele effects outside declared targets")
    species = [(f"P{s + 1}", "plant") for s in range(len(c.plant_sizes))]
    species += [(f"H{s + 1}", "insect") for s in range(len(c.insect_sizes))]
    if set(g["founder_frequencies"]) != {s for s, _ in species}:
        raise ValueError("founder frequency catalog must include every species slot")
    for species_name, kind in species:
        frequencies = g["founder_frequencies"][species_name]
        if set(frequencies) != {key for key, x in loci.items() if x["organism"] == kind}:
            raise ValueError("founder frequencies must cover every locus of the organism")
        for key, entries in frequencies.items():
            ids, probs = [], []
            for entry in entries:
                if set(entry) != {"allele_id", "frequency"}:
                    raise ValueError("allele_id/frequency records required")
                idx, p = entry["allele_id"], entry["frequency"]
                if (
                    type(idx) is not int
                    or not 0 <= idx < len(alleles)
                    or alleles[idx]["locus_id"] != key
                    or idx in ids
                    or type(p) not in (int, float)
                    or not np.isfinite(p)
                    or p < 0
                ):
                    raise ValueError("invalid allele or founder frequency")
                ids.append(idx)
                probs.append(p)
            if not ids or not np.isclose(sum(probs), 1, rtol=0, atol=1e-12):
                raise ValueError("founder frequencies must sum to one")


@dataclass(frozen=True)
class Allele:
    allele_id: int
    locus_id: str
    parent_allele_id: int | None
    birth_generation: int
    effect: tuple


class GeneticState:
    def __init__(self, c):
        if c.genetics["schema_version"] != "slim-qtl-1":
            raise ValueError(
                "Execution requires slim-qtl-1; legacy results remain readable for audit"
            )
        self.spec = copy.deepcopy(c.genetics)
        self.n = c.n
        self.seed = c.seed
        self.engine_calls = []
        self.loci = {
            kind: [x for x in self.spec["loci"] if x["organism"] == kind] for kind in KINDS
        }
        self.alleles = [
            Allele(a["allele_id"], a["locus_id"], None, 0, tuple(a["effect"]))
            for a in self.spec["alleles"]
        ]
        self.events = []
        self.transmissions = []
        self.rng = {
            (phase, kind): [
                np.random.default_rng(np.random.SeedSequence([c.seed, 140, phase, k, s]))
                for s in range(len(c.plant_sizes if k == 0 else c.insect_sizes))
            ]
            for phase in (0,)
            for k, kind in enumerate(KINDS)
        }

    def founders(self, sizes, kind):
        """F-G4: independent categorical draws of two copies, not a mating model."""
        g = np.empty((sum(sizes), 2, len(self.loci[kind])), dtype=np.int64)
        offset = 0
        prefix = "P" if kind == "plant" else "H"
        for s, size in enumerate(sizes):
            for l, locus in enumerate(self.loci[kind]):
                freq = self.spec["founder_frequencies"][f"{prefix}{s + 1}"][locus["locus_id"]]
                ids = [x["allele_id"] for x in freq]
                p = np.array([x["frequency"] for x in freq], float)
                g[offset : offset + size, :, l] = self.rng[0, kind][s].choice(
                    ids, size=(size, 2), p=p / p.sum()
                )
            offset += size
        return g

    def phenotype(self, g, kind):
        """F-G1: common background + configured scale times allele sum; no species argument."""
        if g.dtype.kind not in "iu" or g.shape[1:] != (2, len(self.loci[kind])):
            raise ValueError("genotype must be N by two copies by organism loci")
        if (g < 0).any() or (g >= len(self.alleles)).any():
            raise ValueError("unknown allele ID")
        for l, locus in enumerate(self.loci[kind]):
            if any(
                self.alleles[int(a)].locus_id != locus["locus_id"] for a in np.unique(g[:, :, l])
            ):
                raise ValueError("allele assigned to wrong locus or organism")
        effects = np.asarray([a.effect for a in self.alleles])
        mapping = self.spec["phenotype_map"]
        raw = np.asarray(self.spec["background"][kind]) + mapping["effect_scale"] * effects[g].sum(
            axis=(1, 2)
        )
        values = np.clip(raw, *mapping["bounds"])
        return {TRAITS[kind][0]: values[:, : self.n], TRAITS[kind][1]: values[:, self.n :]}, raw

    def preflight(self, sizes, stored_values):
        """Conservative allocation ceiling: reserve even the all-loci-mutate case."""
        values = sum(2 * len(self.loci[k]) * sizes[i] for i, k in enumerate(KINDS))
        new_alleles = sum(
            2 * sum(x["mutation_rate"] > 0 for x in self.loci[k]) * sizes[i]
            for i, k in enumerate(KINDS)
            if self.spec["mutation_enabled"][k]
        )
        if stored_values + values > self.spec["max_genotype_values_safety"]:
            raise ComputationalLimit("genotype snapshot value limit; no descendants committed")
        if len(self.alleles) + new_alleles > self.spec["max_alleles_safety"]:
            raise ComputationalLimit(
                "conservative allele allocation limit; no descendants committed"
            )

    def fertilize(self, g, labels, ids, mothers, fathers, kind, next_id, generation):
        """Native SLiM crossing followed by Eidos-defined per-copy QTL mutation."""
        mothers, fathers = np.asarray(mothers, int), np.asarray(fathers, int)
        if mothers.shape != fathers.shape or mothers.ndim != 1:
            raise ValueError("aligned parental indices required")
        if any(((p < 0) | (p >= len(ids))).any() for p in (mothers, fathers)):
            raise ValueError("parent index outside current population")
        if (mothers == fathers).any() or (labels[mothers] != labels[fathers]).any():
            raise ValueError("distinct conspecific parents required")
        if next_id + len(mothers) >= 2**63:
            raise ComputationalLimit("individual ID int64 limit")
        child_labels = labels[mothers]
        child_ids = np.arange(next_id, next_id + len(mothers), dtype=np.int64)
        child = np.empty((len(mothers), 2, len(self.loci[kind])), dtype=np.int64)
        homolog = np.empty_like(child)
        switches = np.zeros_like(child, dtype=bool)
        inherited = np.empty_like(child)
        prefix = "P" if kind == "plant" else "H"
        if len(mothers):
            # Seeds are indexed by generation/taxon, so subprocesses have no rollback state.
            seeds = np.random.SeedSequence(
                [self.seed, 140, 4, KINDS.index(kind), generation]
            ).generate_state(2)
            seeds = [int(x) + 1 for x in seeds]
            loci = copy.deepcopy(self.loci[kind])
            if not self.spec["mutation_enabled"][kind]:
                for locus in loci:
                    locus["mutation_rate"] = 0.0
            material, call = reproduce(
                g,
                mothers,
                fathers,
                loci,
                np.asarray([a.effect for a in self.alleles]),
                self.spec["slim"],
                seeds,
                len(self.alleles),
            )
            new_alleles, new_events = [], []
            for i, h, l in np.ndindex(child.shape):
                row, locus = material[i, h, l], loci[l]
                old, allele, source, changed = (int(x) for x in row[3:7])
                parent = mothers[i] if h == 0 else fathers[i]
                if old != g[parent, source, l]:
                    raise SlimError("SLiM allele does not match recorded parental homolog")
                delta, effect = (
                    row[7 : 7 + len(self.alleles[old].effect)],
                    row[7 + len(self.alleles[old].effect) :],
                )
                expected = np.asarray(self.alleles[old].effect) + delta
                if not np.allclose(effect, expected, rtol=1e-12, atol=1e-14):
                    raise SlimError("SLiM allele effect does not equal inherited effect plus delta")
                deterministic = np.asarray(locus["mutation_sd"]) == 0
                if changed:
                    if locus["mutation_rate"] == 0 or allele != len(self.alleles) + len(
                        new_alleles
                    ):
                        raise SlimError("Invalid SLiM mutation identity or disabled mutation")
                    if not np.array_equal(
                        delta[deterministic], np.asarray(locus["mutation_mean"])[deterministic]
                    ):
                        raise SlimError("Invalid deterministic mutation effect")
                    new_alleles.append(
                        Allele(allele, locus["locus_id"], old, generation, tuple(effect))
                    )
                    new_events.append(
                        dict(
                            event_id=allele,
                            role="maternal" if h == 0 else "paternal",
                            chromosome_id=locus["chromosome_id"],
                            generation=generation,
                            kind=kind,
                            species=f"{prefix}{child_labels[i] + 1}",
                            offspring_id=int(child_ids[i]),
                            locus_id=locus["locus_id"],
                            copy=h,
                            parent_allele_id=old,
                            allele_id=allele,
                            delta=delta.tolist(),
                            mutation_stream=f"SLiM/seed{seeds[1]}/{kind}/generation{generation}",
                        )
                    )
                elif allele != old or delta.any() or locus["mutation_rate"] == 1:
                    raise SlimError("Invalid unchanged SLiM allele")
                inherited[i, h, l], child[i, h, l], homolog[i, h, l] = old, allele, source
            self.alleles = self.alleles + new_alleles
            self.events = self.events + new_events
            call.update(kind=kind, generation=generation, offspring_count=len(mothers))
            self.engine_calls = self.engine_calls + [call]
            for indices in validate_genome_map(loci).values():
                for left, right in zip(indices, indices[1:]):
                    switches[:, :, right] = homolog[:, :, left] != homolog[:, :, right]
        values, _ = self.phenotype(child, kind)
        rows = []
        for i, (m, f) in enumerate(zip(mothers, fathers)):
            row = dict(
                generation=generation,
                kind=kind,
                species=f"{prefix}{child_labels[i] + 1}",
                offspring_id=int(child_ids[i]),
                maternal_id=int(ids[m]),
                paternal_id=int(ids[f]),
                maternal_species=f"{prefix}{labels[m] + 1}",
                paternal_species=f"{prefix}{labels[f] + 1}",
                mutations_maternal=int(np.count_nonzero(child[i, 0] != inherited[i, 0])),
                mutations_paternal=int(np.count_nonzero(child[i, 1] != inherited[i, 1])),
            )
            rows.append(row)
            for h, parent in enumerate((m, f)):
                for l, locus in enumerate(self.loci[kind]):
                    old, new = int(inherited[i, h, l]), int(child[i, h, l])
                    self.transmissions.append(
                        dict(
                            generation=generation,
                            kind=kind,
                            offspring_id=int(child_ids[i]),
                            role="maternal" if h == 0 else "paternal",
                            parent_id=int(ids[parent]),
                            chromosome_id=locus["chromosome_id"],
                            locus_id=locus["locus_id"],
                            position_cM=locus["position_cM"],
                            inherited_allele_id=old,
                            child_allele_id=new,
                            transmitted_homolog_before_mutation=int(homolog[i, h, l]),
                            crossover_count=None,
                            switch_before_locus=bool(switches[i, h, l]),
                            mutation_event_id=new if new != old else None,
                        )
                    )
        values[kind + "_genotypes"] = child
        return values, child_labels, child_ids, rows, next_id + len(mothers)

    def summary(self, state, generation, sizes):
        rows = []
        lower, upper = self.spec["phenotype_map"]["bounds"]
        mutation_counts = Counter(
            (event["species"], event["locus_id"])
            for event in self.events
            if event["generation"] == generation
        )
        for kind, count in zip(KINDS, sizes):
            pop = state.plants if kind == "plant" else state.insects
            values, raw = self.phenotype(pop.genotypes, kind)
            prefix = "P" if kind == "plant" else "H"
            for s in range(count):
                mask = pop.species == s
                n = int(mask.sum())
                common = dict(generation=generation, kind=kind, species=f"{prefix}{s + 1}", N=n)
                for l, locus in enumerate(self.loci[kind]):
                    copies = pop.genotypes[mask, :, l]
                    ids, counts = np.unique(copies, return_counts=True)
                    frequencies = dict(zip(ids.tolist(), (counts / (2 * n)).tolist())) if n else {}
                    he = 1 - sum(p * p for p in frequencies.values()) if n else np.nan
                    ho = float(np.mean(copies[:, 0] != copies[:, 1])) if n else np.nan
                    # ponytail: full locus catalogs suit short fixtures; sparse time series needed at large scale.
                    for allele in self.alleles:
                        if allele.locus_id == locus["locus_id"]:
                            rows.append(
                                dict(
                                    common,
                                    record_type="allele",
                                    locus_id=locus["locus_id"],
                                    allele_id=allele.allele_id,
                                    frequency=frequencies.get(allele.allele_id, 0) if n else np.nan,
                                    observed_heterozygosity=ho,
                                    gene_diversity=he,
                                    mutation_count=mutation_counts[
                                        common["species"], locus["locus_id"]
                                    ],
                                )
                            )
                offset = 0
                for key, trait in values.items():
                    for j in range(trait.shape[1]):
                        r = raw[mask, offset + j]
                        rows.append(
                            dict(
                                common,
                                record_type="trait",
                                trait=f"{key}_{j + 1}",
                                mean=float(trait[mask, j].mean()) if n else np.nan,
                                variance=float(trait[mask, j].var()) if n else np.nan,
                                raw_mean=float(r.mean()) if n else np.nan,
                                clipped_fraction=float(np.mean((r < lower) | (r > upper)))
                                if n
                                else np.nan,
                            )
                        )
                    offset += trait.shape[1]
        return rows

    def haplotype_summary(self, state, generation, sizes):
        rows = []
        for kind, count in zip(KINDS, sizes):
            pop = state.plants if kind == "plant" else state.insects
            for indices in validate_genome_map(self.loci[kind]).values():
                for left, right in combinations(indices, 2):
                    for species in range(count):
                        copies = pop.genotypes[pop.species == species][:, :, [left, right]].reshape(
                            -1, 2
                        )
                        pairs, counts = np.unique(copies, axis=0, return_counts=True)
                        for pair, n in zip(pairs, counts):
                            rows.append(
                                dict(
                                    generation=generation,
                                    kind=kind,
                                    species=f"{'P' if kind == 'plant' else 'H'}{species + 1}",
                                    locus_left=self.loci[kind][left]["locus_id"],
                                    locus_right=self.loci[kind][right]["locus_id"],
                                    allele_left=int(pair[0]),
                                    allele_right=int(pair[1]),
                                    frequency=float(n / len(copies)),
                                )
                            )
        return rows

    def records(self):
        def encoded(row):
            return {k: json.dumps(v) if isinstance(v, (tuple, list)) else v for k, v in row.items()}

        return dict(
            genetic_loci=[encoded(x) for x in self.spec["loci"]],
            allele_catalog=[encoded(asdict(x)) for x in self.alleles],
            mutations=[encoded(x) for x in self.events],
            gamete_transmissions=[encoded(x) for x in self.transmissions],
        )
