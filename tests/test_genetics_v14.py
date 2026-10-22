"""G01-G16 implementation checks; no ecological hypothesis tests."""

import copy
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect import Config, simulate
from evo_sim.models.plant_insect.audit import audit_run
from evo_sim.models.plant_insect.batch import save_result
from evo_sim.models.plant_insect.helpers.genetics import GeneticState

ROOT = Path(__file__).resolve().parents[1]


def config(homogeneous=False, rate=None, **changes):
    values = json.loads((ROOT / "tests/fixtures/qtl_development.json").read_text())["base"]
    values.update(generations=3)
    values.update(changes)
    if homogeneous:
        for loci in values["genetics"]["founder_frequencies"].values():
            for entries in loci.values():
                entries[0]["frequency"] = 1.0
                entries[1]["frequency"] = 0.0
    if rate is not None:
        for locus in values["genetics"]["loci"]:
            locus["mutation_rate"] = rate
    return Config(**values)


def checked(c):
    r = simulate(c)
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "run"
        save_result(r, path)
        audit_run(path)
    return r


def assert_same(a, b, full=False):
    names = ["events", "plants", "parents", "populations", "networks", "diagnostics"]
    if full:
        names += ["offspring", "allele_catalog", "mutations", "genetic_summary"]
    for name in names:
        pd.testing.assert_frame_equal(
            pd.DataFrame(a[name]), pd.DataFrame(b[name]), check_exact=True
        )
    assert len(a["snapshots"]) == len(b["snapshots"])
    for x, y in zip(a["snapshots"], b["snapshots"]):
        for key in x.keys() & y.keys():
            if not full and key.endswith("_genotypes"):
                continue
            np.testing.assert_array_equal(x[key], y[key])


def test_reproducibility_legacy_G01_G02_G12_G13():
    c = config(homogeneous=True, rate=0)
    new = checked(c)
    old = copy.deepcopy(c)
    old.genetics = None
    old.evolve_plants = old.evolve_insects = False
    assert_same(new, simulate(old))
    # Neither legacy mutation kernel may be called in QTL mode.
    with (
        patch(
            "evo_sim.models.plant_insect.helpers.traits.mutate",
            side_effect=AssertionError("double mutation"),
        ),
        patch(
            "evo_sim.models.plant_insect.helpers.voc.mutate_Z",
            side_effect=AssertionError("double mutation"),
        ),
    ):
        assert_same(new, simulate(copy.deepcopy(c)), full=True)
    zero = copy.deepcopy(c)
    for locus in zero.genetics["loci"]:
        locus["mutation_rate"] = 1.0
        locus["mutation_sd"] = [0.0] * 4
    r = checked(zero)
    assert len(r["mutations"]) > 0
    assert_same(new, r)
    changing = config()
    a = checked(changing)
    assert_same(a, simulate(copy.deepcopy(changing)), full=True)
    # Observation recording cannot consume genetic/ecological random draws.
    changing.record_observations = False
    b = simulate(changing)
    for key in ("events", "parents", "mutations", "allele_catalog"):
        before, after = pd.DataFrame(a[key]), pd.DataFrame(b[key])
        pd.testing.assert_frame_equal(before[after.columns], after, check_exact=True)


def test_mapping_inheritance_G03_G04_G05_G06_G08_G09():
    c = config(rate=1)
    c.validate()
    g = GeneticState(c)
    parents = g.founders([2, 2], "plant")
    saved = parents.copy()
    values, raw = g.phenotype(parents, "plant")
    manual = np.zeros_like(raw)
    for i in range(len(parents)):
        for h in range(2):
            for l in range(parents.shape[2]):
                manual[i] += np.array(g.alleles[int(parents[i, h, l])].effect) / 2
    np.testing.assert_array_equal(raw, manual)
    equal, _ = g.phenotype(np.repeat(parents[:1], 2, axis=0), "plant")
    np.testing.assert_array_equal(equal["D"][0], equal["D"][1])
    high = copy.deepcopy(c)
    high.genetics["background"]["plant"] = [2.0, -2.0, 2.0, -2.0]
    h = GeneticState(high)
    bounded, unbounded = h.phenotype(parents, "plant")
    np.testing.assert_array_equal(bounded["D"], np.tile([1.0, 0.0], (4, 1)))
    assert (unbounded > 1).any() and (unbounded < 0).any()
    offspring, labels, ids, rows, _ = g.inherit(
        parents, np.array([0, 0, 1, 1]), np.arange(4), np.array([0, 1, 2, 0]), "plant", 4, 1
    )
    np.testing.assert_array_equal(parents, saved)
    assert labels.tolist() == [0, 1, 1] and ids.tolist() == [4, 5, 6]
    assert [r["parent_id"] for r in rows] == [1, 2, 2]
    assert all(r["locus_mutation_count"] == 2 * parents.shape[2] for r in rows)
    assert len({a.allele_id for a in g.alleles}) == len(g.alleles)
    assert not np.array_equal(offspring["plant_genotypes"][1], offspring["plant_genotypes"][2])
    off = config(rate=0)
    off.validate()
    q = GeneticState(off)
    pg = q.founders([100, 0], "plant")
    trait, _ = q.phenotype(pg, "plant")
    assert trait["D"].var(axis=0).max() > 0  # Same species can have distinct phenotypes.
    # With no mutation, unequal reproductive success changes allele frequencies.
    counts = np.zeros(100, int)
    counts[0] = 100
    child, *_ = q.inherit(pg, np.zeros(100, int), np.arange(100), counts, "plant", 100, 1)
    assert not np.array_equal(np.sort(pg.ravel()), np.sort(child["plant_genotypes"].ravel()))
    assert not q.events


def test_mutation_probability_G10():
    c = config(homogeneous=True, rate=0.02)
    c.validate()
    g = GeneticState(c)
    parent = g.founders([1, 0], "plant")
    child, *_ = g.inherit(parent, np.array([0]), np.array([0]), np.array([4000]), "plant", 1, 1)
    observed = (child["plant_genotypes"] != parent).any(axis=(1, 2)).mean()
    expected = 1 - (1 - 0.02) ** (2 * parent.shape[2])
    assert abs(observed - expected) < 0.035


def test_boundaries_G07_G14_G15_G16():
    for changes in (
        dict(plant_sizes=[0, 0]),
        dict(insect_sizes=[0]),
        dict(plant_sizes=[0, 0], insect_sizes=[0]),
        dict(generations=0),
        dict(max_individuals_safety=1),
    ):
        checked(config(**changes))
    c = config()
    c.genetics["max_genotype_values_safety"] = 1
    assert checked(c)["metadata"]["status"]["status"] == "computational_limit"
    c = config()
    c.genetics["max_genotype_values_safety"] = 2 * 8 * 46
    r = checked(c)
    assert r["metadata"]["status"]["status"] == "computational_limit" and not r["offspring"]
    c = config()
    c.genetics["max_alleles_safety"] = len(c.genetics["alleles"])
    assert checked(c)["metadata"]["status"]["status"] == "computational_limit"
    c = config()
    c.insect_sizes = [3, 3]
    c.insect_offspring_lambda["H2"] = c.insect_offspring_lambda["H1"]
    c.max_feeding_demand["H2"] = c.max_feeding_demand["H1"]
    c.parameter_provenance["insect_offspring_lambda"]["H2"] = copy.deepcopy(
        c.parameter_provenance["insect_offspring_lambda"]["H1"]
    )
    for key in ("I", "A"):
        c.founder_traits[key].append([0.6, 0.2])
    c.genetics["founder_frequencies"]["H2"] = copy.deepcopy(c.genetics["founder_frequencies"]["H1"])
    r = checked(c)
    assert {x["species"] for x in r["populations"]} == {"P1", "P2", "H1", "H2"}
    assert "no segregation/recombination/mating/speciation" in r["metadata"]["inheritance"]


def test_corruption_G11():
    c = config()
    r = simulate(c)
    with TemporaryDirectory() as tmp:
        for name in ("allele_catalog.csv", "mutations.csv", "genetic_summary.csv"):
            path = Path(tmp) / name.replace(".csv", "")
            save_result(r, path)
            table = pd.read_csv(path / name, float_precision="round_trip")
            if name == "allele_catalog.csv":
                i = len(c.genetics["alleles"])
                effect = json.loads(table.loc[i, "effect"])
                effect[0] += 0.1
                table.loc[i, "effect"] = json.dumps(effect)
            elif name == "mutations.csv":
                table.loc[0, "parent_allele_id"] = 999999
            else:
                i = table[table.record_type == "trait"].index[0]
                table.loc[i, "mean"] += 0.1
            table.to_csv(path / name, index=False)
            meta = json.loads((path / "metadata.json").read_text())
            meta["output_sha256"][name] = hashlib.sha256((path / name).read_bytes()).hexdigest()
            (path / "metadata.json").write_text(json.dumps(meta))
            try:
                audit_run(path)
            except AssertionError:
                pass
            else:
                raise AssertionError(f"accepted corrupted {name} despite rehashed file")


def test_rectangular_traits_G14():
    c = config(generations=1)
    c.n = 3
    c.B = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    for key in ("D", "I"):
        for row in c.founder_traits[key]:
            row.append(0.4)
    g = c.genetics
    for x in g["background"].values():
        x.insert(2, 0.0)
    for a in g["alleles"]:
        a["effect"].insert(2, 0.0)
    for locus in g["loci"]:
        locus["mutation_sd"].insert(2, 0.0)
        locus["targets"] = [j + 1 if j >= 2 else j for j in locus["targets"]]
    for kind, species in (("plant", ["P1", "P2"]), ("insect", ["H1"])):
        key = kind + "_extra_dimension"
        g["loci"].append(
            dict(
                locus_id=key,
                organism=kind,
                functional_module="extra_QTL",
                targets=[2],
                mutation_rate=0.01,
                mutation_sd=[0.0, 0.0, 0.02, 0.0, 0.0],
            )
        )
        idx = len(g["alleles"])
        g["alleles"].append(dict(allele_id=idx, locus_id=key, effect=[0.0, 0.0, 0.4, 0.0, 0.0]))
        for name in species:
            g["founder_frequencies"][name][key] = [dict(allele_id=idx, frequency=1.0)]
    r = checked(c)
    assert r["snapshots"][0]["D"].shape[1] == 3
    assert r["snapshots"][0]["A"].shape[1] == 2


def test_validation():
    edits = [
        lambda g: g.update(ploidy=1),
        lambda g: g["alleles"][0].update(allele_id=-1),
        lambda g: g["loci"][0].update(mutation_rate=float("nan")),
        lambda g: g["loci"][0].update(mutation_sd=[0.0, 1.0, 0.0, 0.0]),
        lambda g: g["alleles"][0].update(effect=[0.0, 1.0, 0.0, 0.0]),
        lambda g: g["founder_frequencies"]["P1"].pop("plant_1_1"),
        lambda g: g.update(unknown_option=True),
    ]
    for edit in edits:
        c = config()
        edit(c.genetics)
        try:
            c.validate()
        except ValueError:
            pass
        else:
            raise AssertionError("invalid genetic config accepted")


if __name__ == "__main__":
    for name, test in list(globals().items()):
        if name.startswith("test_"):
            test()
            print("PASS", name)
