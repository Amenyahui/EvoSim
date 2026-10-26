"""The production SLiM path, with analytic inheritance and mutation expectations."""

import copy
import json
from pathlib import Path

import numpy as np

from evo_sim.models.plant_insect.helpers.genetics import (
    GeneticState,
    validate_genome_map,
)
from evo_sim.models.plant_insect.helpers.slim_backend import reproduce
from tests.support import config, raises


def test_segregation_and_recombination():
    n = 10000
    loci = json.loads((Path(__file__).parent / "fixtures/linked_map.json").read_text())
    for locus in loci:
        locus.update(mutation_rate=0.0, mutation_mean=[0.0], mutation_sd=[0.0])
    parents = np.array([[[0, 2, 4], [1, 3, 5]]] * 2)
    settings = config().genetics["slim"]

    def draw(seed):
        material, _ = reproduce(
            parents,
            np.zeros(n, int),
            np.ones(n, int),
            loci,
            np.zeros((6, 1)),
            settings,
            [seed, seed + 1],
            6,
        )
        return material[..., 5].astype(int).reshape(-1, 3)

    sources = draw(112)
    assert abs(sources[:, 1].mean() - 0.5) < 0.015
    r = -np.expm1(-0.2) / 2
    for left, right in ((1, 2), (2, 0)):
        observed = np.mean(sources[:, left] != sources[:, right])
        assert abs(observed - r) < 6 * np.sqrt(r * (1 - r) / (2 * n))
        print(f"SLiM 10cM: N={2 * n}, observed={observed:.6f}, expected={r:.6f}")
    a, b = sources[:, 1] != sources[:, 2], sources[:, 2] != sources[:, 0]
    assert abs(np.mean(a & b) - r * r) < 6 * np.sqrt(r * r * (1 - r * r) / (2 * n))
    genotype = sources.reshape(n, 2, 3)[:, :, 0].sum(axis=1)
    assert np.all(np.abs(np.bincount(genotype) / n - [0.25, 0.5, 0.25]) < 0.02)
    for i, locus in enumerate(loci):
        locus["chromosome_id"] = str(i)
    sources = draw(114)
    _, counts = np.unique(sources[:, :2], axis=0, return_counts=True)
    assert np.all(np.abs(counts / (2 * n) - 0.25) < 0.015)
    for bad in (-1, np.nan, np.inf, True, "1"):
        raises(lambda: validate_genome_map([dict(chromosome_id="c", position_cM=bad)]))
    raises(lambda: validate_genome_map([loci[0], loci[0]]))


def test_fertilization_mutation_and_phenotype():
    c = config()
    c.genetics["mutation_enabled"]["plant"] = False
    genetic = GeneticState(c)
    g = genetic.founders([2, 0], "plant")
    for l, locus in enumerate(genetic.loci["plant"]):
        alleles = [a.allele_id for a in genetic.alleles if a.locus_id == locus["locus_id"]]
        g[0, :, l] = alleles[0]
        g[1, :, l] = alleles[1]
    child = genetic.fertilize(g, np.array([0, 0]), np.array([0, 1]), [0], [1], "plant", 2, 1)[0][
        "plant_genotypes"
    ]
    np.testing.assert_array_equal(child[0, 0], g[0, 0])
    np.testing.assert_array_equal(child[0, 1], g[1, 0])
    assert not genetic.events
    c.genetics["mutation_enabled"]["plant"] = True
    for locus in c.genetics["loci"]:
        locus["mutation_rate"] = 1.0
    genetic = GeneticState(c)
    before = copy.deepcopy(genetic.alleles)
    result = genetic.fertilize(g, np.array([0, 0]), np.array([0, 1]), [0], [1], "plant", 2, 1)[0]
    assert len(genetic.events) == 16 and len({e["allele_id"] for e in genetic.events}) == 16
    assert genetic.alleles[: len(before)] == before
    values, raw = genetic.phenotype(result["plant_genotypes"], "plant")
    expected = np.array(c.genetics["background"]["plant"]) + sum(
        np.asarray(genetic.alleles[int(i)].effect) / 2 for i in result["plant_genotypes"].flat
    )
    np.testing.assert_allclose(raw[0], expected)
    np.testing.assert_array_equal(values["D"][0], np.clip(expected, 0, 1)[:2])
    raises(
        lambda: genetic.fertilize(g, np.array([0, 1]), np.array([0, 1]), [0], [1], "plant", 2, 1)
    )


def test_per_offspring_mutation_probability():
    c = config()
    genetic = GeneticState(c)
    g = genetic.founders([2, 0], "plant")
    n = 3000
    result = genetic.fertilize(
        g, np.array([0, 0]), np.array([0, 1]), np.zeros(n, int), np.ones(n, int), "plant", 2, 1
    )
    observed = np.mean([r["mutations_maternal"] + r["mutations_paternal"] > 0 for r in result[3]])
    expected = 1 - 0.99**16
    assert abs(observed - expected) < 6 * np.sqrt(expected * (1 - expected) / n)
    print(f"QTL mutation probability: N={n}, observed={observed:.6f}, expected={expected:.6f}")
    increments = np.array(
        [np.asarray(e["delta"])[np.asarray(e["delta"]) != 0][0] for e in genetic.events]
    )
    assert abs(increments.mean()) < 6 * 0.02 / np.sqrt(len(increments))
    assert abs(increments.var() / 0.02**2 - 1) < 6 * np.sqrt(2 / len(increments))
    c.genetics["mutation_enabled"]["plant"] = False
    control = GeneticState(c)
    control.fertilize(
        g, np.array([0, 0]), np.array([0, 1]), np.zeros(n, int), np.ones(n, int), "plant", 2, 1
    )
    assert [
        (r["inherited_allele_id"], r["transmitted_homolog_before_mutation"])
        for r in genetic.transmissions
    ] == [
        (r["child_allele_id"], r["transmitted_homolog_before_mutation"])
        for r in control.transmissions
    ]
