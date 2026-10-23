"""Fixed seeds; 100,000 gametes, six-standard-error statistical tolerances."""

import copy
import json
from pathlib import Path

import numpy as np

from evo_sim.models.plant_insect.helpers.genetics import (
    GeneticState,
    form_gametes,
    interval_crossovers,
    validate_genome_map,
)
from tests.support import config, raises


def test_segregation_and_recombination():
    n = 100000
    loci = json.loads((Path(__file__).parent / "fixtures/linked_map.json").read_text())
    parents = np.broadcast_to([[0, 0, 0], [1, 1, 1]], (n, 2, 3))
    gametes, sources, counts = form_gametes(parents, loci, np.random.default_rng(112))
    assert np.array_equal(gametes, sources)
    assert abs(gametes[:, 1].mean() - 0.5) < 0.01
    r = -np.expm1(-0.2) / 2
    for left, right in ((1, 2), (2, 0)):
        observed = np.mean(gametes[:, left] != gametes[:, right])
        assert abs(observed - r) < 6 * np.sqrt(r * (1 - r) / n)
        assert np.array_equal(sources[:, right], sources[:, left] ^ (counts[:, right] % 2))
        print(f"10cM: N={n}, r={observed:.6f}, theoretical={r:.6f}")
    assert abs(np.corrcoef(counts[:, 0], counts[:, 2])[0, 1]) < 0.02
    loci[0]["position_cM"] = 10000.0
    gametes, _, _ = form_gametes(parents, loci, np.random.default_rng(113))
    assert abs(np.mean(gametes[:, 0] != gametes[:, 2]) - 0.5) < 0.01
    assert not interval_crossovers(0, n, np.random.default_rng(1)).any()
    for i, locus in enumerate(loci):
        locus["chromosome_id"] = str(i)
    gametes, _, _ = form_gametes(parents, loci, np.random.default_rng(114))
    _, frequency = np.unique(gametes[:, :2], axis=0, return_counts=True)
    assert np.all(np.abs(frequency / n - 0.25) < 0.01)
    homo = np.full((10, 2, 3), 7)
    assert np.all(form_gametes(homo, loci, np.random.default_rng(1))[0] == 7)


def test_parity_and_map_validation():
    class Fixed:
        def integers(self, *args, **kwargs):
            return np.zeros(4, dtype=int)

        def poisson(self, *args, **kwargs):
            return np.arange(4)

    loci = [dict(chromosome_id="c", position_cM=x) for x in (0.0, 10.0)]
    parents = np.broadcast_to([[0, 0], [1, 1]], (4, 2, 2))
    gametes, _, counts = form_gametes(parents, loci, Fixed())
    assert gametes[:, 1].tolist() == [0, 1, 0, 1]
    assert counts[:, 1].tolist() == [0, 1, 2, 3]
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
    # Aa x Aa -> 1:2:1; pure gametes avoid enormous output ledgers.
    hetero = np.repeat(child, 100000, axis=0)
    a, _, _ = form_gametes(hetero, genetic.loci["plant"], np.random.default_rng(18))
    b, _, _ = form_gametes(hetero, genetic.loci["plant"], np.random.default_rng(19))
    genotype = (a[:, 0] == g[0, 0, 0]).astype(int) + (b[:, 0] == g[0, 0, 0])
    assert np.all(np.abs(np.bincount(genotype) / len(genotype) - [0.25, 0.5, 0.25]) < 0.01)
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
