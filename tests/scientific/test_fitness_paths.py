"""Small deterministic causal probes, with native inheritance and fixed event uniforms.

These assays isolate paths; their controlled reproduction draws are not ecological
parameter calibration or evidence that a trait is favored in every environment.
"""

import copy
from types import SimpleNamespace

import numpy as np

from evo_sim.models.plant_insect.biology.agents import PlantPopulation
from evo_sim.models.plant_insect.biology.genetics.qtl import GeneticState
from evo_sim.models.plant_insect.biology.insects.feeding import (
    allocate_resources,
    chemical_survival,
)
from evo_sim.models.plant_insect.biology.insects.host_selection import (
    choose_hosts,
    contact_acceptance,
)
from evo_sim.models.plant_insect.biology.insects.sensory import choice_probabilities
from evo_sim.models.plant_insect.biology.plants.defense import fecundity
from evo_sim.models.plant_insect.biology.plants.demography import plant_resource_update
from evo_sim.models.plant_insect.biology.traits import maternal_fecundity, scores
from evo_sim.models.plant_insect.records.contributions import contributions
from evo_sim.models.plant_insect.simulation.stages.adult_reproduction import sample_parents
from tests.support import config


def genetic_assay(kind, axis):
    """Three parents: low, high, low donor; only one inherited coordinate differs."""
    c = config(plant_sizes=[3, 0], insect_sizes=[3])
    c.genetics["mutation_enabled"] = dict(plant=False, insect=False)
    own = [l for l in c.genetics["loci"] if l["organism"] == kind and axis in l["targets"]]
    alleles = [a for a in c.genetics["alleles"] if a["locus_id"] == own[0]["locus_id"]]
    for entry in c.genetics["alleles"]:
        if entry["locus_id"] in {l["locus_id"] for l in own}:
            entry["effect"][axis] = 0.0
    alleles[1]["effect"][axis] = 1 / (2 * c.genetics["phenotype_map"]["effect_scale"])
    c.validate()
    genetic = GeneticState(c)
    g = genetic.founders([3, 0] if kind == "plant" else [3], kind)
    g[:] = g[0]
    locus = next(i for i, l in enumerate(genetic.loci[kind]) if l["locus_id"] == own[0]["locus_id"])
    g[:, :, locus] = alleles[0]["allele_id"]
    g[1, :, locus] = alleles[1]["allele_id"]
    values = genetic.phenotype(g, kind)[0]
    key = "plant_chemical_traits" if kind == "plant" else "insect_module_traits"
    np.testing.assert_array_equal(values[key][:, axis], [0, 1, 0])
    np.testing.assert_array_equal(np.delete(values[key][0], axis), np.delete(values[key][1], axis))
    return c, genetic, g, values, locus, alleles[1]["allele_id"]


def inherit(genetic, genome, kind, mothers):
    """Native SLiM transmission from two assay mothers to a common low donor."""
    genetic = copy.deepcopy(genetic)
    values, _, _, _, _ = genetic.fertilize(
        genome,
        np.zeros(3, int),
        np.arange(3),
        np.asarray(mothers, int),
        np.full(len(mothers), 2, int),
        kind,
        3,
        1,
    )
    return values


def plant_contribution(c, genetic, genome, values, consumed):
    """Floor Poisson means and zero survival uniforms isolate the maternal path."""

    class Stream:
        def random(self, size):
            return np.zeros(size)

        def poisson(self, mean):
            return np.floor(mean).astype(int)

    chemical = values["plant_chemical_traits"]
    _, _, survival, recovery = plant_resource_update(
        np.ones(3), np.asarray(consumed), np.ones(3), c.plant_demography, chemical[:, 2]
    )
    factor = fecundity(values["D"], chemical[:, 3:5], c, chemical)[0] * recovery
    factor[2] = 0  # Assay donor has no maternal output, but remains available as a father.
    pop = PlantPopulation(np.arange(3), np.zeros(3, int), values["D"], cells=np.zeros((3, 2), int))
    runtime = SimpleNamespace(
        c=c,
        lp=np.array([100, 100]),
        rng={
            2: {0: [Stream(), Stream()]},
            3: {0: [Stream(), Stream()]},
        },
    )
    counts, _, rows, _ = sample_parents(pop, "plant", survival, factor, runtime, 0)
    mothers = np.repeat(np.arange(3), counts)
    children = inherit(genetic, genome, "plant", mothers)
    contributions(rows, mothers, np.full(len(mothers), 2, int), np.ones(len(mothers), bool))
    assert sum(row["realized_genetic_contribution"] for row in rows) == len(mothers)
    return counts, children, rows


def test_psm_genetic_cost_contact_and_larval_survival_paths():
    c, genetic, g, values, locus, allele = genetic_assay("plant", 0)
    chemical, psm = values["plant_chemical_traits"], values["D"]
    female = np.full(10, 0.5)
    female[8] = 0
    acceptance = contact_acceptance(psm[:2], female, c.contact_response)
    assert acceptance[0] > acceptance[1]
    exposure = chemical_survival(
        psm[:2], np.zeros((2, 2)), np.ones(2), np.ones(2), c.chemical_processing
    )
    assert exposure["survival_probability"][0] > exposure["survival_probability"][1]
    cost = fecundity(psm, chemical[:, 3:5], c, chemical)[1]
    assert cost[1] > cost[0]
    counts, children, rows = plant_contribution(c, genetic, g, values, [0, 0, 0])
    assert counts[0] > counts[1] > 0 and rows[0]["W_i"] > rows[1]["W_i"]
    np.testing.assert_allclose(
        (children["plant_genotypes"][:, :, locus] == allele).mean(), counts[1] / (2 * counts.sum())
    )


def test_voc_genetic_signal_changes_host_load_and_plant_contribution():
    c, genetic, g, values, _, _ = genetic_assay("plant", 3)
    female = np.full(10, 0.5)
    female[4:8] = [1, 0, 0, 1]
    volatile = values["plant_chemical_traits"][:2, 3:5]
    probability = choice_probabilities(
        scores(volatile, female, c), c.tau, c.outside_utility, "outside"
    )
    equal = choice_probabilities(
        scores(np.repeat(volatile[[0]], 2, axis=0), female, c), c.tau, c.outside_utility, "outside"
    )
    assert probability[1] > equal[1] and probability[0] < equal[0]
    uniform = (probability[0] + equal[0]) / 2
    selected = [np.searchsorted(p.cumsum(), uniform) for p in (equal, probability)]
    assert selected == [0, 1]  # The same event uniform redirects this egg.
    outcomes = []
    for host in selected:
        allocation = allocate_resources(np.array([host]), np.ones(3), np.array([0.5]))
        outcomes.append(plant_contribution(c, genetic, g, values, allocation["X"]))
    assert outcomes[0][0][1] > outcomes[1][0][1]
    assert (
        outcomes[0][2][1]["realized_genetic_contribution"]
        > outcomes[1][2][1]["realized_genetic_contribution"]
    )


def test_z_genetic_damage_response_changes_newborn_gene_frequency():
    c, genetic, g, values, locus, allele = genetic_assay("plant", 2)
    undamaged, _, _ = plant_contribution(c, genetic, g, values, [0, 0, 0])
    assert undamaged[0] == undamaged[1]
    counts, children, rows = plant_contribution(c, genetic, g, values, [0.5, 0.5, 0])
    assert rows[0]["S_i"] == rows[1]["S_i"] and counts[1] > counts[0]
    assert rows[1]["realized_genetic_contribution"] > rows[0]["realized_genetic_contribution"]
    assert (children["plant_genotypes"][:, :, locus] == allele).mean() > 0.25


def test_i_genetic_transmission_changes_recruitment_without_changing_food():
    c, genetic, g, _, locus, allele = genetic_assay("insect", 0)
    children = inherit(genetic, g, "insect", [0, 1])
    allocation = allocate_resources(np.array([0, 0]), [2], [1, 1])
    psm = np.ones((2, 2))
    result = chemical_survival(
        psm, children["I"], allocation["x"], np.ones(2), c.chemical_processing
    )
    assert result["residual"][0].sum() > result["residual"][1].sum()
    probability = result["survival_probability"]
    assert probability[1] > probability[0]
    survived = np.full(2, probability.mean()) < probability
    np.testing.assert_array_equal(survived, [False, True])
    rows = [{}, {}, {}]
    contributions(rows, np.array([0, 1]), np.array([2, 2]), survived)
    assert (
        rows[0]["realized_genetic_contribution"] == rows[1]["realized_genetic_contribution"] == 0.5
    )
    assert rows[0]["realized_recruitment_contribution"] == 0
    assert rows[1]["realized_recruitment_contribution"] == 0.5
    assert (children["insect_genotypes"][survived, :, locus] == allele).mean() == 0.5
    np.testing.assert_array_equal(allocation["d"], [1, 1])
    np.testing.assert_array_equal(result["dose"][0], result["dose"][1])


def test_sensory_genetic_expression_and_tuning_have_distinct_paths():
    for axis in (4, 6):
        c, genetic, g, values, _, _ = genetic_assay("insect", axis)
        traits = values["insect_module_traits"]
        volatile = np.array([[1, 0], [0, 1]], float)
        probabilities = [
            choice_probabilities(scores(volatile, row, c), c.tau, c.outside_utility, "outside")
            for row in traits[:2]
        ]
        assert not np.allclose(probabilities[0], probabilities[1])
        factor = maternal_fecundity(traits, c)[0]
        if axis == 4:
            assert factor[0] > factor[1]
        else:
            assert factor[0] == factor[1]
        # The same search uniform redirects an egg; the same survival uniform
        # distinguishes a chemically benign host from a toxic one.
        threshold = (probabilities[0][0] + probabilities[1][0]) / 2
        hosts = np.array([np.searchsorted(p.cumsum(), threshold) for p in probabilities])
        mothers = np.flatnonzero(hosts < 2)
        children = inherit(genetic, g, "insect", mothers)
        probability = chemical_survival(
            np.array([[0, 0], [1, 1]])[hosts[mothers]],
            children["I"],
            np.ones(len(mothers)),
            np.ones(len(mothers)),
            c.chemical_processing,
        )["survival_probability"]
        survived = 0.5 < probability
        assert np.count_nonzero(survived) == 1
        rows = [{}, {}, {}]
        contributions(rows, mothers, np.full(len(mothers), 2), survived)
        assert rows[int(mothers[survived][0])]["realized_recruitment_contribution"] == 0.5


def test_g_genetic_female_acceptance_changes_birth_contribution():
    c, genetic, g, values, locus, allele = genetic_assay("insect", 8)
    c.choice_mode = "forced"

    class Stream:
        def choice(self, population, size=None, replace=True, p=None):
            return np.asarray(population)[:size] if size is not None else 0

        def random(self):
            return 0.5

    traits = values["insect_module_traits"]
    found = choose_hosts(
        traits[:2],
        np.zeros(1, int),
        np.zeros(2, int),
        np.zeros((1, 2)),
        c,
        None,
        np.random.default_rng(1),
        [Stream()],
        candidate_pools=[np.array([0])] * 2,
        PSM=np.array([[1, 0]]),
    )
    np.testing.assert_array_equal(found["host_found"], [False, True])
    mothers = np.flatnonzero(found["host_found"])
    children = inherit(genetic, g, "insect", mothers)
    rows = [{}, {}, {}]
    contributions(rows, mothers, np.full(len(mothers), 2), np.ones(len(mothers), bool))
    assert rows[0]["realized_genetic_contribution"] == 0
    assert (
        rows[1]["realized_genetic_contribution"] == rows[2]["realized_genetic_contribution"] == 0.5
    )
    assert (children["insect_genotypes"][:, :, locus] == allele).mean() == 0.5
    np.testing.assert_array_equal(
        maternal_fecundity(traits, c)[0][0], maternal_fecundity(traits, c)[0][1]
    )
