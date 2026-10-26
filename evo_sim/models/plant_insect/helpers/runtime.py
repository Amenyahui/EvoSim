import hashlib
import json
import platform
from copy import deepcopy
from dataclasses import fields
from types import SimpleNamespace

import numpy as np
import pandas as pd

from evo_sim import __version__
from evo_sim.core.helpers import source_hashes

from . import ecology as shared
from . import genetics, spatial, voc
from .agents import CommunityState, InsectPopulation, PlantPopulation
from .config import Config, config_payload, resource_values
from .contributions import mark_uncommitted
from .demography import ComputationalLimit
from .slim_backend import SlimError


def streams(c):
    """Fixed stage/taxon/species IDs, so demographic draws never consume ecological draws."""

    def rng(phase, kind, s):
        return np.random.default_rng(np.random.SeedSequence([c.seed, 130, phase, kind, s]))

    return {
        phase: (
            [rng(phase, 0, s) for s in range(len(c.plant_sizes))],
            [rng(phase, 1, s) for s in range(len(c.insect_sizes))],
        )
        for phase in range(8)
    }


def resource_vector(value, sizes, labels, prefix):
    per = resource_values(value, [1] * len(sizes), prefix)
    return np.broadcast_to(per, (len(sizes),))[labels]


def ecology(state, c, B, w, rng_p, rng_H):
    """F01–F06,F15: same kernel, dynamic labels; empty arrays remain valid."""
    D = state.plants.defense
    I = state.insects.resistance
    A = state.insects.preference
    trait = state.plants.voc_trait
    pl, hl = state.plants.species, state.insects.species
    V_emit = voc.map_voc(D, trait, B, c.rho_voc, c.voc_bounds)
    V_cue = V_emit[rng_p.permutation(len(D))] if c.shuffle_voc else V_emit.copy()
    options = {}
    if c.spatial is not None:
        options = dict(
            candidate_pools=spatial.candidate_pools(state.plants.cells, state.insects.cells, c),
        )
    choice = shared.choose_hosts(A, pl, hl, V_cue, c, w, rng_p, rng_H, **options)
    capacity = resource_vector(c.resource_capacity, c.plant_sizes, pl, "P")
    dmax = resource_vector(c.max_feeding_demand, c.insect_sizes, hl, "H")
    e = shared.settle_feeding(D, I, pl, hl, V_emit, choice["host"], c, capacity, dmax)
    # Explicit Speed beta=0 and common actual resource ledger. No fixed-size reproduction.
    e["wi"] = e["g"] * e["fi"]
    e.update(
        choice,
        V_emit=V_emit,
        V_cue=V_cue,
        R=capacity,
        dmax=dmax,
        S_plant=c.sP0 + (1 - c.sP0) * (1 - e["h"]),
        S_insect=e["g"],
    )
    return e


def mean_or_na(values):
    return float(np.mean(values)) if len(values) else np.nan


def agent_state(community):
    if community is None:
        return {"plants": {}, "insects": {}}
    return {
        name: {f.name: getattr(population, f.name) for f in fields(population)}
        for name, population in (("plants", community.plants), ("insects", community.insects))
    }


def community_state(agents):
    return CommunityState(
        PlantPopulation(**agents["plants"]), InsectPopulation(**agents["insects"])
    )


def initialize_community(config, initial=None):
    runtime = SimpleNamespace(
        c=Config(
            **dict(
                config["parameters"],
                generations=config["simulation_metadata"]["num_steps_per_episode"],
            )
        )
    )
    state = {
        "current_step": 0,
        "current_substep": "0",
        "terminated": False,
        "environment": {"runtime": runtime},
        "agents": {},
        "objects": {},
        "network": {},
    }
    community = None
    runtime.input_sources = deepcopy(config.get("input_sources", {}))
    runtime.B, runtime.w = runtime.c.validate()
    runtime.genetic = genetics.GeneticState(runtime.c)
    if initial is not None:
        raise ValueError("initialization comes only from the founder allele catalog")
    runtime.observe = runtime.c.record_observations
    runtime.reproduction_checkpoint = None
    runtime.lp, runtime.li, runtime.alpha = runtime.c.resolved_demography()
    runtime.rng = streams(runtime.c)
    runtime.result = dict(
        events=[],
        plants=[],
        parents=[],
        offspring=[],
        populations=[],
        networks=[],
        snapshots=[],
        allele_contributions=[],
    )
    if runtime.observe:
        runtime.result["diagnostics"] = []
    runtime.result["genetic_summary"] = []
    runtime.result["haplotype_summary"] = []
    runtime.stored_genotype_values = 0
    runtime.status = dict(status="completed", generation=0, message="")
    runtime.stage = "initialization"
    runtime.extinct = {}
    try:
        if (
            sum(runtime.c.plant_sizes) + sum(runtime.c.insect_sizes)
            > runtime.c.max_individuals_safety
        ):
            raise ComputationalLimit(
                "initial population exceeds safety limit; no individuals allocated"
            )
        initial_values = sum(
            (
                2 * len(runtime.genetic.loci[kind]) * sum(sizes)
                for kind, sizes in zip(
                    genetics.KINDS, (runtime.c.plant_sizes, runtime.c.insect_sizes)
                )
            )
        )
        if initial_values > runtime.c.genetics["max_genotype_values_safety"]:
            raise ComputationalLimit("initial genotype allocation exceeds safety limit")
        gp = runtime.genetic.founders(runtime.c.plant_sizes, "plant")
        gh = runtime.genetic.founders(runtime.c.insect_sizes, "insect")
        pvalues, _ = runtime.genetic.phenotype(gp, "plant")
        hvalues, _ = runtime.genetic.phenotype(gh, "insect")
        initial = (
            pvalues["D"],
            hvalues["I"],
            hvalues["A"],
            pvalues["plant_voc_independent_trait"],
        )
        D, I, A, plant_voc_independent_trait = initial
        pl = np.repeat(np.arange(len(runtime.c.plant_sizes)), runtime.c.plant_sizes)
        hl = np.repeat(np.arange(len(runtime.c.insect_sizes)), runtime.c.insect_sizes)
        sex = np.empty(len(hl), dtype=np.int8)
        for species, rng in enumerate(runtime.rng[6][1]):
            mask = hl == species
            sex[mask] = (
                rng.random(mask.sum())
                >= runtime.c.sexual_reproduction["insect_founder_female_probability"]
            )
        community = CommunityState(
            PlantPopulation(np.arange(len(pl)), pl, D, plant_voc_independent_trait, genotypes=gp),
            InsectPopulation(np.arange(len(hl)), hl, I, A, genotypes=gh, sex=sex),
        )
        if runtime.c.spatial is not None:
            community.plants.cells, community.insects.cells = spatial.initial_positions(
                runtime.c, pl, hl
            )
            runtime.spatial_rng = [
                np.random.default_rng(np.random.SeedSequence([runtime.c.seed, 150, 1, k]))
                for k in range(2)
            ]
        community.validate()
        runtime.next_ids = [len(pl), len(hl)]
    except (ComputationalLimit, MemoryError, ArithmeticError, ValueError) as error:
        fail(state, error)
    state["agents"] = agent_state(community)
    return state


def fail(state, error):
    runtime = state["environment"]["runtime"]
    mark_uncommitted(runtime.result["parents"])
    if runtime.reproduction_checkpoint is not None:
        runtime.rng, runtime.spatial_rng = runtime.reproduction_checkpoint
        runtime.reproduction_checkpoint = None
    runtime.status = dict(
        status="backend_error"
        if isinstance(error, SlimError)
        else "computational_limit"
        if isinstance(error, (ComputationalLimit, MemoryError))
        else "arithmetic_error",
        generation=state["current_step"],
        stage=runtime.stage,
        message=str(error),
    )
    state["terminated"] = True


def finalize(state):
    runtime = state["environment"]["runtime"]
    mark_uncommitted(runtime.result["parents"])
    runtime.result["metadata"] = dict(
        version=runtime.c.model_version,
        config=config_payload(runtime.c),
        input_sources=runtime.input_sources,
        status=runtime.status,
        resolved_B=runtime.B.tolist(),
        resolved_weights=runtime.w.tolist(),
        resolved_alpha=runtime.alpha.tolist(),
        resolved_plant_lambda=runtime.lp.tolist(),
        resolved_insect_lambda=runtime.li.tolist(),
        source_hash=source_hashes(),
        python=platform.python_version(),
        numpy=np.__version__,
        pandas=pd.__version__,
        snapshot_semantics="t=0 initial; each transition t to t+1 draws demography once; final t records ecology and conditional moments without drawing offspring",
        count_semantics="N_next is null when transition is not committed; N_offspring retains all sampled counts at safety stop",
        resource_semantics="edible dry mass per newly recruited plant per feeding window, not regrowth or offspring biomass conservation",
    )
    runtime.result.update(runtime.genetic.records())
    architecture = json.dumps(runtime.c.genetics, sort_keys=True, separators=(",", ":"))
    payload = json.dumps(config_payload(runtime.c), sort_keys=True, separators=(",", ":"))
    runtime.result["metadata"].update(
        software_version=__version__,
        contribution_semantics="parent generation t to offspring generation t+1; realized counts only on committed transitions; missing for observation-only or failed transitions",
        fitness_proxy_semantics="W_i/W_insect/W_plant are S*F proxies, not total fitness; realized_genetic_contribution=(maternal+paternal offspring)/2",
        maternal_expectation_semantics="all plants or all females, including zero expectation for no host/survival/mate; nonmaternal roles are missing; all-individual mean counts those roles as zero",
        allele_contribution_semantics="per-species per-locus parent/offspring allele union; parental contribution includes stochastic reproduction; transmission includes segregation and mutation",
        reference_panel_semantics="reference_p/reference_q snapshot rows align with insect_ids; columns align with reference_plant_species; p is conditional choice without rejection",
        ecology_rule_version=runtime.c.model_version,
        genetics_version=runtime.c.genetics["schema_version"],
        slim_calls=runtime.genetic.engine_calls,
        genetic_architecture_sha256=hashlib.sha256(architecture.encode()).hexdigest(),
        config_sha256=hashlib.sha256(payload.encode()).hexdigest(),
        genetic_rng="NumPy founders: SeedSequence([seed,140,0,taxon,species]); SLiM seeds: SeedSequence([seed,140,4,taxon,generation]).generate_state(2)+1; first meiosis, second Eidos QTL mutation",
        rng="SeedSequence([seed,130,phase,taxon,species]); 0 init, 1 ecology, 2 survival, 3 maternal births, 5 mate choice, 6 founder sex, 7 offspring sex",
        spatial_rng="SeedSequence([seed,150,0]) plants; [seed,150,2] insects; [seed,150,1,taxon] maternal dispersal",
        inheritance="SLiM 5.2 native addCrossed; copy 0 maternal, copy 1 paternal; separate chromosomes; compressed marker boundary r=(1-exp(-2*distance_Morgan))/2",
        mutation_semantics="per transmitted abstract QTL copy; not per nucleotide",
        crossover_semantics="switch_before_locus records marker ancestry changes; crossover_count is missing because compressed intervals do not retain full crossover counts",
        spatial_semantics="periodic Chebyshev distance; independent founder positions; maternal dispersal",
    )
    return runtime.result
