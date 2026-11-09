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
from evo_sim.core.runner import ComputationalLimit
from evo_sim.models.plant_insect.biology import spatial
from evo_sim.models.plant_insect.biology import traits as modules
from evo_sim.models.plant_insect.biology.agents import (
    CommunityState,
    InsectPopulation,
    PlantPopulation,
)
from evo_sim.models.plant_insect.biology.genetics import qtl as genetics
from evo_sim.models.plant_insect.biology.genetics.slim_backend import SlimError
from evo_sim.models.plant_insect.biology.insects import feeding
from evo_sim.models.plant_insect.config.schema import Config, config_payload, resource_values
from evo_sim.models.plant_insect.records.contributions import mark_uncommitted

from .random import streams


def resource_vector(value, sizes, labels, prefix):
    per = resource_values(value, [1] * len(sizes), prefix)
    return np.broadcast_to(per, (len(sizes),))[labels]


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
    runtime.lp, runtime.li = runtime.c.resolved_demography()
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
        oviposition=[],
        larvae=[],
        host_flows=[],
        host_groups=[],
        lagged_networks=[],
        plant_proposals=[],
    )
    if runtime.observe:
        runtime.result["diagnostics"] = []
    runtime.result["genetic_summary"] = []
    runtime.result["haplotype_summary"] = []
    runtime.stored_genotype_values = 0
    runtime.founder_feeding_diagnostics = {"status": "initialization_incomplete"}
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
            pvalues.get("plant_voc_independent_trait"),
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
            PlantPopulation(
                np.arange(len(pl)),
                pl,
                D,
                plant_voc_independent_trait,
                genotypes=gp,
                chemical_traits=pvalues["plant_chemical_traits"],
                resource=resource_vector(
                    runtime.c.resource_capacity, runtime.c.plant_sizes, pl, "P"
                ).copy(),
            ),
            InsectPopulation(
                np.arange(len(hl)),
                hl,
                I,
                A,
                genotypes=gh,
                sex=sex,
                module_traits=hvalues.get("insect_module_traits"),
                birth_host_id=np.full(len(hl), -1, dtype=np.int64),
                birth_host_species=np.full(len(hl), -1, dtype=int),
                adult_cohort=np.zeros(len(hl), dtype=int),
            ),
        )
        if runtime.c.spatial is not None:
            community.plants.cells, community.insects.cells, community.plants.microsite = (
                spatial.initial_positions(runtime.c, pl, hl)
            )
            runtime.spatial_rng = [
                np.random.default_rng(np.random.SeedSequence([runtime.c.seed, 150, 1, k]))
                for k in range(2)
            ]
        community.validate()
        runtime.founder_feeding_diagnostics = feeding.founder_feeding_diagnostics(
            community.plants.defense,
            community.insects.resistance,
            resource_vector(runtime.c.resource_capacity, runtime.c.plant_sizes, pl, "P"),
            resource_vector(runtime.c.max_feeding_demand, runtime.c.insect_sizes, hl, "H"),
            runtime.c.chemical_processing,
        )
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
        founder_feeding_diagnostics=runtime.founder_feeding_diagnostics,
        stored_genotype_values=runtime.stored_genotype_values,
        genotype_storage_semantics="cumulative retained offspring genotype values plus retained snapshot genotype values; not peak live population",
        status=runtime.status,
        resolved_B=None,
        resolved_weights=runtime.w.tolist(),
        resolved_plant_lambda=runtime.lp.tolist(),
        resolved_insect_lambda=runtime.li.tolist(),
        source_hash=source_hashes(),
        python=platform.python_version(),
        numpy=np.__version__,
        pandas=pd.__version__,
        snapshot_semantics="t=0 founder insects/plants; plants retain IDs and microsites across insect generations; round t larvae recruit insect cohort t+1; plant newborns established each round; terminal t is snapshot only",
        lifecycle_schema="dose-tolerance-lifecycle-2",
        count_semantics="plant N_next=N_retained+N_offspring; N_propagules counts pre-establishment proposals; plant offspring means established newborn only; insect N_next counts recruited adults; no terminal/failed event rows",
        resource_semantics="positive capacity K and carried B, 0<=B<=K; actual defoliation=X/B_before, zero if B_before=0; survival=exp(-mu0-muD*defoliation); surviving plants recover at base_rate+Z*tolerance_rate; maternal output=chemical-cost factor*B_recovered/K; founder B=K, newborn B=fraction*K; growth is external production, not a closed carbon budget",
    )
    if modules.enabled(runtime.c):
        runtime.result["metadata"].update(
            functional_schema="expression-tuning-contact-3",
            module_assumptions="fixed chemical features independent of concentration; fixed Gaussian widths; conserved processing capacity with idle channel; expression outside sensory saturation",
            compatibility_semantics="plant species fields track maternal founder-group origin, not current chemical phenotype or mating barriers; insect species fields are fixed ancestral compatibility pools, never dynamically inferred species",
            insect_module_order=[
                "met_E1",
                "met_E2",
                "met_mu1",
                "met_mu2",
                "sense_E1",
                "sense_E2",
                "sense_mu1",
                "sense_mu2",
                "G1_raw",
                "G2_raw",
            ],
            A_semantics="derived unsigned sensory expression profile on fixed chemical catalog, not an odor template",
        )
    runtime.result.update(runtime.genetic.records())
    architecture = json.dumps(runtime.c.genetics, sort_keys=True, separators=(",", ":"))
    payload = json.dumps(config_payload(runtime.c), sort_keys=True, separators=(",", ":"))
    runtime.result["metadata"].update(
        software_version=__version__,
        contribution_semantics="offspring_* and realized_genetic_contribution count all births; recruited_* and realized_recruitment_contribution count recruited adults; parents at round t; children adult_cohort t+1",
        fitness_proxy_semantics="F_i/fp/fi are maternal output multipliers; plant W_i=S*F uses recovered resources; larval W_insect is nutritional survival times toxicity and antinutrition survival; q is full-demand chemical assay survival, never an intake multiplier; established plant births count genetic contribution",
        maternal_expectation_semantics="plants: global surviving nonself donor required; expected_propagules and propagule_variance before shared microsite lottery; expected_maternal_offspring and maternal_offspring_variance missing because lottery couples establishment; insects: potential egg Poisson for females with adult mates, before oviposition rejection; no-mate females zero, males missing; no larval-survival or historical feeding factor in maternal egg expectation",
        allele_contribution_semantics="stage=birth or recruitment; each uses its own parental role-count weights and observed child genotypes; for plant maternal founder-group summaries, weighted parents may include fathers from other groups, so the contribution term includes mixing; recruitment transmission residual also includes larval viability selection",
        insect_lambda_semantics="Poisson lambda_egg * maternal F_m(I_m); no past maternal feeding gain; retained numeric examples are not old-model demographic equivalents",
        host_background_semantics="birth host species; founders unknown; directed parental source to child destination at half each nuclear role; not independent populations or effective migration",
        network_semantics="networks: larval resource and successful recruitment; lagged_networks: parent birth-host actual reproductive contribution; host_groups: reversible ecological categories; offspring: biparental DAG for both plants and insects, no candidate lineage/speciation",
        genetic_summary_mutation_semantics="snapshot allele/trait summaries describe recruited adults/plants; mutation_count counts all births in the cohort, including dead larvae",
        reference_panel_semantics="rows align with insect_ids; columns with reference_plant_species; p is first olfactory approach probability times contact acceptance on equal panel, not eventual acceptance after retries; q is survival at full demand without competition",
        ecology_rule_version=runtime.c.model_version,
        genetics_version=runtime.c.genetics["schema_version"],
        chemical_schema="independent-chemical-2",
        plant_trait_order=["PSM1", "PSM2", "Z", "VOC1", "VOC2"],
        plant_defense_order=["PSM1", "PSM2"],
        odor_trait_order=["VOC1", "VOC2"],
        insect_contact_trait_order=["G1_raw", "G2_raw"],
        contact_weight_mapping="G=2*raw-1; GR-mediated aggregate response, not gene expression",
        chemical_load_semantics="dose=x*PSM*dose_scale; capacity=I*processing_scale; processed=min(dose,capacity); residual=dose-processed; antinutrient=dose; per-window abstract loads",
        plant_resource_semantics="defoliation=X/B_before, zero at B_before=0; Z raises recovery rate; maternal output uses recovered resources after survival; no tolerance fixed cost",
        chemical_observation_window=runtime.c.feeding_window,
        slim_calls=runtime.genetic.engine_calls,
        genetic_architecture_sha256=hashlib.sha256(architecture.encode()).hexdigest(),
        config_sha256=hashlib.sha256(payload.encode()).hexdigest(),
        genetic_rng="NumPy founders: SeedSequence([seed,140,0,taxon,species]); SLiM seeds: SeedSequence([seed,140,4,taxon,generation]).generate_state(2)+1; first meiosis, second Eidos QTL mutation",
        rng="SeedSequence([seed,130,phase,taxon,species]); 0 init, 1 ecology, 2 survival, 3 maternal births, 5 mate choice, 6 founder sex, 7 offspring sex, 8 cue shuffle, 9 shared plant establishment lottery",
        spatial_rng="SeedSequence([seed,150,0]) plants; [seed,150,2] insects; [seed,150,1,taxon] maternal dispersal",
        inheritance="both taxa: SLiM native addCrossed, copy 0 maternal and copy 1 paternal; plants use global random nonself donors across founder groups; insects use local same-pool males; marker boundary r=(1-exp(-2*distance_Morgan))/2",
        mutation_semantics="per transmitted abstract QTL copy per newborn, independently in both parental copies; not per nucleotide or elapsed year; switches and rates remain configuration-defined",
        crossover_semantics="switch_before_locus records marker ancestry changes; crossover_count is missing because compressed intervals do not retain full crossover counts",
        spatial_semantics="periodic Chebyshev distance; adult mating at current coordinates; eggs searched from mother coordinates; recruited insects disperse from birth host once; plants disperse from mother",
    )
    return runtime.result
