"""Deterministic observations on a fixed, equally available ancestral panel."""

import numpy as np

from evo_sim.models.plant_insect.biology.insects.feeding import consumption
from evo_sim.models.plant_insect.biology.insects.sensory import conditional_probabilities
from evo_sim.models.plant_insect.biology.plants import chemistry
from evo_sim.models.plant_insect.biology.plants.chemistry import reference

from .events import host_name


def standardized(state, c, B, w, generation):
    """Return summaries and individual p/q on an equal panel, without ecological RNG."""
    D, V = reference(c, B)
    probability = np.empty((len(state.insects), len(D)))
    ids = c.observation_reference_traits["reference_ids"]
    for i, preference in enumerate(state.insects.module_traits):
        probability[i] = conditional_probabilities(V, preference, c, w)
    q = consumption(
        np.broadcast_to(D, (len(state.insects), *D.shape)),
        state.insects.resistance[:, None, :],
    )
    rows = []
    for h in range(len(c.insect_sizes)):
        mask = state.insects.species == h
        n = int(mask.sum())
        for p in range(len(D)):
            preference, capacity = probability[mask, p], q[mask, p]
            joint = preference * capacity
            mean_p = float(preference.mean()) if n else np.nan
            mean_q = float(capacity.mean()) if n else np.nan
            mean_joint = float(joint.mean()) if n else np.nan
            rows.append(
                dict(
                    generation=generation,
                    insect_species=f"H{h + 1}",
                    reference_plant_species=ids[p],
                    insect_count=n,
                    preference=mean_p,
                    q_unlimited=mean_q,
                    joint_pq=mean_joint,
                    cov_pq=mean_joint - mean_p * mean_q,
                    variance_p=float(preference.var()) if n else np.nan,
                    variance_q=float(capacity.var()) if n else np.nan,
                    variance_pq=float(joint.var()) if n else np.nan,
                    panel="fixed_chemical_reference",
                )
            )
    return rows, probability, q


def observe(runtime, current, generation):
    values = current.plants.genotypes.size + current.insects.genotypes.size
    runtime.genetic.preflight([0, 0], runtime.stored_genotype_values + values)
    snapshot = current.snapshot(generation)
    snapshot["V_emit"] = chemistry.plant_voc(
        current.plants.defense,
        current.plants.voc_trait,
        runtime.c,
        runtime.B,
        current.plants.chemical_traits,
    )
    snapshot["R"] = current.plants.resource.copy()
    runtime.result["snapshots"].append(snapshot)
    runtime.stored_genotype_values += values
    sizes = (len(runtime.c.plant_sizes), len(runtime.c.insect_sizes))
    runtime.result["genetic_summary"].extend(runtime.genetic.summary(current, generation, sizes))
    runtime.result["haplotype_summary"].extend(
        runtime.genetic.haplotype_summary(current, generation, sizes)
    )
    if runtime.observe:
        rows, probability, q = standardized(current, runtime.c, runtime.B, runtime.w, generation)
        runtime.result["diagnostics"].extend(rows)
        snapshot.update(
            reference_p=probability,
            reference_q=q,
            reference_plant_species=np.array(
                runtime.c.observation_reference_traits["reference_ids"]
            ),
        )
    for h in range(sizes[1]):
        for p in range(-1, sizes[0]):
            count = int(
                np.count_nonzero(
                    (current.insects.species == h) & (current.insects.birth_host_species == p)
                )
            )
            runtime.result["host_groups"].append(
                dict(
                    generation=generation,
                    insect_species=f"H{h + 1}",
                    birth_host_species=host_name(p),
                    N=count,
                    status="present" if count else "absent",
                )
            )
    return snapshot
