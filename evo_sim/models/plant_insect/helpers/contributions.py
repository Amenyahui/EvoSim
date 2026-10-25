"""Deterministic ledgers for committed generation transitions; no random draws."""

import numpy as np

ALLELE_COLUMNS = (
    "generation offspring_generation kind species locus_id allele_id N_parent N_offspring "
    "p_parent p_parent_weighted p_offspring delta_parent_contribution delta_transmission"
)
NETWORK_CONTRIBUTION_COLUMNS = (
    "maternal_offspring paternal_offspring genetic_contribution maternal_share "
    "paternal_share genetic_contribution_share contribution_committed"
)


def mark_uncommitted(parents):
    """Sampling is not commitment, including when sampling itself fails."""
    for row in parents:
        row.setdefault("realized_genetic_contribution", None)
        row.setdefault("contribution_committed", False)


def commit_parent_contributions(parent_groups):
    for rows in parent_groups:
        for row in rows:
            row["realized_genetic_contribution"] = (
                row["offspring_count_as_mother"] + row["offspring_count_sired"]
            ) / 2
            row["contribution_committed"] = True


def allele_contributions(current, offspring, parent_groups, genetic, generation, sizes):
    """Parent/child allele union per species/locus; extinct offspring are missing."""
    rows = []
    for k, kind in enumerate(("plant", "insect")):
        before = current.plants if k == 0 else current.insects
        after = offspring.plants if k == 0 else offspring.insects
        lookup = {row["individual_id"]: row for row in parent_groups[k]}
        weights = np.array(
            [lookup[int(i)]["realized_genetic_contribution"] for i in before.ids], float
        )
        for s in range(sizes[k]):
            parent_mask, child_mask = before.species == s, after.species == s
            n, next_n = int(parent_mask.sum()), int(child_mask.sum())
            w = weights[parent_mask]
            if not np.isclose(w.sum(), next_n):
                raise ValueError("biparental contribution does not equal offspring population")
            for l, locus in enumerate(genetic.loci[kind]):
                parental = before.genotypes[parent_mask, :, l]
                child = after.genotypes[child_mask, :, l]
                for allele_id in np.union1d(parental, child):
                    z = (parental == allele_id).sum(axis=1) / 2
                    p = float(z.mean()) if n else np.nan
                    weighted = float(w @ z / next_n) if next_n else np.nan
                    actual = float(np.mean(child == allele_id)) if next_n else np.nan
                    rows.append(
                        dict(
                            generation=generation,
                            offspring_generation=generation + 1,
                            kind=kind,
                            species=f"{'P' if k == 0 else 'H'}{s + 1}",
                            locus_id=locus["locus_id"],
                            allele_id=int(allele_id),
                            N_parent=n,
                            N_offspring=next_n,
                            p_parent=p,
                            p_parent_weighted=weighted,
                            p_offspring=actual,
                            delta_parent_contribution=weighted - p,
                            delta_transmission=actual - weighted,
                        )
                    )
    return rows


def host_contributions(network_rows, community, host, insect_parents):
    """Each parent contributes on its own feeding host, never the child's location."""
    host_species = np.full(len(host), -1, dtype=int)
    found = host >= 0
    host_species[found] = community.plants.species[host[found]]
    lookup = {row["individual_id"]: row for row in insect_parents}
    maternal = np.array(
        [lookup[int(i)]["offspring_count_as_mother"] for i in community.insects.ids], dtype=np.int64
    )
    paternal = np.array(
        [lookup[int(i)]["offspring_count_sired"] for i in community.insects.ids], dtype=np.int64
    )
    if np.any((maternal + paternal)[~found]):
        raise ValueError("unmatched insect cannot have committed offspring")
    result = []
    for original in network_rows:
        row = dict(original)
        h, p = int(row["insect_species"][1:]) - 1, int(row["plant_species"][1:]) - 1
        species = community.insects.species == h
        matched = species & (host_species == p)
        bm, bf = int(maternal[matched].sum()), int(paternal[matched].sum())
        tm, tf = int(maternal[species].sum()), int(paternal[species].sum())
        row.update(
            maternal_offspring=bm,
            paternal_offspring=bf,
            genetic_contribution=(bm + bf) / 2,
            maternal_share=bm / tm if tm else np.nan,
            paternal_share=bf / tf if tf else np.nan,
            genetic_contribution_share=(bm + bf) / (tm + tf) if tm + tf else np.nan,
            contribution_committed=True,
        )
        result.append(row)
    return result
