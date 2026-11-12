"""Exact stage accounting of allele-frequency change from recorded individuals.

No biological outcomes are recalculated. Components describe realized changes,
including demographic sampling; they are not estimates of causal natural selection.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect.analysis.common import sha256

FREQUENCIES = (
    "p_parent p_source_mixed p_parent_weighted p_gamete p_birth p_recruit p_retained p_next"
).split()
COMPONENTS = (
    "delta_source_mixing delta_parent_reweighting delta_meiosis delta_mutation "
    "delta_viability delta_retention"
).split()
COLUMNS = (
    (
        "generation offspring_generation kind scope species locus_id allele_id "
        "N_parent N_birth N_recruit N_retained N_next"
    ).split()
    + FREQUENCIES
    + COMPONENTS
    + ["delta_total", "closure_error"]
)


def decompose_transition(
    before: dict,
    after: dict,
    offspring: pd.DataFrame,
    transmissions: pd.DataFrame,
    loci: dict[str, list[str]],
    species_counts: dict[str, int],
) -> pd.DataFrame:
    """Account for one committed transition, in dimensionless allele frequencies.

    The six deltas sum to p_next-p_parent, weighted by the fractions of recruits
    and retained individuals. Conditional cohort frequencies are NA when empty;
    extinct next populations have undefined changes, never a fabricated zero.
    Source mixing uses realized parental-role source shares and source population
    frequencies before reproduction. Reweighting includes reproductive sampling.
    """
    generation = int(before["generation"])
    if int(after["generation"]) != generation + 1:
        raise ValueError("adjacent committed snapshots required")
    rows = []
    for kind, names in loci.items():
        ids = np.asarray(before[kind + "_ids"])
        next_ids = np.asarray(after[kind + "_ids"])
        labels = np.asarray(before[kind + "_species"])
        next_labels = np.asarray(after[kind + "_species"])
        parent = np.asarray(before[kind + "_genotypes"])
        next_g = np.asarray(after[kind + "_genotypes"])
        born = offspring[
            offspring.generation.eq(generation + 1) & offspring.kind.eq(kind)
        ].sort_values("offspring_id")
        if born.offspring_id.duplicated().any() or len(np.unique(ids)) != len(ids):
            raise ValueError("unique parent and offspring identities required")
        birth_ids = born.offspring_id.to_numpy(dtype=np.int64)
        birth_g = np.asarray([json.loads(g) for g in born.genotype], dtype=np.int64).reshape(
            len(born), 2, len(names)
        )
        lookup = {int(id_): i for i, id_ in enumerate(ids)}
        roles = np.asarray(
            [[lookup[int(m)], lookup[int(f)]] for m, f in zip(born.maternal_id, born.paternal_id)],
            dtype=int,
        ).reshape(len(born), 2)
        birth_labels = labels[roles[:, 0]]
        retained = np.isin(ids, next_ids)
        recruited = np.isin(birth_ids, next_ids)
        if (
            (kind == "insect" and retained.any())
            or len(next_ids) != int(retained.sum() + recruited.sum())
            or np.intersect1d(ids, birth_ids).size
            or len(np.unique(next_ids)) != len(next_ids)
        ):
            raise ValueError("next population must contain exactly retained parents and recruits")
        next_lookup = {int(id_): i for i, id_ in enumerate(next_ids)}
        for source_ids, source_g, source_labels, mask in (
            (ids, parent, labels, retained),
            (birth_ids, birth_g, birth_labels, recruited),
        ):
            positions = np.array([next_lookup[int(id_)] for id_ in source_ids[mask]], dtype=int)
            if not np.array_equal(next_g[positions], source_g[mask]) or not np.array_equal(
                next_labels[positions], source_labels[mask]
            ):
                raise ValueError(
                    "snapshot identity, genotype or source label does not match records"
                )
        tx = transmissions[
            transmissions.generation.eq(generation + 1) & transmissions.kind.eq(kind)
        ].set_index(["offspring_id", "copy", "locus_id"])
        index = pd.MultiIndex.from_product([birth_ids, [0, 1], names], names=tx.index.names)
        if not tx.index.is_unique or len(tx) != len(index):
            raise ValueError("one transmission per offspring, parental copy and locus required")
        tx = tx.reindex(index)
        if tx[["inherited_allele_id", "child_allele_id", "parent_id"]].isna().any().any():
            raise ValueError("missing transmission keys")
        inherited = tx.inherited_allele_id.to_numpy(dtype=np.int64).reshape(birth_g.shape)
        actual = tx.child_allele_id.to_numpy(dtype=np.int64).reshape(birth_g.shape)
        expected_parent_ids = np.broadcast_to(ids[roles, None], (len(born), 2, len(names)))
        if not np.array_equal(actual, birth_g) or not np.array_equal(
            tx.parent_id.to_numpy().reshape(birth_g.shape), expected_parent_ids
        ):
            raise ValueError("transmissions do not match recorded children or parents")
        # A source is a maternal-origin label for plants, a compatibility pool for insects.
        scopes = [("all", "all", None)] + [
            ("type", f"{'P' if kind == 'plant' else 'H'}{s + 1}", s)
            for s in range(species_counts[kind])
        ]
        for scope, species, label in scopes:
            pm = np.ones(len(ids), bool) if label is None else labels == label
            bm = np.ones(len(born), bool) if label is None else birth_labels == label
            nm = np.ones(len(next_ids), bool) if label is None else next_labels == label
            rm, cm = pm & retained, bm & recruited
            n0, nb, nc, nr, nn = map(lambda mask: int(mask.sum()), (pm, bm, cm, rm, nm))
            role_indices = roles[bm].reshape(-1)
            for l, name in enumerate(names):
                alleles = np.unique(
                    np.concatenate(
                        (
                            parent[:, :, l].ravel(),
                            inherited[:, :, l].ravel(),
                            birth_g[:, :, l].ravel(),
                            next_g[:, :, l].ravel(),
                        )
                    )
                )
                for allele in alleles:
                    z = (parent[:, :, l] == allele).mean(axis=1)
                    p0 = float(z[pm].mean()) if n0 else np.nan
                    weighted = float(z[role_indices].mean()) if nb else np.nan
                    mixed = p0 if label is None else np.nan
                    if label is not None and nb:
                        mixed = sum(
                            np.count_nonzero(labels[role_indices] == source)
                            * float(z[labels == source].mean())
                            for source in np.unique(labels[role_indices])
                        ) / (2 * nb)
                    gamete = float((inherited[bm, :, l] == allele).mean()) if nb else np.nan
                    birth = float((birth_g[bm, :, l] == allele).mean()) if nb else np.nan
                    recruit = float((birth_g[cm, :, l] == allele).mean()) if nc else np.nan
                    keep = float(z[rm].mean()) if nr else np.nan
                    nxt = float((next_g[nm, :, l] == allele).mean()) if nn else np.nan
                    deltas = np.full(6, np.nan)
                    if n0 and nn:
                        deltas = np.zeros(6)
                        if nc:
                            deltas[:5] = (
                                nc / nn * np.diff([p0, mixed, weighted, gamete, birth, recruit])
                            )
                        if nr:
                            deltas[5] = nr / nn * (keep - p0)
                        if not np.isclose(deltas.sum(), nxt - p0, rtol=0, atol=1e-12):
                            raise ValueError("allele change decomposition does not close")
                    rows.append(
                        dict(
                            generation=generation,
                            offspring_generation=generation + 1,
                            kind=kind,
                            scope=scope,
                            species=species,
                            locus_id=name,
                            allele_id=int(allele),
                            N_parent=n0,
                            N_birth=nb,
                            N_recruit=nc,
                            N_retained=nr,
                            N_next=nn,
                            **dict(
                                zip(
                                    FREQUENCIES,
                                    [p0, mixed, weighted, gamete, birth, recruit, keep, nxt],
                                )
                            ),
                            **dict(zip(COMPONENTS, deltas)),
                            delta_total=nxt - p0,
                            closure_error=float(deltas.sum() - (nxt - p0)),
                        )
                    )
    return pd.DataFrame(rows, columns=COLUMNS)


def decompose_run(path: Path) -> tuple[pd.DataFrame, dict[str, str]]:
    """Read a previously audited run; verify input hashes and preserve all raw files."""
    path = Path(path)
    meta = json.loads((path / "metadata.json").read_text(encoding="utf-8"))
    if meta["version"] != "0.1.12":
        raise ValueError("current recorded QTL protocol required")
    files = ["offspring.csv", "gamete_transmissions.csv", "genetic_loci.csv"]
    snapshots = sorted((path / "snapshots").glob("g*.npz"))
    snapshot_names = [file.relative_to(path).as_posix() for file in snapshots]
    recorded_names = sorted(
        name
        for name in meta["output_sha256"]
        if name.startswith("snapshots/g") and name.endswith(".npz")
    )
    if snapshot_names != recorded_names:
        raise ValueError("snapshot inventory differs from committed records")
    files += snapshot_names
    hashes = {name: sha256(path / name) for name in files}
    if any(meta["output_sha256"].get(name) != digest for name, digest in hashes.items()):
        raise ValueError("raw analysis inputs do not match recorded hashes")
    hashes["metadata.json"] = sha256(path / "metadata.json")
    offspring = pd.read_csv(path / "offspring.csv")
    tx = pd.read_csv(path / "gamete_transmissions.csv")
    catalog = pd.read_csv(path / "genetic_loci.csv")
    loci = {
        kind: catalog.loc[catalog.organism.eq(kind), "locus_id"].tolist()
        for kind in ("plant", "insect")
    }
    sizes = {kind: len(meta["config"][kind + "_sizes"]) for kind in loci}
    tables = []
    if snapshots:
        with np.load(snapshots[0], allow_pickle=False) as state:
            before = dict(state)
        for file in snapshots[1:]:
            with np.load(file, allow_pickle=False) as state:
                after = dict(state)
            tables.append(decompose_transition(before, after, offspring, tx, loci, sizes))
            before = after
    result = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame(columns=COLUMNS)
    return result, hashes
