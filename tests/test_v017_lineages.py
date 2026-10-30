"""Saved-data counterexamples for operational lineage tracking, without SLiM runs."""

import json
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect.analysis.common import sha256
from evo_sim.models.plant_insect.analysis.summarize_lineages import (
    UNASSIGNED,
    cluster_genotypes,
    read_design,
    summarize,
)

DESIGN = dict(distance_threshold=0.0, minimum_cluster_size=2, persistence_generations=2)


def _genotypes(pairs):
    return np.asarray(pairs, dtype=np.int64).reshape(-1, 2, 1)


def _fixture(root):
    folder = root / "data/case/rep01"
    (folder / "snapshots").mkdir(parents=True)
    (root / "audits").mkdir()
    (root / "settings.json").write_text('{"version": "0.1.7"}', encoding="utf-8")
    (root / "audits/integrity_audit.json").write_text('{"runs": 1}', encoding="utf-8")
    pd.DataFrame([dict(case="case", replicate=1, generation=4, status="completed")]).to_csv(
        root / "data/run_status.csv", index=False
    )
    insects = [
        ([0, 1, 2, 3], [(0, 1)] * 4),
        ([4, 5, 6, 7, 8], [(0, 0)] * 2 + [(1, 1)] * 2 + [(2, 2)]),
        ([10, 11, 12, 13], [(0, 1)] * 4),
        ([14, 15], [(0, 1)] * 2),
        ([], []),
    ]
    offspring, events, parents = [], [], []
    for t, (ids, pairs) in enumerate(insects):
        plants = [2 * t, 2 * t + 1]
        np.savez_compressed(
            folder / f"snapshots/g{t:06d}.npz",
            generation=t,
            plant_ids=np.asarray(plants),
            plant_species=np.zeros(2, dtype=int),
            plant_genotypes=_genotypes([(3, 3)] * 2),
            insect_ids=np.asarray(ids, dtype=int),
            insect_species=np.zeros(len(ids), dtype=int),
            insect_genotypes=_genotypes(pairs),
        )
        if t == 4:
            continue
        for kind, current_ids in (("plant", plants), ("insect", ids)):
            for ident in current_ids:
                candidates = [i for i in current_ids if i != ident]
                parents.append(
                    dict(
                        generation=t,
                        kind=kind,
                        individual_id=ident,
                        eligible_as_mother=True,
                        candidate_mate_ids=json.dumps(candidates),
                        candidate_mate_probabilities=json.dumps(
                            [1 / len(candidates)] * len(candidates)
                        ),
                    )
                )
        for ident in [2 * (t + 1), 2 * (t + 1) + 1]:
            offspring.append(
                dict(
                    generation=t + 1,
                    reproduction_round=t,
                    adult_cohort=t + 1,
                    kind="plant",
                    offspring_id=ident,
                    maternal_id=plants[0],
                    paternal_id=plants[1],
                    survived_to_adult=True,
                    birth_host_id=np.nan,
                )
            )
        next_ids = insects[t + 1][0]
        born = next_ids + ([9] if t == 0 else [])
        for i, ident in enumerate(born):
            if t == 1:
                mother, father = (4, 6) if i % 2 == 0 else (5, 7)
            else:
                mother, father = ids[:2]
            host = plants[i % 2]
            offspring.append(
                dict(
                    generation=t + 1,
                    reproduction_round=t,
                    adult_cohort=t + 1,
                    kind="insect",
                    offspring_id=ident,
                    maternal_id=mother,
                    paternal_id=father,
                    survived_to_adult=ident != 9,
                    birth_host_id=host,
                )
            )
            events.append(
                dict(
                    generation=t,
                    insect_id=ident,
                    host_id=host,
                    x=2.5 if ident == 9 else 0.25,
                )
            )
    for name, rows in (("offspring", offspring), ("events", events), ("parents", parents)):
        pd.DataFrame(rows).to_csv(folder / f"{name}.csv", index=False)
    _rehash(folder)
    return folder


def _rehash(folder):
    hashes = {
        f.relative_to(folder).as_posix(): sha256(f)
        for f in folder.rglob("*")
        if f.is_file() and f.name != "metadata.json"
    }
    (folder / "metadata.json").write_text(
        json.dumps(
            dict(
                version="0.1.7",
                genetics_version="slim-qtl-3",
                output_sha256=hashes,
            )
        ),
        encoding="utf-8",
    )


def test_v017_unphased_matching_and_single_link_bridge():
    design = dict(DESIGN, distance_threshold=0.5)
    # AA -- AB -- BB is a single operational component, despite AA/BB distance 1.
    mapping, groups, _ = cluster_genotypes(
        [7, 3, 9],
        [0, 0, 0],
        _genotypes([(0, 0), (0, 1), (1, 1)]),
        kind="insect",
        generation=0,
        design=design,
    )
    assert len(groups) == 1 and len(set(mapping.values())) == 1
    # Copy order cannot alter distance. AB and BC share one copy, even though sorted slots differ.
    _, groups, _ = cluster_genotypes(
        [0, 1],
        [0, 0],
        _genotypes([(1, 2), (2, 3)]),
        kind="plant",
        generation=0,
        design=design,
    )
    assert len(groups) == 1
    _, groups, _ = cluster_genotypes(
        [0, 1],
        [0, 0],
        _genotypes([(1, 2), (2, 1)]),
        kind="plant",
        generation=0,
        design=DESIGN,
    )
    assert len(groups) == 1
    # Initial compatibility pools remain a biological boundary; host labels are never read.
    mapping, groups, _ = cluster_genotypes(
        [0, 1],
        [0, 1],
        _genotypes([(0, 0), (0, 0)]),
        kind="plant",
        generation=0,
        design=dict(DESIGN, distance_threshold=1),
    )
    assert not groups and set(mapping.values()) == {UNASSIGNED}


def test_v017_saved_split_merge_extinction_and_network_balance():
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        folder = _fixture(root)
        # Current no-flowering protocol; other fixtures retain historical draft protocol 3.
        metadata_path = folder / "metadata.json"
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        metadata["genetics_version"] = "slim-qtl-4"
        metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
        originals = {f: sha256(f) for f in folder.rglob("*") if f.is_file()}
        frames = summarize(root, DESIGN)
        # The documented CLI accepts a relative batch path too.
        original_directory = Path.cwd()
        try:
            os.chdir(root.parent)
            relative = summarize(Path(root.name), DESIGN)
        finally:
            os.chdir(original_directory)
        pd.testing.assert_frame_equal(frames["networks"], relative["networks"])
        insects = frames["clusters"].query("kind == 'insect'")
        assert insects.groupby("generation").size().to_dict() == {0: 1, 1: 2, 2: 1, 3: 1}
        events = frames["lineage_events"].query("kind == 'insect'")
        assert len(events.query("generation == 1 and event == 'split'")) == 2
        assert len(events.query("generation == 2 and event == 'merge'")) == 1
        assert events.query("generation == 4 and event == 'termination'").reason.tolist() == [
            "pool_extinct"
        ]
        assert not insects.query("generation <= 2").persistent_candidate.any()
        assert insects.query("generation == 3").persistent_candidate.all()
        small = frames["memberships"].query("kind == 'insect' and individual_id == 8")
        assert small.cluster_id.tolist() == [UNASSIGNED]
        assert not small.persistent_candidate.any()
        links = frames["transitions"].query("kind == 'insect' and target_generation == 2")
        assert len(links) == 2 and links.recruited_genetic_contribution.sum() == 4
        balance = frames["network_balances"].set_index("reproduction_round")
        assert balance.loc[0, "birth_count"] == 6 and balance.loc[0, "recruitment_count"] == 5
        assert balance.loc[0, "actual_consumption"] == 3.75  # Includes the dead larva.
        assert balance.loc[3, "birth_count"] == balance.loc[3, "actual_consumption"] == 0
        pd.testing.assert_series_equal(
            balance.actual_consumption, balance.maternal_network_consumption, check_names=False
        )
        networks = frames["networks"]
        next_group = networks.query("edge_semantics == 'recruit_cluster'")
        assert (next_group.insect_generation == next_group.plant_generation + 1).all()
        assert next_group.actual_consumption.isna().all()
        actual = frames["mating"].query("kind == 'insect' and reproduction_round == 1")
        assert actual.cross_fraction_of_classified.tolist() == [1.0, 1.0]
        opportunity = frames["mating_opportunities"].query(
            "kind == 'insect' and reproduction_round == 1"
        )
        assert np.isclose(opportunity.probability_weighted_cross_fraction.iloc[0], 2 / 3)
        assert (
            frames["mating"]
            .query("kind == 'insect' and reproduction_round == 3")
            .cross_fraction_of_classified.isna()
            .all()
        )
        assert originals == {f: sha256(f) for f in originals}
        # Threshold changes only observation products, never original saved data.
        merged = summarize(root, dict(DESIGN, distance_threshold=1))
        assert len(merged["clusters"].query("kind == 'insect' and generation == 1")) == 1
        assert originals == {f: sha256(f) for f in originals}


def test_v017_label_permutation_missing_detection_and_input_guards():
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        folder = _fixture(root)
        before = summarize(root, DESIGN)
        # Only ancestor-pool names change; graph and genealogy identities do not.
        for file in (folder / "snapshots").glob("*.npz"):
            with np.load(file) as state:
                values = dict(state)
            values["plant_species"][:] = 6
            values["insect_species"][:] = 4
            np.savez_compressed(file, **values)
        _rehash(folder)
        after = summarize(root, DESIGN)
        for name in ("transitions", "lineage_events", "mating", "networks"):
            pd.testing.assert_frame_equal(before[name], after[name])
        # An undetected small group does not imply true compatibility-pool extinction.
        sparse = summarize(root, dict(DESIGN, minimum_cluster_size=4))
        extinct = sparse["lineage_events"].query("kind == 'insect' and event == 'termination'")
        assert "observation_not_detected" in set(extinct.reason)
        assert not sparse["populations"].query("generation == 1").extinction_observed.any()
        assert (
            sparse["mating"]
            .query("kind == 'insect' and reproduction_round == 1")
            .cross_fraction_of_classified.isna()
            .all()
        )
        # Missing ancestry remains explicit and conserved, rather than becoming a zero-flow group.
        offspring = pd.read_csv(folder / "offspring.csv")
        first_insect = offspring.index[offspring.kind == "insect"][0]
        offspring.loc[first_insect, "maternal_id"] = np.nan
        offspring.to_csv(folder / "offspring.csv", index=False)
        _rehash(folder)
        unknown = summarize(root, DESIGN)
        assert "unknown" in set(unknown["networks"].insect_cluster)
        assert unknown["mating"].unknown_parent_pairs.sum() > 0
        # Audit hashes must not allow silently changed input ledgers.
        with (folder / "events.csv").open("a", encoding="utf-8") as file:
            file.write("\n")
        try:
            summarize(root, DESIGN)
        except ValueError as error:
            assert "hash mismatch" in str(error)
        else:
            raise AssertionError("changed raw input accepted")
    for invalid in (
        dict(DESIGN, distance_threshold=float("nan")),
        dict(DESIGN, distance_threshold=2),
        dict(DESIGN, minimum_cluster_size=1),
        dict(DESIGN, persistence_generations=True),
    ):
        try:
            read_design(invalid)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid observation design accepted")


if __name__ == "__main__":
    test_v017_unphased_matching_and_single_link_bridge()
    test_v017_saved_split_merge_extinction_and_network_balance()
    test_v017_label_permutation_missing_detection_and_input_guards()
    print("PASS v0.1.7 operational lineages")
