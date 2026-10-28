"""Closed-form assay checks, reciprocal effects, contribution identity and censoring."""

from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect.analysis.common import Batch, sha256
from evo_sim.models.plant_insect.analysis.summarize_coevolution import (
    contribution_differentials,
    cross_assay,
    edge_histories,
    read_design,
    summarize,
)
from evo_sim.models.plant_insect.audit import audit
from evo_sim.models.plant_insect.batch import run
from evo_sim.models.plant_insect.helpers import traits, voc
from evo_sim.models.plant_insect.helpers.config import config_payload
from tests.support import config, raises


def test_v014_reciprocal_assay_and_missing_species():
    c = config()
    c.validate()
    with TemporaryDirectory() as tmp:
        design, _ = read_design(Path(tmp))
    plants = dict(
        D=np.array([[0.8, 0.8]]),
        plant_chemical_traits=np.array([[0.8, 0.8, 0.8]]),
        plant_species=np.array([0]),
    )
    insects = dict(I=np.array([[0.2, 0.2]]), A=np.array([[0.8, 0.8]]), insect_species=np.array([0]))
    before = {k: v.copy() for k, v in plants.items()}
    row, absent = cross_assay(plants, insects, c, design)
    a = voc.choice_probabilities(np.array([c.b]), c.tau, c.outside_utility, c.choice_mode)[0]
    assert np.isclose(row["acceptance"], a) and np.isclose(row["q"], 0.25)
    assert np.isclose(row["insect_survival"], a * 0.25)
    assert np.isclose(row["plant_survival"], 1 - (1 - c.sP0) * a * 0.25)
    assert not absent["available"] and np.isnan(absent["q"])
    stronger = dict(insects, I=np.array([[0.4, 0.4]]))
    later = cross_assay(plants, stronger, c, design)[0]
    assert later["plant_survival"] < row["plant_survival"]
    weak = dict(plants, D=np.array([[0.4, 0.4]]), plant_chemical_traits=np.array([[0.4, 0.4, 0.8]]))
    assert cross_assay(weak, insects, c, design)[0]["q"] > row["q"]
    # A resource-limited, forced-acceptance trial: x=min(R,dmax*q), not R*q.
    c.choice_mode = "forced"
    limited = cross_assay(plants, insects, c, dict(design, assay_resource_capacity=0.1))[0]
    assert np.isclose(limited["resource"], 0.1) and np.isclose(limited["plant_damage"], 1)
    c.choice_mode = "outside"
    c.tau = 0
    assert np.isclose(cross_assay(plants, insects, c, design)[0]["acceptance"], 0.5)
    assert all(np.array_equal(plants[k], v) for k, v in before.items())
    # Preserve the joint distribution of A/I; mean traits are not a substitute.
    insects = dict(
        I=np.array([[0.2, 0.2], [0.8, 0.8]]),
        A=np.array([[0.0, 0.0], [0.8, 0.8]]),
        insect_species=np.array([0, 0]),
    )
    c.tau = 1
    row = cross_assay(plants, insects, c, design)[0]
    expected = []
    for i in range(2):
        score = traits.chemical_scores(np.array([[0.8, 0.8]]), insects["A"][i], c, c.weights)
        expected.append(
            voc.choice_probabilities(score, c.tau, c.outside_utility, c.choice_mode)[0]
            * traits.consumption(plants["D"], insects["I"][i])[0]
        )
    assert np.isclose(row["resource"], np.mean(expected))


def test_v014_contribution_identity_and_uncommitted_boundary():
    snap = dict(
        plant_ids=np.array([4, 5]),
        plant_species=np.array([0, 0]),
        D=np.array([[0.2], [0.8]]),
        plant_voc_independent_trait=np.array([[0.1], [0.4]]),
        insect_ids=np.array([], int),
        insect_species=np.array([], int),
        I=np.empty((0, 1)),
        A=np.empty((0, 1)),
    )
    parents = pd.DataFrame(
        dict(
            kind=["plant", "plant"],
            individual_id=[5, 4],
            contribution_committed=[True, True],
            realized_genetic_contribution=[3, 1],
        )
    )
    row = contribution_differentials(snap, snap, parents)[0]
    assert np.isclose(row["mean"], 0.5) and np.isclose(row["contribution_weighted_mean"], 0.65)
    assert np.isclose(row["realized_differential"] + row["transmission_residual"], 0)
    parents["contribution_committed"] = False
    row = contribution_differentials(snap, None, parents)[0]
    assert np.isnan(row["realized_differential"]) and np.isnan(row["next_mean"])


def test_v014_edge_window_and_design_validation():
    with TemporaryDirectory() as tmp:
        path = Path(tmp)
        design, _ = read_design(path)
        (path / "analysis.yaml").write_text("persistence_generations: 0", encoding="utf-8")
        raises(lambda: read_design(path))
    net = pd.DataFrame(
        dict(
            generation=[0, 1, 2, 4],
            insect_species=["H1"] * 4,
            plant_species=["P1"] * 4,
            Q_actual=[1.0] * 4,
            genetic_contribution=[1.0, 1.0, 1.0, np.nan],
            contribution_committed=[True, True, True, False],
        )
    )
    rows = edge_histories(net, design)
    assert rows[2]["persistent"] and rows[3]["consecutive_generations"] == 1
    assert not rows[-1]["observed"] and rows[-1]["present"] is None


def test_v014_saved_batch_analysis_and_tamper_detection():
    c = config(generations=1)
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "batch"
        run(
            dict(base=config_payload(c), cases=[dict(name="check", overrides={})], replicates=1),
            path,
        )
        audit(path)
        raw = {p: sha256(p) for p in (path / "data" / "check").rglob("*") if p.is_file()}
        batch = Batch.open(path)
        result = summarize(batch)
        assert len(result["cross_assay"]) == 8  # 2 periods x 2 periods x 2 plant species.
        frame = pd.DataFrame(result["contribution_differentials"])
        terminal = frame[frame.generation == 1]
        assert not terminal.contribution_committed.any()
        assert terminal.realized_differential.isna().all()
        assert {p: sha256(p) for p in raw} == raw
        saved = pd.read_csv(batch.summary_dir / "host_relations.csv")
        assert saved.plant_abundance.notna().all() and saved.insect_abundance.notna().all()
        assert (batch.summary_dir / "coevolution_metadata.json").is_file()
        events = path / "data/check/rep01/events.csv"
        events.write_bytes(events.read_bytes() + b"\n")
        raises(lambda: summarize(batch))
