"""Independent reconstruction of saved variable-population ledgers (T16)."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from evo_sim.core.helpers import source_hashes

from .batch import write_json
from .helpers.config import Config, config_payload

ROOT = Path(__file__).resolve().parents[3]

if not __debug__:
    raise RuntimeError("Run scientific audits without python -O; assertions are required")


def report_directory(batch):
    return ROOT / "tests" / "reports" / Path(batch).resolve().name


def close(a, b):
    assert np.allclose(
        np.asarray(a, float), np.asarray(b, float), rtol=1e-11, atol=1e-14, equal_nan=True
    )


def flags(series):
    values = series.map({True: True, False: False, "True": True, "False": False})
    assert values.notna().all()
    return values.to_numpy(bool)


def per_species(value, keys):
    if isinstance(value, dict):
        return np.array([value[k] for k in keys], float)
    return np.broadcast_to(value, (len(keys),))


def audit_run(path):
    meta = json.loads((path / "metadata.json").read_text(encoding="utf-8"))
    status = meta["status"]["status"]
    cfg = meta["config"]
    assert json.loads((path / "status.json").read_text(encoding="utf-8")) == meta["status"]
    for name, digest in meta["output_sha256"].items():
        assert hashlib.sha256((path / name).read_bytes()).hexdigest() == digest, name
    if status == "invalid_parameter":
        try:
            Config(**cfg).validate()
        except (ValueError, TypeError, KeyError):
            return dict(
                status=status, states=0, transitions=0, max_overrun=0.0, max_balance_error=0.0
            )
        raise AssertionError("invalid_parameter record has valid configuration")
    c = Config(**cfg)
    c.validate()
    genetic_checks = {}
    if c.genetics is not None and c.genetics.get("mode") == "qtl_additive":
        from .audit_genetics import audit_genetics

        genetic_checks = audit_genetics(path, c, meta)
    tables = {
        name: pd.read_csv(path / f"{name}.csv", float_precision="round_trip", low_memory=False)
        for name in ("events", "plants", "parents", "offspring", "populations", "networks")
    }
    states = []
    for file in sorted((path / "snapshots").glob("*.npz")):
        with np.load(file) as a:
            states.append({k: a[k] for k in a.files})
    if not states:
        from .audit_genetics import genetic_limit

        assert status == "computational_limit" and (
            sum(c.plant_sizes) + sum(c.insect_sizes) > c.max_individuals_safety
            or genetic_limit(c, states, None)
        )
        return dict(
            status=status,
            states=0,
            transitions=0,
            max_overrun=0.0,
            max_balance_error=0.0,
            **genetic_checks,
        )
    ps = [f"P{s + 1}" for s in range(len(c.plant_sizes))]
    hs = [f"H{s + 1}" for s in range(len(c.insect_sizes))]
    alpha = np.array([[c.plant_competition_alpha[s][r] for r in ps] for s in ps])
    lp = np.array([c.plant_offspring_lambda[s] for s in ps])
    li = np.array([c.insect_offspring_lambda[s] for s in hs])
    transitions = 0
    overrun = 0.0
    balance = 0.0
    first_extinct = {}
    seen_ids = [set(), set()]
    assert [int(s["generation"]) for s in states] == list(range(len(states)))
    if c.model_version == "1.6":
        first = states[0]
        for kind, name, stream in (
            ("plant", "initial_plant_cells", 0),
            ("insect", "initial_insect_cells", 2),
        ):
            given = c.spatial[name]
            expected = (
                np.random.default_rng(np.random.SeedSequence([c.seed, 150, stream])).integers(
                    0, c.spatial["grid_shape"], size=(len(first[kind + "_ids"]), 2)
                )
                if given is None
                else np.asarray(given, int)
            )
            assert np.array_equal(first[kind + "_cells"], expected)
    for t, state in enumerate(states):
        pl, hl = state["plant_species"], state["insect_species"]
        D, I, A, Vtrait = (state[k] for k in ("D", "I", "A", "plant_voc_independent_trait"))
        for kind, keys in enumerate((("D", "plant_voc_independent_trait"), ("I", "A"))):
            labels = pl if kind == 0 else hl
            ids = state["plant_ids"] if kind == 0 else state["insect_ids"]
            assert ids.shape == labels.shape and len(np.unique(ids)) == len(ids)
            assert not seen_ids[kind].intersection(map(int, ids))
            seen_ids[kind].update(map(int, ids))
            for key in keys:
                x = state[key]
                assert x.shape[0] == len(ids) and x.ndim == 2
                assert np.isfinite(x).all() and ((x >= 0) & (x <= 1)).all()
            if c.spatial is not None:
                cells = state["plant_cells" if kind == 0 else "insect_cells"]
                assert cells.shape == (len(ids), 2) and cells.dtype.kind in "iu"
                assert ((cells >= 0) & (cells < np.asarray(c.spatial["grid_shape"]))).all()
        if "R" not in state or not (tables["populations"].generation == t).any():
            assert t == len(states) - 1 and status in ("arithmetic_error", "computational_limit")
            continue
        # Explicit masks keep audit selection independent of expression scope.
        ev = tables["events"][tables["events"].generation == t].sort_values("insect_row")
        p = tables["plants"][tables["plants"].generation == t].sort_values("plant_row")
        parent = tables["parents"][tables["parents"].generation == t]
        pop = tables["populations"][tables["populations"].generation == t]
        net = tables["networks"][tables["networks"].generation == t]
        assert len(ev) == len(hl) and len(p) == len(pl)
        assert (
            ev.insect_id.tolist() == state["insect_ids"].tolist()
            and p.plant_id.tolist() == state["plant_ids"].tolist()
        )
        assert ev.species.tolist() == [hs[s] for s in hl] and p.species.tolist() == [
            ps[s] for s in pl
        ]
        R = per_species(c.resource_capacity, ps)[pl]
        dmax = per_species(c.max_feeding_demand, hs)[hl]
        close(state["R"], R)
        close(state["dmax"], dmax)
        close(p.R, R)
        close(ev.dmax, dmax)
        B = np.array(meta["resolved_B"])
        V = np.clip(c.rho_voc * (D @ B.T) + (1 - c.rho_voc) * Vtrait, 0, 1)
        close(V, state["V_emit"])
        assert sorted(map(tuple, state["V_cue"])) == sorted(map(tuple, V))
        if not c.shuffle_voc:
            assert np.array_equal(V, state["V_cue"])
        host = ev.host_row.to_numpy(int)
        accepted = flags(ev.accepted)
        assert np.array_equal(accepted, host >= 0) and ((host >= -1) & (host < len(pl))).all()
        assert ev.loc[~accepted, "q"].isna().all()
        assert ev.loc[~accepted, "host_id"].eq(-1).all()
        assert ev.loc[accepted, "host_id"].tolist() == state["plant_ids"][host[accepted]].tolist()
        q = np.full(len(hl), np.nan)
        if accepted.any():
            ratios = np.divide(
                I[accepted],
                D[host[accepted]],
                out=np.full_like(I[accepted], np.inf),
                where=D[host[accepted]] > 0,
            )
            q[accepted] = np.minimum(1, ratios.min(axis=1))
        close(ev.q, q)
        d = np.zeros(len(hl))
        d[accepted] = dmax[accepted] * q[accepted]
        T = np.bincount(host[accepted], weights=d[accepted], minlength=len(pl))
        gamma = np.ones(len(pl))
        limited = T > R
        gamma[limited] = R[limited] / T[limited]
        x = np.zeros(len(hl))
        x[accepted] = d[accepted] * gamma[host[accepted]]
        X = np.bincount(host[accepted], weights=ev.x.to_numpy()[accepted], minlength=len(pl))
        close(ev.d, d)
        close(ev.x, x)
        close(ev.g, ev.x.to_numpy() / dmax)
        close(p["T"], T)
        close(p.X, X)
        close(p.gamma, gamma)
        close(p.h, X / R)
        close(p.herbivore_load, np.bincount(host[accepted], minlength=len(pl)))
        assert np.all(X <= R) and (ev.g.between(0, 1)).all() and (p.h.between(0, 1)).all()
        close(X, np.minimum(T, R))
        if len(pl):
            overrun = max(overrun, float(np.maximum(0, X - R).max()))
            balance = max(balance, float((np.abs(X - np.minimum(T, R)) / R).max()))
        close(p.defense_cost, c.plant_cost * D.sum(axis=1))
        close(p.voc_cost, c.voc_cost * V.mean(axis=1))
        Sp = c.sP0 + (1 - c.sP0) * (1 - X / R)
        Fp = np.maximum(0, 1 - c.plant_cost * D.sum(axis=1) - c.voc_cost * V.mean(axis=1))
        Si = ev.g.to_numpy()
        Fi = np.maximum(0, 1 - c.insect_cost * I.sum(axis=1))
        close(p.W_plant, Sp * Fp)
        close(ev.W_insect, Si * Fi)
        npcount = np.bincount(pl, minlength=len(ps))
        hicount = np.bincount(hl, minlength=len(hs))
        phi = 1 / (1 + alpha @ npcount)
        individual_phi = phi[pl]
        if c.spatial is not None:
            individual_phi = spatial_phi(state, alpha, c)
            phi = np.array(
                [
                    individual_phi[pl == s].mean() if (pl == s).any() else np.nan
                    for s in range(len(ps))
                ]
            )
            if len(p):
                close(p.density_phi, individual_phi)
            audit_space_events(state, ev, c, meta)
        assert len(pop) == len(ps) + len(hs) and len(net) == len(ps) * len(hs)
        assert net.N_accept.sum() == int(accepted.sum())
        close(net.Q_actual.sum(), ev.x.sum())
        for h, name in enumerate(hs):
            for k, pname in enumerate(ps):
                pair = ev[(ev.species == name) & (ev.plant_species == pname)]
                edge = net[(net.insect_species == name) & (net.plant_species == pname)].iloc[0]
                assert edge.N_accept == len(pair)
                close(edge.Q_actual, pair.x.sum())
                close(edge.L_per_insect, pair.x.sum() / hicount[h] if hicount[h] else np.nan)
        # A runtime numerical failure may leave only a prefix of parent draws; never fabricate it as a transition.
        if len(parent) != len(pl) + len(hl):
            assert status in ("arithmetic_error", "computational_limit") and t == len(states) - 1
            assert not flags(pop.update_applied).any()
            continue
        for k, (kind, labels, S, F, lam, factor, species_keys) in enumerate(
            (
                ("plant", pl, Sp, Fp, lp[pl], individual_phi, ps),
                ("insect", hl, Si, Fi, li[hl], np.ones(len(hl)), hs),
            )
        ):
            ids = state["plant_ids"] if k == 0 else state["insect_ids"]
            par = parent[parent.kind == kind].set_index("individual_id").reindex(ids)
            assert par.species.tolist() == [species_keys[s] for s in labels]
            close(par.S_i, S)
            close(par.F_i, F)
            close(par.W_i, S * F)
            close(par.phi_i, factor)
            close(par.lambda_s, lam)
            mu = lam * F * factor
            close(par.conditional_offspring_mean_mu, mu)
            close(par.expected_offspring, S * mu)
            close(par.offspring_variance, S * mu + S * (1 - S) * mu**2)
            births = np.zeros(len(labels), dtype=np.int64)
            survived = np.zeros(len(labels), bool)
            sampling = flags(par.sampling_applied)
            if sampling.any():
                assert sampling.all()
                survived = flags(par.survived_to_reproduce)
                births = par.offspring_count.to_numpy(float)
                assert (
                    np.isfinite(births).all()
                    and (births >= 0).all()
                    and (births == np.floor(births)).all()
                )
                assert (
                    np.all(births[~survived] == 0)
                    and not survived[S == 0].any()
                    and survived[S == 1].all()
                )
                assert np.all(births[mu == 0] == 0)
            else:
                assert par.offspring_count.isna().all() and par.survived_to_reproduce.isna().all()
            for s, name in enumerate(species_keys):
                mask = labels == s
                row = pop[pop.species == name].iloc[0]
                assert row.N_start == mask.sum() and row.zero_fitness_count == np.count_nonzero(
                    (S * F)[mask] == 0
                )
                assert row.population_status == ("extant" if mask.any() else "extinct")
                close(row.density_phi, phi[s] if k == 0 else 1.0)
                close(row.mean_g, np.mean(S[mask]) if k == 1 and mask.any() else np.nan)
                close(
                    row.mean_h, np.mean(p.h.to_numpy()[mask]) if k == 0 and mask.any() else np.nan
                )
                close(
                    row.total_resource_transfer,
                    p.X.to_numpy()[mask].sum() if k == 0 else ev.x.to_numpy()[mask].sum(),
                )
                if not mask.any() and name not in first_extinct:
                    first_extinct[name] = t
                if name in first_extinct:
                    assert row.extinct_at_generation == first_extinct[name]
                if sampling.any():
                    assert row.N_survived_to_reproduce == survived[mask].sum()
                    assert row.N_offspring == sum(int(x) for x in births[mask])
                if flags(pd.Series([row.update_applied]))[0]:
                    assert row.N_next == row.N_offspring
                    if row.N_next == 0 and name not in first_extinct:
                        first_extinct[name] = t + 1
                        assert row.extinct_at_generation == t + 1
                else:
                    assert pd.isna(row.N_next)
            if flags(pop.update_applied).any():
                assert sampling.all() and t + 1 < len(states)
                nxt = states[t + 1]
                descendants = tables["offspring"]
                child = descendants[
                    (descendants.generation == t + 1) & (descendants.kind == kind)
                ].sort_values("offspring_id")
                next_labels = nxt[kind + "_species"]
                next_ids = nxt[kind + "_ids"]
                assert len(child) == sum(int(x) for x in births) == len(next_ids)
                assert child.offspring_id.tolist() == next_ids.tolist()
                if c.spatial is not None:
                    audit_dispersal(state, nxt, ev, child, kind, c)
                for species_id, name in enumerate(species_keys):
                    assert (next_labels == species_id).sum() == int(
                        pop[pop.species == name].iloc[0].N_next
                    )
                counts = child.parent_id.value_counts().to_dict()
                assert [int(counts.get(int(i), 0)) for i in ids] == births.astype(int).tolist()
                lookup = {int(v): i for i, v in enumerate(ids)}
                indices = child.parent_id.map(lookup)
                assert indices.notna().all()
                indices = indices.to_numpy(int)
                assert child.species.tolist() == [species_keys[s] for s in next_labels]
                assert np.array_equal(labels[indices], next_labels)
                prefix = "P" if k == 0 else "H"
                assert child.parent_genotype_id.tolist() == [f"{prefix}:{v}" for v in ids[indices]]
                if len(child):
                    for key in ("D", "plant_voc_independent_trait") if k == 0 else ("I", "A"):
                        for col in range(nxt[key].shape[1]):
                            close(child[f"{key}_{col + 1}"], nxt[key][:, col])
                            mutated = flags(child[f"{key}_mutated_{col + 1}"])
                            assert np.array_equal(
                                nxt[key][~mutated, col], state[key][indices[~mutated], col]
                            )
            else:
                descendants = tables["offspring"]
                assert descendants[
                    (descendants.generation == t + 1) & (descendants.kind == kind)
                ].empty
        if flags(pop.update_applied).any():
            assert flags(pop.update_applied).all()
            transitions += 1
    assert states[-1]["generation"] == meta["status"]["generation"]
    if status == "completed":
        assert len(states) == c.generations + 1 and transitions == c.generations
    if status == "all_species_extinct":
        assert len(states[-1]["plant_ids"]) + len(states[-1]["insect_ids"]) == 0
    if status == "computational_limit" and len(tables["parents"]):
        end = tables["populations"][tables["populations"].generation == t]
        if end.N_offspring.notna().any():
            from .audit_genetics import genetic_limit

            if c.genetics is not None and c.genetics.get("mode") == "qtl_additive":
                end.attrs["allele_count"] = len(pd.read_csv(path / "allele_catalog.csv"))
            assert (
                end.N_offspring.sum() > c.max_individuals_safety or genetic_limit(c, states, end)
            ) and not flags(end.update_applied).any()
    if c.model_version in ("1.4", "1.5", "1.6"):
        audit_extensions(path, c, meta, states, tables)
    return dict(
        status=status,
        states=len(states),
        transitions=transitions,
        event_rows=len(tables["events"]),
        plant_rows=len(tables["plants"]),
        offspring_rows=len(tables["offspring"]),
        max_overrun=overrun,
        max_balance_error=balance,
        **genetic_checks,
    )


def wrapped_distance(origin, targets, shape):
    delta = np.abs(np.asarray(targets) - np.asarray(origin))
    wrapped = np.where(delta > np.asarray(shape) / 2, np.asarray(shape) - delta, delta)
    return np.max(wrapped, axis=-1)


def spatial_phi(state, alpha, c):
    positions = state["plant_cells"]
    unique, inverse = np.unique(positions, axis=0, return_inverse=True)
    counts = []
    for cell in unique:
        near = (
            wrapped_distance(cell, positions, c.spatial["grid_shape"])
            <= c.spatial["competition_radius"]
        )
        counts.append(np.bincount(state["plant_species"][near], minlength=len(c.plant_sizes)))
    if not len(positions):
        return np.empty(0)
    return 1 / (1 + (alpha[state["plant_species"]] * np.asarray(counts)[inverse]).sum(axis=1))


def audit_space_events(state, ev, c, meta):
    close(ev.search_row, state["insect_cells"][:, 0]) if len(ev) else None
    close(ev.search_col, state["insect_cells"][:, 1]) if len(ev) else None
    plant_ids = state["plant_ids"]
    lookup = {int(v): i for i, v in enumerate(plant_ids)}
    weights = np.asarray(meta["resolved_weights"])
    for i, event in enumerate(ev.to_dict("records")):
        candidates = json.loads(event["candidate_plant_ids"])
        assert len(candidates) == event["candidate_count"] and len(set(candidates)) == len(
            candidates
        )
        rows = np.array([lookup[k] for k in candidates], int)
        initial = int(state["generation"]) == 0
        origin = state["insect_cells"][i]
        near = (
            wrapped_distance(origin, state["plant_cells"], c.spatial["grid_shape"])
            <= c.spatial["search_radius"]
        )
        assert near[rows].all()
        host = int(event["host_row"])
        if initial and c.model_version == "1.5":
            h = f"H{int(state['insect_species'][i]) + 1}"
            target = c.spatial["initial_host_by_insect_species"][h]
            assert event["plant_species"] == target and len(rows) == 1
            assert (
                event["host_assignment"] == "prescribed_initial"
                and event["reason"] == "initial_host"
            )
            assert np.array_equal(state["plant_cells"][rows[0]], origin)
        else:
            assert len(rows) == min(int(near.sum()), c.candidate_size)
            assert event["host_assignment"] == "local_choice"
            if not len(rows):
                assert event["reason"] == "no_candidate" and host == -1
            elif host >= 0:
                assert event["reason"] == "accepted"
            else:
                assert event["reason"] == "rejected"
        if host >= 0:
            assert host in rows
            close([event["host_cell_row"], event["host_cell_col"]], state["plant_cells"][host])
        else:
            assert event["host_cell_row"] == event["host_cell_col"] == -1
        if len(rows):
            V = state["V_cue"][rows]
            A = state["A"][i]
            dist = ((V - A) ** 2) @ weights
            if c.olfaction_mode == "neutral" or c.tau == 0:
                probabilities = np.full(len(rows), 1 / len(rows))
            else:
                logits = -((dist - dist.min()) / c.sigma) / c.sigma / 2
                mass = np.exp(logits * c.tau)
                probabilities = mass / mass.sum()
            for p in range(len(c.plant_sizes)):
                belongs = state["plant_species"][rows] == p
                if c.record_observations:
                    close(event[f"P{p + 1}_candidate_share"], belongs.mean())
                    close(event[f"P{p + 1}_conditional_probability"], probabilities[belongs].sum())


def audit_dispersal(state, nxt, events, children, kind, c):
    cells = nxt[kind + "_cells"]
    if len(children):
        close(children[["cell_row", "cell_col"]], cells)
    lookup = {int(v): i for i, v in enumerate(state[kind + "_ids"])}
    indices = np.array([lookup[int(v)] for v in children.parent_id], int)
    if kind == "plant":
        origins = state["plant_cells"][indices]
    else:
        host = events.host_row.to_numpy(int)[indices]
        assert (host >= 0).all()
        origins = state["plant_cells"][host]
    distance = wrapped_distance(origins, cells, c.spatial["grid_shape"])
    assert (distance <= c.spatial[f"{kind}_dispersal_radius"]).all()


def audit_extensions(path, c, meta, states, tables):
    """Recompute frozen traits, recorded exposure and reference-panel diagnostics."""
    founders = c.founder_traits
    if founders is None:
        return
    for state in states:
        for kind, evolve, keys in (
            ("plant", c.evolve_plants, ("D", "plant_voc_independent_trait")),
            ("insect", c.evolve_insects, ("I", "A")),
        ):
            if (c.genetics is None or c.genetics.get("mode") != "qtl_additive") and (
                int(state["generation"]) == 0 or not evolve
            ):
                for key in keys:
                    assert np.array_equal(
                        state[key], np.asarray(founders[key])[state[kind + "_species"]]
                    )
            if not evolve:
                offspring = tables["offspring"][tables["offspring"].kind == kind]
                for key in keys:
                    for col in range(np.asarray(founders[key]).shape[1]):
                        if len(offspring):
                            assert not flags(offspring[f"{key}_mutated_{col + 1}"]).any()
    if not c.record_observations:
        return
    diagnostics = pd.read_csv(path / "diagnostics.csv", float_precision="round_trip")
    D = np.asarray(founders["D"], float)
    Z = np.asarray(founders["plant_voc_independent_trait"], float)
    V = np.clip(c.rho_voc * (D @ np.asarray(meta["resolved_B"]).T) + (1 - c.rho_voc) * Z, 0, 1)
    weights = np.asarray(meta["resolved_weights"])
    expected_rows = 0
    for state in states:
        t = int(state["generation"])
        ev = tables["events"][tables["events"].generation == t]
        net = tables["networks"][tables["networks"].generation == t]
        diag = diagnostics[diagnostics.generation == t]
        if net.empty:
            assert diag.empty
            continue
        expected_rows += len(c.plant_sizes) * len(c.insect_sizes)
        assert len(diag) == len(c.plant_sizes) * len(c.insect_sizes)
        assert diag.panel.eq("equal_ancestral_species").all()
        if len(ev):
            exposed = ev.candidate_count.to_numpy() > 0
            shares = ev[[f"P{p + 1}_candidate_share" for p in range(len(D))]].to_numpy()
            probs = ev[[f"P{p + 1}_conditional_probability" for p in range(len(D))]].to_numpy()
            for values in (shares, probs):
                assert np.isnan(values[~exposed]).all()
                assert (
                    np.isfinite(values[exposed]).all()
                    and ((values[exposed] >= -1e-12) & (values[exposed] <= 1 + 1e-12)).all()
                )
                close(values[exposed].sum(axis=1), 1)
            count = shares[exposed] * ev.candidate_count.to_numpy()[exposed, None]
            close(count, np.round(count))
        for h in range(len(c.insect_sizes)):
            name = f"H{h + 1}"
            selected = state["insect_species"] == h
            A, I = state["A"][selected], state["I"][selected]
            distance = ((A[:, None, :] - V[None, :, :]) ** 2) @ weights
            if c.olfaction_mode == "neutral" or c.tau == 0:
                probabilities = np.full(distance.shape, 1 / len(D))
            else:
                logits = -((distance - distance.min(axis=1, keepdims=True)) / c.sigma) / c.sigma / 2
                values = np.exp(logits * c.tau)
                probabilities = values / values.sum(axis=1, keepdims=True)
            ratios = np.divide(
                I[:, None, :],
                D[None, :, :],
                out=np.full((len(I), *D.shape), np.inf),
                where=D[None, :, :] > 0,
            )
            q = np.minimum(1, ratios.min(axis=2))
            hev = ev[ev.species == name]
            for p in range(len(D)):
                pname = f"P{p + 1}"
                row = diag[(diag.insect_species == name) & (diag.reference_plant_species == pname)]
                assert len(row) == 1
                row = row.iloc[0]
                assert row.insect_count == len(I)
                close(row.preference, probabilities[:, p].mean() if len(I) else np.nan)
                close(row.q_unlimited, q[:, p].mean() if len(I) else np.nan)
                edge = net[(net.insect_species == name) & (net.plant_species == pname)].iloc[0]
                exposure = (
                    (hev.candidate_count * hev[f"P{p + 1}_candidate_share"]).sum()
                    if len(hev)
                    else 0
                )
                close(edge.candidate_exposure, exposure)
                if not (state["plant_species"] == p).any() and len(hev):
                    close(hev.loc[hev.candidate_count > 0, f"P{p + 1}_candidate_share"], 0)
    assert len(diagnostics) == expected_rows


def audit(out, report_dir=None):
    out = Path(out)
    settings = json.loads((out / "settings.json").read_text(encoding="utf-8"))
    report_dir = Path(report_dir) if report_dir is not None else report_directory(out)
    report_dir.mkdir(parents=True, exist_ok=True)
    status = pd.read_csv(out / "data" / "run_status.csv")
    rows = []
    assert len(status) == len(settings["spec"]["cases"]) * settings["spec"]["replicates"]
    for case in settings["spec"]["cases"]:
        for rep in range(1, settings["spec"]["replicates"] + 1):
            path = out / "data" / case["name"] / f"rep{rep:02d}"
            meta = json.loads((path / "metadata.json").read_text(encoding="utf-8"))
            assert meta["source_hash"] == settings["source_sha256"]
            expected = dict(settings["spec"]["base"], **case["overrides"])
            expected["seed"] += rep - 1
            if meta["status"]["status"] != "invalid_parameter":
                expected = config_payload(Config(**expected))
            assert meta["config"] == expected
            checked = audit_run(path)
            row = status[(status["case"] == case["name"]) & (status.replicate == rep)].iloc[0]
            assert (
                checked["status"] == row.status and meta["status"]["generation"] == row.generation
            )
            rows.append(dict(case=case["name"], replicate=rep, **checked))
            print("AUDITED", case["name"], rep, checked["status"], flush=True)
    pd.DataFrame(rows).to_csv(report_dir / "development_checks.csv", index=False)
    summary = dict(
        runs=len(rows),
        outcomes=pd.Series([r["status"] for r in rows]).value_counts().to_dict(),
        batch=str(out.resolve()),
        source_matches_current=settings["source_sha256"] == source_hashes(),
        source_verification="comparison with current files only; historical source snapshots are not verified",
        states=sum(r["states"] for r in rows),
        transitions=sum(r["transitions"] for r in rows),
        event_rows=sum(r.get("event_rows", 0) for r in rows),
        plant_rows=sum(r.get("plant_rows", 0) for r in rows),
        offspring_rows=sum(r.get("offspring_rows", 0) for r in rows),
        max_resource_overrun=max(r["max_overrun"] for r in rows),
        max_relative_balance_error=max(r["max_balance_error"] for r in rows),
    )
    if any(r.get("genetic_audit", False) for r in rows):
        summary.update(
            genetically_audited_runs=sum(r.get("genetic_audit", False) for r in rows),
            alleles_checked=sum(r.get("alleles_checked", 0) for r in rows),
            mutations_checked=sum(r.get("mutations_checked", 0) for r in rows),
        )
    write_json(report_dir / "integrity_audit.json", summary)
    print(json.dumps(summary, indent=2))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    audit(parser.parse_args().input)


if __name__ == "__main__":
    main()
