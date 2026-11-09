"""Independent events checks; equations preserved from v0.1.12 audit."""

import hashlib
import json

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect.config.schema import Config, config_payload


def close(actual, expected, context):
    assert np.allclose(
        np.asarray(actual, float),
        np.asarray(expected, float),
        rtol=1e-11,
        atol=1e-14,
        equal_nan=True,
    ), context


def require_columns(frame, columns, name):
    assert set(columns) <= set(frame.columns), (
        name,
        "missing columns",
        set(columns) - set(frame.columns),
    )


def unique_keys(frame, columns, name):
    require_columns(frame, columns, name)
    assert frame[list(columns)].notna().all().all(), (name, "missing key")
    assert not frame.duplicated(list(columns)).any(), (name, "duplicate key")


def flag(value):
    assert value in (True, False, "True", "False"), ("invalid boolean", value)
    return value in (True, "True")


def host_name(value):
    return "unknown" if value < 0 else f"P{value + 1}"


def distance(origin, targets, shape):
    delta = np.abs(np.asarray(targets) - origin)
    return np.minimum(delta, shape - delta).max(axis=-1)


def compare(actual, expected, context):
    if expected is None or (isinstance(expected, float) and np.isnan(expected)):
        assert pd.isna(actual), (context, actual, expected)
    elif isinstance(expected, (bool, np.bool_)):
        assert flag(actual) == expected, (context, actual, expected)
    elif isinstance(expected, str):
        assert actual == expected, (context, actual, expected)
    else:
        close(actual, expected, context)


def check_rows(frame, keys, expected, name):
    """Require exact key coverage; check every declared expected value."""
    unique_keys(frame, keys, name)
    indexed = {tuple(row[k] for k in keys): row for row in frame.to_dict("records")}
    assert set(indexed) == set(expected), (name, "missing/extra keys", set(expected) ^ set(indexed))
    checks = {}
    for key, values in expected.items():
        for column, value in values.items():
            assert column in frame.columns, (name, "missing column", column)
            kind = (
                "missing"
                if value is None or (isinstance(value, float) and np.isnan(value))
                else "bool"
                if isinstance(value, (bool, np.bool_))
                else "str"
                if isinstance(value, str)
                else "numeric"
            )
            actual, wanted = checks.setdefault((column, kind), ([], []))
            actual.append(indexed[key][column])
            wanted.append(value)
    for (column, kind), (actual, wanted) in checks.items():
        context = (name, column, kind)
        if kind == "missing":
            assert pd.isna(actual).all(), context
        elif kind == "bool":
            assert [flag(value) for value in actual] == wanted, context
        elif kind == "str":
            assert actual == wanted, context
        else:
            close(actual, wanted, context)


def read_records(path, meta):
    """Read and verify required raw files, configuration fingerprints and snapshots."""
    assert __debug__, "scientific audits require assertions"
    c = Config(**meta["config"])
    c.validate()
    assert c.model_version == meta["version"] == "0.1.12"
    assert meta["genetics_version"] == c.genetics["schema_version"]
    assert meta["genetics_version"] == "slim-qtl-6"
    assert meta["ecology_rule_version"] == c.model_version
    assert meta["functional_schema"] == "expression-tuning-contact-3"
    assert meta["lifecycle_schema"] == "dose-tolerance-lifecycle-2"
    payload = json.dumps(config_payload(c), sort_keys=True, separators=(",", ":")).encode()
    assert hashlib.sha256(payload).hexdigest() == meta["config_sha256"]
    assert json.loads((path / "status.json").read_text(encoding="utf-8")) == meta["status"]
    for name, digest in meta["output_sha256"].items():
        assert hashlib.sha256((path / name).read_bytes()).hexdigest() == digest, name
    names = (
        "events",
        "plants",
        "plant_proposals",
        "parents",
        "offspring",
        "populations",
        "networks",
        "oviposition",
        "larvae",
        "host_flows",
        "host_groups",
        "lagged_networks",
        "allele_contributions",
        "gamete_transmissions",
        "mutations",
        "allele_catalog",
    )
    tables = {
        name: pd.read_csv(path / f"{name}.csv", float_precision="round_trip", low_memory=False)
        for name in names
    }
    for name in ("offspring", "mutations", "gamete_transmissions"):
        frame = tables[name]
        for key, offset in (("reproduction_round", -1), ("adult_cohort", 0)):
            if meta.get("stage_key_schema") == "reproduction-round-adult-cohort-1":
                assert key in frame, (name, "missing stage key", key)
            if key in frame:
                close(frame[key], frame.generation + offset, (name, key))
    assert {f"{name}.csv" for name in names} <= set(meta["output_sha256"])
    assert {
        "genetic_loci.csv",
        "genetic_summary.csv",
        "haplotype_summary.csv",
        "genotypes.npz",
        "phenotypes.npz",
    } <= set(meta["output_sha256"])
    if c.record_observations:
        assert "diagnostics.csv" in meta["output_sha256"]
    files = sorted((path / "snapshots").glob("*.npz"))
    assert {file.relative_to(path).as_posix() for file in files} == {
        name for name in meta["output_sha256"] if name.startswith("snapshots/")
    }
    states = []
    for file in files:
        with np.load(file) as loaded:
            states.append(dict(loaded))
    assert [int(s["generation"]) for s in states] == list(range(len(states)))
    status = meta["status"]["status"]
    return c, names, tables, states, status
