"""Real SLiM subprocess: native meiosis, explicit QTL mutation, validated return material."""

import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

SCRIPT = Path(__file__).with_name("slim_reproduction.slim")


class SlimError(ValueError):
    """Engine/transport failure; handled by the existing uncommitted-transition path."""


def engine_info(settings):
    requested = os.environ.get("EVOSIM_SLIM", settings["executable"])
    executable = shutil.which(requested)
    if executable is None:
        # Local development install; installed packages should use PATH, YAML or EVOSIM_SLIM.
        candidate = Path(__file__).resolve().parents[4] / "codex/slim_runtime/Library/bin/slim.exe"
        if requested == "slim" and candidate.is_file():
            executable = str(candidate)
    if executable is None:
        raise SlimError(
            f"SLiM executable not found: {requested}; configure genetics.slim.executable"
        )
    try:
        result = subprocess.run(
            [executable, "-v"], capture_output=True, text=True, timeout=settings["timeout_seconds"]
        )
    except (OSError, subprocess.SubprocessError) as error:
        raise SlimError(f"Cannot start SLiM: {error}") from error
    version = (result.stdout + result.stderr).strip()
    if result.returncode or not re.search(r"SLiM version 5\.2(?:\D|$)", version):
        raise SlimError(f"SLiM 5.2 required; executable reports: {version[:500]}")
    return dict(
        executable=str(Path(executable).resolve()),
        version=version,
        binary_sha256=hashlib.sha256(Path(executable).read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),
    )


def _constant(name, value):
    """Only generated numeric literals enter Eidos; user strings never become executable code."""
    array = np.asarray(value).reshape(-1)
    if array.dtype.kind not in "iuf" or not np.isfinite(array).all():
        raise SlimError(f"Invalid numeric SLiM input: {name}")
    values = ",".join(repr(x.item()) for x in array)
    return f'defineConstant("{name}", c({values}));\n'


def reproduce(g, mothers, fathers, loci, effects, settings, seeds, next_allele):
    """Return native SLiM offspring plus homologs and explicit mutation deltas/effects."""
    info = engine_info(settings)
    n, width = len(mothers), effects.shape[1]
    groups = {}
    for l, locus in enumerate(loci):
        groups.setdefault(locus["chromosome_id"], []).append(l)
    chromosomes = np.zeros(len(loci), dtype=int)
    positions = np.zeros(len(loci), dtype=int)
    rates = np.zeros(len(loci))
    for chromosome, indices in enumerate(groups.values(), 1):
        indices.sort(key=lambda l: loci[l]["position_cM"])
        for position, l in enumerate(indices):
            chromosomes[l], positions[l] = chromosome, position
            if position:
                distance = (
                    loci[l]["position_cM"] - loci[indices[position - 1]]["position_cM"]
                ) / 100
                rates[l] = -np.expm1(-2 * distance) / 2
    arrays = dict(
        NP=len(g),
        NL=len(loci),
        WIDTH=width,
        PARENT_G=g,
        EFFECTS=effects,
        MOTHERS=mothers,
        FATHERS=fathers,
        CHROMOSOMES=chromosomes,
        POSITIONS=positions,
        RATES=rates,
        MU=[x["mutation_rate"] for x in loci],
        MEANS=[x["mutation_mean"] for x in loci],
        SDS=[x["mutation_sd"] for x in loci],
        NEXT_ALLELE=next_allele,
        MUTATION_SEED=seeds[1],
    )
    # ponytail: one bounded process per taxon/generation keeps failed transitions disposable;
    # use a persistent worker only if measured serialization cost warrants rollback complexity.
    with TemporaryDirectory(prefix="evosim_slim_") as temporary:
        folder = Path(temporary)
        source, output = folder / "input.eidos", folder / "offspring.csv"
        source.write_text("".join(_constant(k, v) for k, v in arrays.items()), encoding="utf-8")
        command = [
            info["executable"],
            "-s",
            str(seeds[0]),
            "-d",
            "INPUT=" + json.dumps(source.as_posix()),
            "-d",
            "OUTPUT=" + json.dumps(output.as_posix()),
            str(SCRIPT),
        ]
        try:
            completed = subprocess.run(
                command, capture_output=True, text=True, timeout=settings["timeout_seconds"]
            )
        except (OSError, subprocess.SubprocessError) as error:
            raise SlimError(f"SLiM reproduction failed: {error}") from error
        if completed.returncode or not output.is_file():
            detail = (completed.stderr + completed.stdout)[-3000:]
            raise SlimError(f"SLiM reproduction failed ({completed.returncode}): {detail}")
        columns = ["child", "copy", "locus", "inherited", "allele", "homolog", "mutated"]
        columns += [f"delta{j}" for j in range(width)] + [f"effect{j}" for j in range(width)]
        try:
            with output.open(encoding="utf-8", newline="") as stream:
                reader = csv.DictReader(stream)
                if reader.fieldnames != columns:
                    raise ValueError("unexpected output columns")
                rows = list(reader)
            if any(set(row) != set(columns) or None in row.values() for row in rows):
                raise ValueError("incorrect output row width")
            if len(rows) != n * 2 * len(loci):
                raise ValueError("incorrect offspring transmission count")
            data = np.array([[float(row[key]) for key in columns] for row in rows])
            if (
                not np.isfinite(data).all()
                or not np.equal(data[:, :7], np.floor(data[:, :7])).all()
            ):
                raise ValueError("nonfinite or noninteger transmission data")
            if np.abs(data[:, :7]).max(initial=0) >= 2**53:
                raise ValueError("transmission IDs exceed the lossless numeric transport range")
            indices = data[:, :3].astype(np.int64)
            expected = np.indices((n, 2, len(loci))).reshape(3, -1).T
            if not np.array_equal(indices, expected):
                raise ValueError("missing, duplicate or reordered transmission keys")
            material = data.reshape(n, 2, len(loci), -1)
            if not np.isin(material[..., 5], [0, 1]).all():
                raise ValueError("invalid homolog")
            if not np.isin(material[..., 6], [0, 1]).all():
                raise ValueError("invalid mutation flag")
        except (OSError, ValueError, KeyError, TypeError) as error:
            raise SlimError(f"Invalid SLiM offspring material: {error}") from error
    info.update(meiosis_seed=int(seeds[0]), mutation_seed=int(seeds[1]))
    return material, info
