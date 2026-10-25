"""Read editable CSV inputs into the existing model parameter contract."""

import csv
import hashlib
import io
import math
from copy import deepcopy

from evo_sim.core.config import read_config as read_yaml_config


def finite(value):
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("finite number required")
    return number


def optional_number(value):
    return finite(value) if value else None


def nonempty(value):
    if not value:
        raise ValueError("empty field")
    return value


def read_table(path, columns, keys):
    """Check headers, types and duplicate records; retain file row order."""
    data = path.read_bytes()
    reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig")))
    if (
        reader.fieldnames is None
        or len(reader.fieldnames) != len(columns)
        or set(reader.fieldnames) != set(columns)
    ):
        raise ValueError(f"{path}: expected CSV columns {', '.join(columns)}")
    records, seen = [], set()
    for row in reader:
        try:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("row has the wrong number of cells")
            record = {}
            for key, convert in columns.items():
                record[key] = convert(row[key].strip())
            identity = tuple(record[key] for key in keys)
            if identity in seen:
                raise ValueError(f"duplicate record {identity}")
            seen.add(identity)
            records.append(record)
        except (ValueError, TypeError) as error:
            raise ValueError(f"{path}, line {reader.line_num}: {error}") from error
    if not records:
        raise ValueError(f"{path}: empty CSV table")
    return records, dict(path=str(path), sha256=hashlib.sha256(data).hexdigest())


def load_inputs(files, directory, parameters):
    allowed = {
        "species",
        "competition",
        "loci",
        "alleles",
        "founder_frequencies",
        "reference_traits",
    }
    if not isinstance(files, dict) or not files or not set(files) <= allowed:
        raise ValueError(f"input_files requires paths keyed by {sorted(allowed)}")
    p = deepcopy(parameters)
    n, B = p["n"], p["B"]
    if type(n) is not int or n < 1 or not isinstance(B, list) or not B:
        raise ValueError("CSV inputs require explicit positive n and nonempty B")
    k, width = len(B), n + len(B)
    sources = {}

    def table(name, columns, keys):
        value = files[name]
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"input_files.{name} must be a nonempty CSV path")
        path = (directory / value).resolve()
        rows, sources[name] = read_table(path, columns, keys)
        return rows

    if "species" in files:
        rows = table(
            "species",
            dict(
                species_id=nonempty,
                organism=nonempty,
                initial_size=int,
                resource_capacity=optional_number,
                max_feeding_demand=optional_number,
                offspring_lambda=finite,
                source_kind=nonempty,
                source=nonempty,
            ),
            ("species_id",),
        )
        if any(row["organism"] not in ("plant", "insect") for row in rows):
            raise ValueError("species.csv organism must be plant or insect")
        provenance = p.setdefault("parameter_provenance", {})
        for kind, prefix, resource in (
            ("plant", "P", "resource_capacity"),
            ("insect", "H", "max_feeding_demand"),
        ):
            group = [row for row in rows if row["organism"] == kind]
            if not group or {row["species_id"] for row in group} != {
                f"{prefix}{i + 1}" for i in range(len(group))
            }:
                raise ValueError(f"species.csv requires consecutive {prefix}1..{prefix}N IDs")
            # Species slots have always been numeric; table row order must not change their IDs.
            group.sort(key=lambda row: int(row["species_id"][1:]))
            other = "max_feeding_demand" if kind == "plant" else "resource_capacity"
            if any(row[resource] is None or row[other] is not None for row in group):
                raise ValueError(f"species.csv {kind}: fill {resource}, leave {other} empty")
            p[kind + "_sizes"] = [row["initial_size"] for row in group]
            p[resource] = {row["species_id"]: row[resource] for row in group}
            field = kind + "_offspring_lambda"
            p[field] = {row["species_id"]: row["offspring_lambda"] for row in group}
            provenance[field] = {
                row["species_id"]: dict(kind=row["source_kind"], source=row["source"])
                for row in group
            }

    if "competition" in files:
        rows = table(
            "competition",
            dict(
                recipient=nonempty,
                competitor=nonempty,
                alpha=finite,
                source_kind=nonempty,
                source=nonempty,
            ),
            ("recipient", "competitor"),
        )
        matrix, provenance = {}, {}
        for row in rows:
            recipient, competitor = row["recipient"], row["competitor"]
            matrix.setdefault(recipient, {})[competitor] = row["alpha"]
            provenance.setdefault(recipient, {})[competitor] = dict(
                kind=row["source_kind"], source=row["source"]
            )
        p["plant_competition_alpha"] = matrix
        p.setdefault("parameter_provenance", {})["plant_competition_alpha"] = provenance

    genetics = p.setdefault("genetics", {})
    if "loci" in files:
        columns = dict(
            locus_id=nonempty,
            organism=nonempty,
            functional_module=nonempty,
            chromosome_id=nonempty,
            position_cM=finite,
            mutation_rate=finite,
        )
        for field, convert in (("target", int), ("mutation_mean", finite), ("mutation_sd", finite)):
            columns.update({f"{field}_{j}": convert for j in range(width)})
        rows = table("loci", columns, ("locus_id",))
        loci = []
        for row in rows:
            targets = [row.pop(f"target_{j}") for j in range(width)]
            if any(value not in (0, 1) for value in targets):
                raise ValueError(f"loci.csv {row['locus_id']}: targets must be 0 or 1")
            row["targets"] = [j for j, value in enumerate(targets) if value]
            for field in ("mutation_mean", "mutation_sd"):
                row[field] = [row.pop(f"{field}_{j}") for j in range(width)]
            loci.append(
                {
                    key: row[key]
                    for key in (
                        "locus_id",
                        "organism",
                        "functional_module",
                        "targets",
                        "mutation_rate",
                        "mutation_mean",
                        "mutation_sd",
                        "chromosome_id",
                        "position_cM",
                    )
                }
            )
        genetics["loci"] = loci

    if "alleles" in files:
        columns = dict(allele_id=int, locus_id=nonempty)
        columns.update({f"effect_{j}": finite for j in range(width)})
        rows = table("alleles", columns, ("allele_id",))
        for row in rows:
            row["effect"] = [row.pop(f"effect_{j}") for j in range(width)]
        genetics["alleles"] = rows

    if "founder_frequencies" in files:
        rows = table(
            "founder_frequencies",
            dict(species_id=nonempty, locus_id=nonempty, allele_id=int, frequency=finite),
            ("species_id", "locus_id", "allele_id"),
        )
        frequencies = {}
        for row in rows:
            entries = frequencies.setdefault(row["species_id"], {}).setdefault(row["locus_id"], [])
            entries.append(dict(allele_id=row["allele_id"], frequency=row["frequency"]))
        genetics["founder_frequencies"] = frequencies

    if "reference_traits" in files:
        columns = dict(species_id=nonempty)
        columns.update({f"D_{j + 1}": finite for j in range(n)})
        columns.update({f"Z_{j + 1}": finite for j in range(k)})
        rows = table("reference_traits", columns, ("species_id",))
        expected = [f"P{i + 1}" for i in range(len(p["plant_sizes"]))]
        by_species = {row["species_id"]: row for row in rows}
        if set(by_species) != set(expected):
            raise ValueError("reference_traits.csv requires exactly the configured plant species")
        p["observation_reference_traits"] = {
            name: [[by_species[s][f"{prefix}_{j + 1}"] for j in range(size)] for s in expected]
            for name, prefix, size in (("D", "D", n), ("plant_voc_independent_trait", "Z", k))
        }
    return p, sources


def read_config(path):
    """Resolve YAML inheritance and CSV paths relative to their declaring YAML."""
    return read_yaml_config(path, input_loader=load_inputs)
