import csv
import hashlib
import io
from pathlib import Path
from tempfile import TemporaryDirectory

import yaml

from evo_sim.models.plant_insect import Config, get_runner, simulate
from evo_sim.models.plant_insect.audit.run import audit_run
from evo_sim.models.plant_insect.batch import save_result
from evo_sim.models.plant_insect.config.input_loader import load_inputs
from evo_sim.models.plant_insect.config.loader import read_config
from tests.integration.test_runtime import same_results
from tests.support import raises

MODEL = Path(__file__).resolve().parents[2] / "evo_sim/models/plant_insect"


def experiment(folder, names):
    data = folder / "data"
    data.mkdir()
    for name in names:
        (data / f"{name}.csv").write_bytes((MODEL / "inputs/current" / f"{name}.csv").read_bytes())
    path = folder / "experiment.yaml"
    path.write_text(
        yaml.safe_dump(
            dict(
                extends=str(MODEL / "configs/current_base.yaml"),
                input_files={name: f"data/{name}.csv" for name in names},
            )
        ),
        encoding="utf-8",
    )
    return path


def records(path):
    return list(csv.DictReader(io.StringIO(path.read_text(encoding="utf-8-sig"))))


def write_records(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def validate_file(path):
    loaded = read_config(path)
    Config(**dict(loaded["parameters"], generations=1)).validate()


def test_csv_relative_paths_overrides_and_self_contained_results():
    with TemporaryDirectory() as tmp:
        folder = Path(tmp)
        experiment(folder, ("species", "alleles"))
        species = folder / "data/species.csv"
        rows = records(species)
        rows[0]["initial_size"] = "12"
        rows[1]["initial_size"] = "9"
        rows[2]["initial_size"] = "6"
        write_records(species, list(reversed(rows)))  # Numeric IDs retain their species slots.
        alleles = folder / "data/alleles.csv"
        rows = records(alleles)
        rows[0]["effect_0"] = "0.15"
        write_records(alleles, rows)
        nested = folder / "nested"
        nested.mkdir()
        child = nested / "child.yaml"
        child.write_text(
            "extends: ../experiment.yaml\nparameters:\n  plant_sizes: [10, 8]\n",
            encoding="utf-8",
        )
        loaded = read_config(child)
        assert loaded["parameters"]["plant_sizes"] == [10, 8]
        assert loaded["parameters"]["insect_sizes"] == [6]
        assert loaded["parameters"]["genetics"]["alleles"][0]["effect"][0] == 0.15
        expected_hash = hashlib.sha256(alleles.read_bytes()).hexdigest()
        assert loaded["input_sources"]["alleles"] == dict(path=str(alleles), sha256=expected_hash)
        c = Config(**dict(loaded["parameters"], generations=1))
        result = simulate(c, runtime_config=loaded)
        assert result["metadata"]["status"]["status"] == "completed"
        same_results(result, simulate(c))
        assert result["metadata"]["input_sources"] == loaded["input_sources"]
        output = folder / "output"
        save_result(result, output)
        species.unlink()
        alleles.unlink()
        assert audit_run(output)["genetic_audit"]  # Audit uses saved values, not live inputs.


def test_replacing_csv_tables_removes_old_species_instead_of_merging():
    with TemporaryDirectory() as tmp:
        folder = Path(tmp)
        path = experiment(folder, ("species", "founder_frequencies", "reference_traits"))
        for name in ("species", "founder_frequencies", "reference_traits"):
            file = folder / "data" / f"{name}.csv"
            write_records(file, [row for row in records(file) if row["species_id"] != "P2"])
        loaded = read_config(path)
        p = loaded["parameters"]
        assert p["plant_sizes"] == [40]
        assert set(p["plant_offspring_lambda"]) == {"P1"}
        assert set(p["resource_capacity"]) == {"P1"}
        assert set(p["genetics"]["founder_frequencies"]) == {"P1", "H1"}
        assert set(p["parameter_provenance"]["plant_offspring_lambda"]) == {"P1"}
        validate_file(path)


def test_csv_rejects_malformed_records_and_invalid_references():
    with TemporaryDirectory() as tmp:
        folder = Path(tmp)
        path = experiment(folder, ("alleles", "founder_frequencies"))
        alleles, frequencies = folder / "data/alleles.csv", folder / "data/founder_frequencies.csv"
        original, founders = records(alleles), records(frequencies)
        bad_rows = [
            original + [original[0]],
            [dict(original[0], effect_0="nan"), *original[1:]],
            [dict(original[0], effect_0=""), *original[1:]],
            [dict(original[0], allele_id="1.5"), *original[1:]],
            [dict(original[0], locus_id="unknown"), *original[1:]],
        ]
        for rows in bad_rows:
            write_records(alleles, rows)
            raises(lambda: validate_file(path))
        write_records(alleles, original)
        text = alleles.read_text(encoding="utf-8-sig")
        for invalid in (
            text.replace("effect_0", "unexpected", 1),
            text.splitlines()[0] + "\n",  # Header-only file.
            text + "0,plant_1_1\n",  # Too few cells.
            text + "0,plant_1_1,0,0,0,0,extra\n",
        ):
            alleles.write_text(invalid, encoding="utf-8")
            raises(lambda: read_config(path))
        write_records(alleles, original)
        for changes in (dict(frequency="0.7"), dict(species_id="P3"), dict(allele_id="999")):
            write_records(frequencies, [dict(founders[0], **changes), *founders[1:]])
            raises(lambda: validate_file(path))
        write_records(frequencies, founders)
        alleles.unlink()
        raises(lambda: read_config(path), FileNotFoundError)


def test_csv_vector_columns_follow_current_trait_width():
    from tests.support import config

    c = config()
    parameters = dict(n=2, B=None, mapping=c.mapping, genetics=c.genetics)
    locus = next(row["locus_id"] for row in c.genetics["loci"] if row["organism"] == "insect")
    with TemporaryDirectory() as tmp:
        folder = Path(tmp)
        path = folder / "alleles.csv"
        path.write_text(
            "allele_id,locus_id,"
            + ",".join(f"effect_{i}" for i in range(10))
            + f"\n0,{locus},0.1,0,0,0,0,0,0,0,0,0\n"
        )
        loaded, _ = load_inputs(dict(alleles="alleles.csv"), folder, parameters)
        assert loaded["genetics"]["alleles"][0]["effect"] == [0.1] + [0.0] * 9
        path.write_text(path.read_text().replace(",effect_7", "").replace(",0\n", "\n"))
        raises(lambda: load_inputs(dict(alleles="alleles.csv"), folder, parameters))
        raises(lambda: load_inputs(dict(unknown="alleles.csv"), folder, parameters))


def test_builtin_csv_inputs_are_complete_and_resolved():
    loaded = get_runner().config
    assert len(loaded["input_sources"]) == 5 and "input_files" not in loaded
    for source in loaded["input_sources"].values():
        assert hashlib.sha256(Path(source["path"]).read_bytes()).hexdigest() == source["sha256"]
    p = loaded["parameters"]
    assert p["plant_sizes"] == [20] * 10 and p["insect_sizes"] == [200]
    assert len(p["genetics"]["loci"]) == 30 and len(p["genetics"]["alleles"]) == 168
