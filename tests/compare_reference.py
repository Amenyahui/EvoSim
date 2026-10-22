"""Optional migration check against an explicitly supplied original checkout."""

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from evo_sim.models.plant_insect import Config, simulate
from tests.test_runtime import same_results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    args = parser.parse_args()
    package = args.reference / "coevo"
    spec = importlib.util.spec_from_file_location(
        "reference_coevo", package / "__init__.py", submodule_search_locations=[str(package)]
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    cases = []
    for path in sorted((Path(__file__).parent / "fixtures").glob("*.json")):
        fixture = json.loads(path.read_text())
        for case in fixture["cases"]:
            values = dict(fixture["base"], **case["overrides"])
            values["generations"] = min(values["generations"], 4)
            try:
                actual = simulate(Config(**values))
            except (ValueError, TypeError) as error:
                try:
                    module.simulate(module.Config(**values))
                except type(error) as expected:
                    assert str(error) == str(expected)
                else:
                    raise AssertionError("New implementation rejects valid reference input")
            else:
                same_results(actual, module.simulate(module.Config(**values)))
            cases.append(f"{path.stem}/{case['name']}")
    report = Path(__file__).parent / "reports/refactor_equivalence.json"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        json.dumps(
            {
                "reference": str(args.reference.resolve()),
                "exact_cases": cases,
                "comparison": "all record tables and every snapshot array, exact equality; terminal status equality",
            },
            indent=2,
        )
    )
    print(f"PASS {len(cases)} reference cases; every numeric record and snapshot matches exactly")


if __name__ == "__main__":
    main()
