"""Run the project's assertion-based checks without an extra test dependency."""

import importlib
from pathlib import Path

from evo_sim.models.plant_insect.biology.genetics.slim_backend import SlimError, engine_info

if not __debug__:
    raise SystemExit("Run tests without python -O; assertions are required")
try:
    slim = engine_info(dict(executable="slim", timeout_seconds=10))
except SlimError as error:
    raise SystemExit(
        f"SLiM preflight failed: {error}. See README.md installation instructions."
    ) from error
print(f"SLiM preflight: {slim['version'].splitlines()[0]} ({slim['executable']})", flush=True)

count = 0
root = Path(__file__).parent
for path in sorted(root.rglob("test_*.py")):
    module_name = ".".join(path.relative_to(root).with_suffix("").parts)
    module = importlib.import_module(f"tests.{module_name}")
    for name, check in vars(module).items():
        if name.startswith("test_") and callable(check):
            check()
            count += 1
            print(f"PASS {module_name}.{name}", flush=True)
print(f"PASS {count} test groups")
