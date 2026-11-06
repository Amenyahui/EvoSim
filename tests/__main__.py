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
for path in sorted(Path(__file__).parent.glob("test_*.py")):
    module = importlib.import_module(f"tests.{path.stem}")
    for name, check in vars(module).items():
        if name.startswith("test_") and callable(check):
            check()
            count += 1
            print(f"PASS {path.stem}.{name}", flush=True)
print(f"PASS {count} test groups")
