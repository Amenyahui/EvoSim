"""Run the project's assertion-based checks without an extra test dependency."""

import importlib
from pathlib import Path

count = 0
for path in sorted(Path(__file__).parent.glob("test_*.py")):
    module = importlib.import_module(f"tests.{path.stem}")
    for name, check in vars(module).items():
        if name.startswith("test_") and callable(check):
            check()
            count += 1
            print(f"PASS {path.stem}.{name}", flush=True)
print(f"PASS {count} test groups")
