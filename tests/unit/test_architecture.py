"""Enforce dependency boundaries and the retained command/import interfaces."""

import ast
import importlib
import importlib.util
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "evo_sim"
MODEL = "evo_sim.models.plant_insect"


def test_dependency_boundaries_and_independent_audit():
    for path in PACKAGE.rglob("*.py"):
        module = ".".join(path.relative_to(ROOT).with_suffix("").parts)
        package = module.rsplit(".", 1)[0]
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Import):
                imports = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                base = (
                    importlib.util.resolve_name("." * node.level + (node.module or ""), package)
                    if node.level
                    else node.module
                )
                imports = [base] + [f"{base}.{alias.name}" for alias in node.names]
            else:
                continue
            forbidden = []
            if module.startswith("evo_sim.core."):
                forbidden.append("evo_sim.models")
            if module.startswith(f"{MODEL}.biology."):
                forbidden.extend(
                    f"{MODEL}.{layer}" for layer in ("simulation", "analysis", "audit")
                )
            if module.startswith(f"{MODEL}.config."):
                forbidden.append(f"{MODEL}.analysis")
            if module.startswith(f"{MODEL}.audit."):
                forbidden.append(f"{MODEL}.biology")
            for imported in imports:
                assert not any(
                    imported == prefix or imported.startswith(prefix + ".") for prefix in forbidden
                ), (path, node.lineno, imported)


def test_legacy_analysis_commands_and_imports():
    modules = {
        "summarize_lifecycle": ("summaries.lifecycle", "summarize"),
        "summarize_lineages": ("summaries.lineages", "summarize"),
        "summarize_qtl_run": ("summaries.qtl", "summarize"),
        "plot_diagnostics": ("plots.diagnostics", "main"),
    }
    for old, (new, function) in modules.items():
        legacy = importlib.import_module(f"{MODEL}.analysis.{old}")
        current = importlib.import_module(f"{MODEL}.analysis.{new}")
        assert getattr(legacy, function) is getattr(current, function)
    for name in ["audit"] + [f"analysis.{old}" for old in modules]:
        result = subprocess.run(
            [sys.executable, "-m", f"{MODEL}.{name}", "--help"],
            cwd=ROOT,
            capture_output=True,
            timeout=30,
        )
        assert result.returncode == 0 and b"--input" in result.stdout, (name, result.stderr)
