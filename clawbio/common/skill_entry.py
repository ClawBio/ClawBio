"""Which file runs a skill, and which standard flags it takes, read from the skill folder.

Shared by scripts/generate_catalog.py (what the catalog advertises) and `clawbio run`
(what it launches for a skill with no SKILLS entry), so both always agree and neither
needs skills/catalog.json. Stdlib only.
"""

from __future__ import annotations

import ast
from pathlib import Path

# Ordered by precedence, so a skill shipping more than one language resolves
# the same way on every run rather than flapping.
ENTRY_POINT_SUFFIXES = (".py", ".R", ".sh")

# Not entry points: package plumbing, and the thin API wrappers skills expose
# for the orchestrator rather than for a human.
NON_ENTRY_POINT_NAMES = {"__init__.py", "api.py", "conftest.py", "setup.py"}


# Extra inputs required by bundled synthetic demos. The South Asian fixture's
# small AIM panel cannot support automatic ancestry inference.
DEMO_EXTRA_ARGS = {
    "ancestry-risk-profiler": ("--ancestry", "SAS"),
}

# Entry scripts with no --demo: the guessed `<script> --demo` would crash
# (clawpathy needs a workspace and the claude CLI) or exit 0 having run nothing
# (the prioritizer is an importable API; its default panel is not committed).
NO_DEMO = {"clawpathy-autoresearch", "clinical-variant-prioritizer"}


def detect_entry_point(skill_dir: Path, folder_name: str) -> Path | None:
    """Return the skill's runnable entry point, whatever language it is in.

    Looks only at the skill root, never inside `tests/`, and prefers a file
    named after the folder so the choice is stable across runs.
    """
    def usable(f: Path, suffix: str) -> bool:
        if f.name in NON_ENTRY_POINT_NAMES or f.name.startswith("test_"):
            return False
        if f.name.endswith("_test" + suffix):
            return False
        parts = f.relative_to(skill_dir).parts
        # Never inside tests/, at any depth, and never build artefacts.
        return "tests" not in parts and "__pycache__" not in parts

    for suffix in ENTRY_POINT_SUFFIXES:
        # Root first, so a top-level script always wins; then nested, because
        # some skills keep theirs in `scripts/` (xena-tcga-gene-query). The old
        # rglob found those, and a root-only search silently demoted them.
        root = sorted(f for f in skill_dir.glob(f"*{suffix}") if usable(f, suffix))
        nested = sorted(
            f for f in skill_dir.rglob(f"*{suffix}")
            if usable(f, suffix) and f.parent != skill_dir
        )
        candidates = root + nested
        if not candidates:
            continue

        stem = folder_name.replace("-", "_")
        preferred = (
            f"{stem}{suffix}",
            f"{folder_name}{suffix}",
            f"cli{suffix}",
            f"__main__{suffix}",
        )
        by_name = {c.name: c for c in candidates}
        for name in preferred:
            if name in by_name:
                return by_name[name]
        return candidates[0]

    return None


STANDARD_RUN_FLAGS = ("--demo", "--input", "--output")


def run_flags(script: Path | None) -> list[str]:
    """The standard flags (--demo, --input, --output) a Python entry point defines with add_argument.

    `clawbio run` reads this to launch a skill it has no SKILLS entry for, and refuses a
    run that needs a flag the script lacks. Read from the source, never by running it; a
    parser built in another module is not seen, so such a skill is refused, not run wrongly.
    """
    if script is None or script.suffix != ".py":
        return []
    try:
        tree = ast.parse(script.read_text(encoding="utf-8"))
    except (OSError, SyntaxError, UnicodeDecodeError, ValueError):
        return []
    found = {
        arg.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
        and node.func.attr == "add_argument"
        for arg in node.args
        if isinstance(arg, ast.Constant) and arg.value in STANDARD_RUN_FLAGS
    }
    return sorted(found)
