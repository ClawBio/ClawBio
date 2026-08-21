from __future__ import annotations

import os
from pathlib import Path

import pytest


OPT_IN_FLAG = "--include-rna"
RNA_TEST_DIRS = {
    "skills/celltype-specificity-profiler/tests",
    "skills/diff-visualizer/tests",
    "skills/rnaseq-de/tests",
    "skills/scrna-embedding/tests",
    "skills/scrna-orchestrator/tests",
}


GI_KEY_VAR = "GI_API_KEY"
GI_LIVE_FLAG = "--include-gi-live"
GI_TEST_DIRS = {
    "skills/gi-annotation/tests",
    "skills/gi-chromatin/tests",
    "skills/gi-enhancer/tests",
    "skills/gi-expression/tests",
    "skills/gi-promoter/tests",
    "skills/gi-splice/tests",
}


def _demo_key_from_env_example(root: Path) -> str | None:
    """Read the committed GI demo key out of .env.example.

    Minimal parser: no dotenv dependency, and it only ever looks at the one
    variable. Returns None if the file or the assignment is missing.
    """
    example = root / ".env.example"
    try:
        lines = example.read_text().splitlines()
    except OSError:
        return None
    for line in lines:
        line = line.strip()
        if not line.startswith(f"{GI_KEY_VAR}="):
            continue
        value = line.split("=", 1)[1].split("#", 1)[0].strip().strip("'\"")
        return value or None
    return None


def _default_gi_key(root: Path) -> None:
    """Fall back to the committed GI demo key when GI_API_KEY is unset.

    The gi-* skills call the hosted Genomic Intelligence API, so their
    integration tests need a bearer key. A capped demo key is committed to
    .env.example on purpose -- it attributes traffic to ClawBio and gates
    nothing -- and gi_client's own error message tells users to copy it.

    Only applied under `--include-gi-live`. Applying it unconditionally made a
    plain `pytest` transmit: every gi-* skill carries an `integration` test that
    uploads its bundled FASTA, `pytest.ini` deselects nothing, and CI runs
    `pytest skills/$skill/tests` unfiltered, so each PR touching a gi skill
    billed the shared demo key. Every gi SKILL.md promises "Remote inference,
    opt-in required"; the flag is what makes that true.

    An explicit GI_API_KEY always wins. The integration tests shell out to the
    skill scripts, which inherit this environment.
    """
    if os.environ.get(GI_KEY_VAR):
        return
    key = _demo_key_from_env_example(root)
    if key:
        os.environ[GI_KEY_VAR] = key


def pytest_configure(config) -> None:
    if config.getoption(GI_LIVE_FLAG):
        _default_gi_key(Path(str(config.rootpath)))


def pytest_addoption(parser) -> None:
    parser.addoption(
        OPT_IN_FLAG,
        action="store_true",
        default=False,
        help="Run resource-intensive RNA/scRNA skill tests during default collection.",
    )
    parser.addoption(
        GI_LIVE_FLAG,
        action="store_true",
        default=False,
        help=(
            "Run the gi-* skill integration tests, which upload their bundled "
            "FASTA to the hosted Genomic Intelligence API."
        ),
    )


def _in_dirs(path: Path, root: Path, suites: set[str]) -> bool:
    try:
        relative = path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return False
    return any(relative == suite or relative.startswith(f"{suite}/") for suite in suites)


def _is_rna_suite(path: Path, root: Path) -> bool:
    return _in_dirs(path, root, RNA_TEST_DIRS)


def _explicit_targets(config, root: Path) -> list[Path]:
    targets: list[Path] = []
    for arg in config.invocation_params.args:
        if not arg or arg.startswith("-"):
            continue
        target = Path(arg)
        if not target.is_absolute():
            target = root / target
        targets.append(target.resolve())
    return targets


def _is_explicit_target(path: Path, targets: list[Path]) -> bool:
    resolved = path.resolve()
    for target in targets:
        if resolved == target:
            return True
        if resolved.is_relative_to(target):
            return True
        if target.is_relative_to(resolved):
            return True
    return False


def pytest_ignore_collect(collection_path: Path, config) -> bool:
    root = Path(str(config.rootpath)).resolve()
    candidate = Path(collection_path).resolve()

    if not _is_rna_suite(candidate, root):
        return False
    if config.getoption(OPT_IN_FLAG):
        return False
    if _is_explicit_target(candidate, _explicit_targets(config, root)):
        return False
    return True


def pytest_collection_modifyitems(config, items) -> None:
    """Skip the gi-* live-API tests unless --include-gi-live is passed.

    Deliberately narrower than the marker: `integration` is also used by
    just-prs-mcp and phylogenetics-builder, which do not call our API and are
    not ours to gate. Scoped to GI_TEST_DIRS for that reason.

    Deliberately without the _is_explicit_target escape hatch that the RNA gate
    has. CI names `skills/gi-expression/tests` explicitly, which is precisely
    the invocation that was transmitting, so an explicit-target exemption would
    gate nothing where it matters.
    """
    if config.getoption(GI_LIVE_FLAG):
        return
    root = Path(str(config.rootpath)).resolve()
    skip = pytest.mark.skip(
        reason=f"gi-* live API test; pass {GI_LIVE_FLAG} to run (uploads to the hosted GI API)"
    )
    for item in items:
        if "integration" not in item.keywords:
            continue
        if _in_dirs(Path(str(item.fspath)), root, GI_TEST_DIRS):
            item.add_marker(skip)
