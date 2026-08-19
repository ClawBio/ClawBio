from __future__ import annotations

import os
from pathlib import Path


OPT_IN_FLAG = "--include-rna"
RNA_TEST_DIRS = {
    "skills/celltype-specificity-profiler/tests",
    "skills/diff-visualizer/tests",
    "skills/rnaseq-de/tests",
    "skills/scrna-embedding/tests",
    "skills/scrna-orchestrator/tests",
}


GI_KEY_VAR = "GI_API_KEY"


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
    nothing -- and gi_client's own error message tells users to copy it. This
    applies that default automatically so the suites also run in CI, where a
    pull request from a fork cannot read repository secrets.

    An explicit GI_API_KEY always wins. The integration tests shell out to the
    skill scripts, which inherit this environment.
    """
    if os.environ.get(GI_KEY_VAR):
        return
    key = _demo_key_from_env_example(root)
    if key:
        os.environ[GI_KEY_VAR] = key


def pytest_configure(config) -> None:
    _default_gi_key(Path(str(config.rootpath)))


def pytest_addoption(parser) -> None:
    parser.addoption(
        OPT_IN_FLAG,
        action="store_true",
        default=False,
        help="Run resource-intensive RNA/scRNA skill tests during default collection.",
    )


def _is_rna_suite(path: Path, root: Path) -> bool:
    try:
        relative = path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return False
    return any(relative == suite or relative.startswith(f"{suite}/") for suite in RNA_TEST_DIRS)


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
