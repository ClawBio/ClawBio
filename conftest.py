from __future__ import annotations

from pathlib import Path


OPT_IN_FLAG = "--include-rna"
RNA_TEST_DIRS = {
    "skills/diff-visualizer/tests",
    "skills/rnaseq-de/tests",
    "skills/scrna-embedding/tests",
    "skills/scrna-orchestrator/tests",
}


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
