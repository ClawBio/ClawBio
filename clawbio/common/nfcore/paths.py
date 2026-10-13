"""Path mechanics with no pipeline-specific validation policy."""

from pathlib import Path


def is_remote_uri(value: str) -> bool:
    return "://" in str(value)


def resolve_path_or_uri(value: str, *, base_dir: Path) -> Path | str:
    if is_remote_uri(value):
        return value
    path = Path(value).expanduser()
    return (path if path.is_absolute() else base_dir / path).resolve()


def is_within(path: Path, anchor: Path) -> bool:
    try:
        path.resolve().relative_to(anchor.resolve())
        return True
    except ValueError:
        return False


def is_under_tmp(path: Path) -> bool:
    """Return True when ``path`` resolves to (or inside) /tmp or /private/tmp.

    Single source of truth for the macOS+Docker /tmp guard, shared by preflight
    (pre-run WARNING) and executor (post-failure hint) so the two cannot diverge.
    On macOS, Colima/Docker does not share /tmp into its VM, so a work-dir under
    /tmp surfaces as a confusing '.command.run: No such file or directory'.
    Resolve-based (not a string prefix) so symlinks and the canonical
    /tmp -> /private/tmp mapping are handled. Mirrors the nfcore-scrnaseq helper.
    """
    try:
        resolved = path.resolve()
    except OSError:
        return False
    tmp = Path("/tmp").resolve()
    private_tmp = Path("/private/tmp").resolve()
    return (
        resolved == tmp
        or resolved == private_tmp
        or tmp in resolved.parents
        or private_tmp in resolved.parents
    )
