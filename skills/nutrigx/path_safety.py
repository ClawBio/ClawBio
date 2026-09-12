"""
path_safety.py — symlink-safe file creation for NutriGx Advisor.

Output filenames are deterministic (nutrigx_report.md, nutrigx_radar.png,
nutrigx_heatmap.png), so anyone able to write into the output directory ahead of
a run can pre-create one of those names as a symbolic link and have the skill
write a genetic report through it, to a location the user never chose.

Python's ordinary write paths follow symlinks: Path.write_text, open(path, "w")
and matplotlib's savefig all do. These helpers do not.

Two components have to be guarded, not one:

  * the final component, via O_NOFOLLOW, so the file itself cannot be a symlink
  * the parent directory, opened with O_DIRECTORY | O_NOFOLLOW and then used as
    a dir_fd, so the write is anchored to a real directory inode and swapping
    the parent for a symlink afterwards cannot move it

O_EXCL is deliberately not used: overwriting a regular file on a re-run is
expected behaviour.

These helpers are deliberately policy-free — they do not decide *where* output
may go, only that the path they are handed is not a symlink. That keeps them
safe to use from api.py, which accepts an arbitrary output_dir from its caller.
"""

import errno
import os
from pathlib import Path


def safe_open_write(path, binary: bool = False):
    """Open ``path`` for writing without following symlinks. Raises ValueError
    if the file or its parent directory is a symbolic link."""
    path = Path(path)
    try:
        dir_fd = os.open(str(path.parent), os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    except OSError as exc:
        if exc.errno in (errno.ELOOP, errno.ENOTDIR):
            raise ValueError(
                f"Refusing to write into '{path.parent}': it is a symbolic link, "
                f"not a real directory."
            ) from exc
        raise

    try:
        fd = os.open(
            path.name,
            os.O_WRONLY | os.O_CREAT | os.O_TRUNC | os.O_NOFOLLOW,
            0o600,
            dir_fd=dir_fd,
        )
    except OSError as exc:
        if exc.errno in (errno.ELOOP, errno.EMLINK):
            raise ValueError(
                f"Refusing to write to '{path}': it is a symbolic link. "
                f"Remove it and re-run."
            ) from exc
        raise
    finally:
        os.close(dir_fd)

    return os.fdopen(fd, "wb" if binary else "w", encoding=None if binary else "utf-8")


def safe_write_text(path, text: str) -> None:
    """Write text to ``path`` without following a symlink."""
    with safe_open_write(path) as fh:
        fh.write(text)
