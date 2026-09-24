"""Cross-platform text writing for reproducibility bundles.

Every text artifact in a reproducibility bundle (commands.sh, params.yaml,
provenance JSONs, environment.yml, checksums, configs) must be byte-identical
regardless of the OS that generated it. Python's ``Path.write_text`` /
``open("w")`` use ``newline=None`` by default, which translates ``\\n`` to
``os.linesep`` — emitting CRLF on Windows. CRLF then:

  * breaks ``commands.sh`` (a ``\\r`` leaks into shell variables / line
    continuations under bash),
  * changes every file's sha256 (so ``manifest.json`` checksums differ between
    a macOS/Linux-generated bundle and a Windows-generated one),
  * makes ``sha256sum -c checksums.sha256`` fail.

``write_text_lf`` writes *bytes* after normalising line endings to LF, so no
OS-level newline translation can reintroduce ``\\r``. This is the single
choke-point all bundle writers route through.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path


def _to_lf_bytes(text: str) -> bytes:
    return text.replace("\r\n", "\n").replace("\r", "\n").encode("utf-8")


def write_text_lf(path: Path | str, text: str) -> Path:
    """Write ``text`` to ``path`` with LF line endings on every OS.

    CRLF and lone-CR sequences in ``text`` are normalised to LF, then the
    result is written as UTF-8 bytes (bypassing text-mode newline translation).
    Returns the written path.
    """
    p = Path(path)
    p.write_bytes(_to_lf_bytes(text))
    return p


def _umask() -> int:
    """Read the process umask. os has no query, so set-and-restore is the only way.

    ponytail: momentarily racy for a threaded writer; use os.umask() directly at
    startup if that ever matters.
    """
    current = os.umask(0o022)
    os.umask(current)
    return current


def write_text_lf_atomic(path: Path | str, text: str) -> Path:
    """Like :func:`write_text_lf`, but replaces ``path`` atomically.

    The content is written to a temporary file in the same directory and then
    ``os.replace``\\ d into place, so a reader that already has ``path`` open keeps
    reading the original (now-unlinked) inode intact instead of seeing a truncated
    file. This matters for ``commands.sh``: an in-place ``--resume`` replay re-invokes
    the wrapper, which regenerates the very script bash is executing — a plain
    truncate-and-rewrite corrupts bash's read mid-run. The existing file's permission
    bits are preserved; a new file lands on the same umask-derived mode a plain write
    would have given it, rather than on mkstemp's 0600. Returns the written path.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    mode = p.stat().st_mode if p.exists() else 0o666 & ~_umask()
    fd, tmp_name = tempfile.mkstemp(dir=str(p.parent), prefix=f".{p.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(_to_lf_bytes(text))
        os.chmod(tmp_name, mode)
        os.replace(tmp_name, p)
    except BaseException:
        try:
            os.unlink(tmp_name)
        except OSError:
            pass
        raise
    return p
