# nfcore_sarek_wrapper / scripts/common.py
"""Small dependency-free stand-ins for the 1.x `clawbio.common.*` shared library.

ClawBio 2.0 has no shared `clawbio.common` package (the 2.0 `clawbio` package is
the CLI itself: __main__.py/commands.py/pipeline.py/state.py/bundle.py — see the
repo root CLAUDE.md). The 1.x wrapper imported exactly one helper from it,
`clawbio.common.textio.write_text_lf`, used by params_builder.py to write
params.yaml with normalized LF line endings. Reimplemented here faithfully
(normalize CRLF/CR -> LF, ensure a trailing newline, UTF-8) so params_builder.py
needs no other change.
"""
from __future__ import annotations

from pathlib import Path


def write_text_lf(path: Path, text: str) -> None:
    """Write ``text`` to ``path`` as UTF-8 with normalized LF line endings.

    Creates parent directories as needed. Ensures the file ends with exactly
    one trailing newline, matching the reproducibility-bundle convention used
    throughout the 1.x wrapper (byte-stable macOS <-> Linux replay).
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    normalized = text.replace("\r\n", "\n").replace("\r", "\n")
    if not normalized.endswith("\n"):
        normalized += "\n"
    with path.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write(normalized)
