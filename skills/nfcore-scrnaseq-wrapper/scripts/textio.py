"""Cross-platform text writing for reproducibility artifacts.

Inlined, dependency-free reimplementation of ClawBio 1.x's
``clawbio.common.textio.write_text_lf`` (the shared-library package does not
exist in ClawBio 2.0). Every text artifact this skill writes (params.yaml,
JSON reports, the samplesheet copy, the composed Nextflow command) should be
byte-identical regardless of the OS that generated it: Python's default
``open("w")`` / ``Path.write_text`` translate ``\\n`` to ``os.linesep``,
which would emit CRLF on Windows and change every output file's sha256
depending on platform. ``write_text_lf`` writes bytes after normalizing line
endings to LF, so no OS-level newline translation can reintroduce ``\\r``.
"""
from __future__ import annotations

from pathlib import Path


def _to_lf_bytes(text: str) -> bytes:
    return text.replace("\r\n", "\n").replace("\r", "\n").encode("utf-8")


def write_text_lf(path: Path | str, text: str) -> Path:
    """Write ``text`` to ``path`` with LF line endings on every OS."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(_to_lf_bytes(text))
    return p
