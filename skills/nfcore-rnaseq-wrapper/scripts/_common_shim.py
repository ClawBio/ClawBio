"""Dependency-free local reimplementation of the two ClawBio 1.x
``clawbio.common`` helpers the ported nf-core/rnaseq wrapper modules need.

1.x's ``preflight.py`` and ``params_builder.py`` import
``clawbio.common.checksums.sha256_file`` and ``clawbio.common.textio.write_text_lf``
respectively. That package (``clawbio/common/``) is part of the 1.x skill-library
runtime and is not vendored into ClawBio 2.0's `skills/` tree, so importing it
here would reach outside this skill's directory. Both functions are tiny and
have no further dependencies, so they are reimplemented verbatim (same
behaviour, same signatures) rather than routing back into 1.x. See
1.x's ``clawbio/common/checksums.py`` and ``clawbio/common/textio.py``.
"""

from __future__ import annotations

import hashlib
from pathlib import Path


def sha256_file(filepath: str | Path) -> str:
    """Compute SHA-256 checksum of a file (full 64-char hex digest)."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def _to_lf_bytes(text: str) -> bytes:
    return text.replace("\r\n", "\n").replace("\r", "\n").encode("utf-8")


def write_text_lf(path: Path | str, text: str) -> Path:
    """Write ``text`` to ``path`` with LF line endings on every OS.

    Byte-identical output regardless of platform, so checksums in the
    reproducibility bundle do not drift between a macOS/Linux and a
    Windows-generated run.
    """
    p = Path(path)
    p.write_bytes(_to_lf_bytes(text))
    return p
