"""Shared security helpers for the RoboTerri bots.

Stdlib-only (no flask / openai / telegram imports) so it is unit-testable in
isolation and importable from every bot adapter. This module encodes the
identity-isolation invariant: every action is attributable to exactly one
authenticated user, with no global / first / any fallback.
"""

from __future__ import annotations

import hashlib
import hmac
import re
import tempfile
from pathlib import Path
from typing import Any, Mapping, Optional, Set


def verify_whatsapp_signature(
    raw_body: bytes,
    signature_header: Optional[str],
    app_secret: str,
) -> bool:
    """Verify a Meta WhatsApp webhook payload signature (``X-Hub-Signature-256``).

    Returns True only if ``app_secret`` is set, ``signature_header`` is a
    well-formed ``sha256=<hex>``, and the HMAC-SHA256 of ``raw_body`` under
    ``app_secret`` matches in constant time. Fails closed (returns False) on any
    missing or malformed input, so an unconfigured secret never silently accepts.
    """
    if not app_secret or not signature_header:
        return False
    if not signature_header.startswith("sha256="):
        return False
    sent = signature_header[len("sha256="):].strip()
    if not sent:
        return False
    expected = hmac.new(app_secret.encode("utf-8"), raw_body, hashlib.sha256).hexdigest()
    return hmac.compare_digest(sent, expected)


def scoped_get(store: Mapping[str, Any], user_id: Optional[str]) -> Optional[Any]:
    """Return the entry for exactly this user, or None. No first / any fallback.

    This is the identity-isolation invariant in code: a user's data is reachable
    only by their own id. Never iterate the store and return an arbitrary entry
    (the ``for _uid, info in store.items(): break`` / ``next(iter(store))``
    patterns the audit flagged are forbidden in user-facing execution paths).
    """
    if not user_id:
        return None
    return store.get(user_id)


def is_sender_allowed(
    user_id: Optional[str],
    allowed: Set[str],
    admin: Optional[str] = None,
    allow_all: bool = False,
) -> bool:
    """True iff the sender may use the bot. Denies by default (fail closed).

    A sender is allowed when they are the configured ``admin`` or appear in the
    explicit ``allowed`` set. Public access (``allow_all``) must be an explicit
    operator choice, never the silent default. An empty / None ``user_id`` is
    always rejected.
    """
    if not user_id:
        return False
    if allow_all:
        return True
    if admin and user_id == admin:
        return True
    return user_id in allowed


# --------------------------------------------------------------------------- #
# Filesystem containment for the chat-facing save_file / write_file tools.
#
# Every adapter routes its model-supplied destination_folder / filename through
# safe_write_path, so the three adapters cannot drift apart.
# --------------------------------------------------------------------------- #


class UnsafePath(ValueError):
    """A model-supplied destination was refused. The message is user-facing."""


# Files the write_file and save_file tools must never overwrite.
# Checked case-insensitively - all entries must be lowercase.
PROTECTED_NAMES = frozenset({
    "soul.md", "claude.md", "agents.md", ".env",
    "roboterri.py", "roboterri_discord.py", "roboterri_whatsapp.py",
    "clawbio.py", "requirements.txt", "contributing.md",
})

ALLOWED_UPLOAD_EXTENSIONS = frozenset({
    ".txt", ".csv", ".vcf", ".fastq", ".fq",   # genetic data (uncompressed)
    ".h5ad",                                     # single-cell AnnData
    ".tif", ".tiff", ".png", ".jpg", ".jpeg", ".heic", ".heif",  # microscopy / photos
    ".tsv",                                      # tab-separated counts
    # .pdf, .html, .md excluded - active content risk / prompt injection
})

# Compound suffixes allowed for gzip-compressed files (e.g. "data.vcf.gz").
# Bare ".gz" is intentionally excluded - it could wrap arbitrary content.
ALLOWED_GZ_STEMS = frozenset({
    ".vcf.gz", ".fastq.gz", ".fq.gz", ".txt.gz", ".tsv.gz", ".csv.gz", ".bed.gz",
})


def is_allowed_extension(filename: str) -> bool:
    """True if the file's extension (or compound .*.gz suffix) is permitted."""
    suffixes = Path(filename).suffixes
    if not suffixes:
        return False
    if "".join(suffixes[-2:]).lower() in ALLOWED_GZ_STEMS:
        return True
    return suffixes[-1].lower() in ALLOWED_UPLOAD_EXTENSIONS


def sanitize_filename(filename: str) -> str:
    """Reduce to a bare basename with no traversal or control characters."""
    filename = Path(filename).name.strip()
    filename = re.sub(r"[\x00-\x1f]", "", filename)
    filename = filename.replace("..", "").replace("/", "").replace("\\", "")
    return filename or "unnamed_file"


def upload_tmp_path(owner_id: Any, filename: str) -> Path:
    """Temp path for an upload, namespaced by the sender or channel that sent it."""
    owner = re.sub(r"[^A-Za-z0-9_-]", "", str(owner_id)) or "unknown"
    return Path(tempfile.gettempdir()) / f"roboterri_{owner}_{sanitize_filename(filename)}"


def safe_write_path(
    folder: str | None,
    filename: str,
    *,
    root: Path,
    require_allowed_extension: bool = False,
) -> Path:
    """Resolve a model-supplied folder + filename to a path safe to write.

    Containment is to ``root``, the user data directory: a destination that
    escapes it falls back to ``root`` rather than being honoured.

    Raises UnsafePath for a protected filename, and for a disallowed extension
    when ``require_allowed_extension`` is set (the save_file path, where the
    model may rename an upload).
    """
    root = Path(root).resolve()

    dest = Path(folder) if folder else root
    if not dest.is_absolute():
        dest = root / dest
    if not dest.resolve().is_relative_to(root):
        dest = root

    filename = sanitize_filename(filename)
    if filename.lower() in PROTECTED_NAMES:
        raise UnsafePath(f"'{filename}' is a protected file and cannot be written.")
    if require_allowed_extension and not is_allowed_extension(filename):
        raise UnsafePath(f"'{filename}' is not an allowed file type.")

    final = dest / filename
    if not final.resolve().is_relative_to(root):
        raise UnsafePath(f"'{filename}' would escape the destination directory.")

    dest.mkdir(parents=True, exist_ok=True)
    return final
