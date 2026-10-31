"""Keep genomic data and skill results out of the language model's context.

Two parts, stdlib only:

- :func:`seal_result` reduces a skill result to what a host needs to deliver the
  files to the user. ``clawbio.py run --sealed`` prints it.
- The genomic-data policy: which files, commands and texts carry genomic data,
  and what an agent may do with them. Agent adapters (integrations/hermes/) map
  their own tools and messages onto :func:`check_path`, :func:`check_command`,
  :func:`replace_uploads` and :func:`scrub_text`.

``dirs`` is the list of ClawBio checkouts (CLAWBIO_DIR) whose ``output/`` folder
holds reports; the adapter decides where it comes from.
"""

from __future__ import annotations

import glob
import gzip
import os
import re
import shlex
import zipfile
from pathlib import Path
from typing import Any, Mapping

# This is an allowlist, so a field a skill adds later stays hidden from the model
# until someone adds it here.
SEALED_FIELDS = ("skill", "success", "exit_code", "duration_seconds", "output_dir", "files")

SEALED_NOTE = "Do not read, open or summarise the output files. The host delivers report.md to the user unopened."


def seal_result(result: Mapping[str, Any]) -> dict:
    """Return the delivery metadata of a skill result, without its contents.

    Drops stdout, stderr, report text and every field not in SEALED_FIELDS.
    ``success`` is True only when the result says exactly ``True``.
    """
    sealed = {key: result[key] for key in SEALED_FIELDS if key in result}
    sealed["success"] = result.get("success") is True
    # A fixed hint, never stderr: a failed run's stderr can carry the input's contents.
    if str(result.get("stderr", "")).startswith("Unknown skill"):
        sealed["error"] = "Unknown skill name. `clawbio.py list` shows the valid names."
    sealed["sealed"] = True
    sealed["note"] = SEALED_NOTE
    return sealed


# --- The genomic-data policy ------------------------------------------------- #

# The genomic entries of bio-orchestrator's EXTENSION_MAP. 23andMe and AncestryDNA exports are
# .txt, caught by content.
GENOMIC_EXTENSIONS = (".vcf", ".vcf.gz", ".fastq", ".fastq.gz", ".fq", ".fq.gz", ".bam", ".cram")

WITHHELD = "[withheld by clawbio-privacy: genomic data is not shown to the model]"

_GENOMIC_EXT_RE = re.compile(r"\S+?(" + "|".join(re.escape(e) for e in GENOMIC_EXTENSIONS) + r")\b", re.I)
_RSID_RE = re.compile(r"\brs\d{3,}\b")
_SEALED_RUN_RE = re.compile(r"^\s*(uv run\s+)?\S*python3?\s+\S*clawbio\.py\s+run\s+[^;&|<>`$\n]*--sealed[^;&|<>`$\n]*$")
_PATHISH_RE = re.compile(r"(?:~|/)[^\s'\";|&<>()]+")
# A VCF header anywhere (BCF puts binary before it); raw-data and SAM headers at a line start, after
# read_file's "12|" line numbers if present.
_GENOMIC_HEADER_RE = re.compile(r"##fileformat=VCF|^(?:\d+\|)?(?:# rsid|rsid\t|@(?:HD|SQ)\t)", re.M | re.I)
_SAVED_AT_RE = re.compile(r"saved at: ([^\]\n]+)\]")
# A relative output/ path, or a ClawBio run folder (lowercase <skill>_YYYYMMDD_HHMMSS, unlike Android's
# IMG_YYYYMMDD_HHMMSS) reached after `cd output`.
# shortcut: blocks every relative output/, not only ClawBio's; the terminal's cwd is not visible here. Upgrade if the host exposes it.
_RELATIVE_OUTPUT_RE = re.compile(r"(?<![\w/~.-])(?:\./)?output/|\b[a-z][a-z0-9_-]*_\d{8}_\d{6}\b")


def is_genomic_text(text: str) -> bool:
    """VCF or raw-data header on any line, or three or more distinct rsIDs."""
    text = text.lstrip("﻿")  # Windows/Excel exports start with a byte-order mark
    return bool(_GENOMIC_HEADER_RE.search(text)) or len(set(_RSID_RE.findall(text))) >= 3


def is_genomic_path(path: str) -> bool:
    """Genomic by extension, or by sniffing the first 4 KB of an existing file."""
    if path.lower().endswith(GENOMIC_EXTENSIONS):
        return True
    try:
        p = Path(path).expanduser()
        if p.is_file():
            return is_genomic_text(_head(p))
    except (OSError, RuntimeError, ValueError, EOFError, zipfile.BadZipFile):  # "~name" with no such user; a NUL in the path; a broken archive
        pass
    return False


def _head(p: Path) -> str:
    """The first 4 KB as text, read through gzip/BGZF, or of every member of a zip."""
    if zipfile.is_zipfile(p):
        with zipfile.ZipFile(p) as z:
            return "\n".join(z.open(n).read(4096).decode(errors="ignore") for n in z.namelist())
    with p.open("rb") as fh:
        data = fh.read(4096)
    if data[:2] == b"\x1f\x8b":
        with gzip.open(p) as fh:
            data = fh.read(4096)
    return data.decode(errors="ignore")


def _paths(command: str) -> list[str]:
    """Absolute paths in a command, with globs expanded."""
    return [x for p in _PATHISH_RE.findall(command) for x in glob.glob(os.path.expanduser(p)) or [p]]


def is_sealed_run(command: str) -> bool:
    """A single `clawbio.py run` command, no shell chaining, with --sealed as its own argument."""
    if not _SEALED_RUN_RE.match(command):
        return False
    try:
        args = shlex.split(command)
    except ValueError:  # unbalanced quotes
        return False
    # argparse takes any prefix of --output, and a report written elsewhere escapes the output/ block.
    return "--sealed" in args and not any(a.startswith("--o") and "--output".startswith(a.split("=")[0]) for a in args)


def is_under_output(path: str, dirs: list[str]) -> bool:
    """A path in a ClawBio output/ folder: under one of ``dirs``, or relative to a cwd this process cannot see."""
    if _RELATIVE_OUTPUT_RE.search(path):
        return True
    try:
        target = Path(path).expanduser().resolve()
    except (OSError, RuntimeError, ValueError):
        return False
    return any(target.is_relative_to(Path(d).expanduser().resolve() / "output") for d in dirs)


def _command_touches_output(command: str, dirs: list[str]) -> bool:
    # ponytail: string matching on shell text; `find /` or an encoded path gets past it, a sandbox does not
    return bool(_RELATIVE_OUTPUT_RE.search(command)) or any(is_under_output(p, dirs) for p in _paths(command))


def _clawbio(dirs: list[str]) -> str:
    """The first of ``dirs`` is a ClawBio checkout with its own .venv."""
    d = next(iter(dirs), None)
    return f"{d}/.venv/bin/python {d}/clawbio.py" if d else "python clawbio.py"


def _sealed_command(path: str, dirs: list[str]) -> str:
    return f"{_clawbio(dirs)} run <skill> --input {path} --sealed"


def _file_kind(content: str) -> str:
    if "##fileformat=vcf" in content[:4096].lower():
        return "VCF"
    return "raw genotype (rsID)" if _RSID_RE.search(content) else "genomic"


def _sealed_note(path: str, dirs: list[str], kind: str) -> str:
    return (
        f"[{kind} file saved at {path}; contents withheld by clawbio-privacy. Do not read it. "
        f"Pick the skill that takes a {kind} file from `{_clawbio(dirs)} list`, then analyse it with: "
        f"{_sealed_command(path, dirs)}. The report is attached automatically to your reply; do not open it]"
    )


def _refusal(dirs: list[str], path: str = "<file>") -> str:
    return (f"{WITHHELD} Genomic files and ClawBio output (reports) are off limits. "
            f"The only allowed command is `{_sealed_command(path, dirs)}`; its report reaches the user directly.")


def check_path(path: str, dirs: list[str]) -> str | None:
    """Why the model may not read or send this file, or None if it may."""
    if is_under_output(path, dirs):
        return _refusal(dirs)
    if is_genomic_path(path):
        return _refusal(dirs, path)
    return None


def check_command(command: str, dirs: list[str]) -> str | None:
    """Why the model may not run this shell command, or None if it may. A sealed run is the one exception."""
    touches = (_GENOMIC_EXT_RE.search(command) or any(is_genomic_path(p) for p in _paths(command))
               or _command_touches_output(command, dirs))
    return _refusal(dirs) if touches and not is_sealed_run(command) else None


def replace_uploads(text: str, paths: list[str], dirs: list[str]) -> str:
    """Swap each genomic upload's pasted contents for its path and the sealed instruction."""
    for path in paths:
        try:
            content = Path(path).read_text(errors="ignore")
        except OSError:
            continue
        if content and content in text and (path.lower().endswith(GENOMIC_EXTENSIONS) or is_genomic_text(content)):
            text = text.replace(content, _sealed_note(path, dirs, _file_kind(content)))
    return text


def scrub_text(text: str, dirs: list[str]) -> str:
    """Text safe to send to the model: saved uploads swapped for their sealed note, anything else genomic withheld."""
    if not is_genomic_text(text):
        return text
    text = replace_uploads(text, [p.strip() for p in _SAVED_AT_RE.findall(text)], dirs)
    return WITHHELD if is_genomic_text(text) else text
