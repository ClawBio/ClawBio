#!/usr/bin/env python3
"""seqtk-downsample: randomly subsample FASTQ files with `seqtk sample`.

One task only: pick N reads (or a fraction) from one FASTQ, or from a
read-ordered pair, with a fixed seed. Mates are kept in sync by giving both
files the same seed, then verified by comparing read names.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (  # noqa: E402
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)

SKILL_DIR = Path(__file__).resolve().parent
DEMO_R1 = SKILL_DIR / "examples" / "synthetic_demo_R1.fastq.gz"
DEMO_R2 = SKILL_DIR / "examples" / "synthetic_demo_R2.fastq.gz"

DEFAULT_SEED = 12345
DEFAULT_N_READS = 10000
DEMO_N_READS = 100

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)


class DownsampleError(Exception):
    """A user-facing failure (bad input, missing seqtk, mismatched mates)."""


# ──────────────────────────────────────────────
# seqtk + FASTQ helpers
# ──────────────────────────────────────────────

def find_seqtk() -> str:
    """Return the seqtk executable: $SEQTK_BIN if set, otherwise PATH."""
    explicit = os.environ.get("SEQTK_BIN")
    if explicit:
        if not (Path(explicit).is_file() and os.access(explicit, os.X_OK)):
            raise DownsampleError(
                f"seqtk not found at SEQTK_BIN={explicit}. "
                "Install it with: conda install -c bioconda seqtk"
            )
        return explicit
    found = shutil.which("seqtk")
    if not found:
        raise DownsampleError(
            "seqtk not found on PATH. Install it with: "
            "conda install -c bioconda seqtk (or set SEQTK_BIN)."
        )
    return found


def seqtk_version(seqtk: str) -> str:
    """seqtk prints 'Version: 1.5-r133' in its usage text (exit code 1)."""
    proc = subprocess.run([seqtk], capture_output=True, text=True)
    m = re.search(r"Version:\s*(\S+)", proc.stderr + proc.stdout)
    return m.group(1) if m else "unknown"


def normalise_read_name(header: str) -> str:
    """Key shared by both mates: drop '@', any comment, and a trailing /1 or /2."""
    first = header.lstrip("@").split()[0] if header.strip("@ \t") else ""
    return re.sub(r"/[12]$", "", first)


def _open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_keys(path: Path):
    """Yield the normalised read name of every record (header is line 1 of 4)."""
    with _open_text(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                yield normalise_read_name(line.rstrip("\n"))


def count_reads(path: Path) -> int:
    return sum(1 for _ in read_keys(path))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def output_name(src: Path) -> str:
    stem = re.sub(r"\.(gz)$", "", src.name)
    stem = re.sub(r"\.(fastq|fq)$", "", stem)
    return f"{stem}.downsampled.fastq.gz"


def _sample_to_gzip(seqtk: str, src: Path, amount: str, seed: int, dest: Path):
    """Run `seqtk sample -s SEED src AMOUNT`, gzip-compress stdout to dest."""
    proc = subprocess.Popen(
        [seqtk, "sample", "-s", str(seed), str(src), amount],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    # mtime=0 and no embedded filename: same reads -> same bytes, so the
    # checksums in reproducibility/ hold on a re-run (cf. `gzip --no-name`).
    with open(dest, "wb") as raw, gzip.GzipFile(
        filename="", mode="wb", fileobj=raw, mtime=0
    ) as out:
        shutil.copyfileobj(proc.stdout, out)
    stderr = proc.stderr.read().decode(errors="replace")
    if proc.wait() != 0:
        dest.unlink(missing_ok=True)
        raise DownsampleError(f"seqtk sample failed on {src}: {stderr.strip()}")


# ──────────────────────────────────────────────
# Core
# ──────────────────────────────────────────────

def validate_amount(n_reads: int | None, fraction: float | None) -> tuple[str, float | int]:
    """Return (mode, value). Exactly one of n_reads / fraction, else default."""
    if n_reads is not None and fraction is not None:
        raise DownsampleError("use --n-reads or --fraction, not both")
    if fraction is not None:
        # seqtk reads any amount >= 1 as a read COUNT, so 1.0 would silently
        # return one read. NaN/inf/<=0 are nonsense. Hence the open interval.
        if not math.isfinite(fraction) or not (0.0 < fraction < 1.0):
            raise DownsampleError(
                "--fraction must be strictly between 0 and 1 "
                "(seqtk treats 1 or more as a read count)"
            )
        return "fraction", fraction
    n = DEFAULT_N_READS if n_reads is None else n_reads
    if n < 1:
        raise DownsampleError("--n-reads must be a positive integer")
    return "n_reads", n


def downsample(
    inputs: list[Path],
    output_dir: Path,
    n_reads: int | None = None,
    fraction: float | None = None,
    seed: int = DEFAULT_SEED,
) -> dict:
    """Downsample one FASTQ or a pair. Returns the result dict (also written)."""
    if len(inputs) not in (1, 2):
        raise DownsampleError("provide one FASTQ, or two for a read pair")
    mode, value = validate_amount(n_reads, fraction)
    for src in inputs:
        if not src.is_file():
            raise DownsampleError(f"input not found: {src}")
    names = [output_name(src) for src in inputs]
    if len(set(names)) != len(names):
        raise DownsampleError("both inputs have the same file name; rename one")

    seqtk = find_seqtk()
    version = seqtk_version(seqtk)

    out_dir = output_dir / "downsampled"
    targets = [out_dir / name for name in names]
    existing = [p for p in [*targets, output_dir / "result.json"] if p.exists()]
    if existing:
        raise DownsampleError(
            "refusing to overwrite existing output: "
            + ", ".join(str(p) for p in existing)
            + ". Choose a new --output directory."
        )
    out_dir.mkdir(parents=True, exist_ok=True)

    amount = str(value)
    try:
        for src, dest in zip(inputs, targets):
            _sample_to_gzip(seqtk, src, amount, seed, dest)
        if len(inputs) == 2:
            _check_mates_in_sync(targets[0], targets[1])
    except DownsampleError:
        for dest in targets:
            dest.unlink(missing_ok=True)
        raise

    files = []
    for src, dest in zip(inputs, targets):
        files.append({
            "input": str(src),
            "output": dest.relative_to(output_dir).as_posix(),
            "reads_out": count_reads(dest),
            "sha256": sha256_file(dest),
        })
    reads_out = files[0]["reads_out"]
    result = {
        "skill": "seqtk-downsample",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "seqtk_version": version,
        "mode": mode,
        "requested": value,
        "seed": seed,
        "paired": len(inputs) == 2,
        "fewer_than_requested": mode == "n_reads" and reads_out < value,
        "files": files,
        "disclaimer": DISCLAIMER,
    }
    return result


def _check_mates_in_sync(r1: Path, r2: Path) -> None:
    for i, (a, b) in enumerate(zip(read_keys(r1), read_keys(r2), strict=False)):
        if a != b:
            raise DownsampleError(
                f"read order mismatch at output read {i + 1}: {a!r} vs {b!r}. "
                "The two files are not in the same read order, so a shared "
                "seed cannot keep the mates paired."
            )
    if count_reads(r1) != count_reads(r2):
        raise DownsampleError(
            "read order mismatch: the downsampled mates have different read counts"
        )


# ──────────────────────────────────────────────
# Report + reproducibility
# ──────────────────────────────────────────────

def write_report(output_dir: Path, result: dict) -> Path:
    amount = (
        f"{result['requested']} reads"
        if result["mode"] == "n_reads"
        else f"fraction {result['requested']}"
    )
    rows = "\n".join(
        f"| `{Path(f['input']).name}` | `{f['output']}` | {f['reads_out']} |"
        for f in result["files"]
    )
    note = ""
    if result["fewer_than_requested"]:
        note = (
            f"\n> **Note:** the input held fewer reads than the {result['requested']} "
            f"requested, so all {result['files'][0]['reads_out']} reads were kept.\n"
        )
    text = f"""# seqtk-downsample Report

**Date**: {result['generated_at']}
**Tool**: seqtk {result['seqtk_version']} (`seqtk sample`)
**Amount**: {amount}
**Seed**: {result['seed']}
**Layout**: {'paired-end (mates checked in sync)' if result['paired'] else 'single file'}

| Input | Output | Reads kept |
|-------|--------|-----------:|
{rows}
{note}
Reads are a uniform random sample without replacement; the same seed and the
same input reproduce the same reads.

*{DISCLAIMER}*
"""
    path = output_dir / "report.md"
    path.write_text(text)
    return path


def write_bundle(output_dir: Path, result: dict, argv_inputs: list[Path]) -> None:
    flag = (
        f"--n-reads {result['requested']}"
        if result["mode"] == "n_reads"
        else f"--fraction {result['requested']}"
    )
    in_flags = f"--input {shlex.quote(str(argv_inputs[0]))}"
    if len(argv_inputs) == 2:
        in_flags += f" --input2 {shlex.quote(str(argv_inputs[1]))}"
    lines = [
        "# ClawBio",
        f"python skills/seqtk-downsample/seqtk_downsample.py {in_flags} "
        f"--output <new_output_dir> {flag} --seed {result['seed']}",
        "# Equivalent raw seqtk (per input file):",
    ]
    amount = result["requested"]
    for f in result["files"]:
        lines.append(
            f"# seqtk sample -s {result['seed']} {shlex.quote(f['input'])} "
            f"{amount} | gzip > {Path(f['output']).name}"
        )
    write_commands_sh(output_dir, "\n".join(lines))
    major_minor = result["seqtk_version"].split("-")[0]
    write_environment_yml(
        output_dir,
        "clawbio-seqtk-downsample",
        pip_deps=[],
        conda_deps=[f"seqtk={major_minor}"],
        channels=["conda-forge", "bioconda", "nodefaults"],
    )
    write_checksums(
        [output_dir / f["output"] for f in result["files"]]
        + [output_dir / "report.md", output_dir / "result.json"],
        output_dir,
        anchor=output_dir,
    )


# ──────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Downsample FASTQ reads with seqtk sample (fixed seed, "
        "paired-end aware).")
    p.add_argument("--input", type=Path, help="FASTQ (.fastq/.fq, optionally .gz); read 1 if paired")
    p.add_argument("--input2", type=Path, help="Read 2 FASTQ, in the same read order as --input")
    p.add_argument("--output", type=Path, help="Output directory (must not already hold results)")
    p.add_argument("--demo", action="store_true", help="Run on bundled synthetic paired FASTQs")
    p.add_argument("--n-reads", type=int, default=None,
                   help=f"Number of reads to keep (default {DEFAULT_N_READS})")
    p.add_argument("--fraction", type=float, default=None,
                   help="Fraction of reads to keep, strictly between 0 and 1")
    p.add_argument("--seed", type=int, default=DEFAULT_SEED,
                   help=f"Random seed (default {DEFAULT_SEED})")
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.output is None:
        parser.error("--output is required")
    if args.demo:
        inputs = [DEMO_R1, DEMO_R2]
        if args.n_reads is None and args.fraction is None:
            args.n_reads = DEMO_N_READS
    else:
        if args.input is None:
            parser.error("--input is required (or use --demo)")
        inputs = [args.input] + ([args.input2] if args.input2 else [])
    if args.n_reads is not None and args.fraction is not None:
        parser.error("use --n-reads or --fraction, not both")

    try:
        result = downsample(inputs, args.output, args.n_reads, args.fraction, args.seed)
        (args.output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
        write_report(args.output, result)
        write_bundle(args.output, result, inputs)
    except DownsampleError as exc:
        print(f"seqtk-downsample: error: {exc}", file=sys.stderr)
        return 1

    kept = result["files"][0]["reads_out"]
    print(f"Kept {kept} reads per file -> {args.output / 'downsampled'}")
    print(f"Report: {args.output / 'report.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
