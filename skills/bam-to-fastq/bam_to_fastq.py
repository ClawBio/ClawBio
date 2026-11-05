#!/usr/bin/env python3
"""bam-to-fastq: convert a BAM back to FASTQ.

One task only. Paired-end uses `bedtools bamtofastq`; single-end uses
`samtools fastq` (bedtools wrote single-end reads twice in testing). Every
run counts primary records before and after each step and fails on any
mismatch.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)

SKILL_DIR = Path(__file__).resolve().parent
DEMO_BAM = SKILL_DIR / "examples" / "synthetic_demo.bam"
DEFAULT_CHECK_N = 10000

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)


class BamToFastqError(Exception):
    """A user-facing failure (bad input, missing tool, failed check)."""


# ──────────────────────────────────────────────
# Tool helpers
# ──────────────────────────────────────────────

def require_tools(paired: bool) -> dict[str, str]:
    needed = ["samtools"] + (["bedtools"] if paired else [])
    paths = {}
    for tool in needed:
        found = shutil.which(tool)
        if not found:
            raise BamToFastqError(
                f"{tool} not found on PATH. Install it with: "
                "conda install -c bioconda bedtools samtools"
            )
        paths[tool] = found
    return paths


def tool_version(tool: str) -> str:
    proc = subprocess.run([tool, "--version"], capture_output=True, text=True, check=False)
    m = re.search(r"(\d+\.\d+(?:\.\d+)*)", proc.stdout + proc.stderr)
    return m.group(1) if m else "unknown"


def _run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False, **kw)
    if proc.returncode != 0:
        raise BamToFastqError(f"{cmd[0]} failed: {proc.stderr.strip()}")
    return proc


# ──────────────────────────────────────────────
# Checks
# ──────────────────────────────────────────────

def count_bam(bam: Path) -> int:
    """Primary records (excludes secondary/supplementary)."""
    return int(_run(["samtools", "view", "-c", "-F", "0x900", str(bam)]).stdout)


def _open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def fastq_names(path: Path):
    """Yield each record's read name, without '@', comment, or /1 /2."""
    with _open_text(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                first = line.split()[0].lstrip("@") if line.strip("@ \t\n") else ""
                yield re.sub(r"/[12]$", "", first)


def count_fastq(path: Path) -> int:
    return sum(1 for _ in fastq_names(path))


def check_name_grouped(bam: Path, n: int = DEFAULT_CHECK_N) -> tuple[bool, str, int, str | None]:
    """Sample the first n records; check each read name is one contiguous block.

    Mates must be adjacent for paired mode. Returns
    (ok, header_sort_order, reads_sampled, first_offending_name).
    """
    hdr = _run(["samtools", "view", "-H", str(bam)]).stdout
    order = next(
        (f for line in hdr.splitlines() if line.startswith("@HD")
         for f in line.split("\t") if f.startswith(("SO:", "GO:"))),
        "none",
    )
    proc = subprocess.Popen(["samtools", "view", str(bam)],
                            stdout=subprocess.PIPE, text=True)
    seen: set[str] = set()
    prev, sampled, bad = None, 0, None
    for line in proc.stdout:
        name = line.split("\t", 1)[0]
        if name != prev:
            if name in seen:
                bad = name
                break
            seen.add(name)
            prev = name
        sampled += 1
        if sampled >= n:
            break
    proc.stdout.close()
    proc.terminate()
    proc.wait()
    return bad is None, order, sampled, bad


def check_pair_names(r1: Path, r2: Path) -> None:
    n = 0
    for n, (a, b) in enumerate(zip(fastq_names(r1), fastq_names(r2)), 1):
        if a != b:
            raise BamToFastqError(f"R1/R2 names differ at record {n}: {a!r} vs {b!r}")
    if count_fastq(r1) != count_fastq(r2):
        raise BamToFastqError("R1 and R2 have different read counts")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def gzip_deterministic(src: Path, dest: Path) -> None:
    """mtime=0 and no filename: same reads -> same bytes, so checksums hold."""
    with open(src, "rb") as fin, open(dest, "wb") as raw, gzip.GzipFile(
        filename="", mode="wb", fileobj=raw, mtime=0
    ) as out:
        shutil.copyfileobj(fin, out)


# ──────────────────────────────────────────────
# Core
# ──────────────────────────────────────────────

def convert(
    bam: Path,
    output_dir: Path,
    paired: bool = False,
    prefix: str | None = None,
    tags: bool = False,
    name_sort: bool = False,
    skip_check: bool = False,
    check_n: int = DEFAULT_CHECK_N,
) -> dict:
    """Convert one BAM. Returns the result dict (the caller writes the files)."""
    if not bam.is_file():
        raise BamToFastqError(f"input not found: {bam}")
    require_tools(paired)
    prefix = prefix or re.sub(r"\.bam$", "", bam.name)
    out_dir = output_dir / "fastq"
    names = [f"{prefix}_R1.fastq.gz", f"{prefix}_R2.fastq.gz"] if paired else [f"{prefix}.fastq.gz"]
    targets = [out_dir / n for n in names]
    existing = [p for p in [*targets, output_dir / "result.json"] if p.exists()]
    if existing:
        raise BamToFastqError(
            "refusing to overwrite existing output: "
            + ", ".join(str(p) for p in existing)
            + ". Choose a new --output directory."
        )

    counts: dict = {}
    counts["input_bam"] = count_bam(bam)
    if counts["input_bam"] == 0:
        raise BamToFastqError("input BAM has no primary records")
    if paired and counts["input_bam"] % 2:
        raise BamToFastqError(
            f"odd primary record count ({counts['input_bam']}) in paired mode: "
            "the BAM has orphan reads or is not paired. Drop --paired for single-end."
        )

    sort_check = None
    if paired and not name_sort and not skip_check:
        ok, order, sampled, bad = check_name_grouped(bam, check_n)
        sort_check = {"header": order, "reads_sampled": sampled, "grouped": ok}
        if not ok:
            raise BamToFastqError(
                f"read {bad!r} reappears after other reads: the BAM is not "
                "name-sorted/collated. Re-run with --name-sort (or --skip-check)."
            )

    output_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    commands: list[list[str]] = []
    try:
        # scratch next to the output (never /tmp: the OS overlay is small)
        with tempfile.TemporaryDirectory(dir=output_dir, prefix=".tmp-") as tmp:
            tmp = Path(tmp)
            src = bam
            if name_sort:
                src = tmp / "namesorted.bam"
                cmd = ["samtools", "sort", "-n", "-o", str(src), str(bam)]
                _run(cmd)
                commands.append(cmd)
                counts["name_sorted_bam"] = count_bam(src)
                if counts["name_sorted_bam"] != counts["input_bam"]:
                    raise BamToFastqError("record count changed during name sort")
            if paired:
                raw = [tmp / "R1.fastq", tmp / "R2.fastq"]
                cmd = ["bedtools", "bamtofastq", "-i", str(src),
                       "-fq", str(raw[0]), "-fq2", str(raw[1])]
                if tags:
                    cmd.append("-tags")
            else:
                raw = [tmp / "out.fastq"]
                cmd = ["samtools", "fastq", "-F", "0x900", "-0", str(raw[0])]
                if tags:
                    cmd += ["-T", "*"]
                cmd.append(str(src))
            # samtools fastq sends reads it cannot place to stdout; discard them
            # (the count check below catches any loss).
            _run(cmd)
            commands.append(cmd)
            for r, t in zip(raw, targets):
                gzip_deterministic(r, t)
    except BamToFastqError:
        for t in targets:
            t.unlink(missing_ok=True)
        raise

    counts["fastq"] = [count_fastq(t) for t in targets]
    expected = [counts["input_bam"] // 2] * 2 if paired else [counts["input_bam"]]
    if counts["fastq"] != expected:
        for t in targets:
            t.unlink(missing_ok=True)
        raise BamToFastqError(
            f"count mismatch: expected {expected} FASTQ records, got {counts['fastq']} "
            f"(input primary records: {counts['input_bam']}). "
            + ("" if paired else "A paired BAM needs --paired.")
        )
    if paired:
        try:
            check_pair_names(*targets)
        except BamToFastqError:
            for t in targets:
                t.unlink(missing_ok=True)
            raise

    versions = {"samtools": tool_version("samtools")}
    if paired:
        versions["bedtools"] = tool_version("bedtools")
    return {
        "skill": "bam-to-fastq",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "input": str(bam),
        "paired": paired,
        "tool": "bedtools bamtofastq" if paired else "samtools fastq",
        "tool_versions": versions,
        "name_sorted_first": name_sort,
        "sort_check": sort_check,
        "counts": counts,
        "commands": [shlex.join(c) for c in commands],
        "files": [
            {"output": t.relative_to(output_dir).as_posix(),
             "reads": c, "sha256": sha256_file(t)}
            for t, c in zip(targets, counts["fastq"])
        ],
        "disclaimer": DISCLAIMER,
    }


# ──────────────────────────────────────────────
# Report + reproducibility
# ──────────────────────────────────────────────

def write_report(output_dir: Path, result: dict) -> Path:
    c = result["counts"]
    steps = [("Input BAM (primary records)", c["input_bam"])]
    if "name_sorted_bam" in c:
        steps.append(("After name sort", c["name_sorted_bam"]))
    steps.append(("FASTQ records written", " + ".join(str(n) for n in c["fastq"])))
    step_rows = "\n".join(f"| {k} | {v} |" for k, v in steps)
    file_rows = "\n".join(
        f"| `{f['output']}` | {f['reads']} |" for f in result["files"]
    )
    sc = result["sort_check"]
    sort_line = (
        f"sampled {sc['reads_sampled']} reads (header `{sc['header']}`), names grouped"
        if sc else ("name-sorted first with `samtools sort -n`" if result["name_sorted_first"]
                    else "not applicable (single-end or skipped)")
    )
    pair_line = "R1/R2 read names match for every record" if result["paired"] else "n/a (single-end)"
    text = f"""# bam-to-fastq Report

**Date**: {result['generated_at']}
**Input**: `{result['input']}`
**Tool**: {result['tool']} ({', '.join(f'{k} {v}' for k, v in result['tool_versions'].items())})
**Layout**: {'paired-end' if result['paired'] else 'single-end'}

## Read-count checks

| Step | Count |
|------|------:|
{step_rows}

Result: **PASS** (the counts above agree).

- Sort check: {sort_line}
- Pair-name check: {pair_line}

## Output

| File | Reads |
|------|------:|
{file_rows}

*{DISCLAIMER}*
"""
    path = output_dir / "report.md"
    path.write_text(text)
    return path


def write_bundle(output_dir: Path, result: dict, argv_flags: str) -> None:
    lines = [
        "# ClawBio",
        f"python skills/bam-to-fastq/bam_to_fastq.py --input {shlex.quote(result['input'])} "
        f"--output <new_output_dir>{argv_flags}",
        "# Equivalent raw commands (temp paths shown as run):",
        *[f"# {c}" for c in result["commands"]],
    ]
    write_commands_sh(output_dir, "\n".join(lines))
    conda = [f"samtools={result['tool_versions']['samtools']}"]
    if "bedtools" in result["tool_versions"]:
        conda.append(f"bedtools={result['tool_versions']['bedtools']}")
    write_environment_yml(
        output_dir,
        "clawbio-bam-to-fastq",
        pip_deps=[],
        conda_deps=conda,
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
        description="Convert a BAM to FASTQ (paired: bedtools bamtofastq; "
        "single-end: samtools fastq), with read-count checks.")
    p.add_argument("--input", type=Path, help="Input BAM")
    p.add_argument("--output", type=Path, help="Output directory (must not already hold results)")
    p.add_argument("--demo", action="store_true", help="Run on the bundled synthetic paired BAM")
    p.add_argument("--paired", action="store_true",
                   help="Paired-end: write _R1/_R2 (needs a name-sorted BAM, or --name-sort)")
    p.add_argument("--name-sort", action="store_true",
                   help="Run `samtools sort -n` first (needed for coordinate-sorted BAMs)")
    p.add_argument("--prefix", help="Output file prefix (default: BAM name)")
    p.add_argument("--tags", action="store_true", help="Keep SAM tags in FASTQ headers")
    p.add_argument("--skip-check", action="store_true",
                   help="Skip the pre-run name-sort check (paired-end only)")
    p.add_argument("--check-n", type=int, default=DEFAULT_CHECK_N,
                   help=f"Reads to sample for the sort check (default {DEFAULT_CHECK_N})")
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.output is None:
        parser.error("--output is required")
    if args.demo:
        args.input, args.paired, args.name_sort = DEMO_BAM, True, True
    elif args.input is None:
        parser.error("--input is required (or use --demo)")

    try:
        result = convert(args.input, args.output, args.paired, args.prefix,
                         args.tags, args.name_sort, args.skip_check, args.check_n)
        (args.output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
        write_report(args.output, result)
        flags = "".join([
            " --paired" if args.paired else "",
            " --name-sort" if args.name_sort else "",
            f" --prefix {shlex.quote(args.prefix)}" if args.prefix else "",
            " --tags" if args.tags else "",
        ])
        write_bundle(args.output, result, flags)
    except BamToFastqError as exc:
        print(f"bam-to-fastq: error: {exc}", file=sys.stderr)
        return 1

    print(f"Count checks passed: {result['counts']}")
    print(f"FASTQ: {args.output / 'fastq'}")
    print(f"Report: {args.output / 'report.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
