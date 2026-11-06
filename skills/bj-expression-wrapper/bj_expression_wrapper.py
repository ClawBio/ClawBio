#!/usr/bin/env python3
"""bj-expression-wrapper: run BioSkryb's BJ-Expression RNA-seq pipeline from ClawBio.

BJ-Expression takes already-demultiplexed FASTQ pairs (one pair per biosample)
and runs seqtk downsampling, cutadapt/fastp trimming, Salmon, STAR + HTSeq,
Qualimap and MultiQC.

This wrapper validates the input sheet, checks the tools, builds one fixed
`nextflow run` command, runs it, and summarises what the pipeline published.
It does NOT contain or download any BioSkryb code: the user supplies their own
checkout of the pinned tag, and must accept BioSkryb's terms to run it.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import zlib
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
DEMO_CSV = SKILL_DIR / "examples" / "demo_samplesheet.csv"

PROG = "bj-expression-wrapper"
PIPELINE_REPO = "https://github.com/BioSkryb/bj-expression.git"
DEFAULT_PIPELINE_VERSION = "2.0.0"
# bj-expression hard-codes `seqtk sample -s123` for read 2 (read 1 uses seqtk_sample_seed), so any
# other seed would put the two mates of each pair out of sync.
REQUIRED_SEED = 123
DEMO_READS_PER_FILE = 10
# Nextflow 26.x cannot parse bj-expression 2.0.0's config: it includes two
# profile files (batch_dev_gpu.config, batch_prod_gpu.config) that are not in
# the public repo. 25.04.6 and 24.10.5 both parse it (config parse only).
DEFAULT_NEXTFLOW_VERSION = "25.04.6"
LICENSE_URL = "https://www.bioskryb.com/basejumper/local-license/"
# Fixed so the published directory is predictable: <publish_dir>_<timestamp>.
RUN_TIMESTAMP = "clawbio"
GENOMES = ("GRCh38", "GRCm39")  # the two the BJ-Expression README documents
PROFILES = ("docker", "singularity")
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
REMOTE_PREFIXES = ("s3://", "http://", "https://", "gs://", "az://", "ftp://")

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)

# Published by workflow/rnaseq.nf under <publish_dir>_<timestamp>/ (required
# ones must exist for a run to count as a success).
# Names follow the pipeline's own CODE and tests/main.nf.test, not its README: the
# README still lists pipeline_metrics_summary.csv, which 2.0.0 does not write
# (the module CREATE_QC_REPORT emits qualimap_stats_mqc.csv and the _percents file).
# Single-nucleus version: Salmon is not run. MultiQC, cell typing and the HTSeq matrix in the
# main results folder come from the INTRON-INCLUSIVE counts (HTSeq -t gene).
REQUIRED_OUTPUTS = {
    "multiqc_report": "multiqc/multiqc_report.html",
    "qualimap_stats": "secondary_analyses/secondary_metrics/qualimap_stats_mqc.csv",
    "pipeline_metrics_summary_percents": "secondary_analyses/secondary_metrics/pipeline_metrics_summary_percents.csv",
    "htseq_gene_counts": "secondary_analyses/quantification_htseq/df_gene_counts_starhtseq.tsv",
    "cell_typing": "tertiary_analyses/classification_cell_typing/df_cell_typing_summary_singler_hpca_gtex_tcga.tsv",
}
# The duplicate built from EXON counts (HTSeq -t exon), published to <results>_exon.
EXON_DUPLICATE_OUTPUTS = {
    "multiqc_report": "multiqc/multiqc_report.html",
    "htseq_gene_counts": "secondary_analyses/quantification_htseq/df_gene_counts_starhtseq.tsv",
    "cell_typing": "tertiary_analyses/classification_cell_typing/df_cell_typing_summary_singler_hpca_gtex_tcga.tsv",
}
OPTIONAL_OUTPUTS = {
    "dynamic_range_expression": "secondary_analyses/secondary_metrics/df_dynamicrange_expression.tsv",
}
ALIGNMENT_DIR = "secondary_analyses/alignment_htseq"

REQUIRED_COLUMNS = ("biosampleName", "read1", "read2")


class WrapperError(Exception):
    """A user-facing failure; printed as one line and exits 1."""


# ──────────────────────────────────────────────
# Input sheet
# ──────────────────────────────────────────────

def is_safe_sample_name(name: str) -> bool:
    """Sample names become file names and Nextflow tags: keep them plain."""
    return bool(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", name)) and ".." not in name


def _resolve_fastq(value: str, base_dir: Path) -> tuple[str, bool]:
    """Return (path-or-url, is_remote); local paths become absolute."""
    if value.lower().startswith(REMOTE_PREFIXES):
        return value, True
    p = Path(value).expanduser()
    if not p.is_absolute():
        p = base_dir / p
    return str(p.resolve()), False


def validate_input_csv(path: Path, require_groups: bool,
                       base_dir: Path | None = None) -> list[dict]:
    """Parse and validate the sample sheet; return normalised rows.

    Relative FASTQ paths resolve against `base_dir` (default: the current
    directory, as Nextflow does). Remote URLs are accepted without a check.
    """
    path = Path(path)
    if not path.is_file():
        raise WrapperError(f"input sheet not found: {path}")
    base = Path(base_dir) if base_dir else Path.cwd()

    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        header = [h.strip() for h in (reader.fieldnames or [])]
        missing = [c for c in REQUIRED_COLUMNS if c not in header]
        if missing:
            raise WrapperError(
                f"input sheet is missing required column(s): {', '.join(missing)} "
                f"(need biosampleName,read1,read2[,groups])")
        raw = [{(k or "").strip(): (v or "").strip() for k, v in r.items()}
               for r in reader]
    if not raw:
        raise WrapperError("input sheet has no samples (header only)")

    rows, seen = [], set()
    for i, r in enumerate(raw, start=2):  # line numbers, header is line 1
        name = r["biosampleName"]
        if not is_safe_sample_name(name):
            raise WrapperError(
                f"line {i}: biosampleName {name!r} is not allowed; use letters, "
                "digits, '.', '_' and '-' only, starting with a letter or digit")
        if name in seen:
            raise WrapperError(f"line {i}: duplicate biosampleName {name!r}")
        seen.add(name)

        fq, remote = [], False
        for col in ("read1", "read2"):
            value = r[col]
            if not value:
                raise WrapperError(f"line {i}: {col} is empty for {name}")
            if not value.lower().endswith(FASTQ_SUFFIXES):
                raise WrapperError(
                    f"line {i}: {col} {value!r} is not a FASTQ "
                    f"({', '.join(FASTQ_SUFFIXES)})")
            resolved, is_remote = _resolve_fastq(value, base)
            if not is_remote and not Path(resolved).is_file():
                raise WrapperError(f"line {i}: {col} file not found: {resolved}")
            remote = remote or is_remote
            fq.append(resolved)
        if fq[0] == fq[1]:
            raise WrapperError(
                f"line {i}: read1 and read2 are the same file for {name}")

        rows.append({"biosampleName": name, "read1": fq[0], "read2": fq[1],
                     "groups": r.get("groups", ""), "_remote": remote})

    if require_groups and any(not r["groups"] for r in rows):
        raise WrapperError(
            "the 'groups' column is required for every sample when QC plots are "
            "on (per the BJ-Expression README). Add it, or pass --skip-qc-plots.")
    return rows


def write_normalised_csv(rows: list[dict], dest: Path) -> Path:
    """Write the sheet the pipeline will read: absolute paths, known columns."""
    with_groups = any(r["groups"] for r in rows)
    cols = [*REQUIRED_COLUMNS, *(["groups"] if with_groups else [])]
    with open(dest, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow([r[c] for c in cols])
    return dest


# ──────────────────────────────────────────────
# Pipeline source, tools
# ──────────────────────────────────────────────

def _git(path: Path, *args: str) -> str:
    try:
        proc = subprocess.run(["git", *args], cwd=str(path), capture_output=True,
                              text=True, timeout=20)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise WrapperError(f"could not run git in {path}: {exc}") from exc
    return proc.stdout.strip() if proc.returncode == 0 else ""


def check_pipeline_source(pipeline_dir: Path | None, version: str,
                          allow_dirty: bool) -> dict:
    """The checkout must be the pinned tag, unmodified (unless allowed)."""
    how_to = (
        f"Get the pipeline yourself, from BioSkryb, at the pinned tag:\n"
        f"  git clone --branch {version} --depth 1 {PIPELINE_REPO}\n"
        f"then pass --pipeline-dir <that folder>. Using it is subject to "
        f"BioSkryb's terms: {LICENSE_URL}")
    if pipeline_dir is None:
        raise WrapperError(f"--pipeline-dir is required.\n{how_to}")
    d = Path(pipeline_dir).expanduser().resolve()
    if not (d / "main.nf").is_file() or not (d / "nextflow.config").is_file():
        raise WrapperError(f"no BJ-Expression checkout at {d} (main.nf or "
                           f"nextflow.config missing).\n{how_to}")
    tag = _git(d, "describe", "--tags", "--exact-match")
    commit = _git(d, "rev-parse", "HEAD")
    if tag != version:
        raise WrapperError(
            f"{d} is at {tag or 'an untagged commit'}, not the pinned {version}. "
            f"Run: git -C {d} checkout {version}")
    dirty = bool(_git(d, "status", "--porcelain"))
    if dirty and not allow_dirty:
        raise WrapperError(
            f"{d} has uncommitted changes, so the run would not match tag "
            f"{version}. Revert them, or pass --allow-dirty-pipeline "
            "(recorded in the result).")
    return {"repo": PIPELINE_REPO, "path": str(d), "version": version,
            "commit": commit, "dirty": dirty, "license_url": LICENSE_URL}


def nextflow_env(nxf_version: str) -> dict:
    """Environment for every nextflow call: pin NXF_VER unless told 'system'."""
    env = dict(os.environ)
    if nxf_version != "system":
        env["NXF_VER"] = nxf_version
    return env


def _version_line(cmd: list[str], pattern: str, env: dict | None = None) -> str:
    try:
        # the first call with a pinned NXF_VER downloads that Nextflow release
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600,
                              env=env)
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"
    m = re.search(pattern, proc.stdout + proc.stderr)
    return m.group(1) if m else "unknown"


def preflight_tools(profile: str, nxf_version: str = DEFAULT_NEXTFLOW_VERSION) -> dict:
    """Check Java, Nextflow and the container backend; report every problem."""
    problems, tools = [], {}
    if shutil.which("java"):
        tools["java"] = _version_line(["java", "-version"], r'version "([^"]+)"')
    else:
        problems.append(
            "java not found on PATH. Nextflow needs a JDK "
            "(for example: conda install -c conda-forge openjdk).")
    if shutil.which("nextflow"):
        tools["nextflow"] = _version_line(["nextflow", "-version"],
                                          r"version\s+(\S+)",
                                          env=nextflow_env(nxf_version))
    else:
        problems.append("nextflow not found on PATH. Install it from "
                        "https://www.nextflow.io/docs/latest/install.html")

    if profile == "docker":
        if not shutil.which("docker"):
            problems.append(
                "docker not found on PATH. Install Docker (needs a Carolina "
                "Cloud VM, not a container), or use --profile singularity.")
        else:
            ok = subprocess.run(["docker", "info"], capture_output=True).returncode == 0
            if ok:
                tools["docker"] = "available"
            else:
                problems.append("docker is installed but its daemon is not "
                                "reachable (`docker info` failed).")
    else:
        sif = shutil.which("apptainer") or shutil.which("singularity")
        if sif:
            tools["singularity"] = Path(sif).name
        else:
            problems.append("--profile singularity needs apptainer or "
                            "singularity on PATH; neither was found.")
    if problems:
        raise WrapperError("preflight failed:\n  - " + "\n  - ".join(problems))
    return tools


# ──────────────────────────────────────────────
# Reads available, and the explicit per-file choice
# ──────────────────────────────────────────────
# BioSkryb's SEQTK_SAMPLE divides n_reads by the reads of R1 AND R2 together, so its default
# n_reads = 100000 meant ~50,000 reads per file. It trims to read_length only when it samples
# (otherwise it copies the files unchanged), and it samples by fraction, so counts are approximate.

def count_fastq_reads(path: Path) -> int:
    """Number of records in a (gzipped) FASTQ; clear errors for unreadable or truncated files."""
    path = Path(path)
    lines, last = 0, b"\n"
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(path, "rb") as fh:
            while chunk := fh.read(1 << 22):
                lines += chunk.count(b"\n")
                last = chunk[-1:]
    except (OSError, EOFError, zlib.error) as exc:
        raise WrapperError(f"cannot read FASTQ {path.name}: {exc}") from exc
    if last != b"\n":
        lines += 1  # final record without a trailing newline
    if lines % 4:
        raise WrapperError(f"{path.name} is not a valid FASTQ: {lines} lines is not a multiple of 4")
    return lines // 4


def count_sample_reads(rows: list[dict]) -> dict[str, dict]:
    """{sample: {"r1": n, "r2": n}} (None for remote files, which are not counted)."""
    jobs = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, min(8, 2 * len(rows)))) as pool:
        for r in rows:
            for col in ("read1", "read2"):
                jobs[(r["biosampleName"], col)] = (
                    None if r["_remote"] else pool.submit(count_fastq_reads, Path(r[col])))
        out: dict[str, dict] = {}
        for (name, col), fut in jobs.items():
            out.setdefault(name, {})["r1" if col == "read1" else "r2"] = fut.result() if fut else None
    for name, c in out.items():
        if c["r1"] is not None and c["r2"] is not None and c["r1"] != c["r2"]:
            raise WrapperError(f"{name}: R1 has {c['r1']} reads but R2 has {c['r2']}; "
                               "the two files of a pair must match")
    return out


def format_read_table(counts: dict[str, dict]) -> str:
    w = max([len(n) for n in counts] + [6])
    lines = [f"  {'sample':<{w}}  {'R1':>12}  {'R2':>12}  {'pairs':>12}"]
    for name, c in counts.items():
        f = lambda v: "unknown" if v is None else f"{v:,}"  # noqa: E731
        lines.append(f"  {name:<{w}}  {f(c['r1']):>12}  {f(c['r2']):>12}  {f(c['r1']):>12}")
    return "\n".join(lines)


def resolve_read_choice(args) -> tuple[str, int | None] | None:
    rpf = args.reads_per_file
    if rpf is not None and rpf != "all" and args.skip_subsampling:
        raise WrapperError("--reads-per-file N and --skip-subsampling contradict each other; "
                           "use one number, or --reads-per-file all")
    if rpf == "all" or (rpf is None and args.skip_subsampling):
        return ("all", None)
    return None if rpf is None else ("reads_per_file", rpf)


def pipeline_default(pipeline_dir: Path | None, key: str) -> int | None:
    """A default read from the pinned checkout's nextflow.config (not from the README)."""
    if pipeline_dir is None:
        return None
    try:
        text = (Path(pipeline_dir) / "nextflow.config").read_text()
    except OSError:
        return None
    m = re.search(rf"^\s*{key}\s*=\s*(\d+)", text, re.M)
    return int(m.group(1)) if m else None


def build_subsampling(choice, counts, read_length, min_reads) -> dict:
    mode, n = choice
    warnings: list[str] = []
    samples = {}
    for name, c in counts.items():
        pairs = c["r1"]
        active = None if pairs is None else (mode == "reads_per_file" and n < pairs)
        samples[name] = {"r1": c["r1"], "r2": c["r2"], "pairs": pairs, "sampling_active": active,
                         "expected_reads_per_file": (None if pairs is None else
                                                     (min(n, pairs) if mode == "reads_per_file" else pairs))}
        if mode == "reads_per_file" and pairs is not None and not active:
            warnings.append(f"{name}: {n:,} is not below its {pairs:,} reads per file, so the pipeline "
                            f"copies the files unchanged and they are not trimmed to "
                            f"{read_length if read_length else 'read_length'} bp (it trims only when it samples)")
        if pairs is not None and min_reads is not None and 2 * pairs < min_reads:
            warnings.append(f"{name}: R1+R2 = {2 * pairs:,} reads is below min_reads ({min_reads:,}); "
                            "the pipeline drops this sample")
    if mode == "all":
        warnings.insert(0, "no subsampling: the SEQTK step is skipped, so reads are not trimmed to "
                           f"{read_length if read_length else 'read_length'} bp either")
    return {"mode": mode, "reads_per_file": n, "pipeline_n_reads": None if n is None else 2 * n,
            "seed": REQUIRED_SEED, "trim_to_bp": (read_length if mode == "reads_per_file" else None),
            "min_reads": min_reads, "samples": samples, "warnings": warnings}


# ──────────────────────────────────────────────
# Params and command
# ──────────────────────────────────────────────

def build_params(args: argparse.Namespace, csv_path: Path, out_dir: Path,
                 choice: tuple[str, int | None] | None = None) -> dict:
    """Required values plus only the overrides the user asked for.

    Anything not listed keeps the pipeline's own default; passing our own
    defaults would silently change BioSkryb's method.
    """
    params: dict = {
        "input_csv": str(csv_path),
        "publish_dir": str(out_dir / "results"),
        "timestamp": RUN_TIMESTAMP,
        "genome": args.genome,
    }
    for attr, key in (("read_length", "read_length"),
                      ("min_reads", "min_reads"), ("seed", "seqtk_sample_seed"),
                      ("instrument", "instrument"), ("max_cpus", "max_cpus"),
                      ("max_memory", "max_memory"), ("genomes_base", "genomes_base")):
        value = getattr(args, attr)
        if value is not None:
            params[key] = value
    if choice and choice[0] == "reads_per_file":
        params["n_reads"] = 2 * choice[1]  # the pipeline divides n_reads by R1+R2 together
    elif choice and choice[0] == "all":
        params["skip_subsampling"] = True
    for flag in ("skip_cutadapt", "skip_qc_plots"):
        if getattr(args, flag):
            params[flag] = True
    return params


def new_trace_path(run_dir: Path, now: datetime | None = None) -> Path:
    """A trace file name that cannot already exist: Nextflow refuses to overwrite one.

    Microseconds in the name, plus a numeric suffix if that name is somehow taken
    (two invocations in the same instant)."""
    now = now or datetime.now(timezone.utc)
    stem = f"trace_{now:%Y%m%dT%H%M%S}.{now.microsecond:06d}Z"
    path, i = run_dir / f"{stem}.tsv", 1
    while path.exists():
        path = run_dir / f"{stem}-{i}.tsv"
        i += 1
    return path


def build_command(pipeline_dir: str, params_file: Path, work_dir: Path,
                  profile: str, resume: bool,
                  extra_config: Path | None = None,
                  trace: Path | None = None) -> list[str]:
    cmd = ["nextflow", "run", str(Path(pipeline_dir) / "main.nf"),
           "-params-file", str(params_file), "-work-dir", str(work_dir),
           "-ansi-log", "false"]
    if trace is not None:
        cmd += ["-with-trace", str(trace)]
    if extra_config is not None:
        cmd += ["-c", str(extra_config)]
    if profile != "docker":  # docker is enabled by the pipeline's own config
        cmd += ["-profile", profile]
    if resume:
        cmd.append("-resume")
    return cmd


# ──────────────────────────────────────────────
# Execution and outputs
# ──────────────────────────────────────────────

def run_nextflow(cmd: list[str], run_dir: Path,
                 env: dict | None = None) -> tuple[int, str]:
    """Run in run_dir (Nextflow litters its cwd); logs go to files, tail returned."""
    run_dir.mkdir(parents=True, exist_ok=True)
    out_log, err_log = run_dir / "nextflow.stdout.log", run_dir / "nextflow.stderr.log"
    with open(out_log, "w") as so, open(err_log, "w") as se:
        proc = subprocess.run(cmd, cwd=run_dir, stdout=so, stderr=se, env=env)
    tail = "\n".join(err_log.read_text(errors="replace").splitlines()[-40:])
    return proc.returncode, tail


FAILED_STATES = ("FAILED", "ABORTED")


def parse_trace(path: Path) -> dict | None:
    """Summarise a Nextflow trace by each task's FINAL state.

    BJ-Expression sets errorStrategy to retry twice and then *ignore*, so a run
    can exit 0 with steps missing. A task name whose last row is FAILED/ABORTED
    was ignored; one that failed and later completed was only retried.
    Returns None if the trace is missing or has no name/status columns.
    """
    path = Path(path)
    if not path.is_file():
        return None
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or not {"name", "status"} <= set(reader.fieldnames):
            return None
        history: dict[str, list[str]] = {}
        for row in reader:
            if row.get("name"):
                history.setdefault(row["name"], []).append((row.get("status") or "").strip())
    ignored = sorted(n for n, h in history.items() if h[-1] in FAILED_STATES)
    retried = sorted(n for n, h in history.items()
                     if h[-1] not in FAILED_STATES and any(x in FAILED_STATES for x in h))
    return {"tasks": len(history), "ignored_failures": ignored, "retried_tasks": retried}


def collect_exon_duplicate(exon_dir: Path) -> dict:
    """Outputs of the exon-count duplicate run (cell typing and MultiQC)."""
    found, missing = {}, []
    for key, rel in EXON_DUPLICATE_OUTPUTS.items():
        if (exon_dir / rel).is_file():
            found[key] = str(exon_dir / rel)
        else:
            missing.append(rel)
    return {"results_dir": str(exon_dir), "outputs": found, "missing_outputs": missing}


def collect_outputs(results_dir: Path, samples: list[str]) -> dict:
    found, missing = {"results_dir": str(results_dir)}, []
    for key, rel in REQUIRED_OUTPUTS.items():
        if (results_dir / rel).is_file():
            found[key] = str(results_dir / rel)
        else:
            missing.append(rel)
    for key, rel in OPTIONAL_OUTPUTS.items():
        if (results_dir / rel).is_file():
            found[key] = str(results_dir / rel)
    bams = sorted(p.name[:-4] for p in (results_dir / ALIGNMENT_DIR).glob("*.bam")) \
        if (results_dir / ALIGNMENT_DIR).is_dir() else []
    found["aligned_samples"] = bams
    return {"outputs": found, "missing_outputs": missing,
            "samples_missing_from_outputs": [s for s in samples if s not in bams]}


# ──────────────────────────────────────────────
# Real files instead of symlinks
# ──────────────────────────────────────────────

def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class _CopyMismatch(Exception):
    pass


def _replace_link_with_file(link: Path, target: Path) -> int:
    tmp = link.with_name(link.name + ".materialise-tmp")
    try:
        shutil.copy2(target, tmp)
        if _sha256(tmp) != _sha256(target):
            raise _CopyMismatch(str(link))
        size = tmp.stat().st_size
        os.replace(tmp, link)  # atomically swaps the symlink for the real file
        return size
    finally:
        if tmp.exists() or tmp.is_symlink():
            tmp.unlink()


def _replace_link_with_dir(link: Path, target: Path) -> tuple[int, int]:
    tmp = link.with_name(link.name + ".materialise-tmp")
    try:
        shutil.copytree(target, tmp, symlinks=False)  # follows links inside, too
        n = size = 0
        for root, _dirs, files in os.walk(target, followlinks=True):
            for name in files:
                src = Path(root) / name
                dst = tmp / src.relative_to(target)
                if not dst.is_file() or _sha256(dst) != _sha256(src):
                    raise _CopyMismatch(str(link))
                n += 1
                size += dst.stat().st_size
        link.unlink()
        os.rename(tmp, link)
        return n, size
    finally:
        if tmp.exists():
            shutil.rmtree(tmp, ignore_errors=True)


def materialise_results(results_dir: Path) -> dict:
    """Replace every symlink under results_dir with a verified real copy.

    BioSkryb's publishDir directives set no mode, so Nextflow publishes symlinks
    into work/, and a config override does not change that. Each copy is checked
    by SHA-256 before it replaces its link; the files in work/ are never modified.
    Broken links and failed copies are reported, not fatal. Safe to repeat.
    """
    results_dir = Path(results_dir)
    if not results_dir.is_dir():
        raise WrapperError(f"results directory not found: {results_dir}")
    report: dict = {"files": 0, "dirs": 0, "bytes": 0, "already_real": 0,
                    "broken": [], "failed": []}
    links: list[Path] = []
    for root, dirs, files in os.walk(results_dir):  # does not descend into links
        for name in dirs + files:
            p = Path(root) / name
            if p.is_symlink():
                links.append(p)
            elif p.is_file():
                report["already_real"] += 1
    for link in links:
        rel = link.relative_to(results_dir).as_posix()
        target = Path(os.path.realpath(link))
        if not target.exists():
            report["broken"].append(rel)
            continue
        try:
            if target.is_dir():
                n, size = _replace_link_with_dir(link, target)
                report["dirs"] += 1
                report["files"] += n
                report["bytes"] += size
            else:
                report["bytes"] += _replace_link_with_file(link, target)
                report["files"] += 1
        except (_CopyMismatch, OSError):
            report["failed"].append(rel)
    return report


# ──────────────────────────────────────────────
# Report and bundle
# ──────────────────────────────────────────────

def write_report(out_dir: Path, result: dict) -> Path:
    status = result["status"]
    lines = [
        "# bj-expression-wrapper Report", "",
        f"**Date**: {result['generated_at']}",
        f"**Status**: {status}",
        f"**Genome**: {result['params']['genome']}",
        f"**Samples**: {len(result['samples'])} ({', '.join(result['samples'])})",
    ]
    pipe = result.get("pipeline")
    if pipe:
        lines.append(f"**Pipeline**: BioSkryb bj-expression {pipe['version']} "
                     f"(commit {pipe['commit'][:12] or 'unknown'}"
                     f"{', modified' if pipe['dirty'] else ''})")
    if result.get("tools"):
        lines.append("**Tools**: " + ", ".join(f"{k} {v}" for k, v in result["tools"].items()))
    lines.append("")

    if status == "demo_plan_only":
        lines += ["> **Plan only.** This demo validated a synthetic sample sheet and "
                  "wrote the plan. Nothing was run: no pipeline, Nextflow or "
                  "container was used.", ""]
    elif status == "check_passed":
        lines += ["> **Check only.** The inputs, pipeline checkout and tools passed; "
                  "the pipeline was not run.", ""]

    outs = result.get("outputs") or {}
    if status in ("success", "incomplete") and outs:
        lines += ["## Outputs", "", "| Output | Path |", "|---|---|"]
        for key in (*REQUIRED_OUTPUTS, *OPTIONAL_OUTPUTS):
            if key in outs:
                lines.append(f"| {key} | `{outs[key]}` |")
        lines.append("")
    if result.get("missing_outputs"):
        lines += ["**Expected but not produced:** "
                  + ", ".join(f"`{m}`" for m in result["missing_outputs"]), ""]
    sub = result.get("subsampling")
    if sub:
        lines += ["## Subsampling", ""]
        if sub["mode"] == "all":
            lines += ["**All reads used** (no subsampling). The SEQTK step is skipped, so the reads are "
                      "**not trimmed** to `read_length` either.", ""]
        else:
            trim = sub["trim_to_bp"] or "read_length"
            lines += [f"**{sub['reads_per_file']:,} reads from each file.** The pipeline's `n_reads` was set "
                      f"to {sub['pipeline_n_reads']:,} because it counts R1 and R2 together. It samples by "
                      f"fraction with seed {sub['seed']}, so the numbers are approximate. Reads are trimmed to "
                      f"{trim} bp only where sampling happens.", ""]
        lines += ["| Sample | Reads per file available | Sampled? | Expected per file |", "|---|--:|---|--:|"]
        for name, c in sub["samples"].items():
            av = "unknown (remote)" if c["pairs"] is None else f"{c['pairs']:,}"
            sam = "n/a" if c["sampling_active"] is None else ("yes" if c["sampling_active"] else "no, copied unchanged")
            ex = "unknown" if c["expected_reads_per_file"] is None else f"{c['expected_reads_per_file']:,}"
            lines.append(f"| {name} | {av} | {sam} | {ex} |")
        lines.append("")
        for w in sub["warnings"]:
            lines += [f"> **Warning:** {w}", ""]
    ed = result.get("exon_duplicate")
    if ed:
        lines += ["## Counts used for MultiQC and cell typing", "",
                  "The main results use **intron-inclusive** counts (HTSeq `-t gene`) for MultiQC "
                  "and cell typing. Duplicates built from **exon** counts (HTSeq `-t exon`) are in "
                  f"`{Path(ed['results_dir']).name}/`: their own MultiQC report and cell-typing tables.", ""]
        if ed["missing_outputs"]:
            lines += ["> **Warning:** missing from the exon duplicate: "
                      + ", ".join(f"`{m}`" for m in ed["missing_outputs"]), ""]
    mat = result.get("results_materialised")
    if mat:
        lines += [f"**Results are real files:** {mat['files']} files ({mat['bytes'] / 1e6:.1f} MB) "
                  "were copied out of the work directory into `results_clawbio/`, because Nextflow "
                  "publishes symlinks. The work directory is no longer needed for the results; it "
                  "is still needed for `--resume`.", ""]
        if mat["broken"] or mat["failed"]:
            lines += [f"> **Warning:** {len(mat['broken'])} broken link(s) and {len(mat['failed'])} "
                      "copy failure(s) were left as symlinks: "
                      + ", ".join(f"`{x}`" for x in (mat["broken"] + mat["failed"])[:8]), ""]
    dropped = result.get("samples_missing_from_outputs") or []
    if dropped and status in ("success", "incomplete"):
        lines += [f"> **Warning:** no alignment was published for: {', '.join(dropped)}. "
                  "BJ-Expression filters out samples with fewer reads than `min_reads` "
                  "before quantification, so they will be missing from the count "
                  "tables. Check the pipeline log before assuming a sample was analysed.",
                  ""]
    if status == "failed":
        lines += ["## Nextflow stderr (tail)", "", "```", result["stderr_tail"], "```", ""]

    ignored_tasks = result.get("ignored_failures")
    if ignored_tasks:
        lines += ["## Tasks that failed and were ignored", "",
                  "BJ-Expression retries a failed task twice and then ignores it, so the "
                  "run continues and can exit 0. These tasks did not finish:", ""]
        lines += [f"- `{t}`" for t in ignored_tasks]
        lines += ["", "Anything downstream of them may be missing or partial. Check "
                  "the trace file named in `result.json` (`trace.path`) and the task work directories.", ""]
    elif ignored_tasks is None and status in ("success", "success_with_ignored_failures",
                                              "incomplete"):
        lines += ["> **No Nextflow trace was available**, so ignored task failures could "
                  "not be checked. Do not read this status as proof every step finished.", ""]

    overridden = sorted(set(result["params"]) - {"input_csv", "publish_dir",
                                                 "timestamp", "genome"})
    lines += ["## Parameters", "",
              "Only these were set by the wrapper beyond the required ones: "
              + (", ".join(overridden) if overridden else "none")
              + ". Everything else is BioSkryb's default for the pinned tag.", "",
              "## Pipeline licence", "",
              f"BJ-Expression is BioSkryb software, used under BioSkryb's terms "
              f"({LICENSE_URL}). This wrapper contains none of its code.", "",
              f"*{DISCLAIMER}*", ""]
    path = out_dir / "report.md"
    path.write_text("\n".join(lines))
    return path


def write_bundle(out_dir: Path, command: str | None,
                 extra_config: Path | None = None) -> None:
    write_commands_sh(out_dir, command or "# plan only: no command was run")
    write_environment_yml(
        out_dir, "clawbio-bj-expression-wrapper", pip_deps=[],
        conda_deps=["openjdk", "nextflow"],
        channels=["conda-forge", "bioconda", "nodefaults"])
    repro = out_dir / "reproducibility"
    write_checksums(
        [out_dir / "report.md", out_dir / "result.json", repro / "commands.sh",
         repro / "params.json", repro / "samplesheet.normalised.csv",
         *([extra_config] if extra_config else [])],
        out_dir, anchor=out_dir)


# ──────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────

def _int(minimum: int, what: str):
    def conv(value: str) -> int:
        try:
            n = int(value)
        except ValueError:
            raise argparse.ArgumentTypeError(f"{what} must be an integer")
        if n < minimum:
            raise argparse.ArgumentTypeError(f"{what} must be >= {minimum}")
        return n
    return conv


def _nxf_version(value: str) -> str:
    if value == "system" or re.fullmatch(r"\d+\.\d+\.\d+", value):
        return value
    raise argparse.ArgumentTypeError(
        "--nextflow-version must look like 25.04.6 (or 'system' to use the "
        "installed Nextflow unchanged)")


def _reads_per_file(value: str):
    if value == "all":
        return "all"
    return _int(1, "--reads-per-file")(value)


def _old_n_reads(value: str):
    raise argparse.ArgumentTypeError(
        "--n-reads was replaced by --reads-per-file. The pipeline's n_reads counts R1 and R2 "
        "together, so the old default of 100000 gave about 50,000 reads per file. Give the "
        "reads you want from each file: --reads-per-file N (or 'all').")


def _seed(value: str) -> int:
    n = _int(0, "--seed")(value)
    if n != REQUIRED_SEED:
        raise argparse.ArgumentTypeError(
            f"--seed must be {REQUIRED_SEED}: the pipeline hard-codes -s{REQUIRED_SEED} for read 2 and "
            "applies the seed parameter only to read 1, so any other seed would put the two mates "
            "of each pair out of sync.")
    return n


def _memory(value: str) -> str:
    if not re.fullmatch(r"\d+(?:\.\d+)?[. ]?(?:KB|MB|GB|TB)", value):
        raise argparse.ArgumentTypeError("--max-memory looks like 50.GB or 512.MB")
    return value


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog=PROG,
        description="Run BioSkryb BJ-Expression (RNA-seq from demultiplexed FASTQ pairs).")
    p.add_argument("--input", type=Path,
                   help="CSV: biosampleName,read1,read2[,groups]; relative paths "
                        "resolve against the current directory")
    p.add_argument("--output", type=Path, help="Output directory")
    p.add_argument("--demo", action="store_true",
                   help="Validate a bundled synthetic sheet and write the plan; runs nothing")
    p.add_argument("--check", action="store_true",
                   help="Validate inputs, pipeline checkout and tools; do not run")
    p.add_argument("--pipeline-dir", type=Path,
                   help="Your checkout of BioSkryb/bj-expression at the pinned tag")
    p.add_argument("--pipeline-version", default=DEFAULT_PIPELINE_VERSION,
                   help=f"Tag the checkout must be at (default {DEFAULT_PIPELINE_VERSION})")
    p.add_argument("--accept-bioskryb-terms", action="store_true",
                   help=f"Confirm you accept BioSkryb's terms ({LICENSE_URL}); needed to run")
    p.add_argument("--allow-dirty-pipeline", action="store_true",
                   help="Allow uncommitted changes in the checkout (recorded)")
    p.add_argument("--genome", choices=GENOMES, default="GRCh38")
    p.add_argument("--profile", choices=PROFILES, default="docker",
                   help="Container backend (default docker)")
    p.add_argument("--reads-per-file", type=_reads_per_file, default=None, metavar="N|all",
                   help="REQUIRED for a run or --check: reads to keep from EACH file of a pair (a "
                        "number), or 'all' for no subsampling. The wrapper shows the reads "
                        "available first (--count-reads shows them without running).")
    p.add_argument("--count-reads", action="store_true",
                   help="Count the reads in every FASTQ of --input, print the table, and stop")
    p.add_argument("--n-reads", type=_old_n_reads, default=None, help=argparse.SUPPRESS)
    p.add_argument("--read-length", type=_int(1, "--read-length"),
                   help="Read length to trim to (pipeline default if omitted)")
    p.add_argument("--min-reads", type=_int(0, "--min-reads"),
                   help="Minimum reads for a sample to be analysed (pipeline default if omitted)")
    p.add_argument("--seed", type=_seed, default=None,
                   help="seqtk seed; only 123 is allowed (the pipeline hard-codes it for read 2)")
    p.add_argument("--skip-subsampling", action="store_true")
    p.add_argument("--skip-cutadapt", action="store_true")
    p.add_argument("--skip-qc-plots", action="store_true")
    p.add_argument("--instrument", help="Sequencer, e.g. NovaSeq (selects two-colour chemistry)")
    p.add_argument("--max-cpus", type=_int(1, "--max-cpus"))
    p.add_argument("--max-memory", type=_memory)
    p.add_argument("--genomes-base", help="Mirror of the reference bundle (default: BioSkryb's public S3)")
    p.add_argument("--resume", action="store_true", help="Pass -resume to Nextflow")
    p.add_argument("--strict", action="store_true",
                   help="Exit 1 if any task failed and was ignored by the pipeline")
    p.add_argument("--keep-symlinks", action="store_true",
                   help="Leave Nextflow's symlinks in results_clawbio/ (default: replace "
                        "them with verified real copies after a successful run)")
    p.add_argument("--materialise-only", action="store_true",
                   help="Do not run anything: replace the symlinks in an existing run's "
                        "results_clawbio/ (use with --output <that run's folder>)")
    p.add_argument("--nextflow-version", type=_nxf_version,
                   default=DEFAULT_NEXTFLOW_VERSION,
                   help=f"Nextflow release to run (default {DEFAULT_NEXTFLOW_VERSION}; "
                        "'system' = whatever is installed). Newer 26.x releases "
                        "cannot parse this pipeline's config.")
    p.add_argument("--extra-config", type=Path,
                   help="Extra Nextflow config passed with -c (for example container "
                        "settings for hosts without Docker); copied into the bundle")
    return p


def run(args: argparse.Namespace) -> int:
    if args.count_reads:
        if args.input is None:
            raise WrapperError("--input is required with --count-reads")
        counts = count_sample_reads(validate_input_csv(args.input, require_groups=False))
        print("Reads available per sample (R1 / R2 / pairs):\n" + format_read_table(counts))
        return 0
    if args.output is None:
        raise WrapperError("--output is required")
    if args.materialise_only:
        out = args.output.expanduser().resolve()
        rep = materialise_results(out / f"results_{RUN_TIMESTAMP}")
        (out / "materialise_report.json").write_text(json.dumps(rep, indent=2) + "\n")
        print(f"{PROG}: materialised {rep['files']} files ({rep['bytes'] / 1e6:.1f} MB) as real "
              f"copies; {rep['already_real']} were already real, {len(rep['broken'])} broken, "
              f"{len(rep['failed'])} failed -> {out / 'materialise_report.json'}")
        return 1 if rep["broken"] or rep["failed"] else 0
    if args.demo:
        mode = "demo"
    elif args.check:
        mode = "check"
    else:
        mode = "run"
        if args.input is None:
            raise WrapperError("--input is required (or use --demo / --check)")
        if not args.accept_bioskryb_terms:
            raise WrapperError(
                "BJ-Expression is BioSkryb software under BioSkryb's terms "
                f"({LICENSE_URL}). Read them, then re-run with "
                "--accept-bioskryb-terms to confirm.")
    if mode == "check" and args.input is None:
        raise WrapperError("--input is required with --check")

    # Nextflow runs from <output>/nextflow_run, so every path must be absolute.
    out_dir: Path = args.output.expanduser().resolve()
    if (out_dir / "result.json").exists() and not (mode == "run" and args.resume):
        raise WrapperError(
            f"{out_dir} already holds results; choose a new --output "
            "(or pass --resume to continue a Nextflow run).")

    if mode == "demo":
        csv_in, base = DEMO_CSV, DEMO_CSV.parent
    else:
        csv_in, base = args.input, None
    rows = validate_input_csv(csv_in, require_groups=not args.skip_qc_plots,
                              base_dir=base)
    samples = [r["biosampleName"] for r in rows]

    choice = resolve_read_choice(args)
    counts = count_sample_reads(rows)
    if mode == "demo" and choice is None:
        choice = ("reads_per_file", DEMO_READS_PER_FILE)
    if choice is None:
        raise WrapperError(
            "choose how many reads to use from each file before the run starts.\n"
            "Reads available per sample (R1 / R2 / pairs):\n" + format_read_table(counts) + "\n"
            "Give --reads-per-file N (reads from EACH file of a pair; the pipeline samples by "
            "fraction, so the result is approximate) or --reads-per-file all to use every read. "
            "The pipeline's own default (100000 reads counted across R1 and R2 together, about "
            "50,000 per file) is no longer applied silently.")

    pipeline, tools = None, {}
    if mode != "demo":
        pipeline = check_pipeline_source(args.pipeline_dir, args.pipeline_version,
                                         args.allow_dirty_pipeline)
        tools = preflight_tools(args.profile, args.nextflow_version)

    pdir = pipeline["path"] if pipeline else None
    sub = build_subsampling(
        choice, counts,
        args.read_length if args.read_length is not None else pipeline_default(pdir, "read_length"),
        args.min_reads if args.min_reads is not None else pipeline_default(pdir, "min_reads"))
    print(f"{PROG}: reads available per sample (R1 / R2 / pairs):\n{format_read_table(counts)}")
    print(f"{PROG}: using " + ("ALL reads (no subsampling)" if sub["mode"] == "all" else
          f"{sub['reads_per_file']:,} reads from each file (pipeline n_reads = {sub['pipeline_n_reads']:,})"))
    for w in sub["warnings"]:
        print(f"{PROG}: WARNING: {w}")

    if args.extra_config is not None and not args.extra_config.is_file():
        raise WrapperError(f"--extra-config file not found: {args.extra_config}")

    out_dir.mkdir(parents=True, exist_ok=True)
    repro = out_dir / "reproducibility"
    repro.mkdir(exist_ok=True)
    kept_config = None
    if args.extra_config is not None:
        kept_config = repro / "extra.config"
        shutil.copyfile(args.extra_config, kept_config)
    csv_path = write_normalised_csv(rows, repro / "samplesheet.normalised.csv")
    params = build_params(args, csv_path, out_dir, choice)
    (repro / "params.json").write_text(json.dumps(params, indent=2) + "\n")

    run_dir = out_dir / "nextflow_run"
    # One trace file per invocation: Nextflow will not overwrite an existing trace
    # file, so a fixed name made a --resume report the PREVIOUS run's trace.
    trace_path = new_trace_path(run_dir)
    cmd = build_command(pipeline["path"] if pipeline else "<path-to-bj-expression-checkout>",
                        repro / "params.json", out_dir / "work", args.profile, args.resume,
                        kept_config, trace_path)
    pin = "" if args.nextflow_version == "system" else f"NXF_VER={args.nextflow_version} "
    command_text = f"cd {shlex.quote(str(run_dir))} && {pin}{shlex.join(cmd)}"

    result = {
        "skill": PROG,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": {"demo": "demo_plan_only", "check": "check_passed", "run": "running"}[mode],
        "terms_accepted": bool(args.accept_bioskryb_terms),
        "pipeline": pipeline,
        "tools": tools,
        "profile": args.profile,
        "results_materialised": None,
        "exon_materialised": None,
        "exon_duplicate": None,
        "subsampling": sub,
        "nextflow_version_requested": args.nextflow_version,
        "extra_config": str(kept_config) if kept_config else None,
        "params": params,
        "samples": samples,
        "command": command_text,
        "disclaimer": DISCLAIMER,
    }
    exit_code = 0
    if mode == "run":
        code, tail = run_nextflow(cmd, run_dir, nextflow_env(args.nextflow_version))
        result["nextflow_exit_code"] = code
        result["stderr_tail"] = tail
        found = collect_outputs(Path(params["publish_dir"] + "_" + RUN_TIMESTAMP), samples)
        exon_dir = Path(params["publish_dir"] + "_exon_" + RUN_TIMESTAMP)
        result["exon_duplicate"] = collect_exon_duplicate(exon_dir)
        found["missing_outputs"] += [f"exon duplicate: {m}" for m in result["exon_duplicate"]["missing_outputs"]]
        result.update(found)
        results_dir = Path(params["publish_dir"] + "_" + RUN_TIMESTAMP)
        if code == 0 and not args.keep_symlinks and results_dir.is_dir():
            result["results_materialised"] = materialise_results(results_dir)
            if exon_dir.is_dir():
                result["exon_materialised"] = materialise_results(exon_dir)
        trace = parse_trace(trace_path)
        result["trace"] = ({"available": True, "path": str(trace_path), "tasks": trace["tasks"],
                            "retried_tasks": trace["retried_tasks"]} if trace
                           else {"available": False, "path": str(trace_path), "tasks": 0,
                                 "retried_tasks": []})
        result["ignored_failures"] = trace["ignored_failures"] if trace else None
        if code != 0:
            result["status"], exit_code = "failed", 1
        elif found["missing_outputs"]:
            result["status"], exit_code = "incomplete", 1
        elif result["ignored_failures"]:
            result["status"] = "success_with_ignored_failures"
            exit_code = 1 if args.strict else 0
        else:
            result["status"] = "success"

    (out_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    write_report(out_dir, result)
    write_bundle(out_dir, None if mode == "demo" else command_text, kept_config)

    if result["status"] == "failed":
        print(f"{PROG}: error: Nextflow exited {result['nextflow_exit_code']}\n"
              f"{result['stderr_tail']}", file=sys.stderr)
    elif result["status"] == "incomplete":
        print(f"{PROG}: error: Nextflow exited 0 but did not publish: "
              + ", ".join(result["missing_outputs"]), file=sys.stderr)
    else:
        print(f"{PROG}: {result['status']} -> {out_dir / 'report.md'}")
    if result["status"] == "success_with_ignored_failures":
        names = ", ".join(result["ignored_failures"])
        msg = f"tasks failed and were ignored by the pipeline: {names}"
        if exit_code:
            print(f"{PROG}: error: {msg} (--strict)", file=sys.stderr)
        else:
            print(f"{PROG}: WARNING: {msg}")
    return exit_code


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return run(args)
    except WrapperError as exc:
        print(f"{PROG}: error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
