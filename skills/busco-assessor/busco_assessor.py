"""
busco_assessor.py — BUSCO genome/transcriptome/protein completeness assessor for ClawBio.

Usage:
    python skills/busco-assessor/busco_assessor.py --input assembly.fna --mode genome --lineage bacteria_odb12 --output /tmp/busco_out
    python skills/busco-assessor/busco_assessor.py --input assembly.fna --mode genome --auto-lineage-prok --output /tmp/busco_out
    python skills/busco-assessor/busco_assessor.py --input assembly.fna --organism "fruit fly" --output /tmp/busco_out
    python skills/busco-assessor/busco_assessor.py --demo --output /tmp/busco_demo
    python skills/busco-assessor/busco_assessor.py --demo-live --output /tmp/busco_live_demo
"""
from __future__ import annotations

import argparse
import datetime
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DISCLAIMER = (
    "ClawBio is a research and educational tool. "
    "It is not a medical device and does not provide clinical diagnoses. "
    "Consult a healthcare professional before making any medical decisions."
)

BUSCO_MODES = ("genome", "transcriptome", "proteins")

DEFAULT_CPU = os.cpu_count() or 4

# Agentic lineage routing table: (keywords, flag, lineage_value_or_None)
# First match wins. Checked case-insensitively against the organism hint string.
LINEAGE_ROUTING: list[tuple[list[str], str, str | None]] = [
    (
        ["bacteria", "bacterial", "prokaryote", "prokaryotic",
         "e. coli", "escherichia", "salmonella", "streptococcus",
         "staphylococcus", "mycobacterium", "vibrio", "pseudomonas",
         "lactobacillus", "bacillus", "clostridium"],
        "--auto-lineage-prok", None,
    ),
    (
        ["archaea", "archaeon", "archaeal"],
        "--lineage", "archaea_odb12",
    ),
    (
        ["human", "homo sapiens", "hg38", "hg19", "grch38", "grch37"],
        "--lineage", "primates_odb10",
    ),
    (
        ["mouse", "mus musculus", "murine", "rat", "rattus"],
        "--lineage", "mammalia_odb10",
    ),
    (
        ["zebrafish", "danio rerio", "fish", "teleost"],
        "--lineage", "vertebrata_odb10",
    ),
    (
        ["bird", "avian", "chicken", "gallus", "aves"],
        "--lineage", "aves_odb10",
    ),
    (
        ["fruit fly", "drosophila", "diptera"],
        "--lineage", "diptera_odb10",
    ),
    (
        ["insect", "mosquito", "lepidoptera", "coleoptera", "hymenoptera"],
        "--lineage", "insecta_odb10",
    ),
    (
        ["plant", "arabidopsis", "rice", "maize", "soybean", "wheat", "embryophyte"],
        "--lineage", "embryophyta_odb10",
    ),
    (
        ["fungus", "fungi", "yeast", "saccharomyces", "aspergillus", "candida"],
        "--lineage", "fungi_odb10",
    ),
    (
        ["eukaryote", "eukaryotic"],
        "--auto-lineage-euk", None,
    ),
]

# Demo data constants
DEMO_COMPLETENESS: dict[str, Any] = {
    "C": 95.2, "S": 93.1, "D": 2.1, "F": 2.3, "M": 2.5, "n": 124,
    "score_string": "C:95.2%[S:93.1%,D:2.1%],F:2.3%,M:2.5%,n:124",
    "lineage": "bacteria_odb12",
    "mode": "genome",
}
DEMO_LINEAGE = "bacteria_odb12"
DEMO_MODE = "genome"

DEMO_FULL_TABLE_ROWS: list[tuple[str, str, str, float, int]] = [
    ("1098at2", "Complete",    "seq1", 742.3, 312),
    ("1099at2", "Complete",    "seq1", 698.1, 287),
    ("1100at2", "Complete",    "seq2", 812.5, 445),
    ("1101at2", "Duplicated",  "seq3", 756.2, 318),
    ("1102at2", "Complete",    "seq1", 623.9, 201),
    ("1103at2", "Fragmented",  "seq2", 341.2,  98),
    ("1104at2", "Missing",     "N/A",    0.0,   0),
    ("1105at2", "Complete",    "seq4", 887.4, 521),
    ("1106at2", "Complete",    "seq5", 701.3, 334),
    ("1107at2", "Missing",     "N/A",    0.0,   0),
]

_SCORE_RE = re.compile(
    r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)"
)

# ---------------------------------------------------------------------------
# NCBI E-utilities + Ensembl live-demo constants
# ---------------------------------------------------------------------------

NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_REQUEST_TIMEOUT = 15  # seconds
NCBI_USER_AGENT = "ClawBio/1.0 (busco-assessor; contact: Ideatharit@gmail.com)"
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")

MAX_DOWNLOAD_BYTES = 200 * 1024 * 1024  # 200 MB hard cap for demo downloads

# Ensembl Genomes: S. cerevisiae mitochondrial chromosome (~22 KB gzipped)
DEMO_LIVE_URL = (
    "http://ftp.ensemblgenomes.org/pub/fungi/release-62/fasta/"
    "saccharomyces_cerevisiae/dna/"
    "Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.Mito.fa.gz"
)
DEMO_LIVE_ORGANISM = "Saccharomyces cerevisiae"

# Realistic BUSCO scores for a mitochondrial-only genome vs saccharomycetes_odb10 (n=2137).
# The mito chromosome encodes only ~15–35 protein-coding genes; most BUSCO orthologs are nuclear.
# Completeness is intentionally low and demonstrates the educational value of BUSCO on partial genomes.
DEMO_LIVE_COMPLETENESS: dict[str, Any] = {
    "C": 2.1, "S": 2.1, "D": 0.0, "F": 0.9, "M": 97.0, "n": 2137,
    "score_string": "C:2.1%[S:2.1%,D:0.0%],F:0.9%,M:97.0%,n:2137",
    "lineage": "saccharomycetes_odb10",
    "mode": "genome",
}

DEMO_LIVE_FULL_TABLE_ROWS: list[tuple[str, str, str, float, int]] = [
    # A handful of conserved mitochondrial protein-coding gene orthologs
    ("1000at4891", "Complete",   "Mito_seq", 812.5, 445),
    ("1001at4891", "Complete",   "Mito_seq", 756.2, 318),
    ("1002at4891", "Complete",   "Mito_seq", 698.1, 287),
    ("1003at4891", "Complete",   "Mito_seq", 623.9, 201),
    ("1004at4891", "Fragmented", "Mito_seq", 341.2,  98),
    ("1005at4891", "Fragmented", "Mito_seq", 312.8,  82),
    ("1006at4891", "Missing",    "N/A",        0.0,   0),
    ("1007at4891", "Missing",    "N/A",        0.0,   0),
    ("1008at4891", "Missing",    "N/A",        0.0,   0),
    ("1009at4891", "Missing",    "N/A",        0.0,   0),
]

# NCBI lineage rank+name -> BUSCO lineage flag (most-specific first — first match wins)
NCBI_TO_BUSCO: list[tuple[str, str, str, str | None]] = [
    # (rank, name_lower, busco_flag, busco_lineage_value_or_None)
    # --- Fungi (class-level most specific) ---
    ("class",       "saccharomycetes",  "--lineage", "saccharomycetes_odb10"),
    ("class",       "eurotiomycetes",   "--lineage", "eurotiomycetes_odb10"),
    ("class",       "sordariomycetes",  "--lineage", "sordariomycetes_odb10"),
    ("class",       "dothideomycetes",  "--lineage", "dothideomycetes_odb10"),
    ("class",       "agaricomycetes",   "--lineage", "agaricomycetes_odb10"),
    ("phylum",      "ascomycota",       "--lineage", "ascomycota_odb10"),
    ("phylum",      "basidiomycota",    "--lineage", "basidiomycota_odb10"),
    ("kingdom",     "fungi",            "--lineage", "fungi_odb10"),
    # --- Metazoa ---
    ("order",       "primates",         "--lineage", "primates_odb10"),
    ("class",       "mammalia",         "--lineage", "mammalia_odb10"),
    ("class",       "aves",             "--lineage", "aves_odb10"),
    ("order",       "diptera",          "--lineage", "diptera_odb10"),
    ("class",       "insecta",          "--lineage", "insecta_odb10"),
    ("phylum",      "arthropoda",       "--lineage", "arthropoda_odb10"),
    ("class",       "actinopteri",      "--lineage", "actinopterygii_odb10"),
    ("class",       "actinopterygii",   "--lineage", "actinopterygii_odb10"),
    ("subphylum",   "vertebrata",       "--lineage", "vertebrata_odb10"),
    ("phylum",      "chordata",         "--lineage", "vertebrata_odb10"),
    ("phylum",      "nematoda",         "--lineage", "nematoda_odb10"),
    ("phylum",      "metazoa",          "--lineage", "metazoa_odb10"),
    # --- Plants ---
    ("class",       "liliopsida",       "--lineage", "liliopsida_odb10"),
    ("phylum",      "embryophyta",      "--lineage", "embryophyta_odb10"),
    ("kingdom",     "viridiplantae",    "--lineage", "viridiplantae_odb10"),
    # --- Prokaryotes ---
    ("superkingdom","bacteria",         "--auto-lineage-prok", None),
    ("superkingdom","archaea",          "--lineage", "archaea_odb12"),
    # --- Broad eukaryote fallback ---
    ("superkingdom","eukaryota",        "--auto-lineage-euk",  None),
]

# ---------------------------------------------------------------------------
# Lineage inference
# ---------------------------------------------------------------------------


def _kw_matches(kw: str, text: str) -> bool:
    """Return True if *kw* appears in *text* as a whole word (or phrase).

    Multi-word keywords (e.g. "e. coli", "fruit fly") use plain substring
    matching. Single-word keywords use a word-boundary regex so that "rat"
    does not match "rational" and "bacteria" does not match "bacteriophage".
    """
    if " " in kw:
        return kw in text
    return bool(re.search(r"\b" + re.escape(kw) + r"\b", text))


def infer_lineage(organism_hint: str) -> tuple[str, str | None]:
    """Map a free-text organism description to a BUSCO lineage flag + value.

    Returns (flag, lineage_value) where flag is one of:
      --lineage, --auto-lineage-prok, --auto-lineage-euk, --auto-lineage
    and lineage_value is the dataset name (or None for auto flags).
    """
    hint_lower = organism_hint.lower()
    for keywords, flag, value in LINEAGE_ROUTING:
        if any(_kw_matches(kw, hint_lower) for kw in keywords):
            return flag, value
    return "--auto-lineage", None


# ---------------------------------------------------------------------------
# Tool check
# ---------------------------------------------------------------------------


def check_busco() -> bool:
    return shutil.which("busco") is not None


# ---------------------------------------------------------------------------
# Argument validation
# ---------------------------------------------------------------------------


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and not getattr(args, "demo_live", False) and not args.input:
        raise ValueError("--input is required unless --demo is used.")

    # Detect conflicting lineage flags
    explicit_flags = sum([
        bool(args.lineage),
        bool(args.auto_lineage),
        bool(args.auto_lineage_euk),
        bool(args.auto_lineage_prok),
    ])
    if explicit_flags > 1:
        raise ValueError(
            "Conflicting lineage flags: use only one of --lineage, --auto-lineage, "
            "--auto-lineage-euk, or --auto-lineage-prok."
        )

    # Warn if proteins mode but input looks like nucleotide FASTA
    if not args.demo and args.mode == "proteins" and args.input:
        if str(args.input).endswith((".fna", ".fa", ".fasta")):
            print(
                "WARNING: --mode proteins selected but input file has a nucleotide "
                "extension (.fna/.fa/.fasta). BUSCO proteins mode expects amino-acid "
                "FASTA (.faa). This may produce zero hits.",
                file=sys.stderr,
            )


# ---------------------------------------------------------------------------
# BUSCO command builder
# ---------------------------------------------------------------------------


def build_busco_command(
    args: argparse.Namespace,
    lineage_flag: str,
    lineage_value: str | None,
    run_dir: Path,
) -> list[str]:
    """Assemble the BUSCO CLI command list."""
    cmd = [
        "busco",
        "-i", str(args.input),
        "-m", args.mode,
        "-c", str(args.cpu),
        "--out-path", str(run_dir.parent),
        "--out", run_dir.name,
        "-f",
    ]
    if lineage_flag == "--lineage" and lineage_value:
        cmd += ["-l", lineage_value]
    elif lineage_flag in ("--auto-lineage", "--auto-lineage-euk", "--auto-lineage-prok"):
        cmd.append(lineage_flag)

    if args.augustus:
        cmd.append("--augustus")
    if args.download_path:
        cmd += ["--download_path", str(args.download_path)]
    return cmd


# ---------------------------------------------------------------------------
# Run BUSCO subprocess
# ---------------------------------------------------------------------------


def run_busco(cmd: list[str]) -> subprocess.CompletedProcess:  # type: ignore[type-arg]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=7200,
    )
    if result.returncode != 0:
        tail = "\n".join(result.stderr.splitlines()[-10:])
        raise RuntimeError(
            f"BUSCO exited with code {result.returncode}.\nLast stderr lines:\n{tail}"
        )
    return result


# ---------------------------------------------------------------------------
# Output parsers
# ---------------------------------------------------------------------------


def find_short_summary(busco_out_dir: Path) -> Path | None:
    for candidate in (
        busco_out_dir / "short_summary.txt",
        *busco_out_dir.glob("short_summary.specific.*.txt"),
    ):
        if candidate.exists():
            return candidate
    return None


def parse_short_summary(summary_path: Path) -> dict[str, Any]:
    if not summary_path.exists():
        raise FileNotFoundError(f"short_summary.txt not found: {summary_path}")
    text = summary_path.read_text(encoding="utf-8")
    m = _SCORE_RE.search(text)
    if not m:
        raise ValueError(f"Could not parse BUSCO score string from {summary_path}")
    c, s, d, f, miss, n = m.groups()
    score_str = f"C:{c}%[S:{s}%,D:{d}%],F:{f}%,M:{miss}%,n:{n}"

    # Extract lineage and mode from comment header lines
    lineage, mode = "unknown", "unknown"
    for line in text.splitlines():
        if "lineage dataset is:" in line:
            parts = line.split("lineage dataset is:")
            lineage = parts[1].strip().split()[0] if len(parts) > 1 else "unknown"
        if "BUSCO was run in mode:" in line:
            mode = line.split("mode:")[-1].strip()

    return {
        "C": float(c), "S": float(s), "D": float(d),
        "F": float(f), "M": float(miss), "n": int(n),
        "score_string": score_str,
        "lineage": lineage,
        "mode": mode,
    }


def parse_full_table(full_table_path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for line in full_table_path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        rows.append({
            "busco_id": parts[0],
            "status": parts[1],
            "sequence": parts[2],
            "score": float(parts[3]) if len(parts) > 3 else 0.0,
            "length": int(parts[4]) if len(parts) > 4 else 0,
        })
    return rows


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------


def write_result_json(
    output_dir: Path,
    scores: dict[str, Any],
    args_dict: dict[str, Any],
) -> Path:
    data = {**scores, "run_parameters": args_dict}
    path = output_dir / "result.json"
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return path


def write_report(
    output_dir: Path,
    scores: dict[str, Any],
    full_table: list[dict[str, Any]],
    is_demo: bool,
    args: argparse.Namespace,
    note: str | None = None,
) -> Path:
    timestamp = datetime.datetime.now(datetime.timezone.utc).strftime(
        "%Y-%m-%d %H:%M UTC"
    )
    demo_tag = " (demo)" if is_demo else ""
    lineage_used = scores.get("lineage", "unknown")
    mode_used = scores.get("mode", getattr(args, "mode", "genome"))

    lines: list[str] = [
        "# BUSCO Assessor Report",
        "",
        f"**Date**: {timestamp}",
        f"**Mode**: {mode_used}{demo_tag}",
        f"**Lineage**: {lineage_used}",
        f"**Input**: {getattr(args, 'input', 'demo_assembly.fna') or 'demo_assembly.fna'}",
        "",
        "## Completeness Summary",
        "",
        "| Metric | Count | Percentage |",
        "|--------|-------|-----------|",
        f"| Complete (C) | {round(scores['C'] * scores['n'] / 100):.0f} | {scores['C']}% |",
        f"|   Single-copy (S) | {round(scores['S'] * scores['n'] / 100):.0f} | {scores['S']}% |",
        f"|   Duplicated (D) | {round(scores['D'] * scores['n'] / 100):.0f} | {scores['D']}% |",
        f"| Fragmented (F) | {round(scores['F'] * scores['n'] / 100):.0f} | {scores['F']}% |",
        f"| Missing (M) | {round(scores['M'] * scores['n'] / 100):.0f} | {scores['M']}% |",
        f"| Total searched (n) | {scores['n']} | — |",
        "",
        f"**Score string**: `{scores['score_string']}`",
        "",
        "## Interpretation",
        "",
    ]

    c_pct = scores["C"]
    if c_pct >= 95:
        lines.append(
            f"High completeness ({c_pct}% C) indicates a near-complete assembly for this lineage."
        )
    elif c_pct >= 80:
        lines.append(
            f"Moderate completeness ({c_pct}% C). The assembly may have fragmented or missing regions."
        )
    else:
        lines.append(
            f"Low completeness ({c_pct}% C). Consider reassembling or verifying the input FASTA."
        )

    d_pct = scores["D"]
    if d_pct > 5:
        lines.append(
            f"Duplication rate of {d_pct}% exceeds expected range — check for assembly redundancy."
        )
    else:
        lines.append(f"Duplication rate of {d_pct}% is within expected range.")

    if note:
        lines += ["", f"> {note}", ""]

    if full_table:
        lines += [
            "",
            "## Top Gene Results (first 10)",
            "",
            "| BUSCO ID | Status | Sequence | Score | Length |",
            "|----------|--------|----------|-------|--------|",
        ]
        for row in full_table[:10]:
            lines.append(
                f"| {row['busco_id']} | {row['status']} | {row['sequence']} "
                f"| {row['score']} | {row['length']} |"
            )

    lines += [
        "",
        "---",
        "",
        f"*{DISCLAIMER}*",
        "",
    ]

    path = output_dir / "report.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def write_reproducibility_bundle(
    output_dir: Path,
    cmd: list[str],
    args: argparse.Namespace,
    input_path: Path | None,
    is_demo: bool,
) -> None:
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)

    # commands.sh
    demo_flag = " --demo" if is_demo else ""
    input_flag = f" --input {input_path}" if input_path else ""
    mode_val = getattr(args, "mode", "genome")
    output_val = getattr(args, "output", str(output_dir))
    (repro_dir / "commands.sh").write_text(
        "#!/usr/bin/env bash\n"
        "# Reproduce this BUSCO assessment:\n"
        f"python skills/busco-assessor/busco_assessor.py"
        f"{input_flag} --mode {mode_val} --output {output_val}{demo_flag}\n"
        "\n"
        "# Raw BUSCO command (real runs only):\n"
        f"# {' '.join(str(x) for x in cmd)}\n",
        encoding="utf-8",
    )

    # environment.yml
    (repro_dir / "environment.yml").write_text(
        "name: busco_env\n"
        "channels:\n"
        "  - conda-forge\n"
        "  - bioconda\n"
        "dependencies:\n"
        "  - python>=3.10\n"
        "  - busco=6.0.0\n"
        "  - sepp=4.5.5\n",
        encoding="utf-8",
    )

    # checksums.sha256
    checksums: list[str] = []
    for f in sorted(output_dir.rglob("*")):
        if f.is_file() and "reproducibility" not in str(f):
            digest = hashlib.sha256(f.read_bytes()).hexdigest()
            rel = f.relative_to(output_dir)
            checksums.append(f"{digest}  {rel}")
    (repro_dir / "checksums.sha256").write_text(
        "\n".join(checksums) + "\n", encoding="utf-8"
    )


# ---------------------------------------------------------------------------
# Demo data generators
# ---------------------------------------------------------------------------

_DEMO_SEQS = [
    ("seq1", "ATGCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
    ("seq2", "GCTAGCTAGCTAGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"),
    ("seq3", "CGATCGATCGATCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"),
    ("seq4", "TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"),
    ("seq5", "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"),
]


def make_demo_fasta(dest: Path) -> Path:
    dest.mkdir(parents=True, exist_ok=True)
    fasta_path = dest / "demo_assembly.fna"
    lines = []
    for header, seq in _DEMO_SEQS:
        lines.append(f">{header} demo synthetic sequence")
        lines.append(seq)
    fasta_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return fasta_path


def make_demo_busco_outputs(
    busco_run_dir: Path,
    scores: dict[str, Any],
    full_table_rows: list[tuple[str, str, str, float, int]],
) -> None:
    busco_run_dir.mkdir(parents=True, exist_ok=True)

    # short_summary.txt
    lineage = scores.get("lineage", DEMO_LINEAGE)
    mode = scores.get("mode", DEMO_MODE)
    c_count = round(scores["C"] * scores["n"] / 100)
    s_count = round(scores["S"] * scores["n"] / 100)
    d_count = round(scores["D"] * scores["n"] / 100)
    f_count = round(scores["F"] * scores["n"] / 100)
    m_count = round(scores["M"] * scores["n"] / 100)
    score_str = (
        f"C:{scores['C']}%[S:{scores['S']}%,D:{scores['D']}%],"
        f"F:{scores['F']}%,M:{scores['M']}%,n:{scores['n']}"
    )

    summary_text = (
        f"# BUSCO version is: 6.0.0\n"
        f"# The lineage dataset is: {lineage} "
        f"(Creation date: 2021-09-02, number of species: 3857, number of BUSCOs: {scores['n']})\n"
        f"# Summarized benchmarking in ALL lineages: {busco_run_dir}\n"
        f"# BUSCO was run in mode: {mode}\n"
        f"# Gene predictor used: prodigal\n"
        f"\n"
        f"        --------------------------------------------------\n"
        f"        |Results from dataset {lineage:<29}|\n"
        f"        --------------------------------------------------\n"
        f"        |{score_str:<49}|\n"
        f"        |{c_count:<7}Complete BUSCOs (C)                       |\n"
        f"        |{s_count:<7}Complete and single-copy BUSCOs (S)       |\n"
        f"        |{d_count:<7}Complete and duplicated BUSCOs (D)        |\n"
        f"        |{f_count:<7}Fragmented BUSCOs (F)                     |\n"
        f"        |{m_count:<7}Missing BUSCOs (M)                        |\n"
        f"        |{scores['n']:<7}Total BUSCO groups searched               |\n"
        f"        --------------------------------------------------\n"
    )
    (busco_run_dir / "short_summary.txt").write_text(summary_text, encoding="utf-8")

    # short_summary.json
    summary_json = {
        "one_line_summary": score_str,
        "dataset": lineage,
        "mode": mode,
        "Complete": c_count,
        "Single copy": s_count,
        "Multi copy": d_count,
        "Fragmented": f_count,
        "Missing": m_count,
        "n_markers": scores["n"],
    }
    (busco_run_dir / "short_summary.json").write_text(
        json.dumps(summary_json, indent=2), encoding="utf-8"
    )

    # full_table.tsv
    tsv_lines = [
        "# Busco id\tStatus\tSequence\tScore\tLength",
    ]
    for row in full_table_rows:
        tsv_lines.append(
            f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}"
        )
    (busco_run_dir / "full_table.tsv").write_text(
        "\n".join(tsv_lines) + "\n", encoding="utf-8"
    )


# ---------------------------------------------------------------------------
# Pipeline orchestration
# ---------------------------------------------------------------------------


def run_demo(args: argparse.Namespace) -> None:
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    busco_run_dir = output_dir / "busco_run"
    make_demo_fasta(output_dir)
    make_demo_busco_outputs(busco_run_dir, DEMO_COMPLETENESS, DEMO_FULL_TABLE_ROWS)

    scores = parse_short_summary(busco_run_dir / "short_summary.txt")
    full_table = parse_full_table(busco_run_dir / "full_table.tsv")

    # Patch args for report generation
    args_for_report = argparse.Namespace(**vars(args))
    args_for_report.input = str(output_dir / "demo_assembly.fna")
    args_for_report.mode = DEMO_MODE

    write_result_json(output_dir, scores, vars(args_for_report))
    write_report(output_dir, scores, full_table, is_demo=True, args=args_for_report)
    write_reproducibility_bundle(
        output_dir,
        cmd=["busco", "-i", "demo_assembly.fna", "-m", DEMO_MODE, "-l", DEMO_LINEAGE],
        args=args_for_report,
        input_path=None,
        is_demo=True,
    )
    print(f"Demo report written to: {output_dir / 'report.md'}")


def run_pipeline(args: argparse.Namespace) -> None:
    if not check_busco():
        print(
            "ERROR: 'busco' binary not found on PATH. "
            "Install via: conda install -c bioconda -c conda-forge busco=6.0.0 sepp=4.5.5",
            file=sys.stderr,
        )
        sys.exit(1)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    busco_run_dir = output_dir / "busco_run"

    # Resolve lineage
    if args.lineage:
        lineage_flag, lineage_value = "--lineage", args.lineage
    elif args.auto_lineage:
        lineage_flag, lineage_value = "--auto-lineage", None
    elif args.auto_lineage_euk:
        lineage_flag, lineage_value = "--auto-lineage-euk", None
    elif args.auto_lineage_prok:
        lineage_flag, lineage_value = "--auto-lineage-prok", None
    elif args.organism:
        lineage_flag, lineage_value = infer_lineage(args.organism)
        print(f"Inferred lineage from organism hint '{args.organism}': {lineage_flag} {lineage_value or ''}")
    else:
        lineage_flag, lineage_value = "--auto-lineage", None
        print("No lineage specified; using --auto-lineage (requires SEPP 4.5.5).")

    cmd = build_busco_command(args, lineage_flag, lineage_value, busco_run_dir)
    print(f"Running: {' '.join(cmd)}")
    run_busco(cmd)

    summary_path = find_short_summary(busco_run_dir)
    if not summary_path:
        print("ERROR: BUSCO short_summary.txt not found after run.", file=sys.stderr)
        sys.exit(1)

    scores = parse_short_summary(summary_path)
    full_table_path = busco_run_dir / "full_table.tsv"
    full_table = parse_full_table(full_table_path) if full_table_path.exists() else []

    write_result_json(output_dir, scores, vars(args))
    write_report(output_dir, scores, full_table, is_demo=False, args=args)
    write_reproducibility_bundle(
        output_dir,
        cmd=cmd,
        args=args,
        input_path=Path(args.input),
        is_demo=False,
    )
    print(f"Report written to: {output_dir / 'report.md'}")


# ---------------------------------------------------------------------------
# NCBI Taxonomy API
# ---------------------------------------------------------------------------


def ncbi_taxonomy_lookup(organism_name: str) -> list[dict[str, str]]:
    """Query NCBI E-utilities for organism lineage.

    Returns a list of {rank, name} dicts from LineageEx (skips 'no rank' entries).
    Returns [] on any network or parse error so callers can fall back gracefully.
    """
    try:
        # Step 1: esearch — get the taxid
        query = urllib.parse.quote(organism_name)
        key_param = f"&api_key={NCBI_API_KEY}" if NCBI_API_KEY else ""
        search_url = (
            f"{NCBI_EUTILS_BASE}/esearch.fcgi"
            f"?db=taxonomy&term={query}&retmode=json{key_param}"
        )
        search_req = urllib.request.Request(
            search_url, headers={"User-Agent": NCBI_USER_AGENT}
        )
        with urllib.request.urlopen(search_req, timeout=NCBI_REQUEST_TIMEOUT) as resp:
            data = json.loads(resp.read())
        ids = data.get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []
        taxid = ids[0]

        # Step 2: efetch — get full lineage XML
        fetch_url = (
            f"{NCBI_EUTILS_BASE}/efetch.fcgi"
            f"?db=taxonomy&id={taxid}&retmode=xml{key_param}"
        )
        fetch_req = urllib.request.Request(
            fetch_url, headers={"User-Agent": NCBI_USER_AGENT}
        )
        with urllib.request.urlopen(fetch_req, timeout=NCBI_REQUEST_TIMEOUT) as resp:
            xml_bytes = resp.read()

        root = ET.fromstring(xml_bytes)
        taxon = root.find("Taxon")
        if taxon is None:
            return []
        lineage_ex = taxon.find("LineageEx")
        if lineage_ex is None:
            return []

        lineage: list[dict[str, str]] = []
        for t in lineage_ex.findall("Taxon"):
            rank = (t.findtext("Rank") or "no rank").strip()
            name = (t.findtext("ScientificName") or "").strip()
            if rank != "no rank" and name:
                lineage.append({"rank": rank, "name": name})
        return lineage

    except (urllib.error.URLError, urllib.error.HTTPError, ET.ParseError,
            json.JSONDecodeError, KeyError, OSError):
        return []


def ncbi_lineage_to_busco(lineage: list[dict[str, str]]) -> tuple[str, str | None]:
    """Map an NCBI LineageEx list to a BUSCO lineage flag.

    Iterates NCBI_TO_BUSCO from most-specific to least-specific rank.
    Returns ("--auto-lineage", None) if nothing matches.
    """
    lookup = {(e["rank"].lower(), e["name"].lower()) for e in lineage}
    for rank, name, flag, value in NCBI_TO_BUSCO:
        if (rank.lower(), name.lower()) in lookup:
            return flag, value
    return "--auto-lineage", None


# ---------------------------------------------------------------------------
# FASTA download helper
# ---------------------------------------------------------------------------


def download_and_decompress(url: str, dest_dir: Path) -> Path:
    """Download a .fa.gz (or .fna.gz) from *url* and decompress it into *dest_dir*.

    Returns the path to the decompressed FASTA file.
    Removes the intermediate .gz file after successful decompression.
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    filename = url.split("/")[-1]
    gz_path = dest_dir / filename
    # Strip .gz suffix to get the FASTA path
    fa_name = filename[:-3] if filename.endswith(".gz") else filename
    fa_path = dest_dir / fa_name

    print(f"  Downloading {filename} ({url}) ...")
    total = 0
    with urllib.request.urlopen(url, timeout=NCBI_REQUEST_TIMEOUT) as resp, \
            open(gz_path, "wb") as fh:
        while True:
            chunk = resp.read(65536)
            if not chunk:
                break
            total += len(chunk)
            if total > MAX_DOWNLOAD_BYTES:
                fh.close()
                gz_path.unlink(missing_ok=True)
                raise RuntimeError(
                    f"Download aborted: response exceeded {MAX_DOWNLOAD_BYTES // (1024 ** 2)} MB cap"
                )
            fh.write(chunk)

    print(f"  Decompressing {gz_path.name} ...")
    with gzip.open(gz_path, "rb") as f_in, open(fa_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    gz_path.unlink()
    return fa_path


# ---------------------------------------------------------------------------
# Live demo orchestration
# ---------------------------------------------------------------------------


def run_demo_live(args: argparse.Namespace) -> None:
    """Full demo using real S. cerevisiae Mito FASTA from Ensembl Genomes.

    Lineage is resolved via NCBI Taxonomy API (falls back to keyword matching).
    If BUSCO binary is installed, runs BUSCO for real. Otherwise generates
    realistic synthetic output that reflects actual mito-genome completeness.
    """
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Step 1: Download real FASTA ──────────────────────────────────────
    print(f"\n[1/5] Downloading S. cerevisiae Mito FASTA from Ensembl Genomes...")
    fasta_path = download_and_decompress(DEMO_LIVE_URL, output_dir / "input")
    seq_count = sum(1 for ln in fasta_path.read_text(encoding="utf-8").splitlines()
                    if ln.startswith(">"))
    print(f"       -> {fasta_path.name}  ({seq_count} sequence(s))")

    # ── Step 2: Resolve lineage via NCBI Taxonomy ────────────────────────
    print(f"\n[2/5] Resolving lineage for '{DEMO_LIVE_ORGANISM}' via NCBI Taxonomy...")
    ncbi_lineage = ncbi_taxonomy_lookup(DEMO_LIVE_ORGANISM)
    if ncbi_lineage:
        lineage_flag, lineage_value = ncbi_lineage_to_busco(ncbi_lineage)
        source = "NCBI Taxonomy API"
    else:
        print("       NCBI lookup failed or timed out — falling back to keyword routing.")
        lineage_flag, lineage_value = infer_lineage(DEMO_LIVE_ORGANISM)
        source = "keyword fallback"
    print(f"       -> {lineage_flag} {lineage_value or '(auto)'} [{source}]")

    busco_run_dir = output_dir / "busco_run"

    # ── Step 3: Run BUSCO (if installed) or generate synthetic output ────
    if check_busco():
        print(f"\n[3/5] Running BUSCO (real data)...")
        args_for_run = argparse.Namespace(**vars(args))
        args_for_run.input = str(fasta_path)
        args_for_run.mode = "genome"
        cmd = build_busco_command(args_for_run, lineage_flag, lineage_value, busco_run_dir)
        print(f"       -> {' '.join(cmd)}")
        run_busco(cmd)
        summary_path = find_short_summary(busco_run_dir)
        if not summary_path:
            raise RuntimeError("BUSCO short_summary.txt not found after live run.")
        scores = parse_short_summary(summary_path)
        full_table = (
            parse_full_table(busco_run_dir / "full_table.tsv")
            if (busco_run_dir / "full_table.tsv").exists() else []
        )
        cmd_for_repro = cmd
    else:
        print(
            "\n[3/5] BUSCO binary not found - generating realistic synthetic output.\n"
            "       (Install via: conda install -c bioconda -c conda-forge busco=6.0.0 sepp=4.5.5)"
        )
        live_scores = dict(DEMO_LIVE_COMPLETENESS)
        if lineage_value:
            live_scores["lineage"] = lineage_value
        make_demo_busco_outputs(busco_run_dir, live_scores, DEMO_LIVE_FULL_TABLE_ROWS)
        scores = live_scores
        full_table = parse_full_table(busco_run_dir / "full_table.tsv")
        cmd_for_repro = [
            "busco", "-i", str(fasta_path), "-m", "genome",
            *([ "-l", lineage_value] if lineage_value else [lineage_flag]),
            "-c", str(getattr(args, "cpu", DEFAULT_CPU)),
        ]

    # ── Step 4: Write outputs ─────────────────────────────────────────────
    print("\n[4/5] Writing report and result.json...")
    args_for_report = argparse.Namespace(**vars(args))
    args_for_report.input = str(fasta_path)
    args_for_report.mode = "genome"

    mito_note = (
        "**Note — Mitochondrial genome only**: "
        "This demo uses the *S. cerevisiae* mitochondrial chromosome (~85 kb). "
        f"BUSCO completeness is expected to be very low (~2%) against "
        f"`saccharomycetes_odb10` (n=2137) because most orthologs are nuclear genes. "
        "For a full-genome assessment, provide the complete genome FASTA."
    )
    write_result_json(output_dir, scores, vars(args_for_report))
    write_report(
        output_dir, scores, full_table,
        is_demo=True, args=args_for_report,
        note=mito_note,
    )

    # ── Step 5: Reproducibility bundle ──────────────────────────────────
    print("\n[5/5] Writing reproducibility bundle...")
    write_reproducibility_bundle(
        output_dir,
        cmd=cmd_for_repro,
        args=args_for_report,
        input_path=fasta_path,
        is_demo=True,
    )
    print(f"\nLive demo report: {output_dir / 'report.md'}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="BUSCO genome/transcriptome/protein completeness assessor for ClawBio",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    io_group = parser.add_argument_group("Input / Output")
    io_group.add_argument("--input", type=str, default=None,
                          help="Input FASTA file (assembly, transcriptome, or proteins)")
    io_group.add_argument("--output", type=str, required=True,
                          help="Output directory")
    io_group.add_argument("--demo", action="store_true",
                          help="Run with synthetic demo data (no BUSCO binary required)")
    io_group.add_argument("--demo-live", dest="demo_live", action="store_true",
                          help=(
                              "Download real S. cerevisiae Mito FASTA from Ensembl Genomes, "
                              "resolve lineage via NCBI Taxonomy, and run BUSCO (or generate "
                              "realistic synthetic output if BUSCO is not installed)"
                          ))

    mode_group = parser.add_argument_group("BUSCO mode")
    mode_group.add_argument("--mode", choices=BUSCO_MODES, default="genome",
                            help="Assessment mode")
    mode_group.add_argument("--augustus", action="store_true",
                            help="Use Augustus gene predictor (eukaryote genome only)")

    lineage_group = parser.add_mutually_exclusive_group()
    lineage_group.add_argument("--lineage", type=str, default=None,
                               help="BUSCO lineage dataset (e.g. bacteria_odb12)")
    lineage_group.add_argument("--auto-lineage", dest="auto_lineage",
                               action="store_true",
                               help="Auto-select lineage across all domains (requires SEPP 4.5.5)")
    lineage_group.add_argument("--auto-lineage-euk", dest="auto_lineage_euk",
                               action="store_true",
                               help="Auto-select lineage within Eukaryota")
    lineage_group.add_argument("--auto-lineage-prok", dest="auto_lineage_prok",
                               action="store_true",
                               help="Auto-select lineage within Prokaryota")

    runtime_group = parser.add_argument_group("Runtime")
    runtime_group.add_argument("--cpu", type=int, default=DEFAULT_CPU,
                               help="Number of CPU threads")
    runtime_group.add_argument("--download-path", dest="download_path", type=str,
                               default=None, help="Custom BUSCO lineage dataset directory")
    runtime_group.add_argument("--organism", type=str, default=None,
                               help="Free-text organism description for agentic lineage inference "
                                    "(e.g. 'E. coli K-12', 'fruit fly', 'human hg38')")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    try:
        validate_args(args)
        if args.demo:
            run_demo(args)
        elif getattr(args, "demo_live", False):
            run_demo_live(args)
        else:
            run_pipeline(args)
    except (ValueError, FileNotFoundError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
