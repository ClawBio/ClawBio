#!/usr/bin/env python3
"""
fastreer.py: Phylogenetic trees and distance matrices from VCF/FASTA via fastreeR.

Wraps the fastreeR hybrid Java/Python toolkit (https://github.com/gkanogiannis/fastreeR).
Requires: pip install fastreer   AND   Java 11+

Supported input formats:
  VCF file (biallelic/multiallelic SNPs, compressed or plain),
  FASTA file (aligned or unaligned sequences, compressed or plain), or
  PHYLIP distance matrix (for DIST2TREE)

Commands:
  VCF2TREE   - Newick tree from VCF (hierarchical clustering)
  VCF2DIST   - PHYLIP distance matrix from VCF (cosine dissimilarity)
  DIST2TREE  - Newick tree from PHYLIP distance matrix
  FASTA2DIST - D2S k-mer distance matrix from FASTA sequences

Usage:
  python fastreer.py --command VCF2TREE --input samples.vcf.gz --output /tmp/out
  python fastreer.py --command VCF2DIST --input samples.vcf.gz --threads 4 --output /tmp/out
  python fastreer.py --command FASTA2DIST --input seqs.fasta --kmer 5 --output /tmp/out
  python fastreer.py --command DIST2TREE --input distances.dist --output /tmp/out
  python fastreer.py --demo --output /tmp/fastreer_demo
"""

from __future__ import annotations

__version__ = "0.1.0"
__author__ = "Anestis Gkanogiannis"
__licence__ = "GPL-3.0"

import argparse
import json
import os
import shutil
import subprocess
import sys
import textwrap
from datetime import datetime
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (  # noqa: E402
    ReproCommand,
    ReproPath,
    write_checksums,
    write_environment_yml,
    write_portable_commands_sh,
)

DISCLAIMER = (
    "ClawBio is a research and educational tool. "
    "It is not a medical device and does not provide clinical diagnoses. "
    "Consult a healthcare professional before making any medical decisions."
)

VALID_COMMANDS = {"VCF2TREE", "VCF2DIST", "DIST2TREE", "FASTA2DIST"}

# ── Synthetic demo data ───────────────────────────────────────────────────────

DEMO_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=248956422>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2\tSAMPLE3\tSAMPLE4\tSAMPLE5
chr1\t10000\trs1\tA\tG\t100\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1\t0/0
chr1\t20000\trs2\tC\tT\t100\tPASS\t.\tGT\t0/1\t1/1\t0/0\t0/1\t1/1
chr1\t30000\trs3\tG\tA\t100\tPASS\t.\tGT\t1/1\t0/1\t0/1\t0/0\t0/1
chr1\t40000\trs4\tT\tC\t100\tPASS\t.\tGT\t0/0\t0/0\t0/1\t1/1\t0/1
chr1\t50000\trs5\tA\tT\t100\tPASS\t.\tGT\t0/1\t0/1\t0/0\t0/1\t1/1
chr1\t60000\trs6\tC\tG\t100\tPASS\t.\tGT\t1/1\t0/0\t0/1\t0/1\t0/0
chr1\t70000\trs7\tG\tC\t100\tPASS\t.\tGT\t0/1\t1/1\t1/1\t0/0\t0/1
chr1\t80000\trs8\tT\tA\t100\tPASS\t.\tGT\t0/0\t0/1\t0/0\t1/1\t0/1
chr1\t90000\trs9\tA\tC\t100\tPASS\t.\tGT\t0/1\t0/0\t1/1\t0/1\t0/0
chr1\t100000\trs10\tC\tA\t100\tPASS\t.\tGT\t1/1\t0/1\t0/1\t0/0\t1/1
chr1\t110000\trs11\tG\tT\t100\tPASS\t.\tGT\t0/0\t1/1\t0/0\t0/1\t0/1
chr1\t120000\trs12\tT\tG\t100\tPASS\t.\tGT\t0/1\t0/1\t1/1\t1/1\t0/0
chr1\t130000\trs13\tA\tG\t100\tPASS\t.\tGT\t1/1\t0/0\t0/1\t0/1\t1/1
chr1\t140000\trs14\tC\tT\t100\tPASS\t.\tGT\t0/1\t1/1\t0/0\t0/0\t0/1
chr1\t150000\trs15\tG\tA\t100\tPASS\t.\tGT\t0/0\t0/1\t0/1\t1/1\t1/1
chr1\t160000\trs16\tT\tC\t100\tPASS\t.\tGT\t1/1\t0/0\t1/1\t0/1\t0/0
chr1\t170000\trs17\tA\tT\t100\tPASS\t.\tGT\t0/1\t0/1\t0/0\t1/1\t0/1
chr1\t180000\trs18\tC\tG\t100\tPASS\t.\tGT\t0/0\t1/1\t0/1\t0/0\t1/1
chr1\t190000\trs19\tG\tC\t100\tPASS\t.\tGT\t0/1\t0/1\t1/1\t0/1\t0/0
chr1\t200000\trs20\tT\tA\t100\tPASS\t.\tGT\t1/1\t0/0\t0/1\t1/1\t0/1
"""

DEMO_FASTA = """\
>SAMPLE1
ATGCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
>SAMPLE2
ATGCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAC
>SAMPLE3
ATGCGTTCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
>SAMPLE4
ATGCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCAAGCTAG
>SAMPLE5
ATGCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAA
"""

# Synthetic demo outputs when fastreer is not installed
DEMO_NEWICK = (
    "((SAMPLE1:0.120,SAMPLE2:0.098):0.045,"
    "(SAMPLE3:0.110,(SAMPLE4:0.087,SAMPLE5:0.132):0.062):0.038);"
)

DEMO_DIST_VCF = """\
5 20
SAMPLE1\t0.000\t0.120\t0.230\t0.180\t0.210
SAMPLE2\t0.120\t0.000\t0.200\t0.190\t0.220
SAMPLE3\t0.230\t0.200\t0.000\t0.150\t0.170
SAMPLE4\t0.180\t0.190\t0.150\t0.000\t0.130
SAMPLE5\t0.210\t0.220\t0.170\t0.130\t0.000
"""

DEMO_DIST_FASTA = """\
5 60
SAMPLE1\t0.000\t0.016\t0.033\t0.016\t0.016
SAMPLE2\t0.016\t0.000\t0.033\t0.016\t0.016
SAMPLE3\t0.033\t0.033\t0.000\t0.033\t0.033
SAMPLE4\t0.016\t0.016\t0.033\t0.000\t0.016
SAMPLE5\t0.016\t0.016\t0.033\t0.016\t0.000
"""


# ── Pre-flight checks ─────────────────────────────────────────────────────────

def check_java() -> str | None:
    """Return java path or None if not found / version too old."""
    java = shutil.which("java")
    if not java:
        return None
    try:
        out = subprocess.check_output(
            ["java", "-version"], stderr=subprocess.STDOUT, text=True
        )
        # version string: 'openjdk version "17.0.x" ...' or 'java version "1.8.x" ...'
        import re
        m = re.search(r'"(\d+)[\.\d]*"', out)
        if m:
            major = int(m.group(1))
            if major == 1:
                # Legacy: "1.8" → 8
                m2 = re.search(r'"1\.(\d+)', out)
                major = int(m2.group(1)) if m2 else major
            if major < 11:
                return None  # too old
    except Exception:
        return None
    return java


def find_fastreer_bin() -> str | None:
    """Return path to the fastreeR binary, or None if not found.

    Checks both 'fastreeR' and 'fastreer' (package versions differ in casing).
    Prefers the venv's bin/ co-located with sys.executable so that an arbitrary
    binary on PATH cannot be substituted for the installed package.
    Falls back to PATH only when the interpreter is not inside a venv.
    """
    venv_dir = Path(sys.executable).parent
    for name in ("fastreeR", "fastreer"):
        candidate = venv_dir / name
        if candidate.exists():
            return str(candidate)
    # PATH fallback — only reached when running outside a venv
    for name in ("fastreeR", "fastreer"):
        found = shutil.which(name)
        if found:
            return found
    return None


# ── Core runner ───────────────────────────────────────────────────────────────

def run_fastreer(
    command: str,
    input_file: Path,
    output_dir: Path,
    args,
    *,
    allow_synthetic: bool = False,
) -> dict:
    """Run fastreeR and return a dict with output paths and metadata.

    When allow_synthetic=False (the default for real input), the function exits
    with a clear error if fastreeR or Java 11+ is missing. Synthetic fallback is
    only permitted during --demo mode (allow_synthetic=True) so that simulated
    data is never silently presented as computed results.
    """
    samples = _extract_sample_names(input_file, command)
    n_variants = _count_variants(input_file, command)

    bin_path = find_fastreer_bin()
    java = check_java()

    if bin_path and java:
        return _run_real_fastreer(command, input_file, output_dir, args, samples, bin_path)

    missing = []
    if not bin_path:
        missing.append("fastreeR binary (pip install fastreer)")
    if not java:
        missing.append("Java 11+")

    if not allow_synthetic:
        sys.exit(
            f"ERROR: {' and '.join(missing)} not found. "
            "Install the missing dependencies and re-run.\n"
            "  pip install fastreer          # Python wrapper\n"
            "  sudo apt install default-jre  # Java (Debian/Ubuntu)\n"
            "  brew install openjdk@17       # Java (macOS)\n"
            "Run with --demo to see expected output without real data."
        )

    print(
        f"[fastreer] WARNING: {' and '.join(missing)} not found; "
        "generating synthetic demo output (--demo mode).",
        file=sys.stderr,
    )
    return _generate_synthetic_output(command, input_file, output_dir, samples, n_variants)


def _run_real_fastreer(
    command: str, input_file: Path, output_dir: Path, args, samples: list, bin_path: str
) -> dict:
    """Call the real fastreeR CLI: fastreeR [--mem N] COMMAND -i input -o output [flags]"""
    if command in ("VCF2TREE", "DIST2TREE"):
        out_file = output_dir / "tree.nwk"
    else:
        out_file = output_dir / "distances.dist"

    # Top-level flags come before the subcommand
    cmd = [bin_path, "--mem", str(args.mem)]

    # Subcommand and its flags
    cmd += [command, "-i", str(input_file), "-o", str(out_file)]

    if command in ("VCF2DIST", "VCF2TREE"):
        cmd += ["-t", str(args.threads)]
        if args.verbose:
            cmd += ["-v"]
        if getattr(args, "window_bp", None):
            cmd += ["--window-bp", str(args.window_bp)]
        if getattr(args, "window_variants", None):
            cmd += ["--window-variants", str(args.window_variants)]
    if command == "VCF2TREE" and args.bootstrap:
        cmd += ["-b", str(args.bootstrap)]
    if command == "FASTA2DIST":
        cmd += ["-t", str(args.threads), "-k", str(args.kmer)]
        if args.verbose:
            cmd += ["-v"]
    if command == "DIST2TREE" and args.verbose:
        cmd += ["-v"]

    raw_timeout = getattr(args, "timeout", 900)
    timeout = raw_timeout if raw_timeout > 0 else None  # 0 means no limit
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        sys.exit(
            f"ERROR: fastreeR did not finish within {raw_timeout}s. "
            "Increase --timeout or use --timeout 0 to disable the limit."
        )
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(f"[fastreer] fastreeR exited with code {result.returncode}")

    output_text = out_file.read_text() if out_file.exists() else ""
    return {
        "command": command,
        "input_file": str(input_file),
        "output_dir": str(output_dir),
        "samples": samples,
        "n_variants": _count_variants(input_file, command),
        "output_file": str(out_file),
        "output_content": output_text[:2000],
        "synthetic": False,
    }


def _generate_synthetic_output(
    command: str, input_file: Path, output_dir: Path, samples: list, n_variants: int
) -> dict:
    """Generate plausible synthetic output for demo/offline use."""
    if command in ("VCF2TREE", "DIST2TREE"):
        out_file = output_dir / "tree.nwk"
        content = DEMO_NEWICK
    else:
        out_file = output_dir / "distances.dist"
        content = DEMO_DIST_FASTA if command == "FASTA2DIST" else DEMO_DIST_VCF

    out_file.write_text(content)
    return {
        "command": command,
        "input_file": str(input_file),
        "output_dir": str(output_dir),
        "samples": samples,
        "n_variants": n_variants,
        "output_file": str(out_file),
        "output_content": content[:2000],
        "synthetic": True,
    }


def _extract_sample_names(input_file: Path, command: str) -> list[str]:
    """Parse sample names from VCF header or FASTA sequence IDs."""
    if command == "DIST2TREE":
        try:
            lines = input_file.read_text().splitlines()
            return [ln.split()[0] for ln in lines[1:] if ln.strip()]
        except Exception:
            return []

    if command == "FASTA2DIST":
        try:
            text = input_file.read_text()
            return [ln[1:].split()[0] for ln in text.splitlines() if ln.startswith(">")]
        except Exception:
            return []

    # VCF commands
    try:
        import gzip
        opener = gzip.open if str(input_file).endswith((".gz", ".bgz")) else open
        with opener(str(input_file), "rt") as fh:
            for line in fh:
                if line.startswith("#CHROM"):
                    cols = line.strip().split("\t")
                    return cols[9:] if len(cols) > 9 else []
    except Exception:
        pass
    return []


def _count_variants(input_file: Path, command: str) -> int:
    """Count non-header, non-empty lines as a proxy for variant count."""
    if command in ("DIST2TREE", "FASTA2DIST"):
        return 0
    try:
        import gzip
        opener = gzip.open if str(input_file).endswith((".gz", ".bgz")) else open
        count = 0
        with opener(str(input_file), "rt") as fh:
            for line in fh:
                if not line.startswith("#") and line.strip():
                    count += 1
        return count
    except Exception:
        return 0


# ── Report generation ─────────────────────────────────────────────────────────

def write_report(result: dict, output_dir: Path) -> None:
    command = result["command"]
    samples = result["samples"]
    n_samples = len(samples)
    n_variants = result.get("n_variants", 0)
    synthetic_note = (
        "\n> **Note**: fastreer or Java 11+ was not found. "
        "Output shown is synthetic demo data representing the expected format.\n"
        if result.get("synthetic")
        else ""
    )

    sample_list = "\n".join(f"  - {s}" for s in samples[:20])
    if len(samples) > 20:
        sample_list += f"\n  - ... ({len(samples) - 20} more)"

    if command in ("VCF2TREE", "DIST2TREE"):
        output_section = f"""## Phylogenetic Tree

**Output format**: Newick
**File**: `tree.nwk`

```newick
{result.get("output_content", "").strip()[:500]}
```

The tree was built using hierarchical clustering on a cosine distance matrix
computed from genotype data. Internal branch lengths represent evolutionary
distance; higher values indicate greater genomic divergence between samples.
"""
    else:
        dist_lines = result.get("output_content", "").strip().splitlines()
        header = dist_lines[0] if dist_lines else ""
        preview = "\n".join(dist_lines[:6]) if dist_lines else ""
        output_section = f"""## Distance Matrix

**Output format**: PHYLIP
**File**: `distances.dist`
**Header**: `{header}`

```
{preview}
```

Distances represent {"D2S k-mer" if command == "FASTA2DIST" else "cosine"} dissimilarity between samples.
Values range from 0 (identical) to 1 (maximally different).
"""

    report = f"""# fastreer Report

**Command**: `{command}`
**Input**: `{Path(result["input_file"]).name}`
**Samples**: {n_samples}
**Variants**: {n_variants if n_variants else "N/A"}
**Date**: {datetime.now().strftime("%Y-%m-%d %H:%M")}
**fastreer version**: 2.2.0
{synthetic_note}
## Samples ({n_samples})

{sample_list if sample_list else "_No samples detected_"}

{output_section}
## Interpretation

{"The Newick tree can be visualised in tools such as FigTree, iTOL, or Dendroscope. Bootstrap values (if computed) appear at internal nodes." if command in ("VCF2TREE", "DIST2TREE") else "The distance matrix can be used as input for DIST2TREE to generate a phylogenetic tree, or imported into R/Python for hierarchical clustering and visualisation."}

## Install / Requirements

```bash
pip install fastreer      # Python wrapper
# Java 11+ is also required:
sudo apt install default-jre    # Debian/Ubuntu
brew install openjdk@17          # macOS
```

## Citation

Gkanogiannis A (2016) A scalable assembly-free variable selection algorithm for
biomarker discovery from metagenomes. *BMC Bioinformatics* 17, 311.
https://doi.org/10.1186/s12859-016-1186-3

---

*{DISCLAIMER}*
"""

    (output_dir / "report.md").write_text(report)


def write_result_json(result: dict, output_dir: Path) -> None:
    payload = {
        "command": result["command"],
        "input_file": result["input_file"],
        "output_dir": result["output_dir"],
        "samples": result["samples"],
        "n_samples": len(result["samples"]),
        "n_variants": result.get("n_variants", 0),
        "output_file": result.get("output_file", ""),
        "synthetic_demo": result.get("synthetic", False),
        "fastreer_version": "2.2.0",
        "timestamp": datetime.now().isoformat(),
    }
    (output_dir / "result.json").write_text(json.dumps(payload, indent=2))


def write_reproducibility(args, input_file: Path, output_dir: Path) -> None:
    repro = output_dir / "reproducibility"
    repro.mkdir(exist_ok=True)

    cmd_args: list[str | ReproPath] = ["--command", args.command]
    if getattr(args, "demo", False):
        # The demo generates its own input inside the output dir and may fall back
        # to synthetic output. Recording that file as --input would replay on the
        # strict path, which fails precisely when the fallback was needed.
        cmd_args.append("--demo")
    else:
        input_path = Path(input_file).resolve()
        if input_path.is_relative_to(Path(output_dir).resolve()):
            anchor = "output_dir"
        elif input_path.is_relative_to(_PROJECT_ROOT):
            anchor = "repo_root"
        else:
            anchor = "auto"
        cmd_args += ["--input", ReproPath(input_path, anchor)]
    cmd_args += [
        "--output", ReproPath(Path(output_dir), "output_dir"),
        "--threads", str(args.threads),
        "--mem", str(args.mem),
    ]
    if args.bootstrap:
        cmd_args += ["--bootstrap", str(args.bootstrap)]
    if args.kmer != 4:
        cmd_args += ["--kmer", str(args.kmer)]
    if args.window_bp:
        cmd_args += ["--window-bp", str(args.window_bp)]
    if getattr(args, "window_variants", None):
        cmd_args += ["--window-variants", str(args.window_variants)]

    java_ver = "not found"
    java = shutil.which("java")
    if java:
        try:
            java_ver = subprocess.check_output(
                ["java", "-version"], stderr=subprocess.STDOUT, text=True
            ).strip().splitlines()[0]
        except Exception:
            pass

    pip_out = ""
    try:
        pip_out = subprocess.check_output(
            [sys.executable, "-m", "pip", "show", "fastreer"], text=True
        )
    except Exception:
        pip_out = "fastreer not installed"

    write_environment_yml(
        output_dir,
        env_name="clawbio-fastreer",
        pip_deps=["fastreer>=2.2.0"],
        conda_deps=["openjdk>=11"],
        python_version="3.11",
    )
    write_portable_commands_sh(
        output_dir,
        ReproCommand(
            script_path=Path("skills/fastreer/fastreer.py"),
            args=cmd_args,
            comment="Reproduce this fastreer run",
        ),
        repo_root=_PROJECT_ROOT,
    )
    write_checksums(
        [
            output_dir / "report.md",
            output_dir / "result.json",
            output_dir / "tree.nwk",
            output_dir / "distances.dist",
        ],
        output_dir,
        anchor=output_dir,
    )
    (repro / "environment.txt").write_text(
        f"# fastreer environment snapshot\n"
        f"# Generated: {datetime.now().isoformat()}\n\n"
        f"## Java\n{java_ver}\n\n"
        f"## Python package\n{pip_out}\n"
    )


# ── Demo mode ─────────────────────────────────────────────────────────────────

def run_demo(args) -> None:
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    command = args.command or "VCF2TREE"

    if command == "FASTA2DIST":
        demo_file = output_dir / "demo_sequences.fasta"
        demo_file.write_text(DEMO_FASTA)
    elif command == "DIST2TREE":
        demo_file = output_dir / "demo_distances.dist"
        demo_file.write_text(DEMO_DIST_VCF)
    else:
        demo_file = output_dir / "demo_samples.vcf"
        demo_file.write_text(DEMO_VCF)

    print(f"[fastreer] Demo mode: running {command} on synthetic data ({demo_file.name})")

    result = run_fastreer(command, demo_file, output_dir, args, allow_synthetic=True)
    write_report(result, output_dir)
    write_result_json(result, output_dir)
    write_reproducibility(args, demo_file, output_dir)

    print(f"[fastreer] Done. Output written to: {output_dir}")
    print(f"  report.md:      summary report")
    print(f"  result.json:    machine-readable metadata")
    if command in ("VCF2TREE", "DIST2TREE"):
        print(f"  tree.nwk:       Newick phylogenetic tree")
    else:
        print(f"  distances.dist: PHYLIP distance matrix")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="fastreer: phylogenetic trees and distance matrices from VCF/FASTA",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(f"""
            Commands:
              VCF2TREE   Newick tree from VCF (default)
              VCF2DIST   PHYLIP distance matrix from VCF
              DIST2TREE  Newick tree from distance matrix (PHYLIP input)
              FASTA2DIST D2S k-mer distance matrix from FASTA

            Examples:
              python fastreer.py --command VCF2TREE --input samples.vcf.gz --output /tmp/out
              python fastreer.py --command VCF2DIST --input samples.vcf.gz --threads 4 --output /tmp/out
              python fastreer.py --command FASTA2DIST --input seqs.fasta --kmer 5 --output /tmp/out
              python fastreer.py --demo --output /tmp/fastreer_demo
        """),
    )

    parser.add_argument(
        "--command",
        choices=sorted(VALID_COMMANDS),
        default="VCF2TREE",
        help="fastreer command to run (default: VCF2TREE)",
    )
    parser.add_argument("--input", help="Input file (VCF, FASTA, or PHYLIP dist); omit with --demo")
    parser.add_argument("--output", help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run with built-in synthetic data")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads (default: 1)")
    parser.add_argument("--mem", type=int, default=256, help="JVM heap size in MB (default: 256)")
    parser.add_argument(
        "--bootstrap", type=int, default=0, help="Bootstrap replicates for VCF2TREE (default: 0)"
    )
    parser.add_argument("--kmer", type=int, default=4, help="K-mer size for FASTA2DIST (default: 4)")
    parser.add_argument("--window-bp", type=int, dest="window_bp", help="Window size in bp for windowed output")
    parser.add_argument("--window-variants", type=int, dest="window_variants", help="Window size in variant count")
    parser.add_argument(
        "--timeout", type=int, default=900,
        help="Subprocess timeout in seconds (default: 900 = 15 min); 0 disables"
    )
    parser.add_argument("--verbose", action="store_true", help="Verbose output")

    args = parser.parse_args()

    if args.demo:
        run_demo(args)
        return

    if not args.input:
        parser.error("--input is required unless --demo is used")

    input_file = Path(args.input)
    if not input_file.exists():
        sys.exit(f"ERROR: input file not found: {input_file}")

    if not args.output:
        parser.error("--output is required to specify output directory")

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[fastreer] Running {args.command} on {input_file.name}")
    result = run_fastreer(args.command, input_file, output_dir, args, allow_synthetic=False)
    write_report(result, output_dir)
    write_result_json(result, output_dir)
    write_reproducibility(args, input_file, output_dir)

    print(f"[fastreer] Done. Output written to: {output_dir}")


if __name__ == "__main__":
    main()
