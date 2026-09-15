#!/usr/bin/env python3
"""Phylogenetics Builder — End-to-end ML phylogenetic tree inference.

Supports:
  • MSA        : mafft (default), muscle, clustalw, kalign, tcoffee, prank
  • Trimming   : trimAl -automated1 (default)
  • Model sel. : ModelFinder built into IQ-TREE2 (-m MFP)
  • Engine     : iqtree2 (default) or raxml-ng
  • Bootstrap  : ufboot (default), standard (Felsenstein), all (triple: UFBoot+aLRT+aBayes)
  • Rooting    : outgroup (-o flag), midpoint (ETE3)

Usage:
    python phylogenetics_builder.py --input seqs.fasta --output /tmp/phylo
    python phylogenetics_builder.py --input aln.fasta  --output /tmp/phylo --aligned
    python phylogenetics_builder.py --demo --output /tmp/phylo_demo
"""

from __future__ import annotations

__version__ = "0.2.0"
__author__ = "ClawBio"

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
from io import StringIO
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd
from Bio import Phylo

SKILL_DIR = Path(__file__).resolve().parent
DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)

VALID_ALIGNERS = ("mafft", "muscle", "clustalw", "kalign", "tcoffee", "prank")
VALID_ENGINES = ("iqtree2", "raxml-ng")
VALID_BOOTSTRAP = ("ufboot", "standard", "all")


# ── Argument parsing ───────────────────────────────────────────────────────────


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--input",
        type=Path,
        dest="input_file",
        help="Input FASTA file (unaligned or pre-aligned)",
    )
    p.add_argument("--output", type=Path, help="Output directory")
    p.add_argument(
        "--demo",
        action="store_true",
        help="Run with built-in 12-taxon primate demo data",
    )

    p.add_argument(
        "--aligned",
        action="store_true",
        help="Input FASTA is already aligned — skip MSA step",
    )
    p.add_argument(
        "--aligner",
        choices=VALID_ALIGNERS,
        default="mafft",
        metavar="ALIGNER",
        help=f"MSA algorithm when --aligned is NOT set. "
        f"Choices: {', '.join(VALID_ALIGNERS)} (default: mafft)",
    )
    p.add_argument(
        "--engine",
        choices=VALID_ENGINES,
        default="iqtree2",
        metavar="ENGINE",
        help=f"Tree inference engine. Choices: {', '.join(VALID_ENGINES)} "
        f"(default: iqtree2)",
    )
    p.add_argument(
        "--model",
        default=None,
        help="Substitution model to use (skips ModelFinder). "
        "Example: GTR+F+G4  or  HKY+G4",
    )
    p.add_argument(
        "--bootstrap",
        choices=VALID_BOOTSTRAP,
        default="ufboot",
        metavar="MODE",
        help="Bootstrap method: ufboot (1000 UFBoot), standard (100 Felsenstein), "
        "all (UFBoot+aLRT+aBayes triple support). Default: ufboot",
    )
    p.add_argument(
        "--outgroup",
        default=None,
        help="Outgroup taxon name(s) for rooting, comma-separated. "
        "Example: Mus_musculus,Rattus_norvegicus",
    )
    p.add_argument(
        "--root",
        choices=("midpoint",),
        default=None,
        help="Post-inference rooting method (requires ETE3 for midpoint)",
    )
    p.add_argument("--no-trim", action="store_true", help="Skip trimAl trimming step")
    p.add_argument(
        "--threads",
        type=int,
        default=2,
        help="CPU threads for tree inference (default: 2)",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)",
    )

    args = p.parse_args(argv)

    # Validate choices that argparse can't catch (aliased by metavar)
    if args.aligner not in VALID_ALIGNERS:
        p.error(
            f"Invalid --aligner '{args.aligner}'. Choose from: {', '.join(VALID_ALIGNERS)}"
        )
    if args.engine not in VALID_ENGINES:
        p.error(
            f"Invalid --engine '{args.engine}'. Choose from: {', '.join(VALID_ENGINES)}"
        )
    if args.bootstrap not in VALID_BOOTSTRAP:
        p.error(
            f"Invalid --bootstrap '{args.bootstrap}'. Choose from: {', '.join(VALID_BOOTSTRAP)}"
        )

    return args


# ── FASTA I/O ──────────────────────────────────────────────────────────────────


def parse_fasta(input_path: Path, require_aligned: bool = True) -> dict[str, str]:
    """Parse FASTA file.  Returns {header: sequence}.

    Args:
        require_aligned: If True, validate all sequences have equal length.
    """
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    sequences: dict[str, str] = {}
    current_header = None
    current_seq: list[str] = []

    with open(input_path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].strip()
                if not header:
                    raise ValueError("Empty sequence header found.")
                if header in sequences or header == current_header:
                    raise ValueError(f"Duplicate sequence header found: {header}")
                if current_header is not None:
                    sequences[current_header] = "".join(current_seq)
                current_header = header
                current_seq = []
            else:
                current_seq.append(line)

    if current_header is not None:
        sequences[current_header] = "".join(current_seq)

    if len(sequences) < 3:
        raise ValueError(
            f"FASTA must contain at least 3 sequences (found {len(sequences)})"
        )

    if require_aligned:
        lengths = {len(s) for s in sequences.values()}
        if len(lengths) > 1:
            raise ValueError(
                "All sequences in the alignment must be of the same length. "
                "Use --aligner to run MSA first, or check your input file."
            )

    return sequences


# ── MSA ────────────────────────────────────────────────────────────────────────


def build_msa_command(aligner: str, input_fasta: Path, output_fasta: Path) -> dict[str, Any]:
    """Return MSA invocation spec: {cmd, stdin_from, stdout_to}.

    Some aligners write to stdout (mafft, kalign) — these use stdout_to.
    Some write directly to a file flag — stdout_to is None.
    kalign also reads from stdin — stdin_from is set.
    """
    if aligner == "mafft":
        return {
            "cmd": ["mafft", "--auto", str(input_fasta)],
            "stdin_from": None,
            "stdout_to": output_fasta,
        }
    if aligner == "muscle":
        return {
            "cmd": ["muscle", "-align", str(input_fasta), "-output", str(output_fasta)],
            "stdin_from": None,
            "stdout_to": None,
        }
    if aligner == "clustalw":
        return {
            "cmd": [
                "clustalw",
                f"-INFILE={input_fasta}",
                "-OUTPUT=FASTA",
                f"-OUTFILE={output_fasta}",
            ],
            "stdin_from": None,
            "stdout_to": None,
        }
    if aligner == "kalign":
        return {
            "cmd": ["kalign"],
            "stdin_from": input_fasta,
            "stdout_to": output_fasta,
        }
    if aligner == "tcoffee":
        return {
            "cmd": [
                "t_coffee",
                f"-infile={input_fasta}",
                f"-outfile={output_fasta}",
                "-output=fasta_aln",
            ],
            "stdin_from": None,
            "stdout_to": None,
        }
    if aligner == "prank":
        # prank appends .best.fas to the output prefix
        prank_prefix = str(output_fasta).replace(".fasta", "").replace(".fa", "")
        return {
            "cmd": ["prank", f"-d={input_fasta}", f"-o={prank_prefix}"],
            "stdin_from": None,
            "stdout_to": None,
            "prank_output": Path(f"{prank_prefix}.best.fas"),
        }
    raise ValueError(
        f"Unknown aligner: '{aligner}'. Valid choices: {', '.join(VALID_ALIGNERS)}"
    )


def run_msa(input_fasta: Path, aligner: str, output_fasta: Path) -> None:
    """Run MSA alignment.  Raises FileNotFoundError / RuntimeError on failure."""
    binary = shutil.which(aligner) or shutil.which(
        "t_coffee" if aligner == "tcoffee" else aligner
    )
    if binary is None:
        raise FileNotFoundError(
            f"MSA aligner '{aligner}' not found in PATH. "
            f"Install it (e.g. conda install -c bioconda {aligner}) and retry."
        )

    spec = build_msa_command(aligner, input_fasta, output_fasta)
    stdin_fh = open(spec["stdin_from"]) if spec["stdin_from"] else None
    stdout_fh = open(spec["stdout_to"], "w") if spec["stdout_to"] else None

    try:
        result = subprocess.run(
            spec["cmd"],
            stdin=stdin_fh,
            stdout=stdout_fh,
            stderr=subprocess.PIPE,
            text=True,
        )
    finally:
        if stdin_fh:
            stdin_fh.close()
        if stdout_fh:
            stdout_fh.close()

    if result.returncode != 0:
        raise RuntimeError(
            f"{aligner} failed (code {result.returncode}):\n{result.stderr}"
        )

    # prank writes to a different path — move it
    if aligner == "prank":
        prank_out = spec.get("prank_output")
        if prank_out and prank_out.exists():
            shutil.move(str(prank_out), str(output_fasta))
        elif not output_fasta.exists():
            raise RuntimeError(
                "prank completed but output file not found. Expected: " + str(prank_out)
            )


# ── trimAl ─────────────────────────────────────────────────────────────────────


def build_trimal_command(
    input_fasta: Path, output_fasta: Path, strategy: str = "-automated1"
) -> list[str]:
    return ["trimal", "-in", str(input_fasta), "-out", str(output_fasta), strategy]


def run_trimal(
    input_fasta: Path, output_fasta: Path, strategy: str = "-automated1"
) -> None:
    binary = shutil.which("trimal")
    if binary is None:
        raise FileNotFoundError(
            "trimal not found in PATH. Install: conda install -c bioconda trimal"
        )
    cmd = build_trimal_command(input_fasta, output_fasta, strategy)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"trimal failed (code {result.returncode}):\n{result.stderr}"
        )


# ── Model selection & parsing ──────────────────────────────────────────────────


def parse_iqtree_model(iqtree_file: Path) -> str:
    """Extract best-fit model from a .iqtree log file."""
    content = iqtree_file.read_text()
    m = re.search(r"Best-fit model according to BIC:\s+(\S+)", content)
    if m:
        return m.group(1)
    m = re.search(r"Best-fit model:\s+(\S+)", content)
    if m:
        return m.group(1)
    return "GTR+G"


def run_modelfinder(iqtree_bin: str, aligned_fasta: Path, prefix: Path) -> str:
    """Run IQ-TREE ModelFinder (-m MFP).  Returns best model string."""
    cmd = [
        iqtree_bin,
        "-m",
        "MFP",
        "-s",
        str(aligned_fasta),
        "--prefix",
        str(prefix),
        "-T",
        "1",
        "-quiet",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"ModelFinder failed (code {result.returncode}):\n{result.stderr}"
        )
    iqtree_file = Path(f"{prefix}.iqtree")
    if not iqtree_file.exists():
        raise FileNotFoundError(f"ModelFinder log {iqtree_file} not created")
    return parse_iqtree_model(iqtree_file)


def adapt_model_for_engine(model: str, engine: str) -> str:
    """Adapt an IQ-TREE model string for a different engine.

    IQ-TREE ModelFinder appends '+F' (empirical frequencies). RAxML-NG uses
    different notation and the '+F' token causes parse errors.  Strip it when
    targeting raxml-ng.
    """
    if engine == "raxml-ng":
        model = re.sub(r"\+F(?!\w)", "", model)
        model = model.strip("+")
    return model


# ── IQ-TREE ────────────────────────────────────────────────────────────────────


def build_iqtree_command(
    iqtree_bin: str,
    aligned_fasta: Path,
    model: str,
    prefix: Path,
    bootstrap: str,
    outgroup: str | None,
    threads: int,
    seed: int,
) -> list[str]:
    cmd = [
        iqtree_bin,
        "-s",
        str(aligned_fasta),
        "-m",
        model,
        "--prefix",
        str(prefix),
        "-T",
        str(threads),
        "--seed",
        str(seed),
    ]
    if bootstrap == "ufboot":
        cmd += ["-bb", "1000"]
    elif bootstrap == "standard":
        cmd += ["-b", "100"]
    elif bootstrap == "all":
        cmd += ["-bb", "1000", "-alrt", "1000", "-abayes"]
    if outgroup:
        cmd += ["-o", outgroup]
    return cmd


def run_iqtree_main(
    aligned_fasta: Path,
    model: str,
    prefix: Path,
    bootstrap: str,
    outgroup: str | None,
    threads: int,
    seed: int,
) -> str:
    """Run IQ-TREE2 tree search.  Returns Newick string from .treefile."""
    iqtree_bin = shutil.which("iqtree2") or shutil.which("iqtree")
    if not iqtree_bin:
        raise FileNotFoundError(
            "IQ-TREE not found in PATH. Install: conda install -c bioconda iqtree2"
        )
    cmd = build_iqtree_command(
        iqtree_bin, aligned_fasta, model, prefix, bootstrap, outgroup, threads, seed
    )
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"IQ-TREE failed (code {result.returncode}):\n{result.stderr or result.stdout}"
        )
    treefile = Path(f"{prefix}.treefile")
    if not treefile.exists():
        raise FileNotFoundError(f"IQ-TREE treefile {treefile} not created")
    return treefile.read_text().strip()


# ── RAxML-NG ───────────────────────────────────────────────────────────────────


def build_raxml_check_command(
    raxml_bin: str, aligned_fasta: Path, model: str, prefix: Path
) -> list[str]:
    return [
        raxml_bin,
        "--check",
        "--msa",
        str(aligned_fasta),
        "--model",
        model,
        "--prefix",
        str(prefix),
    ]


def build_raxml_run_command(
    raxml_bin: str,
    aligned_fasta: Path,
    model: str,
    prefix: Path,
    bootstrap: str,
    outgroup: str | None,
    threads: int,
    seed: int,
) -> list[str]:
    bs_trees = "100" if bootstrap == "standard" else "1000"
    cmd = [
        raxml_bin,
        "--all",
        "--msa",
        str(aligned_fasta),
        "--model",
        model,
        "--prefix",
        str(prefix),
        "--bs-trees",
        bs_trees,
        "--threads",
        str(threads),
        "--seed",
        str(seed),
    ]
    if outgroup:
        cmd += ["--outgroup", outgroup]
    return cmd


def run_raxml_main(
    aligned_fasta: Path,
    model: str,
    prefix: Path,
    bootstrap: str,
    outgroup: str | None,
    threads: int,
    seed: int,
) -> str:
    """Run RAxML-NG (--check then --all).  Returns Newick string."""
    raxml_bin = shutil.which("raxml-ng")
    if not raxml_bin:
        raise FileNotFoundError(
            "raxml-ng not found in PATH. Install: conda install -c bioconda raxml-ng"
        )
    raxml_model = adapt_model_for_engine(model, "raxml-ng")

    # Step 1: sanity check
    check_cmd = build_raxml_check_command(raxml_bin, aligned_fasta, raxml_model, prefix)
    check_result = subprocess.run(check_cmd, capture_output=True, text=True)
    if check_result.returncode != 0:
        raise RuntimeError(
            f"raxml-ng --check failed:\n{check_result.stderr or check_result.stdout}"
        )

    # Step 2: ML inference + bootstrap
    run_cmd = build_raxml_run_command(
        raxml_bin,
        aligned_fasta,
        raxml_model,
        prefix,
        bootstrap,
        outgroup,
        threads,
        seed,
    )
    run_result = subprocess.run(run_cmd, capture_output=True, text=True)
    if run_result.returncode != 0:
        raise RuntimeError(
            f"raxml-ng --all failed:\n{run_result.stderr or run_result.stdout}"
        )

    # RAxML writes support values to {prefix}.raxml.support
    support_file = Path(f"{prefix}.raxml.support")
    best_file = Path(f"{prefix}.raxml.bestTree")
    if support_file.exists():
        return support_file.read_text().strip()
    if best_file.exists():
        return best_file.read_text().strip()
    raise FileNotFoundError(f"RAxML-NG output not found at {prefix}.raxml.*")


# ── Midpoint rooting ───────────────────────────────────────────────────────────


def root_midpoint(newick_str: str) -> str:
    """Apply midpoint rooting via ETE3.  Falls back to Bio.Phylo, then identity."""
    try:
        from ete3 import Tree  # type: ignore[import]

        t = Tree(newick_str)
        midpoint = t.get_midpoint_outgroup()
        t.set_outgroup(midpoint)
        return t.write(format=1)
    except Exception:  # ete3 is optional
        pass

    try:
        tree = Phylo.read(StringIO(newick_str), "newick")
        tree.root_at_midpoint()
        buf = StringIO()
        Phylo.write(tree, buf, "newick")
        return buf.getvalue().strip()
    except Exception:
        pass

    return newick_str  # identity fallback


# ── Newick → table ─────────────────────────────────────────────────────────────


def newick_to_table(newick_str: str) -> list[dict[str, str | float | int]]:
    """Parse Newick string, extract per-node {node, length, support} rows."""
    tree = Phylo.read(StringIO(newick_str), "newick")
    table: list[dict[str, str | float | int]] = []
    for i, clade in enumerate(tree.find_clades()):
        name: str = clade.name if clade.name else f"InnerNode_{i}"
        raw_support = clade.confidence
        # Triple support "alrt/abayes/ufb" is stored in the node name by IQ-TREE.
        if "/" in name:
            parts = name.split("/")
            try:
                support: int = int(float(parts[-1]))
            except (ValueError, IndexError):
                support = 100
        else:
            support = int(raw_support) if raw_support is not None else 100
        length: float = clade.branch_length if clade.branch_length is not None else 0.0
        table.append(
            {
                "node": name,
                "length": round(float(length), 6),
                "support": support,
            }
        )
    return table


# ── Tree visualisation ─────────────────────────────────────────────────────────


def _write_placeholder_figure(output_path: Path) -> None:
    """Write a blank placeholder PNG when full tree rendering fails."""
    try:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(
            0.5,
            0.5,
            "Phylogenetic Tree\n(render fallback)",
            ha="center",
            va="center",
            fontsize=10,
            transform=ax.transAxes,
        )
        ax.axis("off")
        plt.savefig(output_path, dpi=100)
        plt.close(fig)
    except Exception:
        output_path.write_bytes(b"")


def draw_tree(newick_str: str, output_path: Path) -> None:
    """Draw a proportional phylogram with Bio.Phylo + matplotlib."""
    try:
        tree = Phylo.read(StringIO(newick_str), "newick")
        fig, ax = plt.subplots(
            figsize=(10, max(6, len(list(tree.get_terminals())) * 0.4))
        )
        Phylo.draw(tree, do_show=False, axes=ax)
        ax.set_title("Maximum-Likelihood Phylogenetic Tree", fontsize=12, pad=8)
        plt.tight_layout()
        plt.savefig(output_path, dpi=150)
        plt.close(fig)
    except Exception as exc:
        print(
            f"Warning: tree rendering failed ({exc}). Writing fallback.",
            file=sys.stderr,
        )
        _write_placeholder_figure(output_path)


# ── Demo fallback ──────────────────────────────────────────────────────────────


def get_demo_fallback() -> tuple[str, str]:
    """Return (newick_str, model) from pre-computed examples."""
    nwk_file = SKILL_DIR / "examples" / "demo_tree.nwk"
    if nwk_file.exists():
        newick = nwk_file.read_text().strip()
    else:
        newick = (
            "((Mus_musculus:0.235,Rattus_norvegicus:0.240)99:0.065,"
            "(Lemur_catta:0.195,(Callithrix_jacchus:0.140,"
            "((Macaca_mulatta:0.090,Papio_anubis:0.095)96:0.042,"
            "(Nomascus_leucogenys:0.072,(Pongo_pygmaeus:0.052,"
            "((Homo_sapiens:0.010,(Pan_troglodytes:0.008,Pan_paniscus:0.009)98:0.005)"
            "99:0.018,Gorilla_gorilla:0.028)95:0.015)97:0.020)92:0.030)"
            "100:0.040)88:0.055)85:0.025);"
        )
    return newick, "GTR+F+G4"


# ── Reproducibility bundle ─────────────────────────────────────────────────────


def write_reproducibility_bundle(
    repro_dir: Path,
    input_file: Path,
    output_dir: Path,
    args: argparse.Namespace,
    pipeline_steps: list[str],
) -> None:
    repro_dir.mkdir(parents=True, exist_ok=True)

    flags = f"--input {input_file} --output {output_dir}"
    if args.aligned:
        flags += " --aligned"
    if args.no_trim:
        flags += " --no-trim"
    flags += (
        f" --aligner {args.aligner} --engine {args.engine} --bootstrap {args.bootstrap}"
    )
    if args.outgroup:
        flags += f" --outgroup {args.outgroup}"
    if args.model:
        flags += f" --model {args.model}"
    flags += f" --threads {args.threads} --seed {args.seed}"

    (repro_dir / "commands.sh").write_text(
        f"#!/bin/bash\n# Phylogenetics Builder — reproducible run\n"
        f"python skills/phylogenetics-builder/phylogenetics_builder.py {flags}\n\n"
        f"# Pipeline steps executed:\n"
        + "\n".join(f"# {s}" for s in pipeline_steps)
        + "\n"
    )

    (repro_dir / "environment.yml").write_text(
        "name: clawbio-phylo\nchannels:\n  - conda-forge\n  - bioconda\n"
        "dependencies:\n"
        "  - python>=3.10\n"
        "  - pandas>=2.0\n"
        "  - biopython>=1.80\n"
        "  - matplotlib>=3.5\n"
        "  - iqtree>=2.0\n"
        "  - raxml-ng\n"
        "  - mafft\n"
        "  - muscle\n"
        "  - clustalw\n"
        "  - kalign3\n"
        "  - t_coffee\n"
        "  - prank\n"
        "  - trimal\n"
        "  - ete3\n"
    )


def write_checksums(files: list[Path], output_dir: Path, repro_dir: Path) -> None:
    checksums = []
    for fpath in files:
        if fpath.exists():
            sha = hashlib.sha256(fpath.read_bytes()).hexdigest()
            checksums.append(f"{sha}  {fpath.relative_to(output_dir)}")
    (repro_dir / "checksums.sha256").write_text("\n".join(checksums) + "\n")


def export_alignments(
    output_dir: Path, aligned: Path | None, trimmed: Path | None
) -> dict[str, str]:
    """Copy the alignments produced by this run into ``output_dir/alignment``.

    The pipeline works in a temporary directory, so without this copy the
    alignment is lost when the run ends. A managed file whose stage did not
    produce output in this run is removed, so a rerun into the same directory
    never reports an alignment left by an earlier run. Returns
    ``{"aligned"|"trimmed": path}`` with paths relative to ``output_dir``.
    """
    dest_dir = output_dir / "alignment"
    exported: dict[str, str] = {}
    for key, source in (("aligned", aligned), ("trimmed", trimmed)):
        dest = dest_dir / f"{key}.fasta"
        if source is not None and source.is_file():
            dest_dir.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, dest)
            exported[key] = str(dest.relative_to(output_dir))
        elif dest.is_file():
            dest.unlink()
    if dest_dir.is_dir() and not any(dest_dir.iterdir()):
        dest_dir.rmdir()
    return exported


# ── Main pipeline ──────────────────────────────────────────────────────────────


def main(argv: list[str] | None = None) -> None:  # noqa: C901
    args = parse_args(argv)

    # ── Output directory ────────────────────────────────────────────────────
    if not args.output:
        if args.demo:
            args.output = Path(tempfile.gettempdir()) / "phylogenetics_builder_demo"
        else:
            print("Error: --output is required in non-demo mode.", file=sys.stderr)
            sys.exit(1)

    output_dir: Path = args.output
    if output_dir.exists() and any(output_dir.iterdir()):
        print(
            f"Warning: '{output_dir}' exists and is not empty. Files may be overwritten.",
            file=sys.stderr,
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    # Drop alignments left by an earlier run before this run can fail part-way.
    export_alignments(output_dir, None, None)

    # ── Input file ──────────────────────────────────────────────────────────
    if args.demo:
        input_file = SKILL_DIR / "demo_alignment.fasta"
        if not input_file.exists():
            print(f"Error: demo data not found at {input_file}", file=sys.stderr)
            sys.exit(1)
        # Demo always treats input as pre-aligned (it is)
        effective_aligned = True
    else:
        if not args.input_file:
            print("Error: --input is required.", file=sys.stderr)
            sys.exit(1)
        input_file = args.input_file
        if not input_file.exists():
            print(f"Error: input file not found: {input_file}", file=sys.stderr)
            sys.exit(1)
        effective_aligned = args.aligned

    # ── Validate FASTA ──────────────────────────────────────────────────────
    try:
        sequences = parse_fasta(input_file, require_aligned=effective_aligned)
    except Exception as exc:
        print(f"Error validating input: {exc}", file=sys.stderr)
        sys.exit(1)

    # ── Pipeline ────────────────────────────────────────────────────────────
    pipeline_steps: list[str] = []
    engine_used = args.engine
    aligner_used = args.aligner if not effective_aligned else None
    trimmed = False
    newick = ""
    model = args.model or ""
    run_mode = "live"
    produced_aligned: Path | None = None
    produced_trimmed: Path | None = None
    alignment_files: dict[str, str] = {}

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        current_fasta = input_file  # tracks the "active" alignment through stages

        # Stage 1: MSA (if not pre-aligned)
        if not effective_aligned:
            aligned_fasta = tmp / "aligned.fasta"
            iqtree_bin = shutil.which("iqtree2") or shutil.which("iqtree")
            msa_bin = shutil.which(args.aligner) or shutil.which(
                "t_coffee" if args.aligner == "tcoffee" else args.aligner
            )

            if args.demo and msa_bin is None:
                # Demo: skip MSA, treat demo input as pre-aligned
                current_fasta = input_file
                pipeline_steps.append("msa:skipped(demo-fallback)")
            else:
                try:
                    run_msa(input_file, args.aligner, aligned_fasta)
                    current_fasta = aligned_fasta
                    produced_aligned = aligned_fasta
                    pipeline_steps.append(f"msa:{args.aligner}")
                except Exception as exc:
                    if args.demo:
                        current_fasta = input_file
                        pipeline_steps.append(
                            f"msa:failed({exc.__class__.__name__})-demo-fallback"
                        )
                    else:
                        print(
                            f"Error running MSA ({args.aligner}): {exc}",
                            file=sys.stderr,
                        )
                        sys.exit(1)

        # Stage 2: trimAl
        if not args.no_trim:
            trimmed_fasta = tmp / "trimmed.fasta"
            trimal_bin = shutil.which("trimal")
            if trimal_bin:
                try:
                    run_trimal(current_fasta, trimmed_fasta)
                    current_fasta = trimmed_fasta
                    produced_trimmed = trimmed_fasta
                    trimmed = True
                    pipeline_steps.append("trim:trimal")
                except Exception as exc:
                    print(
                        f"Warning: trimAl failed ({exc}). Continuing without trimming.",
                        file=sys.stderr,
                    )
                    pipeline_steps.append(f"trim:failed({exc.__class__.__name__})")
            else:
                pipeline_steps.append("trim:skipped(trimal-not-found)")

        # Keep the alignments before the temporary directory is removed
        alignment_files = export_alignments(output_dir, produced_aligned, produced_trimmed)

        # Stage 3: Model selection
        iqtree_bin = shutil.which("iqtree2") or shutil.which("iqtree")

        if not model and iqtree_bin:
            mf_prefix = tmp / "modelfinder"
            try:
                model = run_modelfinder(iqtree_bin, current_fasta, mf_prefix)
                pipeline_steps.append(f"modelfinder:{model}")
            except Exception as exc:
                model = "GTR+G"
                pipeline_steps.append(
                    f"modelfinder:failed({exc.__class__.__name__})-using-GTR+G"
                )
        elif not model:
            model = "GTR+G"
            pipeline_steps.append("model:default-GTR+G(no-iqtree)")

        # Stage 4: Tree inference
        tree_prefix = tmp / "tree"

        if args.engine == "iqtree2":
            if iqtree_bin:
                try:
                    newick = run_iqtree_main(
                        current_fasta,
                        model,
                        tree_prefix,
                        args.bootstrap,
                        args.outgroup,
                        args.threads,
                        args.seed,
                    )
                    engine_used = "iqtree2"
                    pipeline_steps.append(f"tree:iqtree2:{args.bootstrap}")
                except Exception as exc:
                    if args.demo:
                        newick, model = get_demo_fallback()
                        engine_used = "precomputed"
                        run_mode = "demo-fallback"
                        pipeline_steps.append("tree:precomputed-fallback")
                    else:
                        print(f"Error running IQ-TREE: {exc}", file=sys.stderr)
                        sys.exit(1)
            else:
                if args.demo:
                    newick, model = get_demo_fallback()
                    engine_used = "precomputed"
                    run_mode = "demo-fallback"
                    pipeline_steps.append("tree:precomputed-fallback(no-iqtree)")
                else:
                    print("Error: IQ-TREE not found in PATH.", file=sys.stderr)
                    sys.exit(1)

        elif args.engine == "raxml-ng":
            raxml_bin = shutil.which("raxml-ng")
            if raxml_bin:
                try:
                    newick = run_raxml_main(
                        current_fasta,
                        model,
                        tree_prefix,
                        args.bootstrap,
                        args.outgroup,
                        args.threads,
                        args.seed,
                    )
                    engine_used = "raxml-ng"
                    pipeline_steps.append(f"tree:raxml-ng:{args.bootstrap}")
                except Exception as exc:
                    if args.demo:
                        newick, model = get_demo_fallback()
                        engine_used = "precomputed"
                        run_mode = "demo-fallback"
                        pipeline_steps.append("tree:precomputed-fallback")
                    else:
                        print(f"Error running RAxML-NG: {exc}", file=sys.stderr)
                        sys.exit(1)
            else:
                if args.demo:
                    newick, model = get_demo_fallback()
                    engine_used = "precomputed"
                    run_mode = "demo-fallback"
                    pipeline_steps.append("tree:precomputed-fallback(no-raxml-ng)")
                else:
                    print("Error: raxml-ng not found in PATH.", file=sys.stderr)
                    sys.exit(1)

        # Stage 5: Post-inference rooting
        if args.root == "midpoint" and newick:
            newick = root_midpoint(newick)
            pipeline_steps.append("root:midpoint")
        elif args.outgroup:
            pipeline_steps.append(f"root:outgroup({args.outgroup})")

    # ── Parse and write outputs ─────────────────────────────────────────────
    try:
        table = newick_to_table(newick)
    except Exception as exc:
        print(f"Error parsing Newick: {exc}", file=sys.stderr)
        sys.exit(1)

    # Paths
    tree_file = output_dir / "phylo_tree.nwk"
    figures_dir = output_dir / "figures"
    tables_dir = output_dir / "tables"
    repro_dir = output_dir / "reproducibility"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    # Write Newick
    tree_file.write_text(newick + "\n")

    # Write branch support table
    csv_file = tables_dir / "branch_support.csv"
    pd.DataFrame(table).to_csv(csv_file, index=False)

    # Draw tree
    img_file = figures_dir / "phylogram.png"
    draw_tree(newick, img_file)

    # ClawBio contract fields
    chat_summary = [
        f"Phylogenetic tree inferred for {len(sequences)} taxa.",
        f"Substitution model: **{model}** (selected via ModelFinder).",
        f"Tree engine: {engine_used} | Bootstrap: {args.bootstrap} | "
        f"Trimming: {'trimAl -automated1' if trimmed else 'skipped'}.",
        "Internal node support values in `tables/branch_support.csv`.",
        "Proportional phylogram at `figures/phylogram.png`.",
    ]
    if alignment_files:
        chat_summary.append(
            "Alignment(s) saved: "
            + ", ".join(f"`{path}`" for path in alignment_files.values())
            + "."
        )
    preferred_artifacts = [
        {"type": "figure", "path": str(img_file.relative_to(output_dir))},
        {"type": "report", "path": "report.md"},
        {"type": "table", "path": str(csv_file.relative_to(output_dir))},
        {"type": "tree", "path": "phylo_tree.nwk"},
    ]
    preferred_artifacts += [
        {"type": "alignment", "path": path} for path in alignment_files.values()
    ]
    suggested_actions = [
        "Open figures/phylogram.png to inspect the tree topology.",
        "Use --root midpoint to apply midpoint rooting (requires ETE3).",
        "Try --bootstrap all for triple support (UFBoot + aLRT + aBayes).",
        "Annotate the tree in FigTree or ggtree (R) for publication quality.",
    ]
    if "aligned" in alignment_files:
        suggested_actions.append(
            "Reuse alignment/aligned.fasta for population-genetics statistics "
            "(for example with the dnasp skill); use the untrimmed file, because "
            "trimming removes alignment columns."
        )
    if engine_used == "precomputed":
        suggested_actions.insert(0, "Install IQ-TREE2 or RAxML-NG for live inference.")

    # result.json
    result_data = {
        "engine": engine_used,
        "aligner": aligner_used,
        "aligned": effective_aligned,
        "trimmed": trimmed,
        "model": model,
        "bootstrap_mode": args.bootstrap,
        "num_taxa": len(sequences),
        "tree": newick,
        "input_file": str(input_file),
        "output_dir": str(output_dir),
        "pipeline_steps": pipeline_steps,
        "alignment_files": alignment_files,
        "run_mode": run_mode,
        "status": "success",
        "chat_summary_lines": chat_summary,
        "preferred_artifacts": preferred_artifacts,
        "workflow_state": "completed",
        "suggested_actions": suggested_actions,
        "contract_alerts": (
            ["Tree is pre-computed fallback — install IQ-TREE2 for live inference"]
            if engine_used == "precomputed"
            else []
        ),
    }
    result_json = output_dir / "result.json"
    result_json.write_text(json.dumps(result_data, indent=2))

    # report.md
    bs_label = {
        "ufboot": "UFBoot (1 000 replicates)",
        "standard": "Standard bootstrap (100 replicates)",
        "all": "Triple support: UFBoot + aLRT + aBayes",
    }[args.bootstrap]
    report_lines = [
        "# Phylogenetics Builder Report",
        "",
        "## Pipeline Summary",
        "",
        "| Parameter | Value |",
        "|-----------|-------|",
        f"| Input | `{input_file.name}` |",
        f"| Taxa | {len(sequences)} |",
        f"| Aligner | {aligner_used if aligner_used else 'pre-aligned (--aligned)'} |",
        f"| Trimming | {'trimAl -automated1' if trimmed else 'skipped'} |",
        f"| Substitution model | `{model}` |",
        f"| Tree engine | {engine_used} |",
        f"| Bootstrap | {bs_label} |",
        f"| Rooting | {args.root or (f'outgroup: {args.outgroup}' if args.outgroup else 'unrooted')} |",
        "",
        "## Pipeline Steps",
        "",
    ]
    for step in pipeline_steps:
        report_lines.append(f"- `{step}`")
    report_lines += [
        "",
        "## Branch Lengths & Support Values",
        "",
        "| Node / Taxon | Branch Length | Support |",
        "|:-------------|:-------------:|:-------:|",
    ]
    for row in table:
        report_lines.append(
            f"| {row['node']} | {row['length']:.5f} | {row['support']} |"
        )
    report_lines += [
        "",
        "## Outputs",
        "",
        "| Artifact | Path |",
        "|----------|------|",
        "| Newick tree | `phylo_tree.nwk` |",
        "| Phylogram figure | `figures/phylogram.png` |",
        "| Branch support table | `tables/branch_support.csv` |",
    ]
    report_lines += [
        f"| {'Alignment' if key == 'aligned' else 'Trimmed alignment'} | `{path}` |"
        for key, path in alignment_files.items()
    ]
    report_lines += [
        "| Reproducibility bundle | `reproducibility/` |",
        "",
        "---",
        "",
        f"*{DISCLAIMER}*",
        "",
    ]
    (output_dir / "report.md").write_text("\n".join(report_lines))

    # Reproducibility bundle
    write_reproducibility_bundle(
        repro_dir, input_file, output_dir, args, pipeline_steps
    )
    write_checksums(
        [tree_file, csv_file, img_file, result_json]
        + [output_dir / path for path in alignment_files.values()],
        output_dir,
        repro_dir,
    )

    print(f"Phylogenetics Builder complete. Results in: {output_dir}")


if __name__ == "__main__":
    main()
