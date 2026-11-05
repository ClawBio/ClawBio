"""Shared constants and pipeline logic for the phylogenetics_builder stage scripts.

Ported faithfully from ClawBio 1.x's phylogenetics_builder.py
(https://github.com/ClawBio/ClawBio/tree/main/skills/phylogenetics-builder).
The science (MSA invocation, trimAl trimming, IQ-TREE2 ModelFinder, IQ-TREE2/
RAxML-NG tree inference, ETE3/Bio.Phylo midpoint rooting, Newick parsing,
phylogram rendering) is unchanged. What changed is the CLI surface: the 1.x
script was one monolithic `main()` that threaded a mutable `pipeline_steps`
list and a temp directory through six inline stages; here the same logic is
split across seven stage scripts (align.py, trim.py, select_model.py,
build_tree.py, root.py, render.py, report.py) so each maps onto one clawbio
pipeline.yaml step, with intermediate results (alignment FASTA, model
string, Newick tree) handed off as files under `{out}` instead of local
variables, and each stage's provenance recorded in a small `*_meta.json`
sidecar that report.py fans in on at the end (mirroring claw_metagenomics's
stage-script split). The internal reproducibility-bundle writer from 1.x
(`write_reproducibility_bundle`/`write_checksums`) is dropped — ClawBio
2.0's own `clawbio bundle` verb now owns that job.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
from io import StringIO
from pathlib import Path
from typing import Any, Optional

from Bio import Phylo

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)

VALID_ALIGNERS = ("mafft", "muscle", "clustalw", "kalign", "tcoffee", "prank")
VALID_ENGINES = ("iqtree2", "raxml-ng")
VALID_BOOTSTRAP = ("ufboot", "standard", "all")

# Sentinel values threaded through pipeline.yaml as plain-string inputs
# (pipeline.yaml has no conditional/optional-value mechanism, so every
# referenced input is `required: true` and these strings stand in for
# "nothing" — see SKILL.md).
SENTINEL_AUTO = "auto"      # --model auto -> run ModelFinder
SENTINEL_NONE = "none"      # --outgroup none / --root none -> not set
TRUE_STRINGS = {"1", "true", "yes", "on"}


def as_bool(value: str) -> bool:
    """Parse a pipeline.yaml string input ('true'/'false') as bool."""
    return str(value).strip().lower() in TRUE_STRINGS


def none_if_sentinel(value: Optional[str]) -> Optional[str]:
    """Return None if value is empty or the literal sentinel 'none'."""
    if value is None:
        return None
    v = value.strip()
    if not v or v.lower() == SENTINEL_NONE:
        return None
    return v


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
                "Use an aligner to run MSA first, or check your input file."
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


def msa_binary(aligner: str) -> Optional[str]:
    return shutil.which(aligner) or (shutil.which("t_coffee") if aligner == "tcoffee" else None)


def run_msa(input_fasta: Path, aligner: str, output_fasta: Path) -> None:
    """Run MSA alignment.  Raises FileNotFoundError / RuntimeError on failure."""
    binary = msa_binary(aligner)
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


def iqtree_binary() -> Optional[str]:
    return shutil.which("iqtree2") or shutil.which("iqtree")


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
    outgroup: Optional[str],
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
    outgroup: Optional[str],
    threads: int,
    seed: int,
) -> str:
    """Run IQ-TREE2 tree search.  Returns Newick string from .treefile."""
    iqtree_bin = iqtree_binary()
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
    outgroup: Optional[str],
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
    outgroup: Optional[str],
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

    check_cmd = build_raxml_check_command(raxml_bin, aligned_fasta, raxml_model, prefix)
    check_result = subprocess.run(check_cmd, capture_output=True, text=True)
    if check_result.returncode != 0:
        raise RuntimeError(
            f"raxml-ng --check failed:\n{check_result.stderr or check_result.stdout}"
        )

    run_cmd = build_raxml_run_command(
        raxml_bin, aligned_fasta, raxml_model, prefix, bootstrap, outgroup, threads, seed,
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
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(
            0.5, 0.5, "Phylogenetic Tree\n(render fallback)",
            ha="center", va="center", fontsize=10, transform=ax.transAxes,
        )
        ax.axis("off")
        plt.savefig(output_path, dpi=100)
        plt.close(fig)
    except Exception:
        output_path.write_bytes(b"")


def draw_tree(newick_str: str, output_path: Path) -> None:
    """Draw a proportional phylogram with Bio.Phylo + matplotlib."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

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


def get_demo_fallback(skill_dir: Path) -> tuple[str, str]:
    """Return (newick_str, model) from the bundled pre-computed demo tree."""
    nwk_file = Path(skill_dir) / "examples" / "demo_tree.nwk"
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
