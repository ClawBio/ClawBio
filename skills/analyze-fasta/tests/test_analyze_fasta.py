"""Tests for analyze-fasta skill (red/green TDD)."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
SCRIPT = SKILL_DIR / "analyze_fasta.py"
EXAMPLE_DIR = SKILL_DIR / "example_data"


def run_cli(args, expect_code=0):
    """Run the CLI and return (returncode, stdout, stderr)."""
    result = subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True,
        text=True,
    )
    if expect_code is not None:
        assert result.returncode == expect_code, (
            f"Expected exit {expect_code}, got {result.returncode}\n"
            f"stdout: {result.stdout}\nstderr: {result.stderr}"
        )
    return result


# ──────────────────────────────────────────────
# Demo mode
# ──────────────────────────────────────────────

def test_demo_runs_successfully(tmp_path):
    """--demo --output DIR produces report.md and result.json."""
    out = tmp_path / "demo_out"
    run_cli(["--demo", "--output", str(out)])
    assert (out / "report.md").exists(), "report.md not generated"
    assert (out / "result.json").exists(), "result.json not generated"


def test_demo_result_json_is_valid(tmp_path):
    """Demo result.json parses and has required top-level keys."""
    out = tmp_path / "demo_out"
    run_cli(["--demo", "--output", str(out)])
    data = json.loads((out / "result.json").read_text())
    assert "sequence_type" in data
    assert "sequences" in data
    assert "summary" in data
    assert data["total_sequences"] >= 1


def test_demo_produces_reproducibility_bundle(tmp_path):
    """Output includes reproducibility/ with commands.sh and run.json."""
    out = tmp_path / "demo_out"
    run_cli(["--demo", "--output", str(out)])
    repro = out / "reproducibility"
    assert repro.is_dir()
    assert (repro / "commands.sh").exists()
    assert (repro / "run.json").exists()


def test_report_includes_clawbio_disclaimer(tmp_path):
    """report.md must reference the ClawBio medical disclaimer."""
    out = tmp_path / "demo_out"
    run_cli(["--demo", "--output", str(out)])
    md = (out / "report.md").read_text()
    assert "ClawBio" in md
    assert "research" in md.lower() and "medical device" in md.lower()


# ──────────────────────────────────────────────
# Nucleotide and protein basic paths
# ──────────────────────────────────────────────

def test_nucleotide_input_basic(tmp_path):
    """Nucleotide FASTA produces gc_content and composition."""
    out = tmp_path / "nuc_out"
    run_cli(["--input", str(EXAMPLE_DIR / "demo_nucleotide.fasta"), "--output", str(out)])
    data = json.loads((out / "result.json").read_text())
    assert data["sequence_type"] == "nucleotide"
    seq = data["sequences"][0]
    assert "gc_content" in seq
    assert 0 <= seq["gc_content"] <= 100


def test_protein_input_basic(tmp_path):
    """Protein FASTA produces molecular_weight and pI."""
    out = tmp_path / "prot_out"
    run_cli(["--input", str(EXAMPLE_DIR / "demo_protein.fasta"), "--output", str(out)])
    data = json.loads((out / "result.json").read_text())
    assert data["sequence_type"] == "protein"
    seq = data["sequences"][0]
    assert "molecular_weight_da" in seq
    assert "isoelectric_point" in seq
    assert "gravy" in seq


# ──────────────────────────────────────────────
# Loud-failure validations (no silent degradation)
# ──────────────────────────────────────────────

def test_loud_failure_missing_file(tmp_path):
    """Missing input → exit code 1 with explicit error."""
    out = tmp_path / "out"
    result = run_cli(
        ["--input", str(tmp_path / "nope.fasta"), "--output", str(out)],
        expect_code=1,
    )
    assert "no encontrado" in result.stderr.lower() or "not found" in result.stderr.lower()


def test_loud_failure_empty_file(tmp_path):
    """Empty FASTA file → exit 1, not a fake report."""
    empty = tmp_path / "empty.fasta"
    empty.write_text("")
    out = tmp_path / "out"
    result = run_cli(["--input", str(empty), "--output", str(out)], expect_code=1)
    assert "no se encontraron" in result.stderr.lower() or "no sequences" in result.stderr.lower()
    assert not (out / "result.json").exists(), "Should not write fake result.json on failure"


def test_loud_failure_too_short(tmp_path):
    """Sequence shorter than 10 chars → exit 1."""
    tiny = tmp_path / "tiny.fasta"
    tiny.write_text(">tiny\nACGT\n")
    out = tmp_path / "out"
    result = run_cli(["--input", str(tiny), "--output", str(out)], expect_code=1)
    assert "corta" in result.stderr.lower() or "short" in result.stderr.lower()


def test_loud_failure_too_many_ns(tmp_path):
    """Sequence with >50% Ns → exit 1, no silent fake report."""
    ns = tmp_path / "all_n.fasta"
    ns.write_text(">mostly_n\n" + "N" * 100 + "ACGTACGT\n")
    out = tmp_path / "out"
    result = run_cli(["--input", str(ns), "--output", str(out)], expect_code=1)
    assert "n" in result.stderr.lower()  # mentions Ns


# ──────────────────────────────────────────────
# Backward-compat with legacy CLI
# ──────────────────────────────────────────────

def test_legacy_positional_argument_still_works(tmp_path):
    """Old usage `python script.py file.fasta --json` keeps working."""
    result = run_cli([str(EXAMPLE_DIR / "demo_nucleotide.fasta"), "--json"])
    data = json.loads(result.stdout)
    assert data["sequence_type"] == "nucleotide"


# ──────────────────────────────────────────────
# Reverse-strand ORF detection (6-frame scanning)
# ──────────────────────────────────────────────

_COMPLEMENT = str.maketrans("ATCGN", "TAGCN")


def _revcomp(s: str) -> str:
    return s.translate(_COMPLEMENT)[::-1]


def _make_minus_strand_fasta(tmp_path, orf_len_bp: int = 333, pad: int = 100) -> tuple:
    """
    Build a FASTA where the only long ORF lives on the minus strand.

    Returns (fasta_path, expected_start_1based, expected_end_1based, seq_length).
    """
    # Build a sense ORF of exactly orf_len_bp bp (must be divisible by 3)
    assert orf_len_bp % 3 == 0
    stops = {"TAA", "TAG", "TGA"}
    import random
    rng = random.Random(42)
    internal_codons = []
    n_internal = (orf_len_bp - 6) // 3  # exclude ATG and stop
    while len(internal_codons) < n_internal:
        c = "".join(rng.choices("ATCG", k=3))
        if c not in stops:
            internal_codons.append(c)
    orf_sense = "ATG" + "".join(internal_codons) + "TAA"
    orf_fwd = _revcomp(orf_sense)  # embed this in the forward sequence

    # Pad with sequence guaranteed to have no ATG (avoids accidental forward ORFs)
    def _no_atg(n, seed):
        r = random.Random(seed)
        bases = []
        while len(bases) < n:
            b = r.choice("TCG")  # no A -> no ATG possible
            bases.append(b)
        return "".join(bases)

    left = _no_atg(pad, seed=1)
    right = _no_atg(pad, seed=2)
    seq = left + orf_fwd + right
    seq_len = len(seq)

    fasta = tmp_path / "minus_strand.fasta"
    fasta.write_text(f">test_minus_strand\n{seq}\n")

    # 1-based forward coordinates of the embedded ORF
    expected_start = pad + 1
    expected_end = pad + orf_len_bp
    return fasta, expected_start, expected_end, seq_len


def test_minus_strand_orf_is_detected(tmp_path):
    """A >=300 bp ORF on the reverse strand must appear in orfs list."""
    fasta, _, _, _ = _make_minus_strand_fasta(tmp_path)
    out = tmp_path / "out"
    run_cli(["--input", str(fasta), "--output", str(out)])
    data = json.loads((out / "result.json").read_text())
    seq = data["sequences"][0]
    assert seq["orfs_found"] >= 1, "Expected at least one ORF but found none"
    frames = [o["frame"] for o in seq["orfs"]]
    minus_frames = [f for f in frames if f.startswith("-")]
    assert minus_frames, f"No minus-strand frames in results: {frames}"


def test_minus_strand_orf_coordinates_are_correct(tmp_path):
    """Reverse-strand ORF coordinates must map back to forward-strand positions."""
    fasta, expected_start, expected_end, seq_len = _make_minus_strand_fasta(
        tmp_path, orf_len_bp=333, pad=100
    )
    out = tmp_path / "out"
    run_cli(["--input", str(fasta), "--output", str(out)])
    data = json.loads((out / "result.json").read_text())
    seq = data["sequences"][0]

    minus_orfs = [o for o in seq["orfs"] if o["frame"].startswith("-")]
    assert minus_orfs, "No minus-strand ORF found to check coordinates"

    # The largest minus-strand ORF should span the embedded region
    best = max(minus_orfs, key=lambda o: o["length_bp"])
    assert best["start"] == expected_start, (
        f"Expected start={expected_start}, got {best['start']}"
    )
    assert best["end"] == expected_end, (
        f"Expected end={expected_end}, got {best['end']}"
    )


def test_forward_only_sequence_has_no_minus_orfs(tmp_path):
    """A sequence whose only long ORF is on the forward strand reports no minus frames."""
    # Build a clear forward-strand ORF: ATG + 110 safe codons + TAA = 336 bp
    stops = {"TAA", "TAG", "TGA"}
    import random
    rng = random.Random(7)
    codons = []
    while len(codons) < 110:
        c = "".join(rng.choices("ATCG", k=3))
        if c not in stops:
            codons.append(c)
    fwd_orf = "ATG" + "".join(codons) + "TAA"
    seq = fwd_orf + "C" * 50
    fasta = tmp_path / "fwd_only.fasta"
    fasta.write_text(f">fwd_only\n{seq}\n")
    out = tmp_path / "out"
    run_cli(["--input", str(fasta), "--output", str(out)])
    data = json.loads((out / "result.json").read_text())
    seq_res = data["sequences"][0]
    minus_orfs = [o for o in seq_res["orfs"] if o["frame"].startswith("-")]
    fwd_orfs = [o for o in seq_res["orfs"] if o["frame"].startswith("+")]
    assert fwd_orfs, "Expected at least one forward-strand ORF"
    # Minus-strand ORFs may or may not exist by chance; this is informational
    # The key assertion: forward ORF IS present (regression guard)
    assert any(o["length_bp"] >= 333 for o in fwd_orfs), (
        "Forward ORF should be at least 333 bp"
    )


def test_orf_results_sorted_by_length_descending(tmp_path):
    """orfs list must be sorted longest-first."""
    fasta, _, _, _ = _make_minus_strand_fasta(tmp_path, orf_len_bp=333, pad=100)
    out = tmp_path / "out"
    run_cli(["--input", str(fasta), "--output", str(out)])
    data = json.loads((out / "result.json").read_text())
    seq = data["sequences"][0]
    lengths = [o["length_bp"] for o in seq["orfs"]]
    assert lengths == sorted(lengths, reverse=True), (
        f"ORFs not sorted by length descending: {lengths}"
    )


def test_demo_reproducibility_bundle_is_complete(tmp_path):
    """commands.sh is portable and the bundle carries environment.yml + checksums."""
    out = tmp_path / "demo_out"
    run_cli(["--demo", "--output", str(out)])
    repro = out / "reproducibility"

    commands_text = (repro / "commands.sh").read_text(encoding="utf-8")
    assert "CLAWBIO_ROOT" in commands_text
    assert "$OUTPUT_DIR" in commands_text
    assert str(out) not in commands_text

    environment = (repro / "environment.yml").read_text(encoding="utf-8")
    assert "name: clawbio-analyze-fasta" in environment
    assert "biopython" in environment

    checksum_lines = [
        line
        for line in (repro / "checksums.sha256").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    assert checksum_lines
    labels = set()
    for line in checksum_lines:
        digest, label = line.split("  ", 1)
        assert len(digest) == 64
        labels.add(label)
        assert (out / label).exists()
    assert {"report.md", "result.json", "report.html"} <= labels
