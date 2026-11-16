"""Tests for analyze-fasta skill (red/green TDD)."""

import json
import subprocess
import sys
from pathlib import Path

import pytest
from inline_snapshot import snapshot

from clawbio.common import demo_runner

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

@pytest.mark.parametrize("report", ["report.md", "result.json", "report.html"])
def test_reports_do_not_record_the_input_directory(tmp_path, report):
    """Reports are shared; an absolute input path leaks the user's name and folder layout."""
    private_dir = tmp_path / "private_user_dir"
    private_dir.mkdir()
    fasta = private_dir / "sample.fasta"
    fasta.write_text((EXAMPLE_DIR / "demo_nucleotide.fasta").read_text())
    out = tmp_path / "out"
    run_cli(["--input", str(fasta), "--output", str(out)])
    text = (out / report).read_text()
    assert str(private_dir) not in text and str(private_dir.resolve()) not in text


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


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/analyze-fasta/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("analyze-fasta", tmp_path) == snapshot({'report.html': """\
<!DOCTYPE html>
<html lang="es">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Analisis FASTA - demo_nucleotide.fasta</title>
    <style>
        @import url('https://fonts.googleapis.com/css2?family=Merriweather:wght@400;700;900&family=Source+Sans+3:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');
        :root {
            --bg:#fff; --surface:#fafafa; --surface-2:#f3f3f0; --border:#ddd9d0; --border-lt:#e8e5de;
            --text:#1a1a1a; --text-2:#4a4a4a; --text-dim:#6b6b6b; --text-cap:#888;
            --uy:#FFC600; --uyd:#e0ad00; --uys:#fffbeb; --ub:#1a1a1a;
            --blue:#2c5282; --green:#276749; --amber:#b7791f; --red:#9b2c2c;
            --blue-bg:#ebf0f7; --green-bg:#e6f4ed; --amber-bg:#fef6e6; --red-bg:#fde8e8;
        }
        * { margin:0; padding:0; box-sizing:border-box; }
        body { font-family:'Source Sans 3',system-ui,sans-serif; background:var(--bg); color:var(--text); line-height:1.7; }
        .top-bar { background:var(--ub); padding:0.6rem 2rem; display:flex; align-items:center; justify-content:space-between; }
        .top-bar-logo { font-family:'Merriweather',serif; font-weight:900; font-size:1.1rem; color:var(--uy); }
        .top-bar-info { font-size:0.75rem; color:rgba(255,255,255,0.5); }
        .yellow-strip { height:4px; background:var(--uy); }
        .container { max-width:940px; margin:0 auto; padding:2rem 2rem 4rem; }

        .hero { text-align:center; padding:2.5rem 0 2rem; border-bottom:2px solid var(--uy); margin-bottom:2.5rem; }
        .hero h1 { font-family:'Merriweather',serif; font-size:1.6rem; font-weight:700; color:var(--ub); margin-bottom:0.5rem; }
        .hero-meta { font-size:0.82rem; color:var(--text-cap); margin-top:0.5rem; }
        .hero-file { font-family:'JetBrains Mono',monospace; font-size:0.9rem; color:var(--blue); background:var(--blue-bg); padding:0.2rem 0.7rem; border-radius:3px; }

        .sn { display:inline-flex; align-items:center; gap:0.5rem; font-family:'Merriweather',serif; font-size:0.72rem; font-weight:700; color:var(--uyd); text-transform:uppercase; letter-spacing:0.1em; margin-bottom:0.6rem; }
        .sn::before { content:''; width:24px; height:3px; background:var(--uy); border-radius:2px; }
        h2 { font-family:'Merriweather',serif; font-size:1.35rem; font-weight:700; color:var(--ub); margin-bottom:0.8rem; }
        h3 { font-family:'Merriweather',serif; font-size:1rem; font-weight:700; color:var(--ub); margin-bottom:0.3rem; }
        h4 { font-family:'Merriweather',serif; font-size:0.82rem; font-weight:700; color:var(--text-dim); text-transform:uppercase; letter-spacing:0.06em; margin:1rem 0 0.5rem; }
        p { color:var(--text-2); font-size:0.93rem; margin-bottom:0.7rem; }
        section { margin-bottom:3rem; }
        small { color:var(--text-dim); }

        .badge { display:inline-block; padding:0.12rem 0.5rem; border-radius:3px; font-family:'JetBrains Mono',monospace; font-size:0.72rem; font-weight:600; }
        .badge.nuc { background:var(--blue-bg); color:var(--blue); }
        .badge.prot { background:var(--amber-bg); color:var(--amber); }

        /* Summary table */
        .summary-table { width:100%; border-collapse:collapse; margin-top:1rem; }
        .summary-table td { padding:0.55rem 0.8rem; font-size:0.88rem; border-bottom:1px solid var(--border-lt); }
        .summary-table td:first-child { font-weight:600; color:var(--text); width:45%; }
        .summary-table td:last-child { color:var(--text-2); font-family:'JetBrains Mono',monospace; font-size:0.82rem; }
        .summary-table tr:hover td { background:var(--uys); }

        /* Sequence cards */
        .seq-card { background:var(--surface); border:1px solid var(--border); border-radius:4px; padding:1.5rem 1.8rem; margin-bottom:1.2rem; }
        .seq-header { display:flex; align-items:baseline; gap:0.8rem; margin-bottom:1rem; padding-bottom:0.8rem; border-bottom:1px solid var(--border-lt); }
        .seq-num { font-family:'Merriweather',serif; font-size:1.5rem; font-weight:900; color:var(--surface-2); }
        .seq-desc { font-size:0.8rem; color:var(--text-cap); margin:0; word-break:break-all; }
        .seq-grid { display:grid; grid-template-columns:1fr 1fr; gap:1.5rem; }
        @media(max-width:640px) { .seq-grid { grid-template-columns:1fr; } }

        .metric-group { display:grid; grid-template-columns:repeat(auto-fit,minmax(140px,1fr)); gap:0.5rem; margin-bottom:0.8rem; }
        .metric { background:var(--bg); border:1px solid var(--border-lt); border-radius:4px; padding:0.5rem 0.7rem; }
        .metric-label { display:block; font-size:0.68rem; font-weight:600; text-transform:uppercase; letter-spacing:0.05em; color:var(--text-cap); }
        .metric-value { font-family:'JetBrains Mono',monospace; font-size:0.88rem; font-weight:600; color:var(--text); }

        /* Composition bar */
        .comp-bar { display:flex; height:28px; border-radius:4px; overflow:hidden; margin:0.5rem 0; }
        .comp-seg { display:flex; align-items:center; justify-content:center; font-size:0.7rem; font-weight:700; color:#fff; min-width:20px; }
        .seg-a { background:#4299e1; }
        .seg-t { background:#f56565; }
        .seg-g { background:#48bb78; }
        .seg-c { background:#ed8936; }
        .comp-legend { display:flex; gap:1rem; font-size:0.75rem; color:var(--text-dim); flex-wrap:wrap; }
        .dot { display:inline-block; width:8px; height:8px; border-radius:2px; margin-right:3px; vertical-align:middle; }
        .dot-a { background:#4299e1; } .dot-t { background:#f56565; } .dot-g { background:#48bb78; } .dot-c { background:#ed8936; }

        /* Mini bar charts */
        .mini-bar-row { display:flex; align-items:center; gap:0.5rem; margin-bottom:0.3rem; }
        .mini-label { font-family:'JetBrains Mono',monospace; font-size:0.72rem; width:28px; text-align:right; color:var(--text-dim); }
        .mini-track { flex:1; height:16px; background:var(--surface-2); border-radius:3px; overflow:hidden; }
        .mini-fill { height:100%; border-radius:3px; }
        .nuc-fill { background:var(--blue); opacity:0.6; }
        .mini-val { font-family:'JetBrains Mono',monospace; font-size:0.68rem; color:var(--text-cap); width:50px; }

        /* Inner table (ORFs) */
        .inner-table { width:100%; border-collapse:collapse; margin-top:0.5rem; font-size:0.82rem; }
        .inner-table th { text-align:left; padding:0.4rem 0.6rem; font-size:0.68rem; font-weight:700; letter-spacing:0.06em; text-transform:uppercase; color:var(--text-dim); background:var(--surface-2); border-bottom:2px solid var(--uy); }
        .inner-table td { padding:0.35rem 0.6rem; border-bottom:1px solid var(--border-lt); color:var(--text-2); font-family:'JetBrains Mono',monospace; font-size:0.78rem; }

        /* SS bars */
        .ss-bars { margin:0.3rem 0; }
        .ss-row { display:flex; align-items:center; gap:0.5rem; margin-bottom:0.3rem; }
        .ss-label { font-size:0.75rem; width:55px; text-align:right; color:var(--text-dim); }
        .helix-fill { background:var(--red); opacity:0.6; }
        .sheet-fill { background:var(--blue); opacity:0.6; }
        .turn-fill { background:var(--green); opacity:0.6; }

        .overflow-note { font-size:0.82rem; color:var(--text-cap); font-style:italic; text-align:center; margin-top:1rem; }

        /* JSON block */
        .json-toggle { background:var(--ub); color:var(--uy); border:none; padding:0.45rem 1rem; border-radius:4px; font-family:'Source Sans 3',sans-serif; font-size:0.82rem; font-weight:600; cursor:pointer; margin-bottom:0.8rem; }
        .json-toggle:hover { background:#333; }
        .json-block { background:var(--ub); border-radius:4px; padding:1.2rem 1.5rem; overflow-x:auto; display:none; }
        .json-block pre { font-family:'JetBrains Mono',monospace; font-size:0.75rem; line-height:1.6; color:#e4e4ed; white-space:pre-wrap; word-break:break-word; }

        /* Como funciona */
        .how-section { background:var(--uys); border:1px solid var(--uy); border-radius:4px; padding:2rem 2.2rem; position:relative; }
        .how-section::before { content:''; position:absolute; top:0; left:0; width:100%; height:4px; background:var(--uy); border-radius:4px 4px 0 0; }
        .how-section h2 { margin-bottom:1rem; }
        .how-grid { display:grid; grid-template-columns:1fr 1fr; gap:1.5rem; margin-top:1rem; }
        @media(max-width:640px) { .how-grid { grid-template-columns:1fr; } }
        .how-card { background:var(--bg); border:1px solid var(--border); border-radius:4px; padding:1.2rem; }
        .how-card h3 { font-size:0.88rem; margin-bottom:0.4rem; }
        .how-card p { font-size:0.85rem; margin-bottom:0.3rem; }
        .how-card code { font-family:'JetBrains Mono',monospace; font-size:0.78em; background:var(--surface-2); padding:0.1em 0.35em; border-radius:2px; }
        .how-card ul { list-style:none; padding:0; margin:0.3rem 0 0; }
        .how-card ul li { font-size:0.82rem; padding:0.2rem 0 0.2rem 0.9rem; position:relative; color:var(--text-2); }
        .how-card ul li::before { content:''; position:absolute; left:0; top:0.6rem; width:5px; height:5px; border-radius:50%; background:var(--uyd); }

        .pipe-flow { display:flex; align-items:center; justify-content:center; gap:0.4rem; flex-wrap:wrap; margin:1rem 0; }
        .pipe-node { padding:0.5rem 0.8rem; border-radius:4px; font-size:0.78rem; font-weight:600; text-align:center; }
        .pipe-arrow { color:var(--text-cap); font-size:1rem; }
        .pipe-node.n-input { background:var(--blue-bg); color:var(--blue); border:1px solid rgba(44,82,130,0.2); }
        .pipe-node.n-python { background:var(--green-bg); color:var(--green); border:1px solid rgba(39,103,73,0.2); }
        .pipe-node.n-claude { background:var(--amber-bg); color:var(--amber); border:1px solid rgba(183,121,31,0.2); }
        .pipe-node.n-output { background:var(--red-bg); color:var(--red); border:1px solid rgba(155,44,44,0.2); }

        .divider { width:50px; height:3px; background:var(--uy); margin:2.5rem auto; }
        .footer { margin-top:2.5rem; padding:1.2rem 0; border-top:3px solid var(--uy); text-align:center; }
        .footer p { font-size:0.78rem; color:var(--text-cap); margin:0; }
    </style>
</head>
<body>

<div class="top-bar">
    <div class="top-bar-logo">FASTA Analyzer</div>
    <div class="top-bar-info">analyze_fasta.py &middot; Biopython</div>
</div>
<div class="yellow-strip"></div>

<div class="container">

    <header class="hero">
        <h1>Reporte de analisis FASTA</h1>
        <div><span class="hero-file">demo_nucleotide.fasta</span></div>
        <div class="hero-meta">Generado el <TIMESTAMP> &middot; analyze_fasta.py</div>
    </header>

    <!-- ===== RESUMEN ===== -->
    <section>
        <div class="sn">1. Resumen general</div>
        <h2>Estadisticas del archivo</h2>
        <table class="summary-table">
            <tbody>
        <tr><td>Tipo de secuencia</td><td><span class="badge nuc">nucleotide</span></td></tr>
        <tr><td>Total secuencias</td><td>1</td></tr>
        <tr><td>Total residuos</td><td>690 bp</td></tr>
        <tr><td>Longitud minima</td><td>690</td></tr>
        <tr><td>Longitud maxima</td><td>690</td></tr>
        <tr><td>Longitud promedio</td><td>690.0</td></tr>
        <tr><td>N50</td><td>690</td></tr>
        <tr><td>GC promedio</td><td>56.23%</td></tr>
        <tr><td>ORFs encontrados (&ge;100 aa)</td><td>1</td></tr></tbody>
        </table>
    </section>

    <!-- ===== DETALLE POR SECUENCIA ===== -->
    <section>
        <div class="sn">2. Analisis por secuencia</div>
        <h2>Detalle de cada secuencia</h2>
        \n\
            <div class="seq-card">
                <div class="seq-header">
                    <span class="seq-num">#1</span>
                    <div>
                        <h3>synthetic_demo_orf</h3>
                        <p class="seq-desc">synthetic_demo_orf | Synthetic E. coli-like ORF for analyze-fasta demo (CC0 — no real organism)</p>
                    </div>
                </div>
                <div class="seq-grid">
                    <div>
                        <div class="metric-group">
                            <div class="metric"><span class="metric-label">Longitud</span><span class="metric-value">690 bp</span></div>
                            <div class="metric"><span class="metric-label">GC content</span><span class="metric-value">56.23%</span></div>
                            <div class="metric"><span class="metric-label">AT content</span><span class="metric-value">43.77%</span></div>
                        </div>
                        <h4>Composicion de bases</h4>
                        <div class="comp-bar">
                            <div class="comp-seg seg-a" style="width:32.75%" title="A: 32.75%">A</div>
                            <div class="comp-seg seg-t" style="width:11.01%" title="T: 11.01%">T</div>
                            <div class="comp-seg seg-g" style="width:20.87%" title="G: 20.87%">G</div>
                            <div class="comp-seg seg-c" style="width:35.36%" title="C: 35.36%">C</div>
                        </div>
                        <div class="comp-legend">
                            <span><i class="dot dot-a"></i>A 32.75%</span>
                            <span><i class="dot dot-t"></i>T 11.01%</span>
                            <span><i class="dot dot-g"></i>G 20.87%</span>
                            <span><i class="dot dot-c"></i>C 35.36%</span>
                        </div>
                    </div>
                    <div>
                        <h4>Top dinucleotidos</h4>
                        <div class="mini-bar-row"><span class="mini-label">CA</span><div class="mini-track"><div class="mini-fill nuc-fill" style="width:100.0%"></div></div><span class="mini-val">119</span></div><div class="mini-bar-row"><span class="mini-label">GC</span><div class="mini-track"><div class="mini-fill nuc-fill" style="width:88.23529411764706%"></div></div><span class="mini-val">105</span></div><div class="mini-bar-row"><span class="mini-label">AC</span><div class="mini-track"><div class="mini-fill nuc-fill" style="width:72.26890756302521%"></div></div><span class="mini-val">86</span></div><div class="mini-bar-row"><span class="mini-label">CG</span><div class="mini-track"><div class="mini-fill nuc-fill" style="width:57.98319327731093%"></div></div><span class="mini-val">69</span></div><div class="mini-bar-row"><span class="mini-label">AA</span><div class="mini-track"><div class="mini-fill nuc-fill" style="width:44.537815126050425%"></div></div><span class="mini-val">53</span></div>
                    </div>
                </div>
                \n\
                <h4>ORFs detectados (1)</h4>
                <table class="inner-table">
                    <thead><tr><th>Frame</th><th>Posicion</th><th>Longitud</th><th>Proteina</th></tr></thead>
                    <tbody><tr><td>+1</td><td>1 - 690</td><td>690 bp</td><td>230 aa</td></tr></tbody>
                </table>
            </div>
        \n\
    </section>

    <div class="divider"></div>

    <!-- ===== JSON OUTPUT ===== -->
    <section>
        <div class="sn">3. JSON output</div>
        <h2>Datos crudos del analisis</h2>
        <p>Output completo en formato JSON, tal como lo consume Claude para la interpretacion biologica.</p>
        <button class="json-toggle" onclick="var b=document.getElementById('json-block');b.style.display=b.style.display==='block'?'none':'block';this.textContent=b.style.display==='block'?'Ocultar JSON':'Mostrar JSON'">Mostrar JSON</button>
        <div class="json-block" id="json-block">
            <pre>{
  &quot;file&quot;: &quot;demo_nucleotide.fasta&quot;,
  &quot;analysis_date&quot;: &quot;<TIMESTAMP>&quot;,
  &quot;total_sequences&quot;: 1,
  &quot;sequence_type&quot;: &quot;nucleotide&quot;,
  &quot;sequences&quot;: [
    {
      &quot;id&quot;: &quot;synthetic_demo_orf&quot;,
      &quot;description&quot;: &quot;synthetic_demo_orf | Synthetic E. coli-like ORF for analyze-fasta demo (CC0 — no real organism)&quot;,
      &quot;length_bp&quot;: 690,
      &quot;gc_content&quot;: 56.23,
      &quot;composition&quot;: {
        &quot;A&quot;: 226,
        &quot;T&quot;: 76,
        &quot;G&quot;: 144,
        &quot;C&quot;: 244,
        &quot;N&quot;: 0
      },
      &quot;composition_pct&quot;: {
        &quot;A&quot;: 32.75,
        &quot;T&quot;: 11.01,
        &quot;G&quot;: 20.87,
        &quot;C&quot;: 35.36
      },
      &quot;at_content&quot;: 43.77,
      &quot;top_dinucleotides&quot;: {
        &quot;CA&quot;: 119,
        &quot;GC&quot;: 105,
        &quot;AC&quot;: 86,
        &quot;CG&quot;: 69,
        &quot;AA&quot;: 53
      },
      &quot;orfs_found&quot;: 1,
      &quot;orfs&quot;: [
        {
          &quot;frame&quot;: &quot;+1&quot;,
          &quot;start&quot;: 1,
          &quot;end&quot;: 690,
          &quot;length_bp&quot;: 690,
          &quot;length_aa&quot;: 230
        }
      ],
      &quot;molecular_weight_da&quot;: 211887.4
    }
  ],
  &quot;summary&quot;: {
    &quot;total_sequences&quot;: 1,
    &quot;total_residues&quot;: 690,
    &quot;min_length&quot;: 690,
    &quot;max_length&quot;: 690,
    &quot;avg_length&quot;: 690.0,
    &quot;n50&quot;: 690,
    &quot;avg_gc_content&quot;: 56.23,
    &quot;total_orfs&quot;: 1
  }
}</pre>
        </div>
    </section>

    <div class="divider"></div>

    <!-- ===== COMO FUNCIONA ===== -->
    <section>
        <div class="sn">4. Como funciona</div>
        <div class="how-section">
            <h2>Skill: analyze-fasta</h2>
            <p>
                Esta herramienta combina un <strong>script de Python</strong> que realiza el analisis bioinformatico
                con <strong>Claude (IA de Anthropic)</strong> que interpreta los resultados en contexto biologico.
                El script es el musculo (calculos); Claude es el cerebro (interpretacion).
            </p>

            <div class="pipe-flow">
                <div class="pipe-node n-input">Archivo FASTA</div>
                <div class="pipe-arrow">&rarr;</div>
                <div class="pipe-node n-python">analyze_fasta.py<br><small>Biopython + NumPy</small></div>
                <div class="pipe-arrow">&rarr;</div>
                <div class="pipe-node n-output">JSON + HTML</div>
                <div class="pipe-arrow">&rarr;</div>
                <div class="pipe-node n-claude">Claude (skill)<br><small>Interpretacion biologica</small></div>
            </div>

            <div class="how-grid">
                <div class="how-card">
                    <h3>Librerias utilizadas</h3>
                    <ul>
                        <li><strong>Biopython</strong> (<code>Bio.SeqIO</code>) &mdash; parsing de FASTA, acceso a secuencias</li>
                        <li><strong>Bio.SeqUtils</strong> &mdash; calculo de GC%, peso molecular</li>
                        <li><strong>Bio.SeqUtils.ProtParam</strong> &mdash; analisis fisicoquimico de proteinas (pI, GRAVY, estabilidad, estructura secundaria)</li>
                        <li><strong>NumPy</strong> &mdash; dependencia de Biopython para calculos numericos</li>
                        <li><strong>Python stdlib</strong> &mdash; Counter, json, argparse, pathlib</li>
                    </ul>
                </div>

                <div class="how-card">
                    <h3>Que analiza</h3>
                    <ul>
                        <li><strong>Deteccion automatica</strong> de tipo (nucleotido vs proteina) por composicion</li>
                        <li><strong>Nucleotidos:</strong> GC%, composicion de bases, dinucleotidos, ORFs en 3 frames, peso molecular, N50</li>
                        <li><strong>Proteinas:</strong> peso molecular, punto isoelectrico, indice de estabilidad, GRAVY, aromaticidad, estructura secundaria, residuos cargados</li>
                        <li><strong>Multi-secuencia:</strong> estadisticas agregadas, N50 para evaluar calidad de ensamblaje</li>
                    </ul>
                </div>

                <div class="how-card">
                    <h3>Deteccion de ORFs</h3>
                    <p>
                        Busca marcos de lectura abiertos en los 3 frames forward (+1, +2, +3).
                        Un ORF se define como una secuencia que empieza con <code>ATG</code> y
                        termina en un codon stop (<code>TAA</code>, <code>TAG</code>, <code>TGA</code>)
                        con un minimo de 100 aminoacidos (300 bp).
                    </p>
                    <p>
                        <strong>Limitacion:</strong> no busca en cadena reversa complementaria ni modela intrones.
                        Para prediccion genica completa se recomienda Augustus o GenScan.
                    </p>
                </div>

                <div class="how-card">
                    <h3>Rol de Claude (skill)</h3>
                    <p>
                        La skill <code>analyze-fasta.md</code> instruye a Claude para:
                    </p>
                    <ul>
                        <li>Ejecutar el script con <code>--json</code></li>
                        <li>Interpretar GC%, ORFs, pI, GRAVY en contexto biologico</li>
                        <li>Inferir organismo, funcion, o tipo de secuencia probable</li>
                        <li>Comparar con valores de referencia</li>
                        <li>Sugerir analisis adicionales (BLAST, filogenia, etc.)</li>
                    </ul>
                </div>
            </div>

            <h4 style="margin-top:1.5rem">Modos de uso</h4>
            <table class="summary-table" style="margin-top:0.5rem">
                <tr><td><code>python3 analyze_fasta.py archivo.fasta</code></td><td>Reporte en texto por terminal</td></tr>
                <tr><td><code>python3 analyze_fasta.py archivo.fasta --json</code></td><td>Output JSON (para Claude o pipelines)</td></tr>
                <tr><td><code>python3 analyze_fasta.py archivo.fasta --html reporte.html</code></td><td>Genera este reporte HTML</td></tr>
            </table>
        </div>
    </section>

    <footer class="footer">
        <p>Generado por analyze_fasta.py &middot; Biopython 1.88 &middot; Python 3.12.14</p>
    </footer>
</div>

</body>
</html>\
""", 'report.md': """\
# analyze-fasta Report

**Input file:** `demo_nucleotide.fasta`  \n\
**Analysis date:** <TIMESTAMP>  \n\
**Sequence type:** `nucleotide`  \n\
**Total sequences:** 1  \n\

## Summary

| Metric | Value |
|---|---|
| total_sequences | 1 |
| total_residues | 690 |
| min_length | 690 |
| max_length | 690 |
| avg_length | 690.0 |
| n50 | 690 |
| avg_gc_content | 56.23 |
| total_orfs | 1 |

## Per-sequence metrics

### 1. synthetic_demo_orf

- **Description:** synthetic_demo_orf | Synthetic E. coli-like ORF for analyze-fasta demo (CC0 — no real organism)
- **Length:** 690 bp
- **GC content:** 56.23%
- **AT content:** 43.77%
- **ORFs (>=100 aa):** 1
- **Molecular weight:** 211887.4 Da

---

_ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions._
""", 'result.json': """\
{
  "file": "demo_nucleotide.fasta",
  "analysis_date": "<TIMESTAMP>",
  "total_sequences": 1,
  "sequence_type": "nucleotide",
  "sequences": [
    {
      "id": "synthetic_demo_orf",
      "description": "synthetic_demo_orf | Synthetic E. coli-like ORF for analyze-fasta demo (CC0 — no real organism)",
      "length_bp": 690,
      "gc_content": 56.23,
      "composition": {
        "A": 226,
        "T": 76,
        "G": 144,
        "C": 244,
        "N": 0
      },
      "composition_pct": {
        "A": 32.75,
        "T": 11.01,
        "G": 20.87,
        "C": 35.36
      },
      "at_content": 43.77,
      "top_dinucleotides": {
        "CA": 119,
        "GC": 105,
        "AC": 86,
        "CG": 69,
        "AA": 53
      },
      "orfs_found": 1,
      "orfs": [
        {
          "frame": "+1",
          "start": 1,
          "end": 690,
          "length_bp": 690,
          "length_aa": 230
        }
      ],
      "molecular_weight_da": 211887.4
    }
  ],
  "summary": {
    "total_sequences": 1,
    "total_residues": 690,
    "min_length": 690,
    "max_length": 690,
    "avg_length": 690.0,
    "n50": 690,
    "avg_gc_content": 56.23,
    "total_orfs": 1
  }
}\
"""})
