"""Length-contract tests for the shared gi-* runner.

These are pure-local: they exercise the gate that keeps a badly sized FASTA
from ever reaching the Genomic Intelligence API, so they need no key and no
network. The numbers mirror the published request schemas (one ``minLength``
per task); if the API moves, these fail first.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from clawbio.gi.gi_runner import (
    REQUEST_MAX_BP,
    TASK_CONTEXT_WINDOW_BP,
    TASK_MIN_BP,
    regime_note,
    validate_sequence_length,
)

REPO_ROOT = Path(__file__).resolve().parents[2]

TASKS = ["promoter", "splice", "enhancer", "chromatin", "annotation", "expression"]


def test_every_task_has_a_floor_and_a_context_entry():
    assert sorted(TASK_MIN_BP) == sorted(TASKS)
    assert sorted(TASK_CONTEXT_WINDOW_BP) == sorted(TASKS)


def test_floors_match_the_published_minlength():
    """Mirrors ``minLength`` on each task's request schema."""
    assert TASK_MIN_BP == {
        "promoter": 300,
        "splice": 100,
        "enhancer": 50,
        "chromatin": 200,
        "annotation": 1000,
        "expression": 9198,
    }
    assert REQUEST_MAX_BP == 500_000


@pytest.mark.parametrize("task", TASKS)
def test_one_bp_under_the_floor_is_rejected(task):
    floor = TASK_MIN_BP[task]
    problem = validate_sequence_length(task, "A" * (floor - 1))
    assert problem is not None
    assert f"at least {floor} bp" in problem


@pytest.mark.parametrize("task", TASKS)
def test_exactly_the_floor_is_accepted(task):
    assert validate_sequence_length(task, "A" * TASK_MIN_BP[task]) is None


@pytest.mark.parametrize("task", TASKS)
def test_over_the_shared_maximum_is_rejected(task):
    problem = validate_sequence_length(task, "A" * (REQUEST_MAX_BP + 1))
    assert problem is not None and "500000 bp" in problem


def test_a_short_promoter_sequence_is_not_rejected_by_the_old_one_bp_rule():
    """An earlier contract published a single 1 bp floor for every task.

    300 bp is the promoter floor, and 299 bp must fail locally rather than
    cost a request.
    """
    assert validate_sequence_length("promoter", "A" * 299) is not None
    assert validate_sequence_length("promoter", "A" * 300) is None


def test_in_regime_note_fires_below_the_context_window_but_does_not_reject():
    """Above the floor, below the context window: accepted, scored on padding."""
    # 100 bp is well above the enhancer floor of 50 and below the 249 bp window.
    assert validate_sequence_length("enhancer", "A" * 100) is None
    note = regime_note("enhancer", 100)
    assert note is not None and "249" in note and "padded" in note


def test_in_regime_note_is_silent_at_or_above_the_context_window():
    assert regime_note("enhancer", 249) is None
    assert regime_note("promoter", 2000) is None


def test_tasks_without_a_sliding_window_never_produce_a_regime_note():
    """annotation and expression report ``context_window_bp: null``."""
    assert TASK_CONTEXT_WINDOW_BP["annotation"] is None
    assert TASK_CONTEXT_WINDOW_BP["expression"] is None
    assert regime_note("annotation", 1_000) is None
    assert regime_note("expression", 9_198) is None


@pytest.mark.parametrize(
    "task,script",
    [(t, REPO_ROOT / "skills" / f"gi-{t}" / f"gi_{t}.py") for t in TASKS],
)
def test_cli_rejects_a_too_short_fasta_without_calling_the_api(task, script, tmp_path, monkeypatch):
    """The gate must fire before the client is built — no key, no network."""
    floor = TASK_MIN_BP[task]
    short = tmp_path / "short.fa"
    short.write_text(">short\n" + "ACGT" * ((floor - 1) // 4) + "\n")
    monkeypatch.delenv("GI_API_KEY", raising=False)
    out = tmp_path / "out"
    result = subprocess.run(
        [sys.executable, str(script), "--input", str(short), "--output", str(out)],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 1, f"stdout={result.stdout} stderr={result.stderr}"
    assert f"at least {floor} bp" in result.stderr
    assert not out.exists()


def test_description_is_not_forwarded_to_a_task_whose_options_forbid_it(tmp_path, monkeypatch):
    """Every ``*Options`` is additionalProperties:false — a stray key is a 422.

    Only expression declares ``description``, so passing it elsewhere must be
    dropped with a warning, not forwarded.
    """
    fasta = tmp_path / "seq.fa"
    fasta.write_text(">seq\n" + "ACGT" * 200 + "\n")  # 800 bp, above the promoter floor
    monkeypatch.delenv("GI_API_KEY", raising=False)
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "skills" / "gi-promoter" / "gi_promoter.py"),
         "--input", str(fasta), "--output", str(tmp_path / "out"),
         "--description", "some cell type"],
        capture_output=True, text=True, timeout=60,
    )
    # No API key, so it cannot get further than building the client — but the
    # warning must already have been emitted.
    assert "--description applies only to gi-expression" in result.stderr


def test_a_long_fasta_header_is_truncated_to_the_schema_cap(tmp_path, monkeypatch):
    """``sequence_name`` is maxLength 128 — a long header must not cause a 422."""
    from clawbio.gi.gi_runner import SEQUENCE_NAME_MAX_CHARS

    assert SEQUENCE_NAME_MAX_CHARS == 128
    fasta = tmp_path / "long_header.fa"
    fasta.write_text(">" + "n" * 400 + "\n" + "ACGT" * 200 + "\n")
    monkeypatch.delenv("GI_API_KEY", raising=False)
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "skills" / "gi-promoter" / "gi_promoter.py"),
         "--input", str(fasta), "--output", str(tmp_path / "out")],
        capture_output=True, text=True, timeout=60,
    )
    assert "truncated to 128 characters" in result.stderr


class TestInputRefusal:
    """The client refuses malformed FASTA instead of silently repairing it.

    Both behaviours shipped: non-ACGTN characters were deleted (shifting every
    downstream coordinate) and multi-record files were concatenated into one
    chimeric sequence. Each produced a confident, well-formed, wrong result.
    """

    def _write(self, tmp_path, content):
        p = tmp_path / "in.fa"
        p.write_text(content)
        return p

    def test_accepts_clean_single_record(self, tmp_path):
        from clawbio.gi.gi_client import read_fasta
        name, seq = read_fasta(self._write(tmp_path, ">chr1 desc\nACGT\nacgt\n"))
        assert name == "chr1" and seq == "ACGTACGT"

    def test_rejects_iupac_ambiguity(self, tmp_path):
        import pytest
        from clawbio.gi.gi_client import FastaError, read_fasta
        with pytest.raises(FastaError) as exc:
            read_fasta(self._write(tmp_path, ">x\nACGTRYKM\n"))
        assert "outside ACGTN" in str(exc.value) and "IUPAC" in str(exc.value)

    def test_rejects_multi_record(self, tmp_path):
        import pytest
        from clawbio.gi.gi_client import FastaError, read_fasta
        with pytest.raises(FastaError) as exc:
            read_fasta(self._write(tmp_path, ">a\nACGT\n>b\nTTTT\n"))
        assert "single FASTA record" in str(exc.value)

    def test_fasta_error_is_value_error(self):
        from clawbio.gi.gi_client import FastaError
        assert issubclass(FastaError, ValueError)


class TestBundledFixtureHeadersAreOneBasedInclusive:
    """Every demo header's coordinate span must equal its own base count.

    The headers are documentation a user reads to derive an offset, and
    --tss-index makes offsets load-bearing: an offset taken from a header that
    is one out lands one base off, and the API answers confidently either way.
    expression_hbb_k562.fa stated chr11:5222472-5231670, which is 9,199 bases
    1-based inclusive for a 9,198 base file, because its start was 0-based
    while the other five were 1-based.

    Settled against GRCh38 rather than by moving whichever end looked wrong:
    the file is the exact reverse complement of Ensembl's
    chr11:5222473-5231670, so the start moved, not the end. This asserts the
    arithmetic for all six, which is the part that can drift again.
    """

    @pytest.mark.parametrize("task", TASKS)
    def test_header_span_matches_the_base_count(self, task):
        import re
        fixtures = sorted((REPO_ROOT / "skills" / f"gi-{task}" / "example_data").glob("*.fa"))
        assert fixtures, f"gi-{task} has no bundled fixture"
        for fa in fixtures:
            lines = fa.read_text().splitlines()
            header, bases = lines[0], "".join(l.strip() for l in lines[1:])
            m = re.search(r":(\d+)-(\d+)", header)
            assert m, f"{fa.name}: header states no coordinate span"
            start, end = int(m.group(1)), int(m.group(2))
            assert end - start + 1 == len(bases), (
                f"{fa.name}: header spans {end - start + 1} bases 1-based inclusive "
                f"but the file holds {len(bases)}"
            )
