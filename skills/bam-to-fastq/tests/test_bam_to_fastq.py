"""Tests for the bam-to-fastq skill.

Unit tests build tiny BAMs in tmp_path with samtools and need samtools and
bedtools on PATH; they skip when either is missing.
"""

from __future__ import annotations

import gzip
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from hypothesis import given
from hypothesis import strategies as st

SKILL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL_DIR))

import bam_to_fastq as b

needs_tools = pytest.mark.skipif(
    not (shutil.which("samtools") and shutil.which("bedtools")),
    reason="samtools and bedtools are required",
)

HEADER = "@HD\tVN:1.6\tSO:{so}\n@SQ\tSN:c1\tLN:1000\n"


def make_bam(path: Path, records: list[str], so: str = "coordinate") -> Path:
    sam = path.with_suffix(".sam")
    sam.write_text(HEADER.format(so=so) + "\n".join(records) + "\n")
    subprocess.run(["samtools", "view", "-b", "-o", str(path), str(sam)], check=True)
    return path


def pair(name: str, p1: int, p2: int) -> list[str]:
    return [
        f"{name}\t99\tc1\t{p1}\t60\t4M\t=\t{p2}\t44\tACGT\tIIII",
        f"{name}\t147\tc1\t{p2}\t60\t4M\t=\t{p1}\t-44\tTTGA\tIIII",
    ]


def coord_sorted_pairs(n: int = 3) -> list[str]:
    recs = []
    for i in range(n):
        recs += pair(f"r{i}", 10 + i * 5, 100 + i * 5)
    return sorted(recs, key=lambda r: int(r.split("\t")[3]))


def read_names(path: Path) -> list[str]:
    with gzip.open(path, "rt") as fh:
        return [line.split()[0] for i, line in enumerate(fh) if i % 4 == 0]


@needs_tools
class TestPaired:
    def test_name_sort_converts_and_counts(self, tmp_path):
        bam = make_bam(tmp_path / "s.bam", coord_sorted_pairs())
        res = b.convert(bam, tmp_path / "out", paired=True, name_sort=True)
        assert res["counts"] == {"input_bam": 6, "name_sorted_bam": 6, "fastq": [3, 3]}
        r1, r2 = (tmp_path / "out" / "fastq" / f"s_R{i}.fastq.gz" for i in (1, 2))
        assert [n.rstrip("/12") for n in read_names(r1)] == [n.rstrip("/12") for n in read_names(r2)]

    def test_coordinate_sorted_without_name_sort_is_rejected(self, tmp_path):
        bam = make_bam(tmp_path / "s.bam", coord_sorted_pairs())
        with pytest.raises(b.BamToFastqError, match="not name-sorted"):
            b.convert(bam, tmp_path / "out", paired=True)
        assert not list((tmp_path / "out" / "fastq").glob("*.fastq.gz"))

    def test_name_sorted_bam_passes_sort_check(self, tmp_path):
        recs = pair("a", 10, 20) + pair("b", 12, 22)
        bam = make_bam(tmp_path / "s.bam", recs, so="queryname")
        res = b.convert(bam, tmp_path / "out", paired=True)
        assert res["sort_check"]["grouped"] is True
        assert res["counts"]["fastq"] == [2, 2]

    def test_odd_record_count_rejected(self, tmp_path):
        bam = make_bam(tmp_path / "s.bam", pair("a", 10, 20) + [pair("b", 12, 22)[0]],
                       so="queryname")
        with pytest.raises(b.BamToFastqError, match="odd primary record count"):
            b.convert(bam, tmp_path / "out", paired=True)

    def test_secondary_alignments_not_counted(self, tmp_path):
        recs = pair("a", 10, 20) + ["a\t355\tc1\t30\t0\t4M\t=\t20\t0\tACGT\tIIII"]
        bam = make_bam(tmp_path / "s.bam", recs, so="queryname")
        assert b.count_bam(bam) == 2


@needs_tools
class TestSingleEnd:
    def records(self):
        return [
            "s1\t0\tc1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII",
            "s2\t16\tc1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII",
        ]

    def test_each_read_written_once(self, tmp_path):
        bam = make_bam(tmp_path / "se.bam", self.records())
        res = b.convert(bam, tmp_path / "out")
        assert res["tool"] == "samtools fastq"
        assert res["counts"]["fastq"] == [2]
        assert len(read_names(tmp_path / "out" / "fastq" / "se.fastq.gz")) == 2

    def test_paired_bam_without_paired_flag_fails_count_check(self, tmp_path):
        bam = make_bam(tmp_path / "s.bam", coord_sorted_pairs(1))
        with pytest.raises(b.BamToFastqError, match="count mismatch"):
            b.convert(bam, tmp_path / "out")
        assert not list((tmp_path / "out" / "fastq").glob("*.fastq.gz"))


@needs_tools
class TestSafety:
    def test_refuses_to_overwrite(self, tmp_path):
        bam = make_bam(tmp_path / "se.bam", TestSingleEnd().records())
        b.convert(bam, tmp_path / "out")
        with pytest.raises(b.BamToFastqError, match="refusing to overwrite"):
            b.convert(bam, tmp_path / "out")

    def test_missing_input(self, tmp_path):
        with pytest.raises(b.BamToFastqError, match="input not found"):
            b.convert(tmp_path / "nope.bam", tmp_path / "out")

    def test_no_scratch_left_behind(self, tmp_path):
        bam = make_bam(tmp_path / "se.bam", TestSingleEnd().records())
        b.convert(bam, tmp_path / "out")
        assert not list((tmp_path / "out").glob(".tmp-*"))

    def test_output_is_byte_reproducible(self, tmp_path):
        bam = make_bam(tmp_path / "s.bam", coord_sorted_pairs())
        a = b.convert(bam, tmp_path / "o1", paired=True, name_sort=True)
        c = b.convert(bam, tmp_path / "o2", paired=True, name_sort=True)
        assert [f["sha256"] for f in a["files"]] == [f["sha256"] for f in c["files"]]


class TestHelpers:
    def test_fastq_names_normalises(self, tmp_path):
        fq = tmp_path / "x.fastq"
        fq.write_text("@r1/1 comment\nA\n+\nI\n@r2\nC\n+\nI\n")
        assert list(b.fastq_names(fq)) == ["r1", "r2"]

    def test_pair_name_mismatch_detected(self, tmp_path):
        r1, r2 = tmp_path / "1.fastq", tmp_path / "2.fastq"
        r1.write_text("@a/1\nA\n+\nI\n")
        r2.write_text("@b/2\nA\n+\nI\n")
        with pytest.raises(b.BamToFastqError, match="names differ"):
            b.check_pair_names(r1, r2)


class TestProperties:
    """Property tests for the FASTQ name parser."""

    name = st.text(
        alphabet=st.characters(min_codepoint=33, max_codepoint=126, blacklist_characters="@/"),
        min_size=1, max_size=40,
    )

    @given(name=name, mate=st.sampled_from(["", "/1", "/2"]),
           comment=st.sampled_from(["", " 1:N:0:ACGT", " BC:Z:AAA"]))
    def test_mate_suffix_and_comment_never_change_the_key(self, tmp_path_factory, name, mate, comment):
        fq = tmp_path_factory.mktemp("fq") / "x.fastq"
        fq.write_text(f"@{name}{mate}{comment}\nA\n+\nI\n")
        assert list(b.fastq_names(fq)) == [name]

    @given(n=st.integers(min_value=0, max_value=50))
    def test_record_count_matches_records_written(self, tmp_path_factory, n):
        fq = tmp_path_factory.mktemp("fq") / "x.fastq"
        fq.write_text("".join(f"@r{i}/1\nA\n+\nI\n" for i in range(n)))
        assert b.count_fastq(fq) == n


@needs_tools
class TestDemoContract:
    """The demo must produce every file SKILL.md's Output Structure promises."""

    def test_demo_outputs(self, tmp_path):
        rc = b.main(["--demo", "--output", str(tmp_path / "demo")])
        assert rc == 0
        out = tmp_path / "demo"
        for rel in ["report.md", "result.json", "fastq/synthetic_demo_R1.fastq.gz",
                    "fastq/synthetic_demo_R2.fastq.gz", "reproducibility/commands.sh",
                    "reproducibility/environment.yml", "reproducibility/checksums.sha256"]:
            assert (out / rel).is_file(), rel
        result = json.loads((out / "result.json").read_text())
        assert result["counts"]["fastq"] == [40, 40]
        assert "not a medical device" in (out / "report.md").read_text()
