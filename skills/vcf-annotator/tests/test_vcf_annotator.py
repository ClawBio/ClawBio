"""
Tests for ClawBio VCF Annotator skill.
Run with: pytest skills/vcf-annotator/tests/test_vcf_annotator.py -v
"""

import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest
from hypothesis import HealthCheck, example, given, settings, strategies as st
from types import SimpleNamespace
from unittest.mock import Mock

from inline_snapshot import snapshot

from clawbio.common import demo_runner

sys.path.insert(0, str(Path(__file__).parent.parent))
import vcf_annotator as annotator
from vcf_annotator import (
    DEMO_ANNOTATIONS,
    DEMO_VCF_CONTENT,
    IMPACT_RANK,
    generate_report,
    parse_vcf,
)

# ── Demo VCF fixture ───────────────────────────────────────────────────────────

SAMPLE_VCF = """##fileformat=VCFv4.2
##reference=GRCh38
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
17\t43063931\trs55770810\tG\tA\t.\tPASS\t.
13\t32316461\trs80359550\tC\tT\t.\tPASS\t.
chr7\t117548628\t.\tCTTT\tC\t.\tPASS\t.
"""


@pytest.fixture
def sample_vcf(tmp_path):
    p = tmp_path / "sample.vcf"
    p.write_text(SAMPLE_VCF)
    return p


# ── VCF parser tests ───────────────────────────────────────────────────────────

class TestParseVCF:
    def test_returns_list(self, sample_vcf):
        variants = parse_vcf(sample_vcf)
        assert isinstance(variants, list)

    def test_correct_count(self, sample_vcf):
        variants = parse_vcf(sample_vcf)
        assert len(variants) == 3

    def test_required_fields_present(self, sample_vcf):
        variants = parse_vcf(sample_vcf)
        for v in variants:
            assert "chrom" in v
            assert "pos" in v
            assert "ref" in v
            assert "alt" in v

    def test_chr_prefix_stripped(self, sample_vcf):
        variants = parse_vcf(sample_vcf)
        for v in variants:
            assert not v["chrom"].startswith("chr")

    def test_dot_id_becomes_empty(self, sample_vcf):
        variants = parse_vcf(sample_vcf)
        # Third variant has no rsID
        assert variants[2]["id"] == ""

    def test_rsid_preserved(self, sample_vcf):
        variants = parse_vcf(sample_vcf)
        assert variants[0]["id"] == "rs55770810"

    def test_skips_header_lines(self, sample_vcf):
        variants = parse_vcf(sample_vcf)
        # Should not include any header lines
        for v in variants:
            assert not v["chrom"].startswith("#")

    def test_empty_vcf_returns_empty_list(self, tmp_path):
        empty = tmp_path / "empty.vcf"
        empty.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        variants = parse_vcf(empty)
        assert variants == []


# ── Impact ranking tests ───────────────────────────────────────────────────────

class TestImpactRank:
    def test_high_ranks_first(self):
        assert IMPACT_RANK["HIGH"] < IMPACT_RANK["MODERATE"]

    def test_moderate_before_low(self):
        assert IMPACT_RANK["MODERATE"] < IMPACT_RANK["LOW"]

    def test_low_before_modifier(self):
        assert IMPACT_RANK["LOW"] < IMPACT_RANK["MODIFIER"]

    def test_unknown_ranks_last(self):
        assert IMPACT_RANK["UNKNOWN"] >= IMPACT_RANK["MODIFIER"]


# ── Report generation tests ────────────────────────────────────────────────────

class TestGenerateReport:
    def test_creates_report_md(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            generate_report(DEMO_ANNOTATIONS, out)
            assert (out / "report.md").exists()

    def test_creates_results_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            generate_report(DEMO_ANNOTATIONS, out)
            assert (out / "results.json").exists()

    def test_creates_csv_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            generate_report(DEMO_ANNOTATIONS, out)
            assert (out / "tables" / "variants.csv").exists()

    def test_creates_reproducibility_bundle(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            generate_report(DEMO_ANNOTATIONS, out)
            assert (out / "reproducibility" / "commands.sh").exists()
            assert (out / "reproducibility" / "environment.yml").exists()
            assert (out / "reproducibility" / "checksums.sha256").exists()

    def test_results_json_is_valid(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            generate_report(DEMO_ANNOTATIONS, out)
            data = json.loads((out / "results.json").read_text())
            assert isinstance(data, list)
            assert len(data) == len(DEMO_ANNOTATIONS)

    def test_report_contains_gene_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            generate_report(DEMO_ANNOTATIONS, out)
            text = (out / "report.md").read_text()
            assert "BRCA1" in text
            assert "BRCA2" in text
            assert "CFTR" in text

    def test_report_contains_disclaimer(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            generate_report(DEMO_ANNOTATIONS, out)
            text = (out / "report.md").read_text()
            assert "research tool" in text.lower()

    def test_report_shows_impact_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            generate_report(DEMO_ANNOTATIONS, out)
            text = (out / "report.md").read_text()
            assert "HIGH impact" in text

    def test_csv_has_correct_headers(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            generate_report(DEMO_ANNOTATIONS, out)
            with open(out / "tables" / "variants.csv") as f:
                reader = csv.DictReader(f)
                headers = reader.fieldnames
            assert "gene" in headers
            assert "impact" in headers
            assert "clinvar_significance" in headers
            assert "gnomad_af" in headers

    def test_checksums_file_not_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            generate_report(DEMO_ANNOTATIONS, out)
            text = (out / "reproducibility" / "checksums.sha256").read_text()
            assert len(text.strip()) > 0

    def test_empty_variants_does_not_crash(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            generate_report([], out)
            assert (out / "report.md").exists()


# ── Demo data integrity tests ──────────────────────────────────────────────────

class TestDemoData:
    def test_demo_has_required_fields(self):
        required = [
            "chrom", "pos", "ref", "alt", "gene", "consequence",
            "impact", "clinvar_significance", "gnomad_af",
        ]
        for v in DEMO_ANNOTATIONS:
            for field in required:
                assert field in v, f"Demo variant missing: {field}"

    def test_demo_impacts_are_valid(self):
        valid = {"HIGH", "MODERATE", "LOW", "MODIFIER", "UNKNOWN"}
        for v in DEMO_ANNOTATIONS:
            assert v["impact"] in valid

    def test_demo_contains_high_impact(self):
        highs = [v for v in DEMO_ANNOTATIONS if v["impact"] == "HIGH"]
        assert len(highs) >= 1

    def test_demo_gnomad_af_is_numeric_or_none(self):
        for v in DEMO_ANNOTATIONS:
            af = v.get("gnomad_af")
            assert af is None or isinstance(af, float)

    def test_demo_sorted_by_impact(self):
        impacts = [v["impact"] for v in DEMO_ANNOTATIONS]
        ranks   = [IMPACT_RANK.get(i, 5) for i in impacts]
        assert ranks == sorted(ranks)

    def test_demo_grch38_coordinates(self):
        expected = {
            "rs55770810": ("17", "43063931"),
            "rs429358": ("19", "44908684"),
            "rs1801133": ("1", "11796321"),
        }
        by_id = {v["id"]: v for v in DEMO_ANNOTATIONS}
        for rsid, (chrom, pos) in expected.items():
            assert by_id[rsid]["chrom"] == chrom
            assert by_id[rsid]["pos"] == pos
        vcf = DEMO_VCF_CONTENT
        assert "19\t44908684\trs429358" in vcf
        assert "1\t11796321\trs1801133" in vcf
        assert "17\t43063931\trs55770810\tG\tA" in vcf

    def test_demo_brca1_is_missense_not_duplication(self):
        brca1 = next(v for v in DEMO_ANNOTATIONS if v["id"] == "rs55770810")
        assert (brca1["chrom"], brca1["pos"]) == ("17", "43063931")
        assert (brca1["ref"], brca1["alt"]) == ("G", "A")
        assert brca1["consequence"] == "missense_variant"
        assert brca1["impact"] == "MODERATE"
        assert brca1["hgvs"] == "NM_007294.4:c.5095C>T (p.Arg1699Trp)"

    def test_demo_vcf_rows_match_annotations(self):
        """Every DEMO_ANNOTATIONS record must agree with its DEMO_VCF_CONTENT
        row. Cross-artifact consistency: the VCF and the annotation table are
        maintained separately, so echoing either one cannot catch drift.
        """
        rows = {}
        for line in DEMO_VCF_CONTENT.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            chrom, pos, rsid, ref, alt = line.split("\t")[:5]
            rows[rsid] = (chrom, pos, ref, alt)
        for v in DEMO_ANNOTATIONS:
            assert v["id"] in rows, f"{v['id']} missing from demo VCF"
            assert rows[v["id"]] == (v["chrom"], v["pos"], v["ref"], v["alt"]), (
                f"{v['id']}: VCF row {rows[v['id']]} disagrees with annotation "
                f"{(v['chrom'], v['pos'], v['ref'], v['alt'])}"
            )

    def test_demo_brca1_hgvs_matches_genomic_alleles(self):
        """Derive, do not echo: for a minus-strand gene the genomic REF/ALT
        must be the complement of the cDNA substitution, and the cDNA position
        must be the first base of the codon named by the protein change.

        External anchors (verified 2026-08-30, live):
          - ClinVar VCV000055396 = NM_007294.4:c.5095C>T (p.Arg1699Trp),
            dbSNP rs55770810, GRCh38 chr17:43063931
          - Ensembl VEP colocates rs55770810 (G/A) at 17:43063931,
            codon Cgg/Tgg (Arg->Trp)
        These checks would have failed on the audited inconsistent record
        (17:43106457 T>A labelled p.Arg1699Trp): complement of c.5095C>T is
        G>A, not T>A.
        """
        import re

        complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
        brca1 = next(v for v in DEMO_ANNOTATIONS if v["gene"] == "BRCA1")
        m = re.search(r"c\.(\d+)([ACGT])>([ACGT])", brca1["hgvs"])
        assert m, f"HGVS lacks a cDNA substitution: {brca1['hgvs']}"
        c_pos, c_ref, c_alt = int(m.group(1)), m.group(2), m.group(3)
        p = re.search(r"p\.\w+?(\d+)\w+?", brca1["hgvs"])
        assert p, f"HGVS lacks a protein residue number: {brca1['hgvs']}"
        residue = int(p.group(1))

        # BRCA1 is on the minus strand: genomic alleles complement the cDNA.
        assert brca1["ref"] == complement[c_ref], (
            f"minus-strand REF should complement c.{c_pos}{c_ref}>{c_alt}: "
            f"expected {complement[c_ref]}, got {brca1['ref']}"
        )
        assert brca1["alt"] == complement[c_alt], (
            f"minus-strand ALT should complement c.{c_pos}{c_ref}>{c_alt}: "
            f"expected {complement[c_alt]}, got {brca1['alt']}"
        )
        # c.5095 is the first base of codon 1699: 3k - 2.
        assert 3 * residue - 2 == c_pos, (
            f"cDNA position {c_pos} is not the first base of codon {residue}"
        )


# ── CLI entry point ────────────────────────────────────────────────────────────
#
# These invoke the script the way a user does, in a subprocess. Importing
# functions directly cannot catch a missing import inside main(), which is how
# `import argparse` stayed absent behind 28 green tests.

SCRIPT = Path(__file__).parent.parent / "vcf_annotator.py"


class TestCLI:
    def test_help_exits_zero(self):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--help"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, (
            f"--help failed with exit {result.returncode}:\n{result.stderr}"
        )
        assert "--input" in result.stdout

    def test_no_args_errors_cleanly(self):
        """Missing --input/--demo must be an argparse error, not a traceback."""
        result = subprocess.run(
            [sys.executable, str(SCRIPT)],
            capture_output=True, text=True,
        )
        assert result.returncode == 2
        assert "Traceback" not in result.stderr

    def test_demo_runs_end_to_end(self, tmp_path):
        out = tmp_path / "report"
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(out)],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, (
            f"--demo failed with exit {result.returncode}:\n{result.stderr}"
        )
        assert (out / "report.md").exists()

    def test_missing_input_file_exits_one(self, tmp_path):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--input", str(tmp_path / "nope.vcf")],
            capture_output=True, text=True,
        )
        assert result.returncode == 1
        assert "not found" in result.stderr.lower()


HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
ROW = "chr1\t1\t.\tAT\tA\t.\tPASS\t.\n"


@pytest.mark.parametrize("content, location", [
    ("", "EOF"),
    ("##fileformat=VCFv4.2\n", "EOF"),
    (HEADER.replace("\t", " "), "line 1"),
    (ROW, "line 1"),
    (HEADER.replace("POS", "POSITION"), "line 1"),
    (HEADER + ROW.replace("\t", " "), "line 2"),
    (HEADER + "1\t1\t.\tA\tT\n", "line 2"),
    *[(HEADER + "\t".join("" if i == column else cell
                          for i, cell in enumerate(ROW.rstrip("\n").split("\t"))) + "\n",
       "line 2") for column in [0, 1, 3, 4, 7]],
    (HEADER + ROW.replace("\t1\t", "\tx\t"), "line 2"),
    (HEADER + ROW.rstrip("\n") + "\tGT\n", "line 2"),
    (HEADER.rstrip("\n") + "\tFORMAT\tsample\n" + ROW, "line 2"),
])
@pytest.mark.parametrize("backend", [False, True])
def test_invalid_structure_precedes_backend(tmp_path, monkeypatch, content, location, backend):
    path = tmp_path / "invalid.vcf"
    path.write_text(content)
    reader = Mock(side_effect=AssertionError("backend must not run"))
    monkeypatch.setattr(annotator, "HAS_CYVCF2", backend)
    monkeypatch.setattr(annotator, "cyvcf2", SimpleNamespace(VCF=reader), raising=False)
    with pytest.raises(ValueError, match=f"not valid VCF.*{location}"):
        parse_vcf(path)
    reader.assert_not_called()


@settings(suppress_health_check=[HealthCheck.function_scoped_fixture], deadline=None)
@example(rows=[("A", "T", 1)], index=0, spaces=False)
@given(rows=st.lists(st.tuples(st.text("ACGT", min_size=1, max_size=3),
                              st.text("ACGT", min_size=1, max_size=3),
                              st.integers(1, 1000)), min_size=1, max_size=5),
       index=st.integers(0, 5), spaces=st.booleans())
def test_corrupt_record_never_returns_partial_results(tmp_path, monkeypatch, rows, index, spaces):
    monkeypatch.setattr(annotator, "HAS_CYVCF2", False)
    path = tmp_path / "property.vcf"
    records = [f"chr1\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\n"
               for ref, alt, pos in rows]
    path.write_text(HEADER + "".join(records))
    assert parse_vcf(path) == [dict(chrom="1", pos=str(pos), id="", ref=ref, alt=alt)
                               for ref, alt, pos in rows]
    index %= len(records) + 1
    bad = ROW.replace("\t", " ") if spaces else "1\t1\t.\n"
    records.insert(index, bad)
    path.write_text(HEADER + "".join(records))
    with pytest.raises(ValueError, match=f"line {index + 2}"):
        parse_vcf(path)


@pytest.mark.parametrize("samples", [False, True])
def test_valid_crlf_and_sample_columns(tmp_path, monkeypatch, samples):
    monkeypatch.setattr(annotator, "HAS_CYVCF2", False)
    suffix = "\tFORMAT\tsample" if samples else ""
    values = "\tGT\t0/1" if samples else ""
    path = tmp_path / "valid.vcf"
    path.write_bytes(("##fileformat=VCFv4.2\r\n" + HEADER.rstrip("\n") + suffix +
                      "\r\n" + ROW.rstrip("\n") + values + "\r\n").encode())
    assert parse_vcf(path) == [dict(chrom="1", pos="1", id="", ref="AT", alt="A")]


@pytest.mark.parametrize("case", ["public", "empty", "header", "first", "middle", "last", "unreadable"])
@pytest.mark.parametrize("existing", [False, True])
def test_cli_rejects_input_without_side_effects(tmp_path, monkeypatch, capsys, case, existing):
    path = tmp_path / "input.vcf"
    if case == "public":
        path = SCRIPT.parents[1] / "pharmgx-reporter" / "demo_patient.txt"
    elif case == "unreadable":
        path.mkdir()
    else:
        rows = [ROW, ROW]
        if case in ["first", "middle", "last"]:
            rows.insert(["first", "middle", "last"].index(case), "bad row\n")
        path.write_text("" if case == "empty" else HEADER if case == "header"
                        else HEADER + "".join(rows))
    out = tmp_path / "output"
    if existing:
        out.mkdir()
        (out / "sentinel").write_text("untouched")
    blocked = Mock(side_effect=AssertionError("must not annotate or write"))
    monkeypatch.setattr(annotator, "HAS_CYVCF2", False)
    for name in ["annotate_variants", "generate_report"]:
        monkeypatch.setattr(annotator, name, blocked)
    monkeypatch.setattr(annotator.urllib.request, "urlopen", blocked)
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), "--input", str(path), "--output", str(out)])
    with pytest.raises(SystemExit) as exc:
        annotator.main()
    assert exc.value.code == 1
    captured = capsys.readouterr()
    diagnostic = ("no variant records in input" if case == "header" else
                  "could not read VCF" if case == "unreadable" else "not valid VCF")
    assert diagnostic in captured.err
    assert "Pathogenic:" not in captured.out
    assert "Report:" not in captured.out
    blocked.assert_not_called()
    if existing:
        assert list(out.iterdir()) == [out / "sentinel"]
        assert (out / "sentinel").read_text() == "untouched"
    else:
        assert not out.exists()


def test_valid_cli_still_generates_report(sample_vcf, tmp_path, monkeypatch):
    out = tmp_path / "report"
    monkeypatch.setattr(annotator, "HAS_CYVCF2", False)
    annotation = Mock(return_value=DEMO_ANNOTATIONS)
    monkeypatch.setattr(annotator, "annotate_variants", annotation)
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), "--input", str(sample_vcf), "--output", str(out)])
    annotator.main()
    annotation.assert_called_once_with(parse_vcf(sample_vcf))
    assert json.loads((out / "results.json").read_text()) == DEMO_ANNOTATIONS


@pytest.mark.parametrize("fails", [False, True])
def test_optional_reader_is_closed(sample_vcf, monkeypatch, fails):
    record = SimpleNamespace(CHROM="chr1", POS=1, ID=None, REF="AT", ALT=["A"])
    def records():
        yield record
        if fails:
            raise ValueError("backend read failed")
    reader = SimpleNamespace(close=Mock())
    class Reader:
        close = reader.close
        def __iter__(self):
            return records()
    monkeypatch.setattr(annotator, "HAS_CYVCF2", True)
    monkeypatch.setattr(annotator, "cyvcf2", SimpleNamespace(VCF=lambda _: Reader()), raising=False)
    if fails:
        with pytest.raises(ValueError):
            parse_vcf(sample_vcf)
    else:
        assert parse_vcf(sample_vcf) == [dict(chrom="1", pos="1", id="", ref="AT", alt="A")]
    reader.close.assert_called_once()


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/vcf-annotator/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("vcf-annotator", tmp_path) == snapshot({'report.md': """\
# 🦖 ClawBio VCF Annotator Report

**Input**: demo_variants.vcf  \n\
**Date**: <TIMESTAMP> UTC  \n\
**Total variants**: 5  \n\
**HIGH impact**: 2 | **MODERATE**: 3 | **LOW**: 0<br>
**ClinVar Pathogenic/Likely Pathogenic**: 3

---

## Summary

Annotated 5 variants against Ensembl VEP, ClinVar, and gnomAD. 2 variant(s) carry HIGH impact consequences. 3 variant(s) are classified as Pathogenic or Likely Pathogenic in ClinVar. Variants are ranked by predicted functional impact below.

---

## Variant Table

| # | Gene | Variant | Consequence | Impact | ClinVar | gnomAD AF |
|---|------|---------|-------------|--------|---------|-----------|
| 1 | BRCA2 | 13:32316461 C>T | stop_gained | **HIGH** | Pathogenic | 0.000004 |
| 2 | CFTR | 7:117548628 CTTT>C | frameshift_variant | **HIGH** | Pathogenic | 0.021000 |
| 3 | BRCA1 | 17:43063931 G>A | missense_variant | **MODERATE** | Pathogenic | 0.000024 |
| 4 | APOE | 19:44908684 T>C | missense_variant | **MODERATE** | Risk factor | 0.147000 |
| 5 | MTHFR | 1:11796321 G>A | missense_variant | **MODERATE** | Benign/Likely benign | 0.312000 |

---

## Detailed Annotations

### 1. BRCA2 — 13:32316461 C>T

| Field | Value |
|-------|-------|
| **rsID** | rs80359550 |
| **HGVS** | NM_000059.4:c.9976A>T |
| **Consequence** | stop_gained |
| **Impact** | HIGH |
| **SIFT** | deleterious |
| **PolyPhen** | probably_damaging |
| **ClinVar** | Pathogenic |
| **ClinVar Condition** | Familial cancer of breast |
| **gnomAD AF (global)** | 0.000004 |
| **gnomAD AF (AFR)** | 0.000000 |
| **gnomAD AF (EUR)** | 0.000006 |

### 2. CFTR — 7:117548628 CTTT>C

| Field | Value |
|-------|-------|
| **rsID** | rs113993960 |
| **HGVS** | NM_000492.4:c.1521_1523delCTT |
| **Consequence** | frameshift_variant |
| **Impact** | HIGH |
| **SIFT** | deleterious |
| **PolyPhen** | N/A |
| **ClinVar** | Pathogenic |
| **ClinVar Condition** | Cystic fibrosis |
| **gnomAD AF (global)** | 0.021000 |
| **gnomAD AF (AFR)** | 0.000300 |
| **gnomAD AF (EUR)** | 0.033000 |

### 3. BRCA1 — 17:43063931 G>A

| Field | Value |
|-------|-------|
| **rsID** | rs55770810 |
| **HGVS** | NM_007294.4:c.5095C>T (p.Arg1699Trp) |
| **Consequence** | missense_variant |
| **Impact** | MODERATE |
| **SIFT** | deleterious |
| **PolyPhen** | probably_damaging |
| **ClinVar** | Pathogenic |
| **ClinVar Condition** | Hereditary breast and ovarian cancer syndrome |
| **gnomAD AF (global)** | 0.000024 |
| **gnomAD AF (AFR)** | N/A |
| **gnomAD AF (EUR)** | N/A |

### 4. APOE — 19:44908684 T>C

| Field | Value |
|-------|-------|
| **rsID** | rs429358 |
| **HGVS** | NM_000041.4:c.388T>C |
| **Consequence** | missense_variant |
| **Impact** | MODERATE |
| **SIFT** | tolerated |
| **PolyPhen** | benign |
| **ClinVar** | Risk factor |
| **ClinVar Condition** | Alzheimer disease |
| **gnomAD AF (global)** | 0.147000 |
| **gnomAD AF (AFR)** | 0.232000 |
| **gnomAD AF (EUR)** | 0.143000 |

### 5. MTHFR — 1:11796321 G>A

| Field | Value |
|-------|-------|
| **rsID** | rs1801133 |
| **HGVS** | NM_005957.5:c.665C>T |
| **Consequence** | missense_variant |
| **Impact** | MODERATE |
| **SIFT** | tolerated |
| **PolyPhen** | benign |
| **ClinVar** | Benign/Likely benign |
| **ClinVar Condition** | Homocystinuria |
| **gnomAD AF (global)** | 0.312000 |
| **gnomAD AF (AFR)** | 0.198000 |
| **gnomAD AF (EUR)** | 0.367000 |

---

---
*ClawBio VCF Annotator is a research tool. Not a clinical diagnostic device. Always consult a qualified geneticist for clinical decisions.*
""", 'results.json': """\
[
  {
    "chrom": "13",
    "pos": "32316461",
    "id": "rs80359550",
    "ref": "C",
    "alt": "T",
    "gene": "BRCA2",
    "consequence": "stop_gained",
    "impact": "HIGH",
    "clinvar_significance": "Pathogenic",
    "clinvar_condition": "Familial cancer of breast",
    "gnomad_af": 4e-06,
    "gnomad_af_afr": 0.0,
    "gnomad_af_eur": 6e-06,
    "hgvs": "NM_000059.4:c.9976A>T",
    "sift": "deleterious",
    "polyphen": "probably_damaging"
  },
  {
    "chrom": "7",
    "pos": "117548628",
    "id": "rs113993960",
    "ref": "CTTT",
    "alt": "C",
    "gene": "CFTR",
    "consequence": "frameshift_variant",
    "impact": "HIGH",
    "clinvar_significance": "Pathogenic",
    "clinvar_condition": "Cystic fibrosis",
    "gnomad_af": 0.021,
    "gnomad_af_afr": 0.0003,
    "gnomad_af_eur": 0.033,
    "hgvs": "NM_000492.4:c.1521_1523delCTT",
    "sift": "deleterious",
    "polyphen": "N/A"
  },
  {
    "chrom": "17",
    "pos": "43063931",
    "id": "rs55770810",
    "ref": "G",
    "alt": "A",
    "gene": "BRCA1",
    "consequence": "missense_variant",
    "impact": "MODERATE",
    "clinvar_significance": "Pathogenic",
    "clinvar_condition": "Hereditary breast and ovarian cancer syndrome",
    "gnomad_af": 2.4e-05,
    "gnomad_af_afr": null,
    "gnomad_af_eur": null,
    "hgvs": "NM_007294.4:c.5095C>T (p.Arg1699Trp)",
    "sift": "deleterious",
    "polyphen": "probably_damaging"
  },
  {
    "chrom": "19",
    "pos": "44908684",
    "id": "rs429358",
    "ref": "T",
    "alt": "C",
    "gene": "APOE",
    "consequence": "missense_variant",
    "impact": "MODERATE",
    "clinvar_significance": "Risk factor",
    "clinvar_condition": "Alzheimer disease",
    "gnomad_af": 0.147,
    "gnomad_af_afr": 0.232,
    "gnomad_af_eur": 0.143,
    "hgvs": "NM_000041.4:c.388T>C",
    "sift": "tolerated",
    "polyphen": "benign"
  },
  {
    "chrom": "1",
    "pos": "11796321",
    "id": "rs1801133",
    "ref": "G",
    "alt": "A",
    "gene": "MTHFR",
    "consequence": "missense_variant",
    "impact": "MODERATE",
    "clinvar_significance": "Benign/Likely benign",
    "clinvar_condition": "Homocystinuria",
    "gnomad_af": 0.312,
    "gnomad_af_afr": 0.198,
    "gnomad_af_eur": 0.367,
    "hgvs": "NM_005957.5:c.665C>T",
    "sift": "tolerated",
    "polyphen": "benign"
  }
]\
""", 'tables/variants.csv': """\
chrom,pos,id,ref,alt,gene,consequence,impact,hgvs,sift,polyphen,clinvar_significance,clinvar_condition,gnomad_af,gnomad_af_afr,gnomad_af_eur
13,32316461,rs80359550,C,T,BRCA2,stop_gained,HIGH,NM_000059.4:c.9976A>T,deleterious,probably_damaging,Pathogenic,Familial cancer of breast,4e-06,0.0,6e-06
7,117548628,rs113993960,CTTT,C,CFTR,frameshift_variant,HIGH,NM_000492.4:c.1521_1523delCTT,deleterious,N/A,Pathogenic,Cystic fibrosis,0.021,0.0003,0.033
17,43063931,rs55770810,G,A,BRCA1,missense_variant,MODERATE,NM_007294.4:c.5095C>T (p.Arg1699Trp),deleterious,probably_damaging,Pathogenic,Hereditary breast and ovarian cancer syndrome,2.4e-05,,
19,44908684,rs429358,T,C,APOE,missense_variant,MODERATE,NM_000041.4:c.388T>C,tolerated,benign,Risk factor,Alzheimer disease,0.147,0.232,0.143
1,11796321,rs1801133,G,A,MTHFR,missense_variant,MODERATE,NM_005957.5:c.665C>T,tolerated,benign,Benign/Likely benign,Homocystinuria,0.312,0.198,0.367
"""})
