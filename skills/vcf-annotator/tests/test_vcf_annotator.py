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

sys.path.insert(0, str(Path(__file__).parent.parent))
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
        empty.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n")
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
