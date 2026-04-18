"""
Tests for ClawBio VCF Annotator skill.
Run with: pytest skills/vcf-annotator/tests/test_vcf_annotator.py -v
"""

import csv
import json
import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from vcf_annotator import (
    DEMO_ANNOTATIONS,
    IMPACT_RANK,
    generate_report,
    parse_vcf,
)

# ── Demo VCF fixture ───────────────────────────────────────────────────────────

SAMPLE_VCF = """##fileformat=VCFv4.2
##reference=GRCh38
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
17\t43044295\trs80357382\tG\tA\t.\tPASS\t.
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
        assert variants[0]["id"] == "rs80357382"

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
