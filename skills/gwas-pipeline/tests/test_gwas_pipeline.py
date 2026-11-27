"""Tests for gwas-pipeline skill.

Validates REGENIE output parsing, lambda GC computation, lead variant
extraction, and end-to-end demo mode against the REGENIE official example
dataset (500 samples, 1000 variants, binary trait Y1).

Integration tests require plink2 and regenie on PATH and are skipped
gracefully if the binaries are not available.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

from gwas_pipeline import (
    GWASResult,
    QCParams,
    _find_binary,
    compute_lambda_gc,
    extract_lead_variants,
    parse_regenie_output,
    run_pipeline,
)

HAS_PLINK2 = _find_binary("plink2") is not None
HAS_REGENIE = _find_binary("regenie") is not None
SKIP_INTEGRATION = not (HAS_PLINK2 and HAS_REGENIE)
SKIP_REASON = "plink2 and/or regenie not on PATH"


# ---------------------------------------------------------------------------
# Unit tests — result parsing
# ---------------------------------------------------------------------------
class TestRegenieParser:
    def test_parse_reference_output(self):
        ref_file = SKILL_DIR / "example_data" / "test_bin_out_firth_Y1.regenie"
        results = parse_regenie_output(ref_file)
        assert len(results) == 1000
        assert all(isinstance(r, GWASResult) for r in results)

    def test_parsed_fields_are_numeric(self):
        ref_file = SKILL_DIR / "example_data" / "test_bin_out_firth_Y1.regenie"
        results = parse_regenie_output(ref_file)
        for r in results[:10]:
            assert isinstance(r.pos, int)
            assert isinstance(r.beta, float)
            assert isinstance(r.se, float)
            assert isinstance(r.log10p, float)
            assert r.log10p > 0

    def test_pvalue_conversion(self):
        r = GWASResult(
            chrom="1", pos=100, snp_id="rs1", allele0="A", allele1="G",
            a1freq=0.3, n=500, beta=0.5, se=0.1, chisq=25.0, log10p=5.0,
        )
        assert abs(r.pvalue - 1e-5) < 1e-10

    def test_pvalue_from_zero_log10p(self):
        r = GWASResult(
            chrom="1", pos=100, snp_id="rs1", allele0="A", allele1="G",
            a1freq=0.3, n=500, beta=0.0, se=0.1, chisq=0.0, log10p=0.0,
        )
        assert r.pvalue == 1.0


# ---------------------------------------------------------------------------
# Unit tests — lambda GC
# ---------------------------------------------------------------------------
class TestLambdaGC:
    def test_lambda_gc_from_reference(self):
        ref_file = SKILL_DIR / "example_data" / "test_bin_out_firth_Y1.regenie"
        results = parse_regenie_output(ref_file)
        lam = compute_lambda_gc(results)
        assert 0.8 < lam < 1.5

    def test_lambda_gc_empty_results(self):
        assert compute_lambda_gc([]) == 1.0

    def test_lambda_gc_uniform_chisq(self):
        results = [
            GWASResult("1", i, f"rs{i}", "A", "G", 0.3, 500, 0.0, 0.1, 0.455, 0.3)
            for i in range(100)
        ]
        lam = compute_lambda_gc(results)
        assert 0.9 < lam < 1.1


# ---------------------------------------------------------------------------
# Unit tests — lead variant extraction
# ---------------------------------------------------------------------------
class TestLeadVariants:
    def test_no_gws_hits_in_reference(self):
        ref_file = SKILL_DIR / "example_data" / "test_bin_out_firth_Y1.regenie"
        results = parse_regenie_output(ref_file)
        leads = extract_lead_variants(results)
        assert len(leads) == 0

    def test_extract_with_low_threshold(self):
        ref_file = SKILL_DIR / "example_data" / "test_bin_out_firth_Y1.regenie"
        results = parse_regenie_output(ref_file)
        leads = extract_lead_variants(results, threshold=1.0)
        assert len(leads) > 0

    def test_leads_sorted_by_significance(self):
        results = [
            GWASResult("1", 100, "rs1", "A", "G", 0.3, 500, 0.5, 0.1, 25.0, 9.0),
            GWASResult("1", 200, "rs2", "A", "G", 0.3, 500, 0.3, 0.1, 9.0, 10.0),
        ]
        leads = extract_lead_variants(results, threshold=1.0)
        assert leads[0].snp_id == "rs2"


# ---------------------------------------------------------------------------
# Unit tests — QC params
# ---------------------------------------------------------------------------
class TestQCParams:
    def test_default_thresholds(self):
        p = QCParams()
        assert p.geno == 0.02
        assert p.mind == 0.02
        assert p.maf == 0.01
        assert p.hwe == 1e-6
        assert p.ld_r2 == 0.9

    def test_custom_thresholds(self):
        p = QCParams(geno=0.05, maf=0.05, hwe=1e-4)
        assert p.geno == 0.05
        assert p.maf == 0.05


# ---------------------------------------------------------------------------
# Integration tests — full demo pipeline
# ---------------------------------------------------------------------------
@pytest.mark.skipif(SKIP_INTEGRATION, reason=SKIP_REASON)
class TestDemoPipeline:
    def test_demo_end_to_end(self, tmp_path):
        data_dir = SKILL_DIR / "example_data"
        summary = run_pipeline(
            bed_prefix=data_dir / "example",
            bgen_file=data_dir / "example.bgen",
            pheno_file=data_dir / "phenotype_bin.txt",
            covar_file=data_dir / "covariates.txt",
            trait_type="bt",
            trait_name="Y1",
            output_dir=tmp_path,
            qc_params=QCParams(),
            demo=True,
        )

        assert summary["total_variants"] == 1000
        assert 0.8 < summary["lambda_gc"] < 1.5
        assert summary["gws_hits"] == 0

        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "result.json").exists()
        assert (tmp_path / "tables" / "gwas_results.tsv").exists()
        assert (tmp_path / "figures" / "manhattan.png").exists()
        assert (tmp_path / "figures" / "qq_plot.png").exists()
        assert (tmp_path / "reproducibility" / "commands.sh").exists()

    def test_demo_report_content(self, tmp_path):
        data_dir = SKILL_DIR / "example_data"
        run_pipeline(
            bed_prefix=data_dir / "example",
            bgen_file=data_dir / "example.bgen",
            pheno_file=data_dir / "phenotype_bin.txt",
            covar_file=data_dir / "covariates.txt",
            trait_type="bt",
            trait_name="Y1",
            output_dir=tmp_path,
            qc_params=QCParams(),
            demo=True,
        )

        report = (tmp_path / "report.md").read_text()
        assert "GWAS Report" in report
        assert "Lambda GC" in report
        assert "ClawBio is a research" in report

        result = json.loads((tmp_path / "result.json").read_text())
        assert result["total_variants_tested"] == 1000
        assert result["method"] == "PLINK2 QC + REGENIE two-step regression"

    def test_demo_results_match_reference(self, tmp_path):
        """Verify our results are comparable to the REGENIE reference output."""
        data_dir = SKILL_DIR / "example_data"
        run_pipeline(
            bed_prefix=data_dir / "example",
            bgen_file=data_dir / "example.bgen",
            pheno_file=data_dir / "phenotype_bin.txt",
            covar_file=data_dir / "covariates.txt",
            trait_type="bt",
            trait_name="Y1",
            output_dir=tmp_path,
            qc_params=QCParams(),
            demo=True,
        )

        our_results = parse_regenie_output(tmp_path / "step2" / "assoc_out_Y1.regenie")
        ref_results = parse_regenie_output(data_dir / "test_bin_out_firth_Y1.regenie")

        assert len(our_results) == len(ref_results) == 1000

        our_top = sorted(our_results, key=lambda r: r.log10p, reverse=True)[:5]
        ref_top = sorted(ref_results, key=lambda r: r.log10p, reverse=True)[:5]
        our_top_ids = {r.snp_id for r in our_top}
        ref_top_ids = {r.snp_id for r in ref_top}
        overlap = our_top_ids & ref_top_ids
        assert len(overlap) >= 3, f"Top-5 overlap: {overlap}"
