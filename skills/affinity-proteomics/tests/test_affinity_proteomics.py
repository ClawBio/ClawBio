"""Tests for affinity-proteomics skill.

Validates Olink NPX parsing, SomaLogic ADAT parsing, differential abundance
testing, and end-to-end demo mode for both platforms.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

from affinity_proteomics import (
    DiffAbundanceResult,
    ProteomicsData,
    differential_abundance,
    generate_report,
    parse_olink_npx,
    run_pipeline,
)

HAS_SOMADATA = True
try:
    import somadata
    SOMA_EXAMPLE = Path(somadata.__file__).parent / "data" / "example_data.adat"
    if not SOMA_EXAMPLE.exists():
        HAS_SOMADATA = False
except ImportError:
    HAS_SOMADATA = False


# ---------------------------------------------------------------------------
# Unit tests — Olink parser
# ---------------------------------------------------------------------------
class TestOlinkParser:
    @pytest.fixture
    def demo_npx(self):
        return SKILL_DIR / "example_data" / "olink_demo_npx.csv"

    @pytest.fixture
    def demo_meta(self):
        return SKILL_DIR / "example_data" / "olink_demo_meta.csv"

    def test_parse_npx_returns_proteomics_data(self, demo_npx, demo_meta):
        data = parse_olink_npx(demo_npx, demo_meta)
        assert isinstance(data, ProteomicsData)
        assert data.platform == "olink"

    def test_parse_npx_sample_count(self, demo_npx, demo_meta):
        data = parse_olink_npx(demo_npx, demo_meta)
        assert data.expression.shape[0] == 78  # 80 - 2 QC warned

    def test_parse_npx_protein_count(self, demo_npx, demo_meta):
        data = parse_olink_npx(demo_npx, demo_meta)
        assert data.expression.shape[1] == 40

    def test_qc_summary_present(self, demo_npx, demo_meta):
        data = parse_olink_npx(demo_npx, demo_meta)
        assert "total_samples" in data.qc_summary
        assert data.qc_summary["qc_warned_samples"] == 2

    def test_sample_info_has_group(self, demo_npx, demo_meta):
        data = parse_olink_npx(demo_npx, demo_meta)
        assert "Group" in data.sample_info.columns

    def test_protein_info_has_gene(self, demo_npx, demo_meta):
        data = parse_olink_npx(demo_npx, demo_meta)
        assert "gene_symbol" in data.protein_info.columns


# ---------------------------------------------------------------------------
# Unit tests — SomaLogic parser
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not HAS_SOMADATA, reason="somadata not installed")
class TestSomaLogicParser:
    def test_parse_adat_returns_proteomics_data(self):
        from affinity_proteomics import parse_somascan_adat
        data = parse_somascan_adat(SOMA_EXAMPLE, group_col="Sex")
        assert isinstance(data, ProteomicsData)
        assert data.platform == "somascan"

    def test_parse_adat_filters_non_samples(self):
        from affinity_proteomics import parse_somascan_adat
        data = parse_somascan_adat(SOMA_EXAMPLE, group_col="Sex")
        assert data.expression.shape[0] < 192  # some are calibrator/QC/buffer

    def test_parse_adat_log10_transform(self):
        from affinity_proteomics import parse_somascan_adat
        import numpy as np
        data = parse_somascan_adat(SOMA_EXAMPLE, group_col="Sex")
        # log10(RFU) values should typically be in range 1-5
        median_val = np.nanmedian(data.expression.values)
        assert 1.0 < median_val < 6.0

    def test_parse_adat_qc_summary(self):
        from affinity_proteomics import parse_somascan_adat
        data = parse_somascan_adat(SOMA_EXAMPLE, group_col="Sex")
        assert "transformation" in data.qc_summary
        assert data.qc_summary["transformation"] == "log10(RFU)"


# ---------------------------------------------------------------------------
# Unit tests — differential abundance
# ---------------------------------------------------------------------------
class TestDifferentialAbundance:
    @pytest.fixture
    def olink_data(self):
        return parse_olink_npx(
            SKILL_DIR / "example_data" / "olink_demo_npx.csv",
            SKILL_DIR / "example_data" / "olink_demo_meta.csv",
        )

    def test_returns_results(self, olink_data):
        results = differential_abundance(olink_data, "Group", ("Case", "Control"))
        assert len(results) > 0
        assert all(isinstance(r, DiffAbundanceResult) for r in results)

    def test_sorted_by_pvalue(self, olink_data):
        results = differential_abundance(olink_data, "Group", ("Case", "Control"))
        pvals = [r.pvalue for r in results]
        assert pvals == sorted(pvals)

    def test_detects_known_de_proteins(self, olink_data):
        results = differential_abundance(olink_data, "Group", ("Case", "Control"))
        sig_ids = {r.protein_id for r in results if r.significant}
        known_de = {"OID00001", "OID00003", "OID00013", "OID00033", "OID00019"}
        recovered = sig_ids & known_de
        assert len(recovered) >= 4, f"Only recovered {recovered} of {known_de}"

    def test_de_direction_is_positive(self, olink_data):
        results = differential_abundance(olink_data, "Group", ("Case", "Control"))
        known_de = {"OID00001", "OID00003", "OID00013", "OID00033", "OID00019"}
        for r in results:
            if r.protein_id in known_de and r.significant:
                assert r.log2fc > 0, f"{r.protein_id} should be up in Case"

    def test_invalid_group_col_raises(self, olink_data):
        with pytest.raises(ValueError, match="not found"):
            differential_abundance(olink_data, "NonExistent", ("A", "B"))


# ---------------------------------------------------------------------------
# Integration tests — Olink demo
# ---------------------------------------------------------------------------
class TestOlinkDemo:
    def test_olink_demo_end_to_end(self, tmp_path):
        summary = run_pipeline(
            platform="olink",
            input_path=SKILL_DIR / "example_data" / "olink_demo_npx.csv",
            meta_path=SKILL_DIR / "example_data" / "olink_demo_meta.csv",
            group_col="Group", contrast=("Case", "Control"),
            output_dir=tmp_path, demo=True,
        )
        assert summary["samples"] == 78
        assert summary["proteins"] == 40
        assert summary["significant"] >= 4

        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "result.json").exists()
        assert (tmp_path / "tables" / "diff_abundance.tsv").exists()
        assert (tmp_path / "figures" / "volcano.png").exists()
        assert (tmp_path / "figures" / "heatmap.png").exists()
        assert (tmp_path / "figures" / "pca.png").exists()

    def test_olink_report_content(self, tmp_path):
        run_pipeline(
            platform="olink", input_path=SKILL_DIR / "example_data" / "olink_demo_npx.csv",
            meta_path=SKILL_DIR / "example_data" / "olink_demo_meta.csv",
            group_col="Group", contrast=("Case", "Control"),
            output_dir=tmp_path, demo=True,
        )
        report = (tmp_path / "report.md").read_text()
        assert "OLINK" in report
        assert "ClawBio is a research" in report

        result = json.loads((tmp_path / "result.json").read_text())
        assert result["platform"] == "olink"
        assert result["total_proteins_tested"] == 40


# ---------------------------------------------------------------------------
# Integration tests — SomaLogic demo
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not HAS_SOMADATA, reason="somadata not installed")
class TestSomaLogicDemo:
    def test_somascan_demo_end_to_end(self, tmp_path):
        summary = run_pipeline(
            platform="somascan",
            input_path=SOMA_EXAMPLE,
            meta_path=None,
            group_col="Sex", contrast=("F", "M"),
            output_dir=tmp_path, demo=True,
        )
        assert summary["samples"] > 100
        assert summary["proteins"] > 4000
        assert summary["significant"] > 0

        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "figures" / "volcano.png").exists()

    def test_somascan_detects_sex_diff_proteins(self, tmp_path):
        summary = run_pipeline(
            platform="somascan", input_path=SOMA_EXAMPLE, meta_path=None,
            group_col="Sex", contrast=("F", "M"),
            output_dir=tmp_path, demo=True,
        )
        result = json.loads((tmp_path / "result.json").read_text())
        top_genes = [h["gene"] for h in result["top_10"]]
        assert "KLK3" in top_genes, f"Expected KLK3 (PSA) in top hits, got: {top_genes}"
