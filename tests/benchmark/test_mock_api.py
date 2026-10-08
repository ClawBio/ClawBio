"""Tests for mock API server and benchmark scorer."""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import pytest

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))

from mock_api_server import MockAPIServer, _ensembl_variation, _ensembl_vep, _gwas_catalog_associations, _clinpgx_gene
from benchmark_scorer import BenchmarkScorer

GROUND_TRUTH = Path(__file__).parent / "ad_ground_truth.json"


# ── Ground Truth Integrity ──────────────────────────────────────────────────


class TestGroundTruth:
    """Validate the AD ground truth file structure and content."""

    @pytest.fixture(autouse=True)
    def load_gt(self):
        with open(GROUND_TRUTH) as f:
            self.gt = json.load(f)

    def test_ground_truth_exists(self):
        assert GROUND_TRUTH.exists()

    def test_has_required_sections(self):
        for key in ["positive_genes", "negative_genes", "lead_variants", "scoring"]:
            assert key in self.gt, f"Missing section: {key}"

    def test_tier1_has_mendelian_genes(self):
        tier1 = {g["gene"] for g in self.gt["positive_genes"]["tier1_causal"]}
        assert {"APP", "PSEN1", "PSEN2", "APOE"} == tier1

    def test_tier2_has_gwas_genes(self):
        tier2 = {g["gene"] for g in self.gt["positive_genes"]["tier2_gwas_replicated"]}
        assert len(tier2) >= 15, f"Expected >=15 tier2 genes, got {len(tier2)}"
        assert "BIN1" in tier2
        assert "TREM2" in tier2

    def test_tier3_has_novel_genes(self):
        tier3 = {g["gene"] for g in self.gt["positive_genes"]["tier3_novel_bellenguez"]}
        assert len(tier3) >= 8

    def test_negative_genes_dont_overlap_positive(self):
        pos = set()
        for tier in ["tier1_causal", "tier2_gwas_replicated", "tier3_novel_bellenguez"]:
            pos |= {g["gene"] for g in self.gt["positive_genes"][tier]}
        neg = {g["gene"] for g in self.gt["negative_genes"]["genes"]}
        overlap = pos & neg
        assert len(overlap) == 0, f"Overlap between positive and negative: {overlap}"

    def test_lead_variants_have_rsids(self):
        variants = self.gt["lead_variants"]["variants"]
        assert len(variants) >= 10
        for v in variants:
            assert v["rsid"].startswith("rs")
            assert "chr" in v
            assert "gene" in v

    def test_scoring_has_minimums(self):
        minimums = self.gt["scoring"]["minimum_acceptable"]
        assert "gene_recovery_rate" in minimums
        assert "precision" in minimums
        assert "f1" in minimums


# ── Mock API Responses ──────────────────────────────────────────────────────


class TestMockAPIResponses:
    """Test that mock API functions return valid deterministic data."""

    def test_ensembl_variation_known_rsid(self):
        resp = _ensembl_variation("rs6733839")
        assert resp["name"] == "rs6733839"
        assert resp["var_class"] == "SNP"
        assert len(resp["mappings"]) == 2
        grch38 = [m for m in resp["mappings"] if m["assembly_name"] == "GRCh38"]
        assert len(grch38) == 1
        assert grch38[0]["seq_region_name"] == "2"  # BIN1 is on chr2

    def test_ensembl_variation_unknown_rsid(self):
        resp = _ensembl_variation("rs999999999")
        assert resp["name"] == "rs999999999"
        assert resp["var_class"] == "SNP"

    def test_ensembl_vep_returns_gene(self):
        resp = _ensembl_vep("rs6733839")
        assert isinstance(resp, list)
        assert len(resp) == 1
        consequences = resp[0]["transcript_consequences"]
        assert consequences[0]["gene_symbol"] == "BIN1"

    def test_gwas_catalog_returns_associations(self):
        resp = _gwas_catalog_associations("rs6733839")
        assocs = resp["_embedded"]["associations"]
        assert len(assocs) == 1
        # GWAS Catalog REST API v2 shape
        assert assocs[0]["efo_traits"][0]["efo_trait"] == "Alzheimer's disease"
        assert assocs[0]["p_value"] == pytest.approx(2.1e-173, rel=0.01)
        assert assocs[0]["snp_effect_allele"] == ["rs6733839-A"]
        assert resp["page"]["totalElements"] == 1

    def test_clinpgx_gene_response(self):
        resp = _clinpgx_gene("CYP2D6")
        assert resp["gene"] == "CYP2D6"
        assert len(resp["guidelines"]) >= 1
        assert resp["guidelines"][0]["source"] == "CPIC"


# ── Mock API Server ─────────────────────────────────────────────────────────


class TestMockAPIServer:
    """Test the HTTP server end-to-end."""

    @pytest.fixture(autouse=True)
    def server(self):
        with MockAPIServer(port=18089) as srv:
            self.base = srv.base_url
            yield

    def _get(self, path: str) -> dict:
        url = f"{self.base}{path}"
        with urllib.request.urlopen(url) as resp:
            return json.loads(resp.read())

    def test_health_endpoint(self):
        resp = self._get("/health")
        assert resp["status"] == "ok"

    def test_ensembl_variation_endpoint(self):
        resp = self._get("/ensembl/variation/human/rs6733839")
        assert resp["name"] == "rs6733839"

    def test_ensembl_vep_endpoint(self):
        resp = self._get("/ensembl/vep/human/id/rs6733839")
        assert isinstance(resp, list)
        assert resp[0]["transcript_consequences"][0]["gene_symbol"] == "BIN1"

    def test_gwas_catalog_endpoint(self):
        resp = self._get("/gwas/rest/api/v2/associations?rs_id=rs6733839")
        assert resp["_embedded"]["associations"][0]["efo_traits"][0]["efo_trait"] == "Alzheimer's disease"

    def test_gwas_catalog_legacy_v1_path_is_gone(self):
        """The real v1 API answers 410; the mock must not keep it alive."""
        try:
            self._get("/gwas/rest/api/singleNucleotidePolymorphisms/rs6733839/associations")
            assert False, "Should have raised"
        except urllib.error.HTTPError as e:
            assert e.code == 410

    def test_clinpgx_gene_endpoint(self):
        resp = self._get("/clinpgx/v1/gene/CYP2D6")
        assert resp["gene"] == "CYP2D6"

    def test_404_for_unknown_path(self):
        try:
            self._get("/nonexistent/path")
            assert False, "Should have raised"
        except urllib.error.HTTPError as e:
            assert e.code == 404


# ── Benchmark Scorer ────────────────────────────────────────────────────────


class TestBenchmarkScorer:
    """Test the benchmark scoring logic."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.scorer = BenchmarkScorer()

    def test_perfect_score(self):
        """All positive genes recovered, no false positives."""
        all_pos = list(self.scorer.all_positive)
        result = self.scorer.score(all_pos)
        assert result["gene_recovery_rate"] == 1.0
        assert result["precision"] == 1.0
        assert result["f1"] == 1.0
        assert result["false_positives"] == 0
        assert result["passes_minimum"] is True

    def test_zero_score(self):
        """No genes at all."""
        result = self.scorer.score([])
        assert result["gene_recovery_rate"] == 0.0
        assert result["f1"] == 0.0
        assert result["passes_minimum"] is False

    def test_only_false_positives(self):
        """Only negative genes submitted."""
        neg = list(self.scorer.negative)
        result = self.scorer.score(neg)
        assert result["true_positives"] == 0
        assert result["false_positives"] == len(neg)
        assert result["precision"] == 0.0
        assert result["passes_minimum"] is False

    def test_mixed_results(self):
        """Some true positives, some false positives."""
        genes = ["APP", "PSEN1", "BIN1", "CLU", "GAPDH", "ACTB"]
        result = self.scorer.score(genes)
        assert result["true_positives"] == 4  # APP, PSEN1, BIN1, CLU
        assert result["false_positives"] == 2  # GAPDH, ACTB
        assert result["precision"] == pytest.approx(4 / 6, rel=0.01)

    def test_tier_weights(self):
        """Tier1 genes contribute more to weighted score than tier3."""
        tier1_only = list(self.scorer.tier1)
        tier3_only = list(self.scorer.tier3)[:4]

        r1 = self.scorer.score(tier1_only)
        r3 = self.scorer.score(tier3_only)
        assert r1["weighted_score"] > r3["weighted_score"]

    def test_unknown_genes_not_penalised(self):
        """Genes not in either set are tracked but don't affect precision."""
        genes = ["APP", "PSEN1", "MYGENE123"]
        result = self.scorer.score(genes)
        assert result["unknown_genes"] == 1
        assert result["false_positives"] == 0
        assert result["precision"] == 1.0

    def test_variant_scoring(self):
        """Test variant-level scoring."""
        variants = [
            {"rsid": "rs6733839"},
            {"rsid": "rs679515"},
            {"rsid": "rs999999"},
        ]
        result = self.scorer.score_variants(variants)
        assert result["lead_variants_recovered"] == 2
        assert "BIN1" in result["recovered_genes"]
        assert "CR1" in result["recovered_genes"]

    def test_markdown_output(self):
        """Summary markdown is generated without errors."""
        result = self.scorer.score(["APP", "BIN1", "GAPDH"])
        md = self.scorer.summary_markdown(result)
        assert "AD Benchmark Results" in md
        assert "Gene recovery rate" in md
        assert "Precision" in md
