"""Tests for bio-orchestrator — skill routing hub."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from orchestrator import (
    EXTENSION_MAP,
    KEYWORD_MAP,
    SKILL_REGISTRY_MAP,
    SKILLS_DIR,
    detect_skill_from_file,
    detect_skill_from_tabular_header,
    detect_skill_from_query,
    detect_skill_with_hint_from_query,
    detect_multiple_skills,
    list_available_skills,
    generate_report_header,
    append_audit_log,
)


# ---------------------------------------------------------------------------
# detect_skill_from_file — extension-based routing
# ---------------------------------------------------------------------------

class TestDetectSkillFromFile:
    def test_vcf(self, tmp_path):
        f = tmp_path / "input.vcf"
        f.write_text("##fileformat=VCFv4.2\n")
        assert detect_skill_from_file(f) == "equity-scorer"

    def test_vcf_gz(self, tmp_path):
        f = tmp_path / "input.vcf.gz"
        f.write_bytes(b"")
        assert detect_skill_from_file(f) == "equity-scorer"

    def test_h5ad(self, tmp_path):
        f = tmp_path / "data.h5ad"
        f.write_bytes(b"")
        assert detect_skill_from_file(f) == "scrna-orchestrator"

    def test_fastq(self, tmp_path):
        f = tmp_path / "reads.fastq"
        f.write_text("@read1\nACGT\n+\nIIII\n")
        assert detect_skill_from_file(f) == "seq-wrangler"

    def test_pdb(self, tmp_path):
        f = tmp_path / "structure.pdb"
        f.write_text("ATOM\n")
        assert detect_skill_from_file(f) == "struct-predictor"

    def test_pkl(self, tmp_path):
        f = tmp_path / "data.pkl"
        f.write_bytes(b"")
        assert detect_skill_from_file(f) == "methylation-clock"

    def test_png(self, tmp_path):
        f = tmp_path / "chart.png"
        f.write_bytes(b"")
        assert detect_skill_from_file(f) == "data-extractor"

    def test_unknown_extension(self, tmp_path):
        f = tmp_path / "data.xyz"
        f.write_text("something")
        assert detect_skill_from_file(f) is None

    def test_directory_returns_none(self, tmp_path):
        d = tmp_path / "somedir"
        d.mkdir()
        assert detect_skill_from_file(d) is None


# ---------------------------------------------------------------------------
# detect_skill_from_tabular_header
# ---------------------------------------------------------------------------

class TestDetectSkillFromTabularHeader:
    def test_de_results_csv(self, tmp_path):
        f = tmp_path / "results.csv"
        f.write_text("gene,log2FoldChange,padj\nTP53,-1.5,0.001\n")
        assert detect_skill_from_tabular_header(f) == "diff-visualizer"

    def test_equity_csv(self, tmp_path):
        f = tmp_path / "pops.csv"
        f.write_text("sample,population,region\nS1,EUR,Europe\n")
        assert detect_skill_from_tabular_header(f) == "equity-scorer"

    def test_rnaseq_metadata(self, tmp_path):
        f = tmp_path / "meta.csv"
        f.write_text("sample_id,condition,batch\nS1,treated,1\n")
        assert detect_skill_from_tabular_header(f) == "rnaseq-de"

    def test_count_matrix(self, tmp_path):
        f = tmp_path / "counts.csv"
        f.write_text("gene,S1,S2,S3,S4\nTP53,100,200,150,300\n")
        assert detect_skill_from_tabular_header(f) == "rnaseq-de"

    def test_empty_file(self, tmp_path):
        f = tmp_path / "empty.csv"
        f.write_text("")
        assert detect_skill_from_tabular_header(f) is None

    def test_methylation_csv(self, tmp_path):
        f = tmp_path / "methyl.csv"
        cg_cols = ",".join([f"cg{i:08d}" for i in range(15)])
        f.write_text(f"sample_id,sex,{cg_cols}\nS1,female,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.1,0.2,0.3,0.4,0.5,0.6\n")
        assert detect_skill_from_tabular_header(f) == "methylation-clock"


# ---------------------------------------------------------------------------
# detect_skill_from_query — keyword-based routing
# ---------------------------------------------------------------------------

class TestDetectSkillFromQuery:
    @pytest.mark.parametrize("query,expected", [
        ("What is the equity score?", "equity-scorer"),
        ("Run polygenic risk score", "gwas-prs"),
        ("Look up rs12345", "gwas-lookup"),
        ("Run scanvi on my h5ad", "scrna-embedding"),
        ("Show me differential expression", "rnaseq-de"),
        ("Predict protein structure", "struct-predictor"),
        ("Search pubmed for papers", "lit-synthesizer"),
        ("methylation clock analysis", "methylation-clock"),
        ("clinpgx gene drug lookup", "clinpgx"),
    ])
    def test_keyword_hits(self, query, expected):
        assert detect_skill_from_query(query) == expected

    def test_no_match(self):
        assert detect_skill_from_query("hello world") is None


# ---------------------------------------------------------------------------
# detect_skill_with_hint_from_query — chain-aware scRNA routing
# ---------------------------------------------------------------------------

class TestDetectSkillWithHint:
    def test_latent_artifact_plus_downstream(self):
        skill, hint = detect_skill_with_hint_from_query(
            "Run clustering on integrated.h5ad after scvi markers"
        )
        assert skill == "scrna-orchestrator"
        assert "downstream" in hint.lower() or "latent" in hint.lower()

    def test_embedding_plus_downstream(self):
        skill, hint = detect_skill_with_hint_from_query(
            "Run scvi batch correction and then cluster annotation"
        )
        assert skill == "scrna-embedding"
        assert "two-step" in hint.lower() or "scrna-embedding" in hint.lower()

    def test_embedding_without_downstream_prefers_scrna_embedding(self):
        skill, hint = detect_skill_with_hint_from_query(
            "Build a latent embedding for this single-cell dataset"
        )
        assert skill == "scrna-embedding"
        assert "embedding-focused" in hint.lower()
        assert "latent space" in hint.lower()

    def test_embedding_with_latent_artifact_prefers_scrna_embedding(self):
        skill, hint = detect_skill_with_hint_from_query(
            "Run embedding on my integrated.h5ad file"
        )
        assert skill == "scrna-embedding"
        assert "embedding-focused" in hint.lower()
        assert "latent artifact" in hint.lower()

    def test_generic_single_cell_query_still_prefers_orchestrator(self):
        skill, hint = detect_skill_with_hint_from_query("single-cell dataset")
        assert skill == "scrna-orchestrator"
        assert hint == ""

    def test_downstream_only_query_still_prefers_orchestrator(self):
        skill, hint = detect_skill_with_hint_from_query("Cluster my single-cell data")
        assert skill == "scrna-orchestrator"
        assert hint == ""


# ---------------------------------------------------------------------------
# detect_multiple_skills
# ---------------------------------------------------------------------------

class TestDetectMultipleSkills:
    def test_single_match(self):
        skills = detect_multiple_skills("run polygenic risk score")
        assert "gwas-prs" in skills

    def test_multiple_matches(self):
        skills = detect_multiple_skills("run equity score and differential expression deseq2")
        assert "equity-scorer" in skills
        assert "rnaseq-de" in skills

    def test_no_matches(self):
        skills = detect_multiple_skills("hello world")
        assert skills == []


# ---------------------------------------------------------------------------
# list_available_skills
# ---------------------------------------------------------------------------

class TestListAvailableSkills:
    def test_returns_non_empty(self):
        skills = list_available_skills()
        assert len(skills) > 0

    def test_all_have_skill_md(self):
        for skill in list_available_skills():
            assert (SKILLS_DIR / skill / "SKILL.md").exists()


# ---------------------------------------------------------------------------
# SKILL_REGISTRY_MAP validation
# ---------------------------------------------------------------------------

class TestSkillRegistryMap:
    def test_all_keys_are_valid_skills(self):
        available = list_available_skills()
        for key in SKILL_REGISTRY_MAP:
            assert key in available, f"Registry key '{key}' is not a valid skill directory"


# ---------------------------------------------------------------------------
# generate_report_header
# ---------------------------------------------------------------------------

class TestGenerateReportHeader:
    def test_contains_title(self, tmp_path):
        f = tmp_path / "input.txt"
        f.write_text("data")
        header = generate_report_header("Test Analysis", ["pharmgx"], [f])
        assert "Test Analysis" in header
        assert "pharmgx" in header

    def test_contains_checksums(self, tmp_path):
        f = tmp_path / "input.txt"
        f.write_text("data")
        header = generate_report_header("T", ["s"], [f])
        assert "input.txt" in header


# ---------------------------------------------------------------------------
# append_audit_log
# ---------------------------------------------------------------------------

class TestAppendAuditLog:
    def test_creates_log(self, tmp_path):
        append_audit_log(tmp_path, "Routed to pharmgx", "input=test.txt")
        log = tmp_path / "analysis_log.md"
        assert log.exists()
        text = log.read_text()
        assert "Routed to pharmgx" in text

    def test_appends_entries(self, tmp_path):
        append_audit_log(tmp_path, "action1")
        append_audit_log(tmp_path, "action2")
        text = (tmp_path / "analysis_log.md").read_text()
        assert "action1" in text
        assert "action2" in text


# ── Stub detection and keyword disambiguation ─────────────────────────────────

class TestStubDetection:
    def test_stub_skill_detected(self):
        """Stub skills (SKILL.md only, no .py) should be flagged."""
        from orchestrator import skill_has_executable
        # repro-enforcer and lit-synthesizer are stubs (no Python executable)
        assert not skill_has_executable("repro-enforcer")
        assert not skill_has_executable("lit-synthesizer")

    def test_executable_skill_detected(self):
        """Skills with Python files should not be flagged as stubs."""
        from orchestrator import skill_has_executable
        assert skill_has_executable("equity-scorer")
        assert skill_has_executable("pharmgx-reporter")


class TestKeywordDisambiguation:
    def test_variant_routes_to_vcf_annotator(self):
        """kw_14: 'variant' alone should route to vcf-annotator."""
        from orchestrator import detect_skill_from_query
        assert detect_skill_from_query("annotate my variants") == "vcf-annotator"

    def test_variant_annotation_routes_to_vcf_annotator(self):
        """Longer phrase 'variant annotation' should beat single keyword matches."""
        from orchestrator import detect_skill_from_query
        assert detect_skill_from_query("variant annotation for diversity panel") == "vcf-annotator"

    def test_diversity_alone_routes_to_equity(self):
        """'diversity' without 'variant' should still route to equity-scorer."""
        from orchestrator import detect_skill_from_query
        assert detect_skill_from_query("compute diversity metrics") == "equity-scorer"
