"""
test_tool_recommender.py — Tests for intelligent Galaxy tool selection
======================================================================
Covers: version dedup, task matching, tool recommendation, EDAM resolution,
workflow suggestions, format detection.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL_DIR))

import tool_recommender as rec  # noqa: E402


# ---------------------------------------------------------------------------
# Sample catalog for testing
# ---------------------------------------------------------------------------

SAMPLE_TOOLS = [
    {
        "id": "toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.63",
        "name": "FastQC",
        "description": "Read Quality reports",
        "version": "0.63",
        "section": "FASTQ Quality Control",
        "edam_topics": ["topic_3168"],
        "edam_operations": ["operation_3218"],
        "inputs": [],
        "outputs": [],
    },
    {
        "id": "toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.74+galaxy0",
        "name": "FastQC",
        "description": "Read Quality reports",
        "version": "0.74+galaxy0",
        "section": "FASTQ Quality Control",
        "edam_topics": ["topic_3168"],
        "edam_operations": ["operation_3218"],
        "inputs": [],
        "outputs": [],
    },
    {
        "id": "toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.74+galaxy1",
        "name": "FastQC",
        "description": "Read Quality reports",
        "version": "0.74+galaxy1",
        "section": "FASTQ Quality Control",
        "edam_topics": ["topic_3168"],
        "edam_operations": ["operation_3218"],
        "inputs": [],
        "outputs": [],
    },
    {
        "id": "toolshed.g2.bx.psu.edu/repos/iuc/kraken2/kraken2/2.1.3",
        "name": "Kraken2",
        "description": "Assign taxonomic labels to short DNA sequences",
        "version": "2.1.3",
        "section": "Metagenomic Analysis",
        "edam_topics": ["topic_3697"],
        "edam_operations": ["operation_3460"],
        "inputs": [],
        "outputs": [],
    },
    {
        "id": "toolshed.g2.bx.psu.edu/repos/iuc/freebayes/freebayes/1.3.6",
        "name": "FreeBayes",
        "description": "Bayesian haplotype-based polymorphism discovery",
        "version": "1.3.6",
        "section": "Variant Calling",
        "edam_topics": ["topic_0199"],
        "edam_operations": ["operation_3227"],
        "inputs": [],
        "outputs": [],
    },
    {
        "id": "toolshed.g2.bx.psu.edu/repos/iuc/snpsift/snpSift_filter/4.3",
        "name": "SnpSift Filter",
        "description": "Filter variants using arbitrary expressions",
        "version": "4.3",
        "section": "VCF/BCF",
        "edam_topics": ["topic_0199"],
        "edam_operations": ["operation_3695"],
        "inputs": [],
        "outputs": [],
    },
    {
        "id": "toolshed.g2.bx.psu.edu/repos/iuc/deseq2/deseq2/1.42.0",
        "name": "DESeq2",
        "description": "Differential gene expression analysis",
        "version": "1.42.0",
        "section": "RNA-seq",
        "edam_topics": ["topic_3170"],
        "edam_operations": ["operation_3223"],
        "inputs": [],
        "outputs": [],
    },
    {
        "id": "toolshed.g2.bx.psu.edu/repos/devteam/samtools_sort/samtools_sort/2.0.5",
        "name": "Samtools sort",
        "description": "Sort alignments by leftmost coordinates or read name",
        "version": "2.0.5",
        "section": "SAM/BAM",
        "edam_topics": [],
        "edam_operations": [],
        "inputs": [],
        "outputs": [],
    },
    {
        "id": "toolshed.g2.bx.psu.edu/repos/iuc/snpeff/snpEff/5.2",
        "name": "SnpEff",
        "description": "Variant annotation and effect prediction",
        "version": "5.2",
        "section": "VCF/BCF",
        "edam_topics": ["topic_0199"],
        "edam_operations": ["operation_3672"],
        "inputs": [],
        "outputs": [],
    },
    {
        "id": "toolshed.g2.bx.psu.edu/repos/iuc/gatk4_haplotypecaller/gatk4_haplotypecaller/4.5",
        "name": "GATK4 HaplotypeCaller",
        "description": "Call germline SNPs and indels via local re-assembly",
        "version": "4.5",
        "section": "Variant Calling",
        "edam_topics": ["topic_0199"],
        "edam_operations": ["operation_3227"],
        "inputs": [],
        "outputs": [],
    },
]

SAMPLE_CATALOG = {"tools": SAMPLE_TOOLS}


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestVersionDedup:
    def test_dedup_reduces_count(self):
        """Multiple FastQC versions collapse to one entry."""
        deduped = rec.build_deduped_index(SAMPLE_TOOLS)
        fastqc_entries = [t for t in deduped.values() if t["name"] == "FastQC"]
        assert len(fastqc_entries) == 1

    def test_dedup_picks_latest_version(self):
        """Dedup picks FastQC 0.74+galaxy1, not 0.63."""
        deduped = rec.build_deduped_index(SAMPLE_TOOLS)
        fastqc = [t for t in deduped.values() if t["name"] == "FastQC"][0]
        assert fastqc["version"] == "0.74+galaxy1"

    def test_dedup_preserves_version_count(self):
        """Deduped entry has version_count = 3 for FastQC."""
        deduped = rec.build_deduped_index(SAMPLE_TOOLS)
        fastqc = [t for t in deduped.values() if t["name"] == "FastQC"][0]
        assert fastqc["version_count"] == 3

    def test_version_sort_key(self):
        """Version sorting: 0.74+galaxy1 > 0.74+galaxy0 > 0.63."""
        k1 = rec._version_sort_key("0.74+galaxy1")
        k0 = rec._version_sort_key("0.74+galaxy0")
        k63 = rec._version_sort_key("0.63")
        assert k1 > k0 > k63

    def test_tool_base_key(self):
        """Base key strips version from ToolShed IDs."""
        tid = "toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.74+galaxy1"
        expected = "toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc"
        assert rec._tool_base_key(tid) == expected

    def test_non_toolshed_id_unchanged(self):
        """Non-ToolShed IDs returned as-is."""
        assert rec._tool_base_key("Summary_Statistics1") == "Summary_Statistics1"


class TestTaskMatching:
    def test_qc_matches(self):
        """'quality control on reads' matches quality_control task."""
        matches = rec._match_tasks("I need to check the quality of my reads")
        task_ids = [m[0] for m in matches]
        assert "quality_control" in task_ids

    def test_metagenomics_matches(self):
        """'classify microbial species' matches metagenomics."""
        matches = rec._match_tasks("classify microbial species in my sample")
        task_ids = [m[0] for m in matches]
        assert "metagenomics" in task_ids

    def test_variant_calling_matches(self):
        """'call variants from BAM' matches variant_calling."""
        matches = rec._match_tasks("call variants from my BAM file")
        task_ids = [m[0] for m in matches]
        assert "variant_calling" in task_ids

    def test_wes_matches(self):
        """'whole exome sequencing' matches wes_wgs_analysis."""
        matches = rec._match_tasks("analyse whole exome sequencing data from Novogene")
        task_ids = [m[0] for m in matches]
        assert "wes_wgs_analysis" in task_ids

    def test_no_match_returns_empty(self):
        """Unrelated query returns no matches."""
        matches = rec._match_tasks("bake a chocolate cake")
        assert matches == []


class TestRecommendTool:
    def test_qc_recommends_fastqc(self):
        """QC task returns FastQC as top recommendation."""
        results = rec.recommend_tool("quality control on reads", SAMPLE_CATALOG)
        assert len(results) > 0
        assert results[0]["name"] == "FastQC"

    def test_metagenomics_recommends_kraken2(self):
        """Metagenomics task returns Kraken2."""
        results = rec.recommend_tool("classify microbial species", SAMPLE_CATALOG)
        names = [r["name"] for r in results]
        assert "Kraken2" in names

    def test_variant_calling_recommends_freebayes(self):
        """Variant calling task returns FreeBayes or GATK."""
        results = rec.recommend_tool("call germline variants", SAMPLE_CATALOG)
        names = [r["name"] for r in results]
        assert "FreeBayes" in names or "GATK4 HaplotypeCaller" in names

    def test_format_filter_bam(self):
        """Providing .bam format boosts variant calling tools."""
        results = rec.recommend_tool("analyse my data", SAMPLE_CATALOG, input_format=".bam")
        names = [r["name"] for r in results[:3]]
        # Should include BAM-compatible tools
        assert any(n in names for n in ["FreeBayes", "GATK4 HaplotypeCaller", "Samtools sort"])

    def test_returns_explanations(self):
        """Each recommendation has a non-empty explanation."""
        results = rec.recommend_tool("quality control", SAMPLE_CATALOG)
        for r in results:
            assert r["explanation"]

    def test_max_results(self):
        """max_results parameter limits output."""
        results = rec.recommend_tool("quality control", SAMPLE_CATALOG, max_results=2)
        assert len(results) <= 2

    def test_wes_novogene_context(self):
        """WES/Novogene context returns variant annotation tools."""
        results = rec.recommend_tool(
            "annotate variants from whole exome sequencing Novogene GATK output",
            SAMPLE_CATALOG,
            input_format=".vcf",
        )
        names = [r["name"] for r in results]
        # Should find SnpEff or SnpSift for VCF annotation
        assert any("Snp" in n for n in names)

    def test_dedup_cache_cleared(self):
        """Clear dedup cache between tests."""
        rec._deduped_cache = None


class TestEdamResolution:
    def test_topic_resolution(self):
        """topic_3168 resolves to 'Sequencing'."""
        assert rec.EDAM_TOPICS["topic_3168"] == "Sequencing"

    def test_operation_resolution(self):
        """operation_3227 resolves to 'Variant calling'."""
        assert rec.EDAM_OPERATIONS["operation_3227"] == "Variant calling"

    def test_resolve_edam_on_tool(self):
        """Tool with EDAM IDs gets resolved labels."""
        tool = {"edam_topics": ["topic_3697"], "edam_operations": ["operation_3460"]}
        labels = rec._resolve_edam(tool)
        assert "Metagenomics" in labels
        assert "Taxonomic classification" in labels


class TestWorkflowSuggestion:
    def test_rnaseq_workflow(self):
        """'RNA-seq' matches rnaseq_de workflow."""
        workflows = rec.suggest_workflow("RNA-seq differential expression")
        names = [w["name"] for w in workflows]
        assert "RNA-seq Differential Expression" in names

    def test_metagenomics_workflow(self):
        """'microbiome' matches metagenomics workflow."""
        workflows = rec.suggest_workflow("microbiome analysis")
        names = [w["name"] for w in workflows]
        assert "Metagenomics Community Profiling" in names

    def test_wes_workflow(self):
        """'exome' matches WES workflows."""
        workflows = rec.suggest_workflow("whole exome sequencing analysis")
        assert len(workflows) > 0
        # Should include WES germline or annotation workflow
        names = [w["name"] for w in workflows]
        assert any("WES" in n or "Exome" in n for n in names)

    def test_workflow_has_steps(self):
        """Each workflow has at least 3 steps."""
        workflows = rec.suggest_workflow("RNA-seq")
        for wf in workflows:
            assert len(wf["steps"]) >= 3
            for step in wf["steps"]:
                assert "tool" in step
                assert "purpose" in step

    def test_no_match_returns_empty(self):
        """Unrelated query returns no workflows."""
        workflows = rec.suggest_workflow("bake a cake")
        assert workflows == []


class TestFormatDetection:
    def test_double_extension_fastq_gz(self):
        """'reads.fastq.gz' detects as .fastq.gz."""
        assert rec.detect_format("reads.fastq.gz") == ".fastq.gz"

    def test_single_extension_bam(self):
        """'alignment.bam' detects as .bam."""
        assert rec.detect_format("alignment.bam") == ".bam"

    def test_vcf_gz(self):
        """'variants.vcf.gz' detects as .vcf.gz."""
        assert rec.detect_format("variants.vcf.gz") == ".vcf.gz"

    def test_simple_vcf(self):
        """'output.vcf' detects as .vcf."""
        assert rec.detect_format("output.vcf") == ".vcf"


class TestToolSlug:
    def test_toolshed_slug(self):
        """Extracts 'fastqc' from full ToolShed ID."""
        tid = "toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.74"
        assert rec._tool_slug(tid) == "fastqc"

    def test_simple_slug(self):
        """Simple IDs return last component."""
        assert rec._tool_slug("upload1") == "upload1"
