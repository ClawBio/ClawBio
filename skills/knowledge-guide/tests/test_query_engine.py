"""Tests for the knowledge-guide query engine."""
from __future__ import annotations
import json
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

# Minimal cache for testing
MINI_CACHE = {
    "topics": [
        {
            "name": "variant-analysis",
            "title": "Variant Analysis",
            "summary": "Genetic variant detection",
            "tutorials": [
                {
                    "title": "Calling very rare variants",
                    "name": "rare-variants",
                    "topic": "variant-analysis",
                    "time_estimation": "3H",
                    "level": "Intermediate",
                    "objectives": ["Process duplex sequencing data",
                                   "Identify rare variants without diploid assumptions"],
                    "tools": ["fastqc", "bwa_mem"],
                    "url": "https://training.galaxyproject.org/.../rare-variants/tutorial.html",
                },
                {
                    "title": "Exome sequencing variant calling",
                    "name": "exome-seq",
                    "topic": "variant-analysis",
                    "time_estimation": "4H",
                    "level": "Intermediate",
                    "objectives": ["Perform variant calling on exome data",
                                   "Annotate variants with functional information"],
                    "tools": ["bwa_mem", "freebayes", "snpeff"],
                    "url": "https://training.galaxyproject.org/.../exome-seq/tutorial.html",
                },
            ],
        },
        {
            "name": "transcriptomics",
            "title": "Transcriptomics",
            "summary": "Gene expression analysis",
            "tutorials": [
                {
                    "title": "Reference-based RNA-Seq data analysis",
                    "name": "ref-based",
                    "topic": "transcriptomics",
                    "time_estimation": "4H",
                    "level": "Intermediate",
                    "objectives": ["Perform RNA-seq alignment",
                                   "Identify differentially expressed genes"],
                    "tools": ["hisat2", "featurecounts", "deseq2"],
                    "url": "https://training.galaxyproject.org/.../ref-based/tutorial.html",
                },
            ],
        },
    ],
    "tool_index": {
        "devteam/fastqc/fastqc": [
            {"path": "variant-analysis/rare-variants",
             "title": "Calling very rare variants",
             "topic": "Variant Analysis",
             "url": "https://training.galaxyproject.org/.../rare-variants/tutorial.html"},
        ],
    },
}


class TestFreeTextQuery:
    """Test free-text query matching."""

    def test_matches_by_title(self):
        from query_engine import search_tutorials
        results = search_tutorials("variant calling", MINI_CACHE)
        assert len(results) > 0
        assert any("variant" in r["title"].lower() for r in results)

    def test_matches_by_objective(self):
        from query_engine import search_tutorials
        results = search_tutorials("differentially expressed genes", MINI_CACHE)
        assert len(results) > 0
        assert results[0]["topic"] == "transcriptomics"

    def test_returns_max_5(self):
        from query_engine import search_tutorials
        results = search_tutorials("analysis", MINI_CACHE)
        assert len(results) <= 5

    def test_empty_query_returns_empty(self):
        from query_engine import search_tutorials
        results = search_tutorials("", MINI_CACHE)
        assert results == []

    def test_no_matches_returns_empty(self):
        from query_engine import search_tutorials
        results = search_tutorials("quantum computing blockchain", MINI_CACHE)
        assert results == []


class TestTopicLookup:
    """Test direct topic lookup."""

    def test_finds_topic(self):
        from query_engine import lookup_topic
        result = lookup_topic("variant-analysis", MINI_CACHE)
        assert result is not None
        assert result["name"] == "variant-analysis"
        assert len(result["tutorials"]) == 2

    def test_unknown_topic_returns_none(self):
        from query_engine import lookup_topic
        result = lookup_topic("nonexistent-topic", MINI_CACHE)
        assert result is None


class TestToolLookup:
    """Test tool→tutorial reverse lookup."""

    def test_finds_tool(self):
        from query_engine import lookup_tool
        results = lookup_tool("fastqc", MINI_CACHE)
        assert len(results) > 0

    def test_unknown_tool_returns_empty(self):
        from query_engine import lookup_tool
        results = lookup_tool("nonexistent_tool_xyz", MINI_CACHE)
        assert results == []


class TestConceptLookup:
    """Test concept→tutorial fuzzy matching."""

    def test_matches_concept(self):
        from query_engine import lookup_concept

        recs = {
            "pharmgx-reporter": {
                "concepts": ["pharmacogenomics", "CYP2D6"],
                "tutorials": [
                    {"id": "variant-analysis/pharmgx", "title": "Pharmacogenomics"}
                ],
            }
        }
        results = lookup_concept("pharmacogenomics", MINI_CACHE, recs)
        assert len(results) > 0
