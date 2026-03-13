"""Tests for the GTN cache builder."""
from __future__ import annotations
import json
import sys
import tempfile
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

# --- Sample GTN API responses for mocking ---

SAMPLE_TOPICS = [
    {"name": "variant-analysis", "title": "Variant Analysis",
     "summary": "Genetic variant detection", "material": []},
    {"name": "transcriptomics", "title": "Transcriptomics",
     "summary": "Gene expression analysis", "material": []},
]

SAMPLE_TOPIC_DETAIL = {
    "variant-analysis": {
        "name": "variant-analysis",
        "title": "Variant Analysis",
        "summary": "Genetic variant detection",
        "materials": [
            {
                "title": "Calling very rare variants",
                "name": "rare-variants",
                "hands_on": True,
                "type": "tutorial",
                "time_estimation": "3H",
                "level": "Intermediate",
                "objectives": ["Process duplex sequencing data", "Identify rare variants"],
                "tools": ["fastqc", "bwa_mem"],
            }
        ],
    },
    "transcriptomics": {
        "name": "transcriptomics",
        "title": "Transcriptomics",
        "summary": "Gene expression analysis",
        "materials": [
            {
                "title": "Reference-based RNA-Seq data analysis",
                "name": "ref-based",
                "hands_on": True,
                "type": "tutorial",
                "time_estimation": "4H",
                "level": "Intermediate",
                "objectives": ["Perform RNA-seq alignment", "Identify differentially expressed genes"],
                "tools": ["hisat2", "featurecounts", "deseq2"],
            }
        ],
    },
}

SAMPLE_TOOL_MAP = {
    "devteam/fastqc/fastqc": {
        "tool_id": [["toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.74", "0.74"]],
        "tutorials": [
            ["variant-analysis/rare-variants", "Calling very rare variants",
             "Variant Analysis", "/topics/variant-analysis/tutorials/rare-variants/tutorial.html"],
        ],
    },
}


class TestBuildGtnCache:
    """Test build_gtn_cache produces a valid cache file."""

    def test_writes_cache_file(self, monkeypatch, tmp_path):
        from gtn_cache_builder import build_gtn_cache
        import gtn_client

        monkeypatch.setattr(gtn_client, "fetch_topics", lambda: SAMPLE_TOPICS)
        monkeypatch.setattr(gtn_client, "fetch_topic_detail",
                            lambda tid: SAMPLE_TOPIC_DETAIL[tid])
        monkeypatch.setattr(gtn_client, "fetch_tool_tutorial_map", lambda: SAMPLE_TOOL_MAP)

        cache_path = tmp_path / "gtn_cache.json"
        build_gtn_cache(output_path=cache_path)

        assert cache_path.exists()
        data = json.loads(cache_path.read_text())
        assert "topics" in data
        assert "tool_index" in data
        assert len(data["topics"]) == 2
        # Check tutorials are nested
        va = next(t for t in data["topics"] if t["name"] == "variant-analysis")
        assert len(va["tutorials"]) == 1
        assert va["tutorials"][0]["title"] == "Calling very rare variants"

    def test_atomic_write_preserves_old_on_failure(self, monkeypatch, tmp_path):
        """If fetching fails mid-build, old cache is preserved."""
        from gtn_cache_builder import build_gtn_cache
        import gtn_client

        # Write an existing cache
        cache_path = tmp_path / "gtn_cache.json"
        cache_path.write_text('{"old": true}')

        # Make fetch_topics fail
        monkeypatch.setattr(gtn_client, "fetch_topics",
                            lambda: (_ for _ in ()).throw(ConnectionError("offline")))

        with pytest.raises(ConnectionError):
            build_gtn_cache(output_path=cache_path)

        # Old cache should be preserved
        assert json.loads(cache_path.read_text()) == {"old": True}


class TestBuildSkillRecommendations:
    """Test build_skill_recommendations produces skill→tutorial mappings."""

    def test_writes_recommendations(self, monkeypatch, tmp_path):
        from gtn_cache_builder import build_skill_recommendations

        # Provide a pre-built cache
        cache_data = {
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
                            "objectives": ["Process duplex sequencing data"],
                            "tools": ["fastqc", "bwa_mem"],
                            "url": "/topics/variant-analysis/tutorials/rare-variants/tutorial.html",
                        }
                    ],
                }
            ],
            "tool_index": {},
        }
        cache_path = tmp_path / "gtn_cache.json"
        cache_path.write_text(json.dumps(cache_data))

        recs_path = tmp_path / "skill_recommendations.json"
        build_skill_recommendations(cache_path=cache_path, output_path=recs_path)

        assert recs_path.exists()
        recs = json.loads(recs_path.read_text())
        # Should have entries for ClawBio skills
        assert isinstance(recs, dict)
        assert len(recs) > 0
        # Each entry should have concepts and tutorials
        for skill_name, entry in recs.items():
            assert "concepts" in entry
            assert "tutorials" in entry
