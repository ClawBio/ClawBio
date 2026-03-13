"""Tests for the knowledge-guide importable API."""
from __future__ import annotations
import json
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))


class TestRun:
    """Test the run() entry point."""

    def test_run_with_query(self):
        from api import run
        result = run(options={"query": "variant calling", "demo": True})
        assert "summary" in result
        assert result["summary"]["matches"] > 0

    def test_run_with_output_dir(self, tmp_path):
        from api import run
        result = run(options={
            "query": "variant calling",
            "demo": True,
            "output_dir": str(tmp_path),
        })
        assert (tmp_path / "result.json").exists()


class TestGetLearnMore:
    """Test the get_learn_more() function for report integration."""

    def test_known_skill(self, tmp_path):
        """Should return tutorials for a registered skill."""
        from api import get_learn_more

        recs = {
            "pharmgx-reporter": {
                "concepts": ["pharmacogenomics", "CYP2D6"],
                "gtn_topics": ["variant-analysis"],
                "tutorials": [
                    {
                        "id": "variant-analysis/pharmgx",
                        "title": "Pharmacogenomics tutorial",
                        "topic": "variant-analysis",
                        "time": "2h",
                        "url": "https://example.com/tutorial",
                        "relevance_score": 10.0,
                    }
                ],
            }
        }
        recs_path = tmp_path / "skill_recommendations.json"
        recs_path.write_text(json.dumps(recs))

        result = get_learn_more("pharmgx-reporter", recommendations_path=recs_path)
        assert result["section_title"] == "Learn More"
        assert len(result["tutorials"]) == 1
        assert "html" in result
        assert result["html"] != ""

    def test_unknown_skill_returns_empty(self, tmp_path):
        """Should return empty section for unregistered skill."""
        from api import get_learn_more

        recs_path = tmp_path / "skill_recommendations.json"
        recs_path.write_text("{}")

        result = get_learn_more("nonexistent-skill", recommendations_path=recs_path)
        assert result["tutorials"] == []
        assert result["html"] == ""
