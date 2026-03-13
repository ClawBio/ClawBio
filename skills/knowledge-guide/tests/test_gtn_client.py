"""Tests for the GTN API client."""
from __future__ import annotations
import json
import pytest
from pathlib import Path
import sys

SKILL_DIR = Path(__file__).resolve().parent.parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))


class TestFetchTopics:
    """Test fetch_topics returns structured topic data."""

    def test_returns_dict(self, monkeypatch):
        from gtn_client import fetch_topics

        sample = [
            {"name": "variant-analysis", "title": "Variant Analysis",
             "summary": "Detecting genetic variants", "material": []}
        ]
        monkeypatch.setattr(
            "gtn_client.requests.get",
            lambda url, **kw: type("R", (), {"json": lambda s: sample, "raise_for_status": lambda s: None})(),
        )
        result = fetch_topics()
        assert isinstance(result, list)
        assert result[0]["name"] == "variant-analysis"

    def test_handles_api_error(self, monkeypatch):
        from gtn_client import fetch_topics
        import requests as req

        def raise_error(url, **kw):
            raise req.ConnectionError("offline")

        monkeypatch.setattr("gtn_client.requests.get", raise_error)
        with pytest.raises(req.ConnectionError):
            fetch_topics()


class TestFetchTopicDetail:
    """Test fetch_topic_detail returns tutorials for a topic."""

    def test_returns_tutorials(self, monkeypatch):
        from gtn_client import fetch_topic_detail

        sample = {
            "name": "variant-analysis",
            "title": "Variant Analysis",
            "materials": [
                {
                    "title": "Calling very rare variants",
                    "name": "rare-variants",
                    "hands_on": True,
                    "type": "tutorial",
                    "time_estimation": "3H",
                    "level": "Intermediate",
                    "objectives": ["Process duplex sequencing data"],
                    "tools": ["fastqc", "bwa_mem"],
                    "url": "/topics/variant-analysis/tutorials/rare-variants/tutorial.html",
                }
            ],
        }
        monkeypatch.setattr(
            "gtn_client.requests.get",
            lambda url, **kw: type("R", (), {"json": lambda s: sample, "raise_for_status": lambda s: None})(),
        )
        result = fetch_topic_detail("variant-analysis")
        assert result["name"] == "variant-analysis"
        assert len(result["materials"]) == 1
        assert result["materials"][0]["title"] == "Calling very rare variants"


class TestFetchToolMap:
    """Test fetch_tool_tutorial_map returns tool→tutorial index."""

    def test_returns_mapping(self, monkeypatch):
        from gtn_client import fetch_tool_tutorial_map

        sample = {
            "devteam/fastqc/fastqc": {
                "tool_id": [["toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.74+galaxy1", "0.74+galaxy1"]],
                "tutorials": [
                    ["assembly/general-introduction", "An Introduction to Genome Assembly", "Assembly",
                     "/topics/assembly/tutorials/general-introduction/tutorial.html"]
                ],
            }
        }
        monkeypatch.setattr(
            "gtn_client.requests.get",
            lambda url, **kw: type("R", (), {"json": lambda s: sample, "raise_for_status": lambda s: None})(),
        )
        result = fetch_tool_tutorial_map()
        assert "devteam/fastqc/fastqc" in result
        assert len(result["devteam/fastqc/fastqc"]["tutorials"]) == 1


class TestFetchTutorialContent:
    """Test fetch_tutorial_content returns full tutorial body."""

    def test_returns_content(self, monkeypatch):
        from gtn_client import fetch_tutorial_content

        sample = {
            "name": "rare-variants",
            "title": "Calling very rare variants",
            "content": "# Introduction\nThis tutorial covers...",
        }
        monkeypatch.setattr(
            "gtn_client.requests.get",
            lambda url, **kw: type("R", (), {"json": lambda s: sample, "raise_for_status": lambda s: None})(),
        )
        result = fetch_tutorial_content("variant-analysis", "rare-variants")
        assert "content" in result
        assert "Introduction" in result["content"]

    def test_returns_none_on_404(self, monkeypatch):
        from gtn_client import fetch_tutorial_content
        import requests as req

        def raise_404(url, **kw):
            resp = req.Response()
            resp.status_code = 404
            raise req.HTTPError(response=resp)

        monkeypatch.setattr("gtn_client.requests.get", raise_404)
        result = fetch_tutorial_content("nonexistent", "nonexistent")
        assert result is None
