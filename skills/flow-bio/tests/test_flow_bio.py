"""Tests for flow_bio.py — Flow.bio API bridge skill."""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
SCRIPT = SKILL_DIR / "flow_bio.py"

# Ensure the skill directory is importable
sys.path.insert(0, str(SKILL_DIR))
import flow_bio  # noqa: E402


# ---------------------------------------------------------------------------
# Demo mode tests
# ---------------------------------------------------------------------------


class TestDemoMode:
    """Tests for demo/overview mode with mocked API."""

    @pytest.fixture(autouse=True)
    def _mock_api(self):
        """Mock FlowClient so tests don't hit the real API."""
        fake_pipelines = [
            {"name": "Genomics", "subcategories": [
                {"name": "RNA", "pipelines": [
                    {"id": "p1", "name": "nf-core/rnaseq", "description": "RNA-seq pipeline"},
                ]}
            ]}
        ]
        fake_organisms = [
            {"id": "1", "name": "Homo sapiens"},
        ]
        fake_sample_types = [
            {"identifier": "rna_seq", "name": "RNA-Seq"},
        ]

        with patch.object(flow_bio.FlowClient, "login", return_value={"token": "t", "user": {"username": "tester"}}), \
             patch.object(flow_bio.FlowClient, "get_pipelines", return_value=fake_pipelines), \
             patch.object(flow_bio.FlowClient, "get_organisms", return_value=fake_organisms), \
             patch.object(flow_bio.FlowClient, "get_sample_types", return_value=fake_sample_types), \
             patch.object(flow_bio.FlowClient, "get_samples_owned", return_value=[
                 {"id": "s1", "name": "Sample1", "sample_type": {"name": "RNA-Seq"}, "organism": {"name": "Homo sapiens"}, "project_name": "Proj1"},
             ]), \
             patch.object(flow_bio.FlowClient, "get_projects_owned", return_value=[
                 {"id": "pr1", "name": "Proj1", "description": "Test project"},
             ]), \
             patch.object(flow_bio.FlowClient, "get_executions_owned", return_value=[
                 {"id": "e1", "status": "COMPLETED", "pipelineVersion": {"pipeline": {"name": "nf-core/rnaseq"}}, "created": "2026-03-19T10:00:00Z"},
             ]):
            yield

    def test_demo_creates_output_files(self, tmp_path):
        result = flow_bio.run_demo(tmp_path, username="u", password="p")
        assert (tmp_path / "result.json").exists()
        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "reproducibility" / "commands.sh").exists()
        assert (tmp_path / "reproducibility" / "environment.yml").exists()

    def test_demo_result_has_expected_keys(self, tmp_path):
        result = flow_bio.run_demo(tmp_path, username="u", password="p")
        assert result["mode"] == "live"
        assert result["authenticated"] is True
        assert len(result["pipelines"]) == 1
        assert len(result["organisms"]) == 1
        assert len(result["sample_types"]) == 1
        assert len(result["samples"]) == 1
        assert len(result["projects"]) == 1
        assert len(result["executions"]) == 1

    def test_demo_result_json_is_valid(self, tmp_path):
        flow_bio.run_demo(tmp_path, username="u", password="p")
        data = json.loads((tmp_path / "result.json").read_text())
        assert data["mode"] == "live"
        assert "pipelines" in data

    def test_demo_report_contains_disclaimer(self, tmp_path):
        flow_bio.run_demo(tmp_path, username="u", password="p")
        report = (tmp_path / "report.md").read_text()
        assert "research and educational tool" in report
        assert "not a medical device" in report

    def test_demo_report_contains_live_data(self, tmp_path):
        flow_bio.run_demo(tmp_path, username="u", password="p")
        report = (tmp_path / "report.md").read_text()
        assert "nf-core/rnaseq" in report
        assert "Sample1" in report
        assert "Proj1" in report

    def test_demo_unauthenticated_still_works(self, tmp_path):
        """Without credentials, public endpoints should still be fetched."""
        result = flow_bio.run_demo(tmp_path)
        assert result["mode"] == "live"
        assert result["authenticated"] is False
        # Public data still present
        assert len(result["pipelines"]) == 1
        assert len(result["organisms"]) == 1
        # Auth-required data empty
        assert len(result["samples"]) == 0
        assert len(result["projects"]) == 0
        assert len(result["executions"]) == 0


# ---------------------------------------------------------------------------
# Flatten pipelines helper
# ---------------------------------------------------------------------------


class TestFlattenPipelines:
    """Tests for _flatten_pipelines helper."""

    def test_flatten_nested_structure(self):
        categories = [
            {"name": "Cat1", "subcategories": [
                {"name": "Sub1", "pipelines": [
                    {"id": "p1", "name": "tool-a", "description": "Desc A"},
                    {"id": "p2", "name": "tool-b", "description": "Desc B"},
                ]},
                {"name": "Sub2", "pipelines": [
                    {"id": "p3", "name": "tool-c", "description": "Desc C"},
                ]},
            ]},
        ]
        flat = flow_bio._flatten_pipelines(categories)
        assert len(flat) == 3
        assert flat[0]["name"] == "tool-a"
        assert flat[0]["category"] == "Cat1"
        assert flat[0]["subcategory"] == "Sub1"
        assert flat[2]["subcategory"] == "Sub2"

    def test_empty_categories(self):
        assert flow_bio._flatten_pipelines([]) == []

    def test_empty_subcategories(self):
        categories = [{"name": "Cat1", "subcategories": []}]
        assert flow_bio._flatten_pipelines(categories) == []


# ---------------------------------------------------------------------------
# FlowClient unit tests (mocked HTTP)
# ---------------------------------------------------------------------------


class TestFlowClient:
    """Unit tests for FlowClient with mocked HTTP responses."""

    @pytest.fixture()
    def mock_requests(self):
        with patch("flow_bio.FlowClient.__init__", return_value=None) as _:
            client = flow_bio.FlowClient.__new__(flow_bio.FlowClient)
            import requests
            client._requests = requests
            client.base_url = "https://app.flow.bio/api"
            client._token = "test_token_123"
            client._last_token_refresh = time.time()
            client._session = MagicMock()
            yield client

    def test_login_sets_token(self, mock_requests):
        client = mock_requests
        mock_resp = MagicMock()
        mock_resp.json.return_value = {"token": "new_token", "user": {"username": "testuser"}}
        mock_resp.raise_for_status = MagicMock()
        client._session.post.return_value = mock_resp

        result = client.login("testuser", "testpass")
        assert client._token == "new_token"
        assert result["user"]["username"] == "testuser"

    def test_is_authenticated(self, mock_requests):
        client = mock_requests
        assert client.is_authenticated is True
        client._token = None
        assert client.is_authenticated is False

    def test_set_token(self, mock_requests):
        client = mock_requests
        client.set_token("abc123")
        assert client._token == "abc123"
        assert client._last_token_refresh is not None

    def test_ensure_token_no_refresh_when_fresh(self, mock_requests):
        """Token should NOT be refreshed if it was obtained recently."""
        client = mock_requests
        client._last_token_refresh = time.time()  # just now
        client._ensure_token()
        # GET /token should not have been called
        client._session.get.assert_not_called()

    def test_ensure_token_refreshes_when_stale(self, mock_requests):
        """Token should be refreshed if it is older than _TOKEN_REFRESH_AFTER."""
        client = mock_requests
        client._last_token_refresh = time.time() - 300  # 5 minutes ago (>240s threshold)
        mock_resp = MagicMock()
        mock_resp.json.return_value = {"token": "refreshed_token"}
        mock_resp.raise_for_status = MagicMock()
        client._session.get.return_value = mock_resp

        client._ensure_token()
        client._session.get.assert_called_once_with(f"{client.base_url}/token")
        assert client._token == "refreshed_token"

    def test_ensure_token_silent_on_failure(self, mock_requests):
        """If refresh fails, _ensure_token should not raise (best-effort)."""
        client = mock_requests
        client._last_token_refresh = time.time() - 300
        client._session.get.side_effect = Exception("network error")

        # Should not raise
        client._ensure_token()
        assert client._token == "test_token_123"  # unchanged

    def test_get_pipelines(self, mock_requests):
        client = mock_requests
        mock_resp = MagicMock()
        mock_resp.json.return_value = [{"id": "p1", "name": "test"}]
        mock_resp.raise_for_status = MagicMock()
        client._session.get.return_value = mock_resp

        result = client.get_pipelines()
        assert len(result) == 1
        assert result[0]["name"] == "test"

    def test_get_sample(self, mock_requests):
        client = mock_requests
        mock_resp = MagicMock()
        mock_resp.json.return_value = {"id": "s1", "name": "Sample1", "sample_type": "RNA-Seq"}
        mock_resp.raise_for_status = MagicMock()
        client._session.get.return_value = mock_resp

        result = client.get_sample("s1")
        assert result["name"] == "Sample1"

    def test_search(self, mock_requests):
        client = mock_requests
        mock_resp = MagicMock()
        mock_resp.json.return_value = {"samples": [{"name": "hit1"}], "projects": []}
        mock_resp.raise_for_status = MagicMock()
        client._session.get.return_value = mock_resp

        result = client.search("RNA-seq")
        assert "samples" in result


# ---------------------------------------------------------------------------
# Report generation tests
# ---------------------------------------------------------------------------


class TestReportGeneration:
    """Tests for report writing functions."""

    def test_write_report_demo(self, tmp_path):
        data = {
            "mode": "live",
            "authenticated": True,
            "user": "tester",
            "pipelines": [{"name": "nf-core/rnaseq", "subcategory": "RNA", "description": "RNA pipeline"}],
            "samples": [{"id": "s1", "name": "S1", "sample_type": {"name": "RNA-Seq"}, "organism": {"name": "Human"}, "project_name": "P1"}],
            "projects": [{"id": "pr1", "name": "P1", "description": "Test"}],
            "executions": [{"id": "e1", "status": "COMPLETED", "pipelineVersion": {"pipeline": {"name": "rnaseq"}}, "created": "2026-01-01"}],
            "organisms": [{"id": "1", "name": "Homo sapiens"}],
            "sample_types": [{"identifier": "rna_seq", "name": "RNA-Seq"}],
        }
        path = flow_bio.write_report(tmp_path, "demo", data)
        assert path.exists()
        content = path.read_text()
        assert "Flow Bio Bridge Report" in content
        assert "nf-core/rnaseq" in content

    def test_write_report_pipelines(self, tmp_path):
        data = {"pipelines": [{"name": "test-pipe", "version": "1.0", "description": "A test"}]}
        path = flow_bio.write_report(tmp_path, "pipelines", data)
        content = path.read_text()
        assert "test-pipe" in content

    def test_write_report_search(self, tmp_path):
        data = {"query": "RNA", "results": {"samples": [{"name": "hit1"}]}}
        path = flow_bio.write_report(tmp_path, "search", data)
        content = path.read_text()
        assert "RNA" in content

    def test_write_result_json(self, tmp_path):
        path = flow_bio.write_result_json(tmp_path, "test", {"key": "value"})
        data = json.loads(path.read_text())
        assert data["skill"] == "flow-bio"
        assert data["summary"]["action"] == "test"
        assert data["data"]["key"] == "value"

    def test_write_reproducibility(self, tmp_path):
        flow_bio.write_reproducibility(tmp_path, "python flow_bio.py --demo", "https://app.flow.bio/api")
        assert (tmp_path / "reproducibility" / "commands.sh").exists()
        assert (tmp_path / "reproducibility" / "environment.yml").exists()
        content = (tmp_path / "reproducibility" / "commands.sh").read_text()
        assert "flow_bio.py --demo" in content


# ---------------------------------------------------------------------------
# CLI integration tests
# ---------------------------------------------------------------------------


class TestCLI:
    """Integration tests for the CLI entry point."""

    def test_no_args_shows_help(self):
        result = subprocess.run(
            [sys.executable, str(SCRIPT)],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "Flow Bio Bridge" in result.stdout or "usage:" in result.stdout.lower()

    def test_help_flag(self):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "--demo" in result.stdout
        assert "--pipelines" in result.stdout
        assert "--samples" in result.stdout
        assert "--search" in result.stdout
        assert "--upload-sample" in result.stdout
