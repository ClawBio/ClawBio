"""Tests for the /jobs endpoints."""

import time
from unittest.mock import patch

from fastapi.testclient import TestClient


# ── Health check ────────────────────────────────────────────────────────────


def test_health(client: TestClient):
    resp = client.get("/health")
    assert resp.status_code == 200
    assert resp.json() == {"status": "ok"}


# ── POST /jobs ──────────────────────────────────────────────────────────────


def _mock_run_skill(**kwargs):
    """Return a canned success result without actually running a subprocess."""
    return {
        "skill": kwargs.get("skill_name", "pharmgx"),
        "success": True,
        "exit_code": 0,
        "output_dir": "/tmp/test_output",
        "files": ["report.md", "figure.png"],
        "stdout": "All done.",
        "stderr": "",
        "duration_seconds": 1.23,
    }


def _mock_run_skill_fail(**kwargs):
    """Return a canned failure result."""
    return {
        "skill": kwargs.get("skill_name", "pharmgx"),
        "success": False,
        "exit_code": 1,
        "output_dir": "/tmp/test_output",
        "files": [],
        "stdout": "",
        "stderr": "Something went wrong",
        "duration_seconds": 0.5,
    }


def test_create_job_returns_202(client: TestClient):
    with patch("api.jobs.run_skill", side_effect=_mock_run_skill):
        resp = client.post("/jobs", json={"skill": "pharmgx", "demo": True})
    assert resp.status_code == 202
    body = resp.json()
    assert "job_id" in body
    assert body["skill"] == "pharmgx"
    assert body["status"] == "pending"


def test_create_job_and_poll_completion(client: TestClient):
    with patch("api.jobs.run_skill", side_effect=_mock_run_skill):
        resp = client.post("/jobs", json={"skill": "pharmgx", "demo": True})
        job_id = resp.json()["job_id"]

        # Background task runs synchronously in TestClient, so it should
        # already be finished by the time we poll.
        detail = client.get(f"/jobs/{job_id}").json()
        assert detail["status"] == "completed"
        assert detail["exit_code"] == 0
        assert "report.md" in detail["files"]
        assert detail["stdout"] == "All done."


# ── GET /jobs ───────────────────────────────────────────────────────────────


def test_list_jobs_empty(client: TestClient):
    resp = client.get("/jobs")
    assert resp.status_code == 200
    assert resp.json() == []


def test_list_jobs_filter_by_skill(client: TestClient):
    with patch("api.jobs.run_skill", side_effect=_mock_run_skill):
        client.post("/jobs", json={"skill": "pharmgx", "demo": True})
        client.post("/jobs", json={"skill": "equity", "demo": True})

    all_jobs = client.get("/jobs").json()
    assert len(all_jobs) == 2

    pharmgx_only = client.get("/jobs", params={"skill": "pharmgx"}).json()
    assert len(pharmgx_only) == 1
    assert pharmgx_only[0]["skill"] == "pharmgx"


# ── GET /jobs/{job_id} ──────────────────────────────────────────────────────


def test_get_nonexistent_job_returns_404(client: TestClient):
    resp = client.get("/jobs/doesnotexist")
    assert resp.status_code == 404


# ── POST /jobs/{job_id}/resume ──────────────────────────────────────────────


def test_resume_failed_job(client: TestClient):
    with patch("api.jobs.run_skill", side_effect=_mock_run_skill_fail):
        resp = client.post("/jobs", json={"skill": "pharmgx", "demo": True})
        original_id = resp.json()["job_id"]

    # Original job should be failed
    detail = client.get(f"/jobs/{original_id}").json()
    assert detail["status"] == "failed"

    # Resume it — should create a new job
    with patch("api.jobs.run_skill", side_effect=_mock_run_skill):
        resp = client.post(f"/jobs/{original_id}/resume")
    assert resp.status_code == 202
    new_id = resp.json()["job_id"]
    assert new_id != original_id

    # New job should complete
    new_detail = client.get(f"/jobs/{new_id}").json()
    assert new_detail["status"] == "completed"


def test_resume_completed_job_returns_409(client: TestClient):
    with patch("api.jobs.run_skill", side_effect=_mock_run_skill):
        resp = client.post("/jobs", json={"skill": "pharmgx", "demo": True})
        job_id = resp.json()["job_id"]

    resp = client.post(f"/jobs/{job_id}/resume")
    assert resp.status_code == 409
    assert "completed" in resp.json()["detail"]


def test_resume_nonexistent_job_returns_404(client: TestClient):
    resp = client.post("/jobs/doesnotexist/resume")
    assert resp.status_code == 404


# ── Validation ──────────────────────────────────────────────────────────────


def test_create_job_missing_skill_returns_422(client: TestClient):
    resp = client.post("/jobs", json={"demo": True})
    assert resp.status_code == 422


def test_create_job_with_input_path(client: TestClient):
    with patch("api.jobs.run_skill", side_effect=_mock_run_skill):
        resp = client.post("/jobs", json={
            "skill": "pharmgx",
            "input_path": "/data/patient.txt",
            "output_dir": "/tmp/out",
        })
    assert resp.status_code == 202
    body = resp.json()
    assert body["request"]["input_path"] == "/data/patient.txt"
    assert body["request"]["output_dir"] == "/tmp/out"
