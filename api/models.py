"""Pydantic models for the ClawBio API."""

from __future__ import annotations

import enum
from datetime import datetime

from pydantic import BaseModel, Field


class JobStatus(str, enum.Enum):
    pending = "pending"
    running = "running"
    completed = "completed"
    failed = "failed"
    timed_out = "timed_out"


# ── Request models ──────────────────────────────────────────────────────────


class JobCreate(BaseModel):
    """Submit a new skill run."""

    skill: str = Field(..., description="Skill name (e.g. pharmgx, equity, nutrigx)")
    input_path: str | None = Field(None, description="Path to input file")
    output_dir: str | None = Field(None, description="Output directory (auto-generated if omitted)")
    demo: bool = Field(False, description="Run with built-in demo data")
    extra_args: list[str] | None = Field(None, description="Extra CLI flags forwarded to the skill")
    timeout: int = Field(300, ge=10, le=3600, description="Timeout in seconds")


# ── Response models ─────────────────────────────────────────────────────────


class JobSummary(BaseModel):
    """Lightweight representation returned in list views."""

    job_id: str
    skill: str
    status: JobStatus
    created_at: datetime
    duration_seconds: float | None = None


class JobDetail(BaseModel):
    """Full job record including outputs."""

    job_id: str
    skill: str
    status: JobStatus
    created_at: datetime
    started_at: datetime | None = None
    finished_at: datetime | None = None
    duration_seconds: float | None = None
    output_dir: str | None = None
    files: list[str] = []
    stdout: str = ""
    stderr: str = ""
    exit_code: int | None = None
    # Original request, so /resume can replay it
    request: JobCreate
