"""Job endpoints — submit, list, check status, and resume skill runs."""

from __future__ import annotations

import asyncio

from fastapi import APIRouter, BackgroundTasks, HTTPException

from api.jobs import store
from api.models import JobCreate, JobDetail, JobStatus, JobSummary

router = APIRouter(prefix="/jobs", tags=["jobs"])


def _to_summary(job: JobDetail) -> JobSummary:
    return JobSummary(
        job_id=job.job_id,
        skill=job.skill,
        status=job.status,
        created_at=job.created_at,
        duration_seconds=job.duration_seconds,
    )


# ── POST /jobs — submit a new skill run ────────────────────────────────────


@router.post("", response_model=JobDetail, status_code=202)
async def create_job(req: JobCreate, background_tasks: BackgroundTasks):
    """Submit a skill run. Returns immediately with a job ID; the skill
    executes in the background.  Poll ``GET /jobs/{job_id}`` to track
    progress."""
    job = await store.create(req)
    background_tasks.add_task(store.run, job.job_id)
    return job


# ── GET /jobs — list jobs ──────────────────────────────────────────────────


@router.get("", response_model=list[JobSummary])
async def list_jobs(skill: str | None = None):
    """List all jobs, newest first.  Optionally filter by skill name."""
    jobs = await store.list_all(skill=skill)
    return [_to_summary(j) for j in jobs]


# ── GET /jobs/{job_id} — get full status / results ────────────────────────


@router.get("/{job_id}", response_model=JobDetail)
async def get_job(job_id: str):
    """Retrieve the full detail of a job.  This is the primary 'resume'
    check — submit a job, close your laptop, come back later, and call
    this endpoint to pick up where you left off."""
    job = await store.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")
    return job


# ── POST /jobs/{job_id}/resume — re-run a failed/timed-out job ────────────


@router.post("/{job_id}/resume", response_model=JobDetail, status_code=202)
async def resume_job(job_id: str, background_tasks: BackgroundTasks):
    """Re-submit a failed or timed-out job using its original parameters.

    Creates a *new* job with a fresh ID so the original record is preserved
    for auditing.
    """
    original = await store.get(job_id)
    if not original:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    if original.status not in (JobStatus.failed, JobStatus.timed_out):
        raise HTTPException(
            status_code=409,
            detail=f"Job {job_id} has status '{original.status.value}' — only failed or timed_out jobs can be resumed",
        )

    # Create a brand-new job from the same request
    new_job = await store.create(original.request)
    background_tasks.add_task(store.run, new_job.job_id)
    return new_job
