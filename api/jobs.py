"""In-memory job store for tracking async skill runs.

Swap this out for Redis/Postgres/SQLite when persistence is needed.
"""

from __future__ import annotations

import asyncio
import time
import uuid
from datetime import datetime, timezone

from api.models import JobCreate, JobDetail, JobStatus

# Import the existing ClawBio skill runner
import sys
from pathlib import Path

CLAWBIO_DIR = Path(__file__).resolve().parent.parent
if str(CLAWBIO_DIR) not in sys.path:
    sys.path.insert(0, str(CLAWBIO_DIR))

from clawbio import run_skill  # noqa: E402


class JobStore:
    """Thread-safe in-memory job store."""

    def __init__(self) -> None:
        self._jobs: dict[str, JobDetail] = {}
        self._lock = asyncio.Lock()

    # ── Queries ─────────────────────────────────────────────────────────

    async def get(self, job_id: str) -> JobDetail | None:
        async with self._lock:
            return self._jobs.get(job_id)

    async def list_all(self, skill: str | None = None) -> list[JobDetail]:
        async with self._lock:
            jobs = list(self._jobs.values())
        if skill:
            jobs = [j for j in jobs if j.skill == skill]
        return sorted(jobs, key=lambda j: j.created_at, reverse=True)

    # ── Mutations ───────────────────────────────────────────────────────

    async def create(self, req: JobCreate) -> JobDetail:
        job_id = uuid.uuid4().hex[:12]
        now = datetime.now(timezone.utc)
        job = JobDetail(
            job_id=job_id,
            skill=req.skill,
            status=JobStatus.pending,
            created_at=now,
            request=req,
        )
        async with self._lock:
            self._jobs[job_id] = job
        return job

    async def _update(self, job_id: str, **fields) -> JobDetail | None:
        async with self._lock:
            job = self._jobs.get(job_id)
            if not job:
                return None
            updated = job.model_copy(update=fields)
            self._jobs[job_id] = updated
            return updated

    async def run(self, job_id: str) -> JobDetail | None:
        """Execute the skill in a thread pool and update the job record."""
        job = await self.get(job_id)
        if not job:
            return None

        now = datetime.now(timezone.utc)
        await self._update(job_id, status=JobStatus.running, started_at=now)

        req = job.request
        loop = asyncio.get_event_loop()
        t0 = time.time()

        try:
            result = await loop.run_in_executor(
                None,
                lambda: run_skill(
                    skill_name=req.skill,
                    input_path=req.input_path,
                    output_dir=req.output_dir,
                    demo=req.demo,
                    extra_args=req.extra_args,
                    timeout=req.timeout,
                ),
            )
        except Exception as exc:
            duration = round(time.time() - t0, 2)
            return await self._update(
                job_id,
                status=JobStatus.failed,
                finished_at=datetime.now(timezone.utc),
                duration_seconds=duration,
                stderr=str(exc),
                exit_code=-1,
            )

        finished = datetime.now(timezone.utc)
        duration = round(time.time() - t0, 2)

        if result.get("stderr", "").endswith("seconds.") and "Timed out" in result.get("stderr", ""):
            status = JobStatus.timed_out
        elif result.get("success"):
            status = JobStatus.completed
        else:
            status = JobStatus.failed

        return await self._update(
            job_id,
            status=status,
            finished_at=finished,
            duration_seconds=duration,
            output_dir=result.get("output_dir"),
            files=result.get("files", []),
            stdout=result.get("stdout", ""),
            stderr=result.get("stderr", ""),
            exit_code=result.get("exit_code"),
        )


# Module-level singleton — shared across the app
store = JobStore()
