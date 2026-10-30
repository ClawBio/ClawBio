"""Run-control options common to the nf-core wrappers' command lines."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from .paths import is_remote_uri


def timeout_hours(value: str) -> float:
    """argparse type for --timeout-hours: a finite, non-negative number of hours.

    0 disables the wall-clock cap. NaN and infinity are usage errors: neither can
    become a wait in seconds, so accepting them only moved the failure past
    preflight, where it surfaced as an internal error.
    """
    try:
        hours = float(value)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError(
            f"--timeout-hours must be a number, got {value!r}"
        ) from None
    if not math.isfinite(hours * 3600):
        raise argparse.ArgumentTypeError(
            f"--timeout-hours must be a finite number of hours, got {value!r}"
        )
    if hours < 0:
        raise argparse.ArgumentTypeError(
            "--timeout-hours must be >= 0 (0 disables the timeout)"
        )
    return hours


def timeout_seconds(hours: float | None, *, default_hours: float) -> int | None:
    """Wall-clock cap in whole seconds; None (no cap) when hours is 0.

    A missing value keeps the default cap. Any positive value waits at least one
    second, because a zero-second wait would terminate Nextflow on launch.
    """
    if hours is None:
        hours = default_hours
    if hours == 0:
        return None
    return max(1, round(hours * 3600))


def resolve_work_dir(raw: str | None, *, output_dir: Path) -> Path | str:
    """Nextflow work directory: <output>/upstream/work unless --work-dir is given.

    Object-store URIs (s3://, gs://, az://) are kept verbatim for cloud executors;
    local paths are resolved so replay does not depend on the launch directory.
    """
    if not raw or not str(raw).strip():
        return output_dir / "upstream" / "work"
    work_dir = str(raw).strip()
    if is_remote_uri(work_dir):
        return work_dir
    return Path(work_dir).expanduser().resolve()


def canonical_work_dir(args, output_dir: Path) -> str:
    """Stable manifest spelling, with the shared resolver owning path handling."""
    raw = getattr(args, "work_dir", None)
    if not raw or not str(raw).strip():
        return "upstream/work"
    work_dir = resolve_work_dir(raw, output_dir=output_dir)
    return work_dir.as_posix() if isinstance(work_dir, Path) else work_dir
