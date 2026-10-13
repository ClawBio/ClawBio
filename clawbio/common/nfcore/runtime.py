"""Runtime probes return facts; adapters retain validation and fallback policy."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from collections.abc import Sequence


def run_probe(
    command: Sequence[str], *, timeout_seconds: float
) -> subprocess.CompletedProcess[str] | None:
    try:
        return subprocess.run(
            command,
            capture_output=True,
            text=True,
            errors="replace",
            timeout=timeout_seconds,
            check=False,
        )
    except (subprocess.TimeoutExpired, OSError):
        return None


def parse_version_tuple(text: str) -> tuple[int, ...]:
    m = re.search(r"\b(\d+)\.(\d+)\.(\d+)\b", text)
    if m:
        return (int(m.group(1)), int(m.group(2)), int(m.group(3)))
    m = re.search(r"\b(\d+)\.(\d+)\b", text)
    if m:
        return (int(m.group(1)), int(m.group(2)))
    m = re.search(r"\b(\d+)\b", text)
    if m:
        return (int(m.group(1)),)
    return ()


def pad_version(t: tuple[int, ...], length: int = 3) -> tuple[int, ...]:
    return t + (0,) * max(0, length - len(t))


def host_memory_gb() -> int | None:
    """Return total physical RAM in whole GB, or None when it cannot be determined."""
    try:
        page_size = os.sysconf("SC_PAGE_SIZE")
        phys_pages = os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError, AttributeError):
        return None
    if page_size <= 0 or phys_pages <= 0:
        return None
    return int((page_size * phys_pages) / (1024**3))


def docker_memory_gb() -> int | None:
    """Return the Docker VM's total memory in whole GB, or None when unknown.

    On macOS, Docker Desktop runs containers in a Linux VM whose RAM is usually
    smaller than the host's. `docker info` reports that VM's MemTotal (bytes), which
    is the true ceiling a container process can use before being OOM-killed
    (audit F-8). Best-effort: any failure (docker absent, daemon down, unparseable
    output) returns None so the caller falls back to the host-derived heuristic.
    """
    if not shutil.which("docker"):
        return None
    try:
        proc = subprocess.run(
            ["docker", "info", "--format", "{{.MemTotal}}"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if proc.returncode != 0:
        return None
    raw = proc.stdout.strip()
    if not raw.isdigit():
        return None
    gb = int(int(raw) / (1024**3))
    return gb or None


def host_memory_bytes() -> int | None:
    """Raw signal; callers decide whether zero/negative means unknown or a floor."""
    try:
        return os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    except (ValueError, OSError, AttributeError):
        return None
