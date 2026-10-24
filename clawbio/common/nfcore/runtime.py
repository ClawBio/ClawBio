"""Runtime probes and environment checks; adapters keep their validation policy."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from collections.abc import Mapping, Sequence

from .errors import SkillError


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


def reported_version(text: str) -> str:
    """The version exactly as a tool prints it (e.g. '26.04.3').

    parse_version_tuple is for comparison only: rebuilding a string from it drops
    zero padding (26.04.3 -> 26.4.3), which is no real Nextflow release and would
    break NXF_VER and conda pins on replay.
    """
    for pattern in (r"\b(\d+\.\d+\.\d+)\b", r"\b(\d+\.\d+)\b", r"\b(\d+)\b"):
        m = re.search(pattern, text)
        if m:
            return m.group(1)
    return ""


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


def nxf_offline(environ: Mapping[str, str] | None = None) -> bool:
    """Whether NXF_OFFLINE asks for no network access (true, 1, yes or on)."""
    value = (os.environ if environ is None else environ).get("NXF_OFFLINE")
    return str(value or "").strip().lower() in {"true", "1", "yes", "on"}


def require_network_for_demo() -> None:
    """Refuse a --demo run under NXF_OFFLINE before Nextflow starts.

    The demos run nf-core's upstream ``-profile test``, whose inputs are remote
    GitHub URLs by design; offline, nf-schema aborts parameter validation with a
    confusing "does not exist". Downloading nf-core's public test data does not
    conflict with local-first, which governs the user's genetic data.
    """
    if nxf_offline():
        raise SkillError(
            stage="preflight",
            error_code="DEMO_REQUIRES_NETWORK",
            message=(
                "--demo runs the upstream nf-core `-profile test`, which downloads its test "
                "datasets and reference files from GitHub, but NXF_OFFLINE is set — Nextflow "
                "would abort during parameter validation because the remote test inputs cannot "
                "be reached."
            ),
            fix=(
                "Unset NXF_OFFLINE to run the demo (it fetches only nf-core's public test data — "
                "no user or genetic data is uploaded), or run a real analysis with your own local "
                "--input samplesheet and references, which is fully offline."
            ),
            details={"NXF_OFFLINE": os.environ.get("NXF_OFFLINE")},
        )


def docker_cpu_count() -> int | None:
    """Best-effort CPU capacity of the active Docker daemon."""
    if not shutil.which("docker"):
        return None
    try:
        result = subprocess.run(
            ["docker", "info", "--format", "{{.NCPU}}"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if result.returncode == 0 and result.stdout.strip().isdigit():
            return int(result.stdout.strip()) or None
    except (OSError, subprocess.SubprocessError):
        # Capacity probing is optional; callers fall back to host capacity.
        pass
    return None
