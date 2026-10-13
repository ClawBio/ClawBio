from __future__ import annotations

from clawbio.common.nfcore.commands import build_nextflow_command as _shared_command

from pathlib import Path
import shlex
import subprocess
import sys


def compose_profile(
    profile: str,
    *,
    demo: bool = False,
    arm: bool = False,
    gpu: bool = False,
    spark: bool = False,
    mutect_profile: bool = False,
) -> str:
    # Nextflow applies profiles left-to-right; later profiles override earlier ones.
    # Order: test (must come first — backend profile then overrides its container settings),
    # then the user-chosen backend (docker/singularity/…), then scientific/runtime modifiers
    # last (gpu, spark, mutect, arm64) so they take final precedence over base settings.
    parts = [part.strip() for part in profile.split(",") if part.strip()]
    if demo:
        parts = [part for part in parts if part != "test"]
        parts.insert(0, "test")
    # Modifiers appended last (Nextflow: later overrides earlier).
    if gpu:
        parts = [part for part in parts if part != "gpu"]
        parts.append("gpu")
    if spark:
        parts = [part for part in parts if part != "spark"]
        parts.append("spark")
    if mutect_profile:
        parts = [part for part in parts if part != "mutect"]
        parts.append("mutect")
    if arm:
        parts = [part for part in parts if part != "arm64"]
        parts.append("arm64")
    return ",".join(dict.fromkeys(parts))


def build_nextflow_command(
    *,
    pipeline_source: dict[str, str | bool],
    profile: str,
    params_path: Path,
    resume: bool,
    work_dir: Path | str | None = None,
    extra_configs: list[Path] | None = None,
    demo: bool = False,
    arm: bool = False,
    gpu: bool = False,
    spark: bool = False,
    mutect_profile: bool = False,
) -> tuple[list[str], str]:
    return _shared_command(pipeline_source=pipeline_source, profile=compose_profile(profile, demo=demo, arm=arm, gpu=gpu, spark=spark, mutect_profile=mutect_profile), params_path=params_path, resume=resume, work_dir=work_dir, extra_configs=extra_configs)
