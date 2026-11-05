from __future__ import annotations

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
    source_kind = str(pipeline_source["source_kind"])
    source_ref = str(pipeline_source["source_ref"])
    resolved_version = str(pipeline_source["resolved_version"])
    composed_profile = compose_profile(
        profile,
        demo=demo,
        arm=arm,
        gpu=gpu,
        spark=spark,
        mutect_profile=mutect_profile,
    )

    command = ["nextflow", "run"]
    if source_kind == "local_checkout":
        command.append(Path(source_ref).as_posix())
    else:
        command.extend([source_ref, "-r", resolved_version])

    command.extend(["-profile", composed_profile, "-params-file", params_path.as_posix()])
    if work_dir is not None:
        # work_dir may be a local Path or a remote object-store URI (str, e.g.
        # s3://bucket/work); preserve URIs verbatim, posix-normalise local paths.
        work_dir_arg = work_dir if isinstance(work_dir, str) else work_dir.as_posix()
        command.extend(["-work-dir", work_dir_arg])
    for cfg in extra_configs or []:
        command.extend(["-c", cfg.as_posix()])
    if resume:
        command.append("-resume")

    if sys.platform == "win32":
        command_str = subprocess.list2cmdline(command)
    else:
        command_str = " ".join(shlex.quote(part) for part in command)

    return command, command_str
