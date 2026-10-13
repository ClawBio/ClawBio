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
    prokaryotic: bool = False,
    rapid_quant: bool = False,
    arm: bool = False,
) -> str:
    # Nextflow applies profiles left-to-right; later profiles override earlier ones.
    # Order: test (must come first — backend profile then overrides its container settings),
    # then the user-chosen backend (docker/singularity/…), then scientific/arch modifiers last
    # (prokaryotic, rapid_quant, arm64) so they take final precedence over base settings.
    parts = [part.strip() for part in profile.split(",") if part.strip()]
    if demo:
        parts = [part for part in parts if part != "test"]
        parts.insert(0, "test")
    if prokaryotic:
        # test_prokaryotic already inherits conf/prokaryotic.config upstream;
        # appending 'prokaryotic' separately would produce test_prokaryotic,prokaryotic.
        if "test_prokaryotic" not in parts:
            parts = [part for part in parts if part != "prokaryotic"]
            parts.append("prokaryotic")
    if rapid_quant:
        parts = [part for part in parts if part != "rapid_quant"]
        parts.append("rapid_quant")
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
    prokaryotic: bool = False,
    rapid_quant: bool = False,
    arm: bool = False,
) -> tuple[list[str], str]:
    return _shared_command(pipeline_source=pipeline_source, profile=compose_profile(profile, demo=demo, arm=arm, prokaryotic=prokaryotic, rapid_quant=rapid_quant), params_path=params_path, resume=resume, work_dir=work_dir, extra_configs=extra_configs)
