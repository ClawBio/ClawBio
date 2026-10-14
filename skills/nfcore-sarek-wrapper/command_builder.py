from __future__ import annotations

import sys
from pathlib import Path

from clawbio.common.nfcore.commands import build_nextflow_command as _shared_command

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))
sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules  # noqa: E402

purge_foreign_bare_modules("schemas")

from schemas import NEXTFLOW_ENV  # noqa: E402

def compose_profile(
    profile: str,
    *,
    demo: bool = False,
    arm: bool = False,
    gpu: bool = False,
    spark: bool = False,
    mutect_profile: bool = False,
) -> str:
    # Canonical ordering: test, backend, then modifiers. Strict config syntax
    # applies profiles in CLI order; legacy syntax uses their declaration order
    # in nextflow.config. Preserve upstream behavior under either parser.
    parts = [part.strip() for part in profile.split(",") if part.strip()]
    if demo:
        parts = [part for part in parts if part != "test"]
        parts.insert(0, "test")
    # Append modifiers last for deterministic strict-parser composition.
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
    return _shared_command(
        pipeline_source=pipeline_source,
        profile=compose_profile(profile, demo=demo, arm=arm, gpu=gpu, spark=spark, mutect_profile=mutect_profile),
        params_path=params_path,
        resume=resume,
        work_dir=work_dir,
        extra_configs=extra_configs,
        env=NEXTFLOW_ENV,
    )
