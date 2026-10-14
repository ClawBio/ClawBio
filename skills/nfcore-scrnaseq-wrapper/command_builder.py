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

def _command_path(value: Path | str) -> str:
    text = str(value)
    if "://" in text:
        return text
    return Path(text).as_posix()


def build_nextflow_command(
    *,
    pipeline_source: dict[str, object],
    profile: str,
    params_path: Path,
    resume: bool,
    work_dir: Path | str | None = None,
    extra_configs: list[Path] | None = None,
) -> tuple[list[str], str]:
    return _shared_command(
        pipeline_source=pipeline_source,
        profile=profile,
        params_path=params_path,
        resume=resume,
        work_dir=_command_path(work_dir) if work_dir is not None else None,
        extra_configs=extra_configs,
        env=NEXTFLOW_ENV,
    )
