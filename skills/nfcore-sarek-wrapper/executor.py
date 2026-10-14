"""nf-core/sarek execution policy on top of the shared Nextflow runner."""

from __future__ import annotations

import sys
from pathlib import Path

from clawbio.common.nfcore import execution as _execution

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))
sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules  # noqa: E402

purge_foreign_bare_modules("schemas")

from schemas import NEXTFLOW_ENV  # noqa: E402

_TIMEOUT_FIX = "Increase resources, inspect the pipeline, or retry with a smaller dataset."


def execute_nextflow(
    command: list[str], *, cwd: Path, output_dir: Path, timeout_seconds: int | None
) -> dict[str, object]:
    # Logs live at <output_dir>/logs/ for all three nf-core wrappers; the `logs`
    # directory is excluded from checksums.sha256.
    return _execution.execute_nextflow(
        command,
        cwd=cwd,
        output_dir=output_dir,
        timeout_seconds=timeout_seconds,
        timeout_fix=_TIMEOUT_FIX,
        env=NEXTFLOW_ENV,
    )
