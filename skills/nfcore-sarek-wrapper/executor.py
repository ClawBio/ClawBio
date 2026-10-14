"""nf-core/sarek execution policy on top of the shared Nextflow runner."""

from __future__ import annotations

from pathlib import Path

from clawbio.common.nfcore import execution as _execution

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
    )
