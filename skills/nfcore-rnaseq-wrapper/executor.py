"""nf-core/rnaseq execution policy on top of the shared Nextflow runner."""

from __future__ import annotations

from pathlib import Path

from clawbio.common.nfcore import execution as _execution

_TIMEOUT_FIX = (
    "Raise --timeout-hours for large cohorts, inspect the pipeline, or retry "
    "with a smaller dataset. Nextflow was terminated, but containers started by "
    "the Docker/Singularity daemon are not in its process group and may still be "
    "running — check for and remove any leftover containers (e.g. `docker ps`)."
)


def execute_nextflow(
    command: list[str], *, cwd: Path, output_dir: Path, timeout_seconds: int | None
) -> dict[str, object]:
    return _execution.execute_nextflow(
        command,
        cwd=cwd,
        output_dir=output_dir,
        timeout_seconds=timeout_seconds,
        timeout_fix=_TIMEOUT_FIX,
    )
