"""nf-core/scrnaseq execution policy on top of the shared Nextflow runner."""

from __future__ import annotations

from pathlib import Path

from clawbio.common.nfcore import execution as _execution

_TIMEOUT_FIX = "Increase resources, inspect the pipeline, or retry with a smaller dataset."


def _cellbender_failure_hint(stdout_path: Path, stderr_path: Path) -> str:
    """Hint when CellBender background removal is the failing process.

    CELLBENDER_REMOVEBACKGROUND estimates ambient RNA from the droplet-count
    distribution and errors on very small or test datasets (a common symptom is
    ``IndexError: index -100 is out of bounds``). CellBender is optional, so a
    re-run with ``--skip-cellbender`` lets the rest of the pipeline produce its
    count matrices. Best-effort: empty string when CellBender is not the failing
    process (or the logs are unreadable).
    """
    blob = (
        _execution.read_log_tail(stdout_path) + "\n" + _execution.read_log_tail(stderr_path)
    ).lower()
    if "cellbender_removebackground" in blob and "error executing process" in blob:
        return (
            " The failing process is CellBender background removal, which estimates "
            "ambient RNA from the droplet distribution and does not work on very small "
            "or test datasets (a common symptom is 'IndexError: index -100 is out of "
            "bounds'). CellBender is optional — re-run with --skip-cellbender to let the "
            "rest of the pipeline produce its count matrices."
        )
    return ""


def execute_nextflow(
    command: list[str], *, cwd: Path, output_dir: Path, timeout_seconds: int | None
) -> dict[str, object]:
    return _execution.execute_nextflow(
        command,
        cwd=cwd,
        output_dir=output_dir,
        timeout_seconds=timeout_seconds,
        timeout_fix=_TIMEOUT_FIX,
        pipeline_failure_hint=_cellbender_failure_hint,
    )
