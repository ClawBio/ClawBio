"""How a wrapper reports a failure: result.json, a box on stdout, JSON on stderr."""

from __future__ import annotations

import json
import sys
import traceback
from pathlib import Path

from .artifacts import write_error_result
from .errors import SkillError


def _indent(text: str, width: int = 4) -> str:
    pad = " " * width
    return "\n".join(pad + line for line in text.splitlines())


def _emit(output_dir: Path, payload: dict[str, object], box: list[str]) -> int:
    write_error_result(output_dir, payload)
    print("\n".join(box), flush=True)
    # Machine-readable on stderr, available even when result.json cannot be written
    # (e.g. the output path is a file or the directory was refused).
    print(json.dumps(payload, indent=2, default=str), file=sys.stderr)
    # 1 for a reported failure; argparse keeps 2 for command-line usage errors.
    return 1


def report_skill_error(output_dir: Path, exc: SkillError, *, verbose: bool) -> int:
    box = [
        "",
        "================ SkillError ================",
        f"  stage:   {exc.stage}",
        f"  code:    {exc.error_code}",
        f"  message: {exc.message}",
        f"  fix:     {exc.fix}",
    ]
    if exc.details and verbose:
        box += ["  details:", _indent(json.dumps(exc.details, indent=2, default=str))]
    box.append("============================================")
    return _emit(output_dir, exc.to_dict(), box)


def report_unexpected_error(output_dir: Path, exc: BaseException, *, verbose: bool) -> int:
    """Report a bug: the traceback goes to result.json, and to stdout with --verbose."""
    trace = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
    payload = {
        "status": "error",
        "ok": False,
        "stage": "internal",
        "error_code": "UNEXPECTED_ERROR",
        "message": str(exc),
        "fix": "Report this as a bug. Include the traceback and command arguments.",
        "details": {"exception_type": type(exc).__name__, "traceback": trace},
    }
    box = [
        "",
        "================ Internal error ================",
        f"  type:    {type(exc).__name__}",
        f"  message: {exc}",
        _indent(trace) if verbose else "  Re-run with --verbose for the full traceback.",
        "================================================",
    ]
    return _emit(output_dir, payload, box)
