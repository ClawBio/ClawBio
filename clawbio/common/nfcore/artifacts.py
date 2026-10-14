"""Artifact mechanics; callers select payloads, names, and file policies."""

import csv
import json
import re
from collections.abc import Callable
from pathlib import Path

from clawbio.common.textio import write_text_lf


def write_json_lf(
    path: Path,
    payload: object,
    *,
    trailing_newline: bool,
    sort_keys: bool = False,
    default: Callable | None = None,
) -> None:
    write_text_lf(
        path,
        json.dumps(payload, indent=2, sort_keys=sort_keys, default=default)
        + ("\n" if trailing_newline else ""),
    )


def write_csv_lf(
    output_path: Path, fieldnames: list[str], rows: list[dict[str, str]]
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def safe_config_basename(name: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._")
    return safe or "nextflow.config"


def write_json_snapshots(directory: Path, payloads: dict[str, object]) -> None:
    for filename, payload in payloads.items():
        write_json_lf(directory / filename, payload, trailing_newline=False)


# Never write into a directory the wrapper refused to use.
_REFUSED_OUTPUT_DIR_CODES = frozenset({"OUTPUT_DIR_NOT_EMPTY", "OUTPUT_DIR_NOT_WRITABLE"})
# A refused --resume leaves the earlier run untouched, so its result.json stands:
# overwriting it would record a failure over outputs that are still on disk.
_REFUSED_RESUME_CODES = frozenset({"INVALID_RESUME_STATE", "RESUME_DRIFT"})


def write_error_result(output_dir: Path, payload: dict[str, object]) -> None:
    """Persist a structured error as <output_dir>/result.json when that is safe.

    Best-effort: the caller always echoes the payload on stderr, so an unwritable
    directory must not mask the original error.
    """
    code = payload.get("error_code")
    if code in _REFUSED_OUTPUT_DIR_CODES:
        return
    result_path = output_dir / "result.json"
    if code in _REFUSED_RESUME_CODES and result_path.exists():
        return
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        write_json_lf(result_path, payload, trailing_newline=False, default=str)
    except OSError:
        return
