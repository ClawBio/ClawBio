"""Artifact mechanics; callers select payloads, names, and file policies."""
import csv
import json
import re
import shutil
from typing import Callable
from pathlib import Path
from clawbio.common.textio import write_text_lf


def write_json_lf(path: Path, payload: object, *, trailing_newline: bool, sort_keys: bool = False, default: Callable | None = None) -> None:
    write_text_lf(path, json.dumps(payload, indent=2, sort_keys=sort_keys, default=default) + ("\n" if trailing_newline else ""))


def write_csv_lf(output_path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def safe_config_basename(name: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._")
    return safe or "nextflow.config"


def stage_config_files(config_paths: list[Path], *, destination_dir: Path, filenames: list[str]) -> list[Path]:
    if len(config_paths) != len(filenames):
        raise ValueError("Each config source needs one destination filename.")
    if config_paths:
        destination_dir.mkdir(parents=True, exist_ok=True)
    copied = []
    for source, filename in zip(config_paths, filenames):
        destination = destination_dir / filename
        shutil.copyfile(source, destination)
        copied.append(destination)
    return copied


def write_json_snapshots(directory: Path, payloads: dict[str, object]) -> None:
    for filename, payload in payloads.items():
        write_json_lf(directory / filename, payload, trailing_newline=False)
