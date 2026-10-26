from __future__ import annotations

import csv
import os
import re
import sys
from collections import Counter
from pathlib import Path

from clawbio.common.nfcore.artifacts import write_csv_lf

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))


sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules

purge_foreign_bare_modules("errors", "schemas")

from errors import ErrorCode, SkillError
from schemas import SUPPORTED_SAMPLE_COLUMNS


_SAMPLE_WHITESPACE_RE = re.compile(r"\s+")
_BASE_OUTPUT_COLUMNS = ("sample", "fastq_1", "fastq_2", "strandedness", "genome_bam", "transcriptome_bam")
_UPSTREAM_OPTIONAL_SAMPLE_COLUMNS = {"seq_platform", "seq_center", "percent_mapped"}
_SUPPORTED_SAMPLE_COLUMNS = set(SUPPORTED_SAMPLE_COLUMNS) | _UPSTREAM_OPTIONAL_SAMPLE_COLUMNS
_FASTQ_COLUMNS = ("fastq_1", "fastq_2")
_BAM_COLUMNS = ("genome_bam", "transcriptome_bam")
_ResolvedPath = Path | str


def validate_and_normalize_samplesheet(
    input_path: Path, output_path: Path, *, skip_alignment: bool = False
) -> dict[str, object]:
    # nf-core/rnaseq 3.26.0 schema: input pattern ^\S+\.csv$ — reject non-.csv files early.
    if input_path.suffix.lower() != ".csv":
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message=f"The samplesheet must be a CSV file (nf-core schema: --input pattern ^\\.csv$); "
                    f"got '{input_path.name}'.",
            fix="Rename the file to end in .csv and ensure it uses comma-separated format.",
            details={"path": str(input_path), "suffix": input_path.suffix},
        )
    fieldnames, rows = _read_samplesheet(input_path)
    unknown_columns = _unknown_columns(fieldnames)
    output_columns = _build_output_columns(fieldnames)
    normalized_rows, sample_names, fastq_paths, bam_paths, strandedness_counts = _normalize_rows(
        rows,
        input_path=input_path,
        output_columns=output_columns,
        skip_alignment=skip_alignment,
    )
    _write_normalized_samplesheet(output_path, output_columns, normalized_rows)

    return {
        "normalized_path": output_path,
        "sample_count": len(sample_names),
        "sample_names": sample_names,
        "fastq_paths": fastq_paths,
        "bam_paths": bam_paths,
        "strandedness_counts": dict(strandedness_counts),
        "unknown_columns": unknown_columns,
    }


def _read_samplesheet(input_path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not input_path.exists():
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.MISSING_INPUT,
            message="Samplesheet was not found.",
            fix="Provide a valid --input path to a samplesheet.csv file.",
            details={"path": str(input_path)},
        )

    try:
        from clawbio.common.nfcore.artifacts import read_delimited_samplesheet

        fieldnames, rows = read_delimited_samplesheet(input_path)
    except (OSError, ValueError, csv.Error) as exc:
        raise SkillError(
            stage="validation", error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="Samplesheet could not be parsed without data loss.",
            fix="Use unique nonempty headers, correctly quoted UTF-8 records, and one value per column.",
            details={"path": str(input_path), "error": str(exc)},
        ) from exc
    if not rows:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="Samplesheet is empty.",
            fix="Add at least one sample row to the input CSV.",
            details={"path": str(input_path)},
        )
    return fieldnames, rows


def _unknown_columns(fieldnames: list[str]) -> list[str]:
    return [name for name in fieldnames if name not in _SUPPORTED_SAMPLE_COLUMNS]


def _build_output_columns(fieldnames: list[str]) -> list[str]:
    extras = [name for name in fieldnames if name not in _BASE_OUTPUT_COLUMNS]
    return list(dict.fromkeys([*_BASE_OUTPUT_COLUMNS, *extras]))


def _normalize_rows(
    rows: list[dict[str, str]],
    *,
    input_path: Path,
    output_columns: list[str],
    skip_alignment: bool = False,
) -> tuple[list[dict[str, str]], list[str], list[_ResolvedPath], list[_ResolvedPath], Counter[str]]:
    normalized_rows: list[dict[str, str]] = []
    sample_names: list[str] = []
    fastq_paths: list[_ResolvedPath] = []
    bam_paths: list[_ResolvedPath] = []
    strandedness_counts: Counter[str] = Counter()
    raw_samples_by_normalized: dict[str, str] = {}
    seen_fastq_rows: set[tuple[str, str, str]] = set()

    for line_number, row in enumerate(rows, start=2):
        normalized, sample, row_fastqs, row_bams = _normalize_row(
            row,
            input_path=input_path,
            output_columns=output_columns,
            line_number=line_number,
            skip_alignment=skip_alignment,
        )
        _reject_sample_name_collision(raw_samples_by_normalized, row, line_number, sample)
        _reject_duplicate_fastq_row(seen_fastq_rows, sample, row_fastqs, line_number)
        normalized_rows.append(normalized)
        if sample not in sample_names:
            sample_names.append(sample)
        fastq_paths.extend(row_fastqs)
        bam_paths.extend(row_bams)
        strandedness_counts[normalized["strandedness"]] += 1

    return normalized_rows, sample_names, fastq_paths, bam_paths, strandedness_counts


def _normalize_row(
    row: dict[str, str],
    *,
    input_path: Path,
    output_columns: list[str],
    line_number: int,
    skip_alignment: bool = False,
) -> tuple[dict[str, str], str, list[_ResolvedPath], list[_ResolvedPath]]:
    # Values nf-core itself validates (assets/schema_input.json: required columns,
    # patterns, the strandedness enum, file existence; checkSamplesAfterGrouping:
    # one strandedness per sample) pass through for nf-schema to name. The wrapper
    # only normalizes them and resolves paths, since Nextflow runs from --output.
    sample = _normalize_sample_name(str(row.get("sample", "")))
    strandedness = str(row.get("strandedness", "")).strip()
    bam_provided = _row_has_bam(row)
    fastq_provided = _row_has_fastq(row)
    if bam_provided and fastq_provided and not skip_alignment:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="A samplesheet row cannot mix FASTQ and BAM columns.",
            fix=(
                "Use either FASTQ columns (fastq_1[, fastq_2]) for fresh runs, OR BAM columns "
                "(genome_bam, transcriptome_bam) with --skip-alignment for reprocessing — never both."
            ),
            details={
                "line": line_number,
                "sample": sample,
                "fastq_columns": [c for c in _FASTQ_COLUMNS if str(row.get(c, "")).strip()],
                "bam_columns": [c for c in _BAM_COLUMNS if str(row.get(c, "")).strip()],
            },
        )
    resolved_fastqs = _resolve_columns(row, _FASTQ_COLUMNS, base_dir=input_path.parent)
    resolved_bams = _resolve_columns(row, _BAM_COLUMNS, base_dir=input_path.parent)
    normalized = {column: str(row.get(column, "")).strip() for column in output_columns}
    normalized.update(
        {
            "sample": sample,
            "fastq_1": _serialize_resolved_path(resolved_fastqs["fastq_1"]) if "fastq_1" in resolved_fastqs else "",
            "fastq_2": _serialize_resolved_path(resolved_fastqs["fastq_2"]) if "fastq_2" in resolved_fastqs else "",
            "strandedness": strandedness,
            "genome_bam": _serialize_resolved_path(resolved_bams["genome_bam"]) if "genome_bam" in resolved_bams else "",
            "transcriptome_bam": (
                _serialize_resolved_path(resolved_bams["transcriptome_bam"])
                if "transcriptome_bam" in resolved_bams
                else ""
            ),
        }
    )
    return normalized, sample, list(resolved_fastqs.values()), list(resolved_bams.values())


def _row_has_fastq(row: dict[str, str]) -> bool:
    return any(str(row.get(column, "")).strip() for column in _FASTQ_COLUMNS)


def _row_has_bam(row: dict[str, str]) -> bool:
    return any(str(row.get(column, "")).strip() for column in _BAM_COLUMNS)


def _normalize_sample_name(raw_sample: str) -> str:
    return _SAMPLE_WHITESPACE_RE.sub("_", raw_sample.strip())


def _resolve_columns(row: dict[str, str], columns: tuple[str, ...], *, base_dir: Path) -> dict[str, _ResolvedPath]:
    resolved: dict[str, _ResolvedPath] = {}
    for column in columns:
        raw_value = str(row.get(column, "")).strip()
        if not raw_value:
            continue
        resolved[column] = _resolve_path(raw_value, base_dir=base_dir)
    return resolved


def _resolve_path(value: str, *, base_dir: Path) -> _ResolvedPath:
    if _is_remote_uri(value):
        return value
    path = Path(os.path.expanduser(value))
    if not path.is_absolute():
        path = base_dir / path
    return path.resolve()


def _is_remote_uri(value: _ResolvedPath) -> bool:
    return "://" in str(value)


def _serialize_resolved_path(value: _ResolvedPath) -> str:
    if _is_remote_uri(value):
        return str(value)
    return Path(value).as_posix()


def _reject_sample_name_collision(
    raw_samples_by_normalized: dict[str, str],
    row: dict[str, str],
    line_number: int,
    sample: str,
) -> None:
    raw_sample = str(row.get("sample", "")).strip()
    previous_raw = raw_samples_by_normalized.setdefault(sample, raw_sample)
    if previous_raw != raw_sample:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="Two sample names collide after whitespace normalization.",
            fix="Rename samples so whitespace-normalized names are unique.",
            details={
                "line": line_number,
                "sample": raw_sample,
                "previous_sample": previous_raw,
                "normalized_sample": sample,
            },
        )


def _reject_duplicate_fastq_row(
    seen_fastq_rows: set[tuple[str, str, str]],
    sample: str,
    fastqs: list[_ResolvedPath],
    line_number: int,
) -> None:
    key = (
        sample,
        _serialize_resolved_path(fastqs[0]) if fastqs else "",
        _serialize_resolved_path(fastqs[1]) if len(fastqs) > 1 else "",
    )
    if key in seen_fastq_rows:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="Duplicate sample/FASTQ row found.",
            fix="Remove duplicate technical replicate rows from the samplesheet.",
            details={"line": line_number, "sample": sample},
        )
    seen_fastq_rows.add(key)


def _write_normalized_samplesheet(
    output_path: Path,
    fieldnames: list[str],
    rows: list[dict[str, str]],
) -> None:
    write_csv_lf(output_path, fieldnames, rows)
