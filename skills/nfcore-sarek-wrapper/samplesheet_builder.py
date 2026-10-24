# ruff: noqa: E402
"""Parse the user's Sarek samplesheet and write the copy Nextflow reads.

nf-core/sarek 3.8.1 validates the samplesheet itself when Nextflow starts, before any
task runs: nf-schema checks ``assets/schema_input.json`` (required ``patient`` and
``sample``, the ``^\\S+$`` patterns, the ``sex``/``status`` enums, ``contamination``'s
type, every file's pattern and existence, ``dependentRequired`` and
``uniqueEntries``), and ``subworkflows/local/samplesheet_to_channel`` refuses rows whose
inputs do not fit ``--step``, duplicate patient-sample-status-lane rows, a sample
under two patients, two normals with exactly one tumour, and the ``--tools`` rules
on tumour/normal samples, sex and varlociraptor's contamination. None of that is
repeated here.

This module resolves paths (Nextflow runs from ``--output``), refuses what Sarek
accepts silently or does not check (rows mixing input modes, BAM/CRAM rows without
an index, recalibrate rows without a table, contamination outside [0, 1]) and
summarises the sheet leniently (patients, tumour/normal pairings) for the wrapper's
own policy and report. A value nf-schema judges is written through unchanged.
"""
from __future__ import annotations

import csv
import math
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, cast

from clawbio.common.nfcore.artifacts import write_csv_lf

import yaml

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))


sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules

purge_foreign_bare_modules("errors", "schemas")

from errors import ErrorCode, SkillError
from schemas import SUPPORTED_SAMPLE_COLUMNS, SUPPORTED_STEPS


_PATH_COLUMNS = ("fastq_1", "fastq_2", "spring_1", "spring_2", "bam", "bai", "cram", "crai", "table", "vcf")
_URI_PREFIXES = ("s3://", "gs://", "az://", "http://", "https://", "ftp://", "ftps://")
# Report keys for the resolved paths, in row order (the remote-input gate reads them).
_SUMMARY_PATHS = {
    "fastq_paths": ("fastq_1", "fastq_2"),
    "bam_paths": ("bam", "bai"),
    "cram_paths": ("cram", "crai"),
    "vcf_paths": ("vcf",),
    "spring_paths": ("spring_1", "spring_2"),
    "tables": ("table",),
}
_ALIGNED_STEPS = ("markduplicates", "prepare_recalibration", "recalibrate", "variant_calling")
# Columns whose raw value nf-schema matches against ^\S+$: written exactly as given.
_RAW_COLUMNS = ("patient", "sample", "lane")

_BASE_OUTPUT_COLUMNS = (
    "patient", "sample", "sex", "status", "lane",
    "fastq_1", "fastq_2", "spring_1", "spring_2",
    "bam", "bai", "cram", "crai", "table", "vcf", "variantcaller",
    "contamination",
)


# --- public API ------------------------------------------------------------

def validate_and_normalize_samplesheet(
    input_path: Path,
    output_path: Path,
    *,
    step: str = "mapping",
) -> dict[str, object]:
    """Write the normalised samplesheet and return the wrapper's summary of it."""
    _validate_step(step)
    fieldnames, rows = _read_samplesheet(input_path)
    output_columns = _build_output_columns(fieldnames)

    normalized_rows: list[dict[str, str]] = []
    paths: dict[str, list[Path | str]] = defaultdict(list)
    sample_names: list[str] = []
    patient_names: list[str] = []
    sex_counts: Counter = Counter()
    status_counts: Counter = Counter()
    rows_by_patient: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for line_number, row in enumerate(rows, start=2):
        normalized, resolved, entry = _normalize_row(
            row, input_path=input_path, output_columns=output_columns, line_number=line_number, step=step,
        )
        normalized_rows.append(normalized)
        for key, columns in _SUMMARY_PATHS.items():
            paths[key].extend(resolved[column] for column in columns if column in resolved)
        if entry["sample"] not in sample_names:
            sample_names.append(entry["sample"])
        if entry["patient"] not in patient_names:
            patient_names.append(entry["patient"])
        sex_counts[entry["sex"]] += 1
        status_counts[entry["status"]] += 1
        rows_by_patient[entry.pop("patient")].append(entry)

    pairings = _build_pairings(rows_by_patient)
    _write_normalized_samplesheet(output_path, output_columns, normalized_rows)
    return {
        "normalized_path": output_path,
        "step": step,
        "sample_count": len(sample_names),
        "patient_count": len(patient_names),
        "sample_names": sample_names,
        "patient_names": patient_names,
        **{key: paths[key] for key in _SUMMARY_PATHS},
        "unknown_columns": _unknown_columns(fieldnames),
        "sex_counts": dict(sex_counts),
        "status_counts": dict(status_counts),
        "pairings": pairings,
        "analysis_mode": _derive_global_analysis_mode(pairings),
        "rows_by_patient": dict(rows_by_patient),
    }


# --- header / file reading -------------------------------------------------

def _validate_step(step: str) -> None:
    if step not in SUPPORTED_STEPS:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_STEP,
            message=f"Unsupported --step value '{step}'.",
            fix=f"Choose one of: {', '.join(sorted(SUPPORTED_STEPS))}.",
            details={"step": step, "supported": sorted(SUPPORTED_STEPS)},
        )


def _read_samplesheet(input_path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not input_path.exists():
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.MISSING_INPUT,
            message="Samplesheet was not found.",
            fix="Provide a valid --input path to a samplesheet.csv file.",
            details={"path": str(input_path)},
        )
    if not input_path.is_file():
        raise SkillError(
            stage="validation", error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="Samplesheet path is not a regular file.",
            fix="Provide a samplesheet file instead of a directory.", details={"path": str(input_path)},
        )
    if input_path.stat().st_size == 0:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="Samplesheet is empty (zero bytes).",
            fix="Add a header row and at least one sample row to the input CSV.",
            details={"path": str(input_path)},
        )

    suffix = input_path.suffix.lower()
    if suffix in {".csv", ".tsv"}:
        delimiter = "\t" if suffix == ".tsv" else ","
        try:
            from clawbio.common.nfcore.artifacts import read_delimited_samplesheet

            fieldnames, rows = read_delimited_samplesheet(input_path, delimiter=delimiter)
        except (OSError, ValueError, csv.Error) as exc:
            raise SkillError(
                stage="validation", error_code=ErrorCode.INVALID_SAMPLESHEET,
                message="Samplesheet could not be parsed without data loss.",
                fix="Use unique nonempty headers, correctly quoted UTF-8 records, and one value per column.",
                details={"path": str(input_path), "error": str(exc)},
            ) from exc

    elif suffix in {".json", ".yaml", ".yml"}:
        try:
            from clawbio.common.nfcore.parameters import load_unique_json, load_unique_yaml

            text = input_path.read_text(encoding="utf-8-sig")
            payload = load_unique_json(text) if suffix == ".json" else load_unique_yaml(text)
        except (OSError, ValueError, yaml.YAMLError) as exc:
            raise SkillError(
                stage="validation", error_code=ErrorCode.INVALID_SAMPLESHEET,
                message="Structured samplesheet could not be parsed.",
                fix="Use valid UTF-8 JSON/YAML with unique mapping keys.",
                details={"path": str(input_path), "error": str(exc)},
            ) from exc
        rows_obj = _extract_structured_rows(payload, input_path=input_path)
        rows = [{str(k): "" if v is None else str(v) for k, v in row.items()} for row in rows_obj]
        fieldnames = list(dict.fromkeys(key for row in rows for key in row))
    else:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message=f"Unsupported samplesheet format '{input_path.suffix}'.",
            fix="Use one of the Sarek-supported formats: .csv, .tsv, .json, .yaml, or .yml.",
            details={"path": str(input_path), "suffix": input_path.suffix},
        )

    if not fieldnames:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="Samplesheet header could not be parsed.",
            fix="Ensure the first line of the CSV is a comma-separated header.",
            details={"path": str(input_path)},
        )
    if not rows:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="Samplesheet is empty (no data rows).",
            fix="Add at least one sample row beneath the header.",
            details={"path": str(input_path)},
        )
    return cast("list[str]", fieldnames), cast("list[dict[str, str]]", rows)


def _extract_structured_rows(payload: Any, *, input_path: Path) -> list[dict[str, Any]]:
    if isinstance(payload, list):
        rows = payload
    elif isinstance(payload, dict):
        rows_val: Any = payload.get("samples") or payload.get("rows") or payload.get("data")
        rows = rows_val
    else:
        rows = None
    if not isinstance(rows, list) or not all(isinstance(row, dict) for row in rows):
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="Structured samplesheet must contain a list of row objects.",
            fix="Use a top-level list, or a mapping with a list under `samples`, `rows`, or `data`.",
            details={"path": str(input_path)},
        )
    return rows


def _unknown_columns(fieldnames: list[str]) -> list[str]:
    return [name for name in fieldnames if name not in SUPPORTED_SAMPLE_COLUMNS]


def _build_output_columns(fieldnames: list[str]) -> list[str]:
    extras = [name for name in fieldnames if name not in _BASE_OUTPUT_COLUMNS]
    return list(dict.fromkeys([*_BASE_OUTPUT_COLUMNS, *extras]))


# --- per-row normalisation -------------------------------------------------

def _normalize_row(
    row: dict[str, str],
    *,
    input_path: Path,
    output_columns: list[str],
    line_number: int,
    step: str,
) -> tuple[dict[str, str], dict[str, Path | str], dict[str, Any]]:
    sample = str(row.get("sample", ""))
    resolved = {
        column: _resolve_path(value, base_dir=input_path.parent)
        for column in _PATH_COLUMNS
        if (value := str(row.get(column, "")).strip())
    }
    _check_inputs_for_step(resolved, step=step, line_number=line_number, sample=sample)
    status = _lenient_status(row)
    contamination = _contamination(row, line_number)
    normalized = {column: str(row.get(column, "")).strip() for column in output_columns}
    normalized.update({column: str(row.get(column, "")) for column in _RAW_COLUMNS})
    normalized["sex"] = normalized["sex"] or "NA"
    normalized["status"] = str(status)
    if isinstance(contamination, float):
        normalized["contamination"] = f"{contamination:g}"
    for column in _PATH_COLUMNS:
        normalized[column] = _path_to_text(resolved[column]) if column in resolved else ""
    entry = {
        "patient": normalized["patient"],
        "sample": sample,
        "status": status,
        "sex": normalized["sex"],
        "contamination": contamination,
        "line": line_number,
    }
    return normalized, resolved, entry


def _lenient_status(row: dict[str, str]) -> int | str:
    """0/1 for the pairing summary; anything else is kept for nf-schema to refuse."""
    raw = str(row.get("status", "")).strip()
    if not raw:
        return 0
    try:
        return int(raw)
    except ValueError:
        return raw


def _contamination(row: dict[str, str], line_number: int) -> float | str | None:
    """A numeric contamination must be a fraction: Sarek's schema bounds only its type."""
    raw = str(row.get("contamination", "")).strip()
    if not raw:
        return None
    try:
        value = float(raw)
    except ValueError:
        return raw  # not a number: nf-schema refuses it by type
    if not (math.isfinite(value) and 0.0 <= value <= 1.0):
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message=f"`contamination` must be in [0, 1]; got {raw}.",
            fix="Provide a fractional contamination value (0 = clean, 1 = fully contaminated).",
            details={"line": line_number, "contamination": raw},
        )
    return value


def _check_inputs_for_step(resolved: dict[str, Path | str], *, step: str, line_number: int, sample: str) -> None:
    """Refuse what nf-core/sarek accepts without a word: it uses the first input mode
    a mixed row matches, passes a null index for a BAM/CRAM row without one unless a
    lane is set, and routes a recalibrate row without a table as if it had one."""
    if step == "mapping":
        modes = [mode for mode, columns in (("fastq", ("fastq_1", "fastq_2")), ("spring", ("spring_1", "spring_2")), ("bam", ("bam",)))
                 if any(column in resolved for column in columns)]
        if len(modes) > 1:
            raise SkillError(
                stage="validation",
                error_code=ErrorCode.INVALID_SAMPLESHEET,
                message="A samplesheet row cannot mix FASTQ, Spring and BAM columns.",
                fix="Use exactly one input mode per row: (fastq_1, fastq_2) OR spring_1[/spring_2] OR bam.",
                details={"line": line_number, "sample": sample, "modes_detected": modes},
            )
        return
    if step not in _ALIGNED_STEPS:
        return
    if "bam" in resolved and "cram" in resolved:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="A samplesheet row cannot mix BAM and CRAM columns.",
            fix="Provide either bam+bai OR cram+crai per row — never both.",
            details={"line": line_number, "sample": sample},
        )
    for data, index in (("bam", "bai"), ("cram", "crai")):
        if data in resolved and index not in resolved:
            raise SkillError(
                stage="validation",
                error_code=ErrorCode.INVALID_SAMPLESHEET,
                message=f"`{index}` is required alongside `{data}` at --step {step}.",
                fix=f"Provide both {data} and {index} for {data.upper()}-mode rows.",
                details={"line": line_number, "sample": sample, "column": index},
            )
    if step == "recalibrate" and "table" not in resolved:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="`table` (BQSR recalibration table) is required at --step recalibrate.",
            fix="Provide a path to the .table file emitted by the prepare_recalibration step.",
            details={"line": line_number, "sample": sample},
        )


# --- low-level path helpers ------------------------------------------------

def _is_uri(value: str) -> bool:
    return any(value.startswith(prefix) for prefix in _URI_PREFIXES)


def _resolve_path(value: str, *, base_dir: Path) -> Path | str:
    if _is_uri(value):
        return value
    path = Path(os.path.expanduser(value))
    if not path.is_absolute():
        path = base_dir / path
    return path.resolve()


def _path_to_text(path: Path | str | None) -> str:
    if path is None:
        return ""
    return path.as_posix() if isinstance(path, Path) else path


# --- tumor / normal pairing ------------------------------------------------

def _build_pairings(rows_by_patient: dict[str, list[dict[str, Any]]]) -> list[dict[str, Any]]:
    """Tumour/normal combinations per patient, as workflows/sarek/main.nf crosses them.

    Sarek itself refuses two normals with exactly one tumour (samplesheet_to_channel).
    """
    pairings: list[dict[str, Any]] = []
    for patient, entries in rows_by_patient.items():
        # Collapse per sample/status pair: upstream deduplicates
        # meta.subMap('sample', 'status'), not sample alone.
        samples_by_status: dict[int | str, list[str]] = defaultdict(list)
        seen: set[tuple[str, int | str]] = set()
        for entry in entries:
            sample = str(entry["sample"])
            status = entry["status"]
            key = (sample, status)
            if key in seen:
                continue
            seen.add(key)
            samples_by_status[status].append(sample)

        normals = samples_by_status.get(0, [])
        tumors = samples_by_status.get(1, [])

        if normals and tumors:
            mode = "somatic_paired"
        elif tumors and not normals:
            mode = "tumor_only"
        else:
            # Germline, or no 0/1 status at all (nf-schema refuses other statuses).
            mode = "germline"

        if mode == "somatic_paired":
            # workflows/sarek/main.nf uses cross() to make an output key for
            # every tumor/normal combination, including recurrent tumors.
            for normal in normals:
                for tumor in tumors:
                    pairings.append({
                        "patient": patient,
                        "normal": normal,
                        "tumor": tumor,
                        "mode": mode,
                    })
        else:
            pairings.append({
                "patient": patient,
                "normal": normals[0] if normals else None,
                "tumors": tumors,
                "mode": mode,
            })
    return pairings


def _derive_global_analysis_mode(pairings: list[dict[str, Any]]) -> str:
    if not pairings:
        return "germline"
    modes = {p["mode"] for p in pairings}
    if len(modes) == 1:
        return next(iter(modes))
    return "mixed"


# --- output writing --------------------------------------------------------

def _write_normalized_samplesheet(
    output_path: Path,
    fieldnames: list[str],
    rows: list[dict[str, str]],
) -> None:
    write_csv_lf(output_path, fieldnames, rows)
