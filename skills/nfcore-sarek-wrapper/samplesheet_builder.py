# ruff: noqa: E402
from __future__ import annotations

import csv
import os
import re
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
from schemas import (
    REQUIRED_SAMPLE_COLUMNS_ANY_STEP,
    SUPPORTED_SAMPLE_COLUMNS,
    SUPPORTED_SEX,
    SUPPORTED_STATUS,
    SUPPORTED_STEPS,
)


# --- regexes and column groups ---------------------------------------------

_SAMPLE_NAME_RE = re.compile(r"^\S+$")
_FASTQ_BASENAME_RE = re.compile(r"^[^\s/]+\.f(?:ast)?q\.gz$")
_SPRING_BASENAME_RE = re.compile(r"^[^\s/]+\.f(?:ast)?q\.gz\.spring$")
_BAM_BASENAME_RE = re.compile(r"^[^\s/]+\.bam$", re.IGNORECASE)
_BAI_BASENAME_RE = re.compile(r"^[^\s/]+\.bai$", re.IGNORECASE)
_CRAM_BASENAME_RE = re.compile(r"^[^\s/]+\.cram$", re.IGNORECASE)
_CRAI_BASENAME_RE = re.compile(r"^[^\s/]+\.crai$", re.IGNORECASE)
# Matches the official schema_input.json pattern `^\S+\.vcf(\.gz)?$`: both plain
# .vcf and bgzipped .vcf.gz are accepted (bgzipped is recommended in the docs).
_VCF_BASENAME_RE = re.compile(r"^[^\s/]+\.vcf(\.gz)?$", re.IGNORECASE)
_TABLE_BASENAME_RE = re.compile(r"^[^\s/]+\.table$", re.IGNORECASE)

_FASTQ_COLUMNS = ("fastq_1", "fastq_2")
_SPRING_COLUMNS = ("spring_1", "spring_2")
_BAM_COLUMNS = ("bam", "bai")
_CRAM_COLUMNS = ("cram", "crai")
_URI_PREFIXES = ("s3://", "gs://", "az://", "http://", "https://", "ftp://", "ftps://")

_BASE_OUTPUT_COLUMNS = (
    "patient", "sample", "sex", "status", "lane",
    "fastq_1", "fastq_2", "spring_1", "spring_2",
    "bam", "bai", "cram", "crai", "table", "vcf", "variantcaller",
    "contamination",
)

# Step-aware required-column logic.  These are the columns the wrapper expects
# to find in the CSV header for each step; per-row content rules are enforced
# in _normalize_row_for_step.
_HEADER_REQUIREMENTS: dict[str, tuple[str, ...]] = {
    "mapping":              ("patient", "sample", "lane"),
    "markduplicates":       ("patient", "sample"),
    "prepare_recalibration":("patient", "sample"),
    "recalibrate":          ("patient", "sample", "table"),
    "variant_calling":      ("patient", "sample"),
    "annotate":             ("patient", "sample", "vcf"),
}


# --- public API ------------------------------------------------------------

def validate_and_normalize_samplesheet(
    input_path: Path,
    output_path: Path,
    *,
    step: str = "mapping",
    tools: list[str] | None = None,
) -> dict[str, object]:
    """Validate the user's Sarek samplesheet and write a normalised copy.

    See the module docstring (and SKILL.md) for the full contract.  The
    function is step-aware: required columns and per-row input modes change
    according to ``step``.
    """
    tools = list(tools or [])
    _validate_step(step)

    fieldnames, rows = _read_samplesheet(input_path)
    _validate_required_columns(fieldnames, step=step)
    unknown_columns = _unknown_columns(fieldnames)
    output_columns = _build_output_columns(fieldnames)

    (
        normalized_rows,
        sample_names,
        patient_names,
        fastq_paths,
        bam_paths,
        cram_paths,
        vcf_paths,
        spring_paths,
        tables,
        sex_counts,
        status_counts,
        rows_by_patient,
    ) = _normalize_rows(
        rows,
        input_path=input_path,
        output_columns=output_columns,
        step=step,
    )

    pairings = _build_pairings(rows_by_patient)
    analysis_mode = _derive_global_analysis_mode(pairings)

    _enforce_varlociraptor_contamination(
        tools=tools,
        rows_by_patient=rows_by_patient,
        analysis_mode=analysis_mode,
    )

    _write_normalized_samplesheet(output_path, output_columns, normalized_rows)

    return {
        "normalized_path": output_path,
        "step": step,
        "sample_count": len(sample_names),
        "patient_count": len(patient_names),
        "sample_names": sample_names,
        "patient_names": patient_names,
        "fastq_paths": fastq_paths,
        "bam_paths": bam_paths,
        "cram_paths": cram_paths,
        "vcf_paths": vcf_paths,
        "spring_paths": spring_paths,
        "tables": tables,
        "unknown_columns": unknown_columns,
        "sex_counts": dict(sex_counts),
        "status_counts": {int(k): v for k, v in status_counts.items()},
        "pairings": pairings,
        "analysis_mode": analysis_mode,
        "rows_by_patient": rows_by_patient,
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


def _validate_required_columns(fieldnames: list[str], *, step: str) -> None:
    missing_any_step = [c for c in REQUIRED_SAMPLE_COLUMNS_ANY_STEP if c not in fieldnames]
    step_required = _HEADER_REQUIREMENTS.get(step, REQUIRED_SAMPLE_COLUMNS_ANY_STEP)
    missing_step = [c for c in step_required if c not in fieldnames]
    missing = list(dict.fromkeys([*missing_any_step, *missing_step]))
    if missing:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message=f"Samplesheet is missing required columns for step '{step}'.",
            fix=(
                f"For --step {step}, the samplesheet must contain: "
                f"{', '.join(step_required)}.  Add the missing columns to the header."
            ),
            details={"step": step, "missing_columns": missing, "have_columns": fieldnames},
        )


def _unknown_columns(fieldnames: list[str]) -> list[str]:
    return [name for name in fieldnames if name not in SUPPORTED_SAMPLE_COLUMNS]


def _build_output_columns(fieldnames: list[str]) -> list[str]:
    extras = [name for name in fieldnames if name not in _BASE_OUTPUT_COLUMNS]
    return list(dict.fromkeys([*_BASE_OUTPUT_COLUMNS, *extras]))


# --- per-row normalisation -------------------------------------------------

def _normalize_rows(
    rows: list[dict[str, str]],
    *,
    input_path: Path,
    output_columns: list[str],
    step: str,
) -> tuple[
    list[dict[str, str]],
    list[str],
    list[str],
    list[Path | str],
    list[Path | str],
    list[Path | str],
    list[Path | str],
    list[Path | str],
    list[Path | str],
    Counter,
    Counter,
    dict[str, list[dict[str, Any]]],
]:
    normalized_rows: list[dict[str, str]] = []
    sample_names: list[str] = []
    patient_names: list[str] = []
    fastq_paths: list[Path | str] = []
    bam_paths: list[Path | str] = []
    cram_paths: list[Path | str] = []
    vcf_paths: list[Path | str] = []
    spring_paths: list[Path | str] = []
    tables: list[Path | str] = []
    sex_counts: Counter = Counter()
    status_counts: Counter = Counter()

    sample_patients: dict[str, str] = {}
    seen_patient_sample_status_lane: set[tuple[str, str, int, str]] = set()
    rows_by_patient: dict[str, list[dict[str, Any]]] = defaultdict(list)

    for line_number, row in enumerate(rows, start=2):
        normalized, summary = _normalize_row_for_step(
            row,
            input_path=input_path,
            output_columns=output_columns,
            line_number=line_number,
            step=step,
        )
        sample = summary["sample"]
        patient = summary["patient"]
        sex = summary["sex"]
        status = summary["status"]
        contamination = summary["contamination"]
        lane = summary["lane"]

        _reject_sample_multiple_patients(sample_patients, patient, sample, line_number)
        # Sarek validates the patient-sample-status-lane uniqueness key BEFORE
        # branching by step (subworkflows/local/samplesheet_to_channel/main.nf),
        # so it applies to every step. For non-mapping steps `lane` is empty,
        # which correctly collapses the key to patient-sample-status.
        _reject_duplicate_patient_sample_status_lane(
            seen_patient_sample_status_lane, patient, sample, status, lane, line_number
        )

        normalized_rows.append(normalized)
        if sample not in sample_names:
            sample_names.append(sample)
        if patient not in patient_names:
            patient_names.append(patient)
        fastq_paths.extend(cast("list[Path | str]", summary["fastq_paths"]))
        spring_paths.extend(cast("list[Path | str]", summary["spring_paths"]))
        bam_paths.extend(cast("list[Path | str]", summary["bam_paths"]))
        cram_paths.extend(cast("list[Path | str]", summary["cram_paths"]))
        vcf_paths.extend(cast("list[Path | str]", summary["vcf_paths"]))
        tables.extend(cast("list[Path | str]", summary["tables"]))
        sex_counts[sex] += 1
        status_counts[status] += 1
        rows_by_patient[patient].append({
            "sample": sample,
            "status": status,
            "sex": sex,
            "contamination": contamination,
            "line": line_number,
        })

    return (
        normalized_rows,
        sample_names,
        patient_names,
        fastq_paths,
        bam_paths,
        cram_paths,
        vcf_paths,
        spring_paths,
        tables,
        sex_counts,
        status_counts,
        dict(rows_by_patient),
    )


def _normalize_row_for_step(
    row: dict[str, str],
    *,
    input_path: Path,
    output_columns: list[str],
    line_number: int,
    step: str,
) -> tuple[dict[str, str], dict[str, Any]]:
    patient = _validate_patient(row, line_number)
    sample = _validate_sample_name(row, line_number)
    sex = _validate_sex(row, line_number)
    status = _validate_status(row, line_number)
    contamination = _validate_contamination(row, line_number)
    lane_raw = str(row.get("lane", ""))
    lane = lane_raw.strip()
    # Schema (assets/schema_input.json) constrains lane to `^\S+$`: when present,
    # the RAW value must contain no whitespace (surrounding or internal).
    if lane and not _SAMPLE_NAME_RE.match(lane_raw):
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message=f"`lane` value '{lane_raw}' must not contain whitespace (schema pattern ^\\S+$).",
            fix="Use a whitespace-free lane identifier (e.g. lane_1).",
            details={"line": line_number, "sample": sample, "lane": lane_raw},
        )

    resolved: dict[str, Path | str] = {}
    fastq_paths: list[Path | str] = []
    spring_paths: list[Path | str] = []
    bam_paths: list[Path | str] = []
    cram_paths: list[Path | str] = []
    vcf_paths: list[Path | str] = []
    tables: list[Path | str] = []

    if step == "mapping":
        if not lane:
            raise SkillError(
                stage="validation",
                error_code=ErrorCode.INVALID_SAMPLESHEET,
                message="`lane` is required for every row at --step mapping.",
                fix="Add a non-empty lane identifier (e.g. lane_1) for each row.",
                details={"line": line_number, "sample": sample},
            )
        _resolve_mapping_inputs(
            row,
            input_path=input_path,
            line_number=line_number,
            sample=sample,
            resolved=resolved,
            fastq_paths=fastq_paths,
            spring_paths=spring_paths,
            bam_paths=bam_paths,
        )

    elif step in ("markduplicates", "prepare_recalibration", "variant_calling"):
        _resolve_aligned_inputs(
            row,
            input_path=input_path,
            line_number=line_number,
            sample=sample,
            resolved=resolved,
            bam_paths=bam_paths,
            cram_paths=cram_paths,
        )

    elif step == "recalibrate":
        _resolve_aligned_inputs(
            row,
            input_path=input_path,
            line_number=line_number,
            sample=sample,
            resolved=resolved,
            bam_paths=bam_paths,
            cram_paths=cram_paths,
        )
        table_path = _resolve_table(row, input_path=input_path, line_number=line_number, sample=sample)
        resolved["table"] = table_path
        tables.append(table_path)

    elif step == "annotate":
        vcf_path = _resolve_vcf(row, input_path=input_path, line_number=line_number, sample=sample)
        resolved["vcf"] = vcf_path
        vcf_paths.append(vcf_path)

    normalized = {column: str(row.get(column, "")).strip() for column in output_columns}
    normalized["patient"] = patient
    normalized["sample"] = sample
    normalized["sex"] = sex
    normalized["status"] = str(status)
    normalized["lane"] = lane
    normalized["contamination"] = "" if contamination is None else f"{contamination:g}"
    for col in ("fastq_1", "fastq_2", "spring_1", "spring_2",
                "bam", "bai", "cram", "crai", "table", "vcf"):
        path = resolved.get(col)
        normalized[col] = _path_to_text(path) if path is not None else ""

    summary = {
        "patient": patient,
        "sample": sample,
        "sex": sex,
        "status": status,
        "contamination": contamination,
        "lane": lane,
        "fastq_paths": fastq_paths,
        "spring_paths": spring_paths,
        "bam_paths": bam_paths,
        "cram_paths": cram_paths,
        "vcf_paths": vcf_paths,
        "tables": tables,
    }
    return normalized, summary


# --- patient / sample / sex / status / contamination -----------------------

def _validate_patient(row: dict[str, str], line_number: int) -> str:
    # Validate the RAW value against the schema pattern ^\S+$ so surrounding
    # whitespace is rejected (not silently stripped), matching nf-core/sarek.
    raw = str(row.get("patient", ""))
    if not raw.strip():
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="`patient` is required and must be non-empty.",
            fix="Provide a patient identifier (e.g. P1) for every row.",
            details={"line": line_number},
        )
    if not _SAMPLE_NAME_RE.match(raw):
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="`patient` must not contain whitespace (schema pattern ^\\S+$).",
            fix="Use a single whitespace-free token (e.g. P_1) for patient identifiers.",
            details={"line": line_number, "patient": raw},
        )
    return raw


def _validate_sample_name(row: dict[str, str], line_number: int) -> str:
    # Validate the RAW value: surrounding whitespace is rejected, not stripped.
    raw_sample = str(row.get("sample", ""))
    if not raw_sample.strip() or not _SAMPLE_NAME_RE.match(raw_sample):
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="Sample names must be non-empty and must not contain whitespace (schema pattern ^\\S+$).",
            fix="Use a single whitespace-free sample token (for example Sample_One).",
            details={"line": line_number, "sample": raw_sample},
        )
    return raw_sample


def _validate_sex(row: dict[str, str], line_number: int) -> str:
    raw = str(row.get("sex", "")).strip()
    if not raw:
        return "NA"
    if raw not in SUPPORTED_SEX:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SEX,
            message=f"Invalid value for `sex`: '{raw}'.",
            fix=f"Use one of: {', '.join(sorted(SUPPORTED_SEX))} (or leave empty for NA).",
            details={"line": line_number, "sex": raw, "supported": sorted(SUPPORTED_SEX)},
        )
    return raw


def _validate_status(row: dict[str, str], line_number: int) -> int:
    raw = str(row.get("status", "")).strip()
    if not raw:
        return 0
    try:
        value = int(raw)
    except ValueError:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_STATUS,
            message=f"Invalid value for `status`: '{raw}'.",
            fix="Use 0 (normal) or 1 (tumor).",
            details={"line": line_number, "status": raw},
        )
    if value not in SUPPORTED_STATUS:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_STATUS,
            message=f"Invalid value for `status`: {value}.",
            fix="Use 0 (normal) or 1 (tumor).",
            details={"line": line_number, "status": value, "supported": sorted(SUPPORTED_STATUS)},
        )
    return value


def _validate_contamination(row: dict[str, str], line_number: int) -> float | None:
    raw = str(row.get("contamination", "")).strip()
    if not raw:
        return None
    try:
        value = float(raw)
    except ValueError:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message=f"Invalid value for `contamination`: '{raw}'.",
            fix="Provide a float between 0 and 1 (e.g. 0.05).",
            details={"line": line_number, "contamination": raw},
        )
    if not (0.0 <= value <= 1.0):
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message=f"`contamination` must be in [0, 1]; got {value}.",
            fix="Provide a fractional contamination value (0 = clean, 1 = fully contaminated).",
            details={"line": line_number, "contamination": value},
        )
    return value


# --- mapping-step input resolution -----------------------------------------

def _resolve_mapping_inputs(
    row: dict[str, str],
    *,
    input_path: Path,
    line_number: int,
    sample: str,
    resolved: dict[str, Path | str],
    fastq_paths: list[Path | str],
    spring_paths: list[Path | str],
    bam_paths: list[Path | str],
) -> None:
    has_fastq = any(str(row.get(c, "")).strip() for c in _FASTQ_COLUMNS)
    has_spring = any(str(row.get(c, "")).strip() for c in _SPRING_COLUMNS)
    has_bam = bool(str(row.get("bam", "")).strip())

    chosen = [name for name, flag in
              (("fastq", has_fastq), ("spring", has_spring), ("bam", has_bam)) if flag]
    if len(chosen) == 0:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="A mapping-step row must provide FASTQ, Spring, or BAM inputs.",
            fix="Populate exactly one of: (fastq_1, fastq_2), spring_1[/spring_2], or bam.",
            details={"line": line_number, "sample": sample},
        )
    if len(chosen) > 1:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="A samplesheet row cannot mix FASTQ, Spring and BAM columns.",
            fix="Use exactly one input mode per row: (fastq_1, fastq_2) OR spring_1[/spring_2] OR bam.",
            details={"line": line_number, "sample": sample, "modes_detected": chosen},
        )

    mode = chosen[0]
    if mode == "fastq":
        fq1_raw = str(row.get("fastq_1", "")).strip()
        if not fq1_raw:
            raise SkillError(
                stage="validation",
                error_code=ErrorCode.INVALID_SAMPLESHEET,
                message="fastq_1 is required when FASTQ inputs are used.",
                fix="Provide a FASTQ path in the fastq_1 column.",
                details={"line": line_number, "sample": sample},
            )
        fq1 = _resolve_path(fq1_raw, base_dir=input_path.parent)
        _validate_fastq_path(fq1, column="fastq_1", line_number=line_number)
        resolved["fastq_1"] = fq1
        fastq_paths.append(fq1)
        fq2_raw = str(row.get("fastq_2", "")).strip()
        if not fq2_raw:
            raise SkillError(
                stage="validation",
                error_code=ErrorCode.INVALID_SAMPLESHEET,
                message="fastq_2 is required when FASTQ inputs are used.",
                fix="Provide paired FASTQ paths in fastq_1 and fastq_2, or use spring_1 for one-file Spring input.",
                details={"line": line_number, "sample": sample},
            )
        fq2 = _resolve_path(fq2_raw, base_dir=input_path.parent)
        _validate_fastq_path(fq2, column="fastq_2", line_number=line_number)
        resolved["fastq_2"] = fq2
        fastq_paths.append(fq2)
    elif mode == "spring":
        # Schema (schema_input.json) declares spring_2 -> spring_1: spring_2 may
        # only appear together with spring_1 (one-file or paired Spring input).
        if str(row.get("spring_2", "")).strip() and not str(row.get("spring_1", "")).strip():
            raise SkillError(
                stage="validation",
                error_code=ErrorCode.INVALID_SAMPLESHEET,
                message="spring_2 requires spring_1 (spring_2 cannot be used on its own).",
                fix="Provide spring_1 (single-file Spring) or both spring_1 and spring_2.",
                details={"line": line_number, "sample": sample},
            )
        for col in _SPRING_COLUMNS:
            raw = str(row.get(col, "")).strip()
            if not raw:
                continue
            path = _resolve_path(raw, base_dir=input_path.parent)
            _validate_spring_path(path, column=col, line_number=line_number)
            resolved[col] = path
            spring_paths.append(path)
    else:  # bam (uBAM)
        bam_raw = str(row.get("bam", "")).strip()
        bam = _resolve_path(bam_raw, base_dir=input_path.parent)
        _validate_bam_path(bam, column="bam", line_number=line_number)
        resolved["bam"] = bam
        bam_paths.append(bam)
        # Sarek propagates the optional bai for a mapping-step BAM ([meta, bam, bai]).
        bai_raw = str(row.get("bai", "")).strip()
        if bai_raw:
            bai = _resolve_path(bai_raw, base_dir=input_path.parent)
            _validate_path_against_regex(
                bai, regex=_BAI_BASENAME_RE, column="bai", line_number=line_number,
                error_code=ErrorCode.INVALID_BAM, expected="BAM index (.bai)",
            )
            resolved["bai"] = bai


# --- aligned-input resolution (BAM+BAI / CRAM+CRAI) ------------------------

def _resolve_aligned_inputs(
    row: dict[str, str],
    *,
    input_path: Path,
    line_number: int,
    sample: str,
    resolved: dict[str, Path | str],
    bam_paths: list[Path | str],
    cram_paths: list[Path | str],
) -> None:
    has_bam = any(str(row.get(c, "")).strip() for c in _BAM_COLUMNS)
    has_cram = any(str(row.get(c, "")).strip() for c in _CRAM_COLUMNS)

    if has_bam and has_cram:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="A samplesheet row cannot mix BAM and CRAM columns.",
            fix="Provide either bam+bai OR cram+crai per row — never both.",
            details={"line": line_number, "sample": sample},
        )
    if not has_bam and not has_cram:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="This step requires aligned-read inputs (BAM+BAI or CRAM+CRAI).",
            fix="Populate either (bam, bai) or (cram, crai) for each row.",
            details={"line": line_number, "sample": sample},
        )

    if has_bam:
        for col, regex in (("bam", _BAM_BASENAME_RE), ("bai", _BAI_BASENAME_RE)):
            raw = str(row.get(col, "")).strip()
            if not raw:
                raise SkillError(
                    stage="validation",
                    error_code=ErrorCode.INVALID_SAMPLESHEET,
                    message=f"`{col}` is required alongside the other BAM column.",
                    fix="Provide both bam and bai for BAM-mode rows.",
                    details={"line": line_number, "sample": sample, "column": col},
                )
            path = _resolve_path(raw, base_dir=input_path.parent)
            _validate_path_against_regex(
                path, regex=regex, column=col, line_number=line_number,
                error_code=ErrorCode.INVALID_BAM,
                expected="BAM (.bam)" if col == "bam" else "BAM index (.bai)",
            )
            resolved[col] = path
            bam_paths.append(path)
    else:
        for col, regex in (("cram", _CRAM_BASENAME_RE), ("crai", _CRAI_BASENAME_RE)):
            raw = str(row.get(col, "")).strip()
            if not raw:
                raise SkillError(
                    stage="validation",
                    error_code=ErrorCode.INVALID_SAMPLESHEET,
                    message=f"`{col}` is required alongside the other CRAM column.",
                    fix="Provide both cram and crai for CRAM-mode rows.",
                    details={"line": line_number, "sample": sample, "column": col},
                )
            path = _resolve_path(raw, base_dir=input_path.parent)
            _validate_path_against_regex(
                path, regex=regex, column=col, line_number=line_number,
                error_code=ErrorCode.INVALID_SAMPLESHEET,
                expected="CRAM (.cram)" if col == "cram" else "CRAM index (.crai)",
            )
            resolved[col] = path
            cram_paths.append(path)


# --- single-file resolvers (table / vcf) -----------------------------------

def _resolve_table(row: dict[str, str], *, input_path: Path, line_number: int, sample: str) -> Path | str:
    raw = str(row.get("table", "")).strip()
    if not raw:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="`table` (BQSR recalibration table) is required at --step recalibrate.",
            fix="Provide a path to the .table file emitted by the prepare_recalibration step.",
            details={"line": line_number, "sample": sample},
        )
    path = _resolve_path(raw, base_dir=input_path.parent)
    _validate_path_against_regex(
        path, regex=_TABLE_BASENAME_RE, column="table", line_number=line_number,
        error_code=ErrorCode.INVALID_SAMPLESHEET, expected="recalibration table (.table)",
    )
    return path


def _resolve_vcf(row: dict[str, str], *, input_path: Path, line_number: int, sample: str) -> Path | str:
    raw = str(row.get("vcf", "")).strip()
    if not raw:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="`vcf` is required at --step annotate.",
            fix="Provide a path to a bgzipped VCF (.vcf.gz) for annotation.",
            details={"line": line_number, "sample": sample},
        )
    path = _resolve_path(raw, base_dir=input_path.parent)
    _validate_path_against_regex(
        path, regex=_VCF_BASENAME_RE, column="vcf", line_number=line_number,
        error_code=ErrorCode.INVALID_SAMPLESHEET, expected="VCF (.vcf or .vcf.gz; bgzipped recommended)",
    )
    return path


# --- low-level path / extension helpers ------------------------------------

def _is_uri(value: str) -> bool:
    return any(value.startswith(prefix) for prefix in _URI_PREFIXES)


def _resolve_path(value: str, *, base_dir: Path) -> Path | str:
    if _is_uri(value):
        return value
    path = Path(os.path.expanduser(value))
    if not path.is_absolute():
        path = base_dir / path
    return path.resolve()


def _path_name(path: Path | str) -> str:
    if isinstance(path, Path):
        return path.name
    return path.rstrip("/").rsplit("/", 1)[-1]


def _path_to_text(path: Path | str | None) -> str:
    if path is None:
        return ""
    return path.as_posix() if isinstance(path, Path) else path


def _path_text_has_whitespace(path: Path | str) -> bool:
    return any(ch.isspace() for ch in _path_to_text(path))


def _path_exists_or_uri(path: Path | str) -> bool:
    if isinstance(path, str) and _is_uri(path):
        return True
    return isinstance(path, Path) and path.exists()


def _validate_fastq_path(path: Path | str, *, column: str, line_number: int) -> None:
    filename = _path_name(path)
    if not _FASTQ_BASENAME_RE.match(filename):
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_FASTQ,
            message="FASTQ basename is not compatible with nf-core/sarek.",
            fix=(
                "FASTQs must be gzipped and end in .fastq.gz or .fq.gz; "
                "basenames cannot contain whitespace."
            ),
            details={"line": line_number, "column": column, "filename": filename},
        )
    if not _path_exists_or_uri(path):
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.MISSING_FASTQ,
            message="FASTQ file was not found.",
            fix="Correct the FASTQ path in the samplesheet.",
            details={"line": line_number, "column": column, "path": _path_to_text(path)},
        )


def _validate_spring_path(path: Path | str, *, column: str, line_number: int) -> None:
    filename = _path_name(path)
    if _path_text_has_whitespace(path):
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_FASTQ,
            message=f"`{column}` path must not contain whitespace (schema pattern ^\\S+).",
            fix="Move or rename the Spring file/path so the complete path contains no spaces.",
            details={"line": line_number, "column": column, "path": _path_to_text(path)},
        )
    if not _SPRING_BASENAME_RE.match(filename):
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_FASTQ,
            message="Spring basename is not compatible with nf-core/sarek.",
            fix="Spring archives must end in .fastq.gz.spring or .fq.gz.spring.",
            details={"line": line_number, "column": column, "filename": filename},
        )
    if not _path_exists_or_uri(path):
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.MISSING_FASTQ,
            message="Spring file was not found.",
            fix="Correct the spring_1/spring_2 path in the samplesheet.",
            details={"line": line_number, "column": column, "path": _path_to_text(path)},
        )


def _validate_bam_path(path: Path | str, *, column: str, line_number: int) -> None:
    filename = _path_name(path)
    if _path_text_has_whitespace(path):
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_BAM,
            message=f"`{column}` path must not contain whitespace (schema pattern ^\\S+).",
            fix="Move or rename the BAM path so the complete path contains no spaces.",
            details={"line": line_number, "column": column, "path": _path_to_text(path)},
        )
    if not _BAM_BASENAME_RE.match(filename):
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_BAM,
            message="BAM path does not have a .bam extension.",
            fix="Ensure the bam column points to a file ending in .bam.",
            details={"line": line_number, "column": column, "filename": filename},
        )
    if not _path_exists_or_uri(path):
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.MISSING_INPUT,
            message="BAM file referenced in samplesheet was not found.",
            fix=f"Correct the {column} path or remove it from the samplesheet.",
            details={"line": line_number, "column": column, "path": _path_to_text(path)},
        )


def _validate_path_against_regex(
    path: Path | str,
    *,
    regex: re.Pattern[str],
    column: str,
    line_number: int,
    error_code: str,
    expected: str,
) -> None:
    filename = _path_name(path)
    if _path_text_has_whitespace(path):
        raise SkillError(
            stage="validation",
            error_code=error_code,
            message=f"`{column}` path must not contain whitespace (schema pattern ^\\S+).",
            fix=f"Move or rename the {column} path so the complete path contains no spaces.",
            details={"line": line_number, "column": column, "path": _path_to_text(path)},
        )
    if not regex.match(filename):
        raise SkillError(
            stage="validation",
            error_code=error_code,
            message=f"`{column}` does not look like a valid {expected} file.",
            fix=f"Ensure the {column} column points to a {expected} file.",
            details={"line": line_number, "column": column, "filename": filename},
        )
    if not _path_exists_or_uri(path):
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.MISSING_INPUT,
            message=f"{expected} referenced in samplesheet was not found.",
            fix=f"Correct the {column} path or remove it from the samplesheet.",
            details={"line": line_number, "column": column, "path": _path_to_text(path)},
        )


# --- cross-row consistency checks ------------------------------------------

def _reject_sample_multiple_patients(
    sample_patients: dict[str, str],
    patient: str,
    sample: str,
    line_number: int,
) -> None:
    previous_patient = sample_patients.setdefault(sample, patient)
    if previous_patient != patient:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message=f"Sample '{sample}' is assigned to multiple patients.",
            fix="Each sample ID must belong to exactly one patient.",
            details={
                "line": line_number,
                "sample": sample,
                "patient": patient,
                "previous_patient": previous_patient,
            },
        )


def _reject_duplicate_patient_sample_status_lane(
    seen_keys: set[tuple[str, str, int, str]],
    patient: str,
    sample: str,
    status: int,
    lane: str,
    line_number: int,
) -> None:
    key = (patient, sample, status, lane)
    if key in seen_keys:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message=f"Duplicate patient/sample/status/lane tuple: {patient}/{sample}/{status}/{lane}.",
            fix="Each row must have a unique patient-sample-status-lane tuple.",
            details={"line": line_number, "patient": patient, "sample": sample, "status": status, "lane": lane},
        )
    seen_keys.add(key)


# --- tumor / normal pairing ------------------------------------------------

def _build_pairings(rows_by_patient: dict[str, list[dict[str, Any]]]) -> list[dict[str, Any]]:
    pairings: list[dict[str, Any]] = []
    for patient, entries in rows_by_patient.items():
        # Collapse per sample/status pair: upstream deduplicates
        # meta.subMap('sample', 'status'), not sample alone.
        samples_by_status: dict[int, list[str]] = defaultdict(list)
        seen: set[tuple[str, int]] = set()
        for entry in entries:
            sample = str(entry["sample"])
            status = int(str(entry["status"]))
            key = (sample, status)
            if key in seen:
                continue
            seen.add(key)
            samples_by_status[status].append(sample)

        normals = samples_by_status.get(0, [])
        tumors = samples_by_status.get(1, [])

        # Mirror Sarek 3.8.1 exactly: its samplesheet_to_channel validation
        # rejects multiple normals only when paired with exactly one tumor.
        # Multiple germline samples are valid, and we must not reject inputs
        # accepted by upstream solely because they share a patient identifier.
        if len(tumors) == 1 and len(normals) > 1:
            raise SkillError(
                stage="validation",
                error_code=ErrorCode.INVALID_PAIRING,
                message=(
                    f"Patient '{patient}' has more than one normal (status=0) sample "
                    "and exactly one tumor (status=1) sample."
                ),
                fix=(
                    "For a single-tumor paired analysis, keep one matched normal; "
                    "otherwise structure patient/sample statuses as accepted by Sarek."
                ),
                details={"patient": patient, "normals": normals, "tumors": tumors},
            )

        if normals and tumors:
            mode = "somatic_paired"
        elif tumors and not normals:
            mode = "tumor_only"
        elif normals and not tumors:
            mode = "germline"
        else:
            # Should not happen — every row carries a status.
            raise SkillError(
                stage="validation",
                error_code=ErrorCode.INVALID_PAIRING,
                message=f"Patient '{patient}' has no usable samples.",
                fix="Ensure each patient has at least one row with a valid status.",
                details={"patient": patient},
            )

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


def _enforce_varlociraptor_contamination(
    *,
    tools: list[str],
    rows_by_patient: dict[str, list[dict[str, Any]]],
    analysis_mode: str,
) -> None:
    if "varlociraptor" not in {t.lower() for t in tools}:
        return
    if analysis_mode == "germline":
        return  # varlociraptor on pure germline doesn't need contamination

    missing: list[dict[str, Any]] = []
    for patient, entries in rows_by_patient.items():
        tumor_entries = [e for e in entries if int(str(e["status"])) == 1]
        if not tumor_entries:
            continue
        for entry in tumor_entries:
            if entry["contamination"] is None:
                missing.append({
                    "patient": patient,
                    "sample": entry["sample"],
                    "line": entry["line"],
                })
    if missing:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.VARLOCIRAPTOR_MISSING_CONTAMINATION,
            message=(
                "Varlociraptor is requested but tumor samples are missing the "
                "`contamination` value."
            ),
            fix=(
                "Add a `contamination` column (float in [0, 1]) to every tumor row "
                "(status=1) when --tools includes varlociraptor."
            ),
            details={"missing": missing, "tools": tools, "analysis_mode": analysis_mode},
        )


# --- output writing --------------------------------------------------------

def _write_normalized_samplesheet(
    output_path: Path,
    fieldnames: list[str],
    rows: list[dict[str, str]],
) -> None:
    write_csv_lf(output_path, fieldnames, rows)
