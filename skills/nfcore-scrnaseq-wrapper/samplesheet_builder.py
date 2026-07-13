from __future__ import annotations

import csv
import re
import sys
from pathlib import Path
from typing import cast

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))
sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules

purge_foreign_bare_modules("errors", "nfcore_4_1_0_contract", "schemas")

from errors import ErrorCode, SkillError
from nfcore_4_1_0_contract import CELLRANGER_FAMILY_PRESETS
from schemas import (
    REQUIRED_SAMPLE_COLUMNS,
    SUPPORTED_FEATURE_TYPE_VALUES,
    SUPPORTED_SAMPLE_COLUMNS,
    SUPPORTED_SAMPLE_TYPE_VALUES,
)

_SAMPLE_NAME_RE = re.compile(r"^\S+$")
_SAMPLE_WHITESPACE_RE = re.compile(r"\s+")
_FASTQ_BASENAME_RE = re.compile(r"^[^\s/]+\.f(ast)?q\.gz$")
_CELLRANGER_READ_MARKER_RE = re.compile(
    r"(?P<prefix>.*)(?P<marker>_R[12])(?P<suffix>(?:_001)?\.f(?:ast)?q\.gz)$"
)
# 10x Illumina naming convention for Cell Ranger ARC. The extension accepts both
# .fastq.gz and .fq.gz to match the upstream samplesheet schema (audit F-4); the
# basename, lane, read-marker and _001 segments stay strict per Cell Ranger.
_TENX_FASTQ_RE = re.compile(
    r"^[^\s/]+_S\d+_L\d{3}_(R1|R2|R3|I1|I2)_001\.f(ast)?q\.gz$"
)
_REMOTE_URI_RE = re.compile(r"^[a-zA-Z][a-zA-Z0-9+.-]*://")
_FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz")
_BASE_OUTPUT_COLUMNS = ("sample", "fastq_1", "fastq_2", "expected_cells", "seq_center")
_REQUIRED_FASTQ_COLUMNS = ("fastq_1", "fastq_2")
_OPTIONAL_FASTQ_COLUMNS = ("fastq_barcode",)

# A resolved FASTQ value is either a normalized local Path or a remote URI kept
# verbatim as a str (s3://, gs://, https://, …). nf-core/scrnaseq supports remote
# inputs and Nextflow stages them in the execution context, so — matching
# nfcore-sarek-wrapper and nfcore-rnaseq-wrapper — the wrapper passes remote URIs
# through unchanged and validates only their basename, deferring existence to
# Nextflow. Local-first (no data exfiltration) is preserved: this only *reads*
# user-specified inputs.
_ResolvedPath = Path | str


def _is_remote_uri(value: str) -> bool:
    return bool(_REMOTE_URI_RE.match(str(value)))


def _resolved_name(value: _ResolvedPath) -> str:
    """Basename of a resolved FASTQ path or remote URI."""
    return Path(str(value)).name


def _resolved_output(value: _ResolvedPath) -> str:
    """Output-samplesheet string: remote URIs verbatim, local paths as posix."""
    return value if isinstance(value, str) else value.as_posix()


def validate_and_normalize_samplesheet(
    input_path: Path,
    output_path: Path,
    *,
    expected_cells_override: int | None = None,
    preset: str | None = None,
) -> dict[str, object]:
    fieldnames, rows = _read_samplesheet(input_path)
    _validate_preset_columns(fieldnames, preset=preset)
    _validate_expected_cells_override_scope(
        rows, expected_cells_override=expected_cells_override
    )
    unknown_columns = _unknown_columns(fieldnames)
    output_columns = _build_output_columns(fieldnames)
    normalized_rows, sample_names, fastq_paths = _normalize_rows(
        rows,
        input_path=input_path,
        output_columns=output_columns,
        expected_cells_override=expected_cells_override,
        preset=preset,
    )
    empty_expected_cells = sum(
        1 for row in normalized_rows if not row.get("expected_cells", "")
    )
    if empty_expected_cells > 0:
        # Human-facing notice → stdout. stderr is reserved for the wrapper's
        # structured error payload (json.dumps on failure); emitting diagnostics
        # there would make that error JSON unparseable for machine consumers.
        print(
            f"WARNING: {empty_expected_cells} sample(s) have an empty 'expected_cells' value. "
            "The upstream pipeline will use auto-estimation for these samples.",
            file=sys.stdout,
        )
    # nf-core/scrnaseq usage docs: "Since cellranger v7, it is not recommended
    # anymore to supply the --expected-cells parameter." Surface that advisory for
    # the whole Cell Ranger family (cellranger/cellrangerarc/cellrangermulti) —
    # cellranger-arc count also has no --expect-cells and auto-estimates, so the
    # guidance applies there too. The value is still passed through (guidance, not
    # a hard error). Uses the centralised family constant so the set cannot drift.
    if preset in CELLRANGER_FAMILY_PRESETS:
        with_expected_cells = sum(
            1 for row in normalized_rows if row.get("expected_cells", "")
        )
        if with_expected_cells > 0:
            print(
                f"WARNING: 'expected_cells' is set for {with_expected_cells} sample(s) under the "
                f"{preset!r} preset. Since Cell Ranger v7 the nf-core/scrnaseq docs no longer "
                "recommend supplying expected_cells; Cell Ranger auto-estimates cell counts.",
                file=sys.stdout,
            )

    _write_normalized_samplesheet(output_path, output_columns, normalized_rows)

    return {
        "normalized_path": output_path,
        # Count SAMPLES, not CSV rows: a lane-split sample occupies several rows but is
        # one sample (nf-core concatenates them). This value is written into the
        # provenance bundle and preferred by report.md, so counting rows overstated the
        # run in its own audit trail. `sample_names` is already deduplicated above.
        "sample_count": len(sample_names),
        "sample_names": sample_names,
        "fastq_paths": fastq_paths,
        "sample_types": sorted(
            {
                row.get("sample_type", "")
                for row in normalized_rows
                if row.get("sample_type", "")
            }
        ),
        "feature_types": sorted(
            {
                row.get("feature_type", "")
                for row in normalized_rows
                if row.get("feature_type", "")
            }
        ),
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

    with input_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        # nf-core/scrnaseq validates the samplesheet by column NAME (assets/
        # schema_input.json), not by position, so we only require presence here.
        # Column order is normalized on write (_build_output_columns leads with
        # sample,fastq_1,fastq_2), so any input order is accepted without loss.
        _validate_required_columns(fieldnames)
        rows = cast("list[dict[str, str]]", list(reader))

    if not rows:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="Samplesheet is empty.",
            fix="Add at least one sample row to the input CSV.",
            details={"path": str(input_path)},
        )
    return fieldnames, rows


def _validate_expected_cells_override_scope(
    rows: list[dict[str, str]],
    *,
    expected_cells_override: int | None,
) -> None:
    if expected_cells_override is None:
        return
    normalized_samples = {
        _normalize_sample_name(str(row.get("sample", "")).strip())
        for row in rows
        if str(row.get("sample", "")).strip()
    }
    if len(normalized_samples) <= 1:
        return
    raise SkillError(
        stage="validation",
        error_code=ErrorCode.INVALID_SAMPLESHEET,
        message="A global expected_cells override is only safe for single-sample samplesheets.",
        fix=(
            "Set expected_cells per row in the samplesheet for multi-sample runs, "
            "or run one sample at a time when using --expected-cells."
        ),
        details={
            "field": "expected_cells_override",
            "sample_count": len(normalized_samples),
            "samples": sorted(normalized_samples),
        },
    )


def _validate_required_columns(fieldnames: list[str]) -> None:
    missing_columns = [
        name for name in REQUIRED_SAMPLE_COLUMNS if name not in fieldnames
    ]
    if missing_columns:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="Samplesheet is missing required columns.",
            fix="Ensure the header contains sample, fastq_1, and fastq_2.",
            details={"missing_columns": missing_columns},
        )


def _validate_preset_columns(fieldnames: list[str], *, preset: str | None) -> None:
    required_by_preset = {
        "cellrangerarc": ("sample_type", "fastq_barcode"),
        "cellrangermulti": ("feature_type",),
    }
    missing_columns = [
        name
        for name in required_by_preset.get(preset or "", ())
        if name not in fieldnames
    ]
    if missing_columns:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="Samplesheet is missing columns required by the selected preset.",
            fix=f"Add the required columns for preset {preset!r}: {', '.join(missing_columns)}.",
            details={"preset": preset, "missing_columns": missing_columns},
        )


def _unknown_columns(fieldnames: list[str]) -> list[str]:
    return [name for name in fieldnames if name not in SUPPORTED_SAMPLE_COLUMNS]


def _build_output_columns(fieldnames: list[str]) -> list[str]:
    # This file is passed directly to nf-core/scrnaseq, so emit only headers in
    # the pinned 4.1.0 assets/schema_input.json contract. Unsupported input
    # metadata remains visible through ``unknown_columns`` without reaching
    # nf-schema, which would report every such header as unidentified.
    extra_columns = [
        name
        for name in fieldnames
        if name not in _BASE_OUTPUT_COLUMNS and name in SUPPORTED_SAMPLE_COLUMNS
    ]
    return list(dict.fromkeys([*_BASE_OUTPUT_COLUMNS, *extra_columns]))


def _normalize_rows(
    rows: list[dict[str, str]],
    *,
    input_path: Path,
    output_columns: list[str],
    expected_cells_override: int | None,
    preset: str | None,
) -> tuple[list[dict[str, str]], list[str], list[Path]]:
    normalized_rows: list[dict[str, str]] = []
    sample_names: list[str] = []
    fastq_paths: list[_ResolvedPath] = []
    seen_rows: set[tuple[str, str, str]] = set()
    raw_samples_by_normalized_sample: dict[str, str] = {}
    metadata_by_sample: dict[str, dict[str, str]] = {}

    for line_number, row in enumerate(rows, start=2):
        normalized, sample, resolved_fastqs = _normalize_samplesheet_row(
            row,
            line_number=line_number,
            input_path=input_path,
            output_columns=output_columns,
            expected_cells_override=expected_cells_override,
            preset=preset,
        )
        _reject_sample_name_collision(
            raw_samples_by_normalized_sample, row, line_number, sample
        )
        _reject_inconsistent_repeated_sample_metadata(
            metadata_by_sample, line_number, sample, normalized
        )
        _reject_duplicate_fastq_row(seen_rows, line_number, sample, resolved_fastqs)
        normalized_rows.append(normalized)
        # One entry per SAMPLE, not per row. nf-core/scrnaseq documents re-sequencing a
        # sample across lanes by repeating the `sample` identifier on several rows (the
        # pipeline concatenates the reads), so rows > samples is normal and every row is
        # still written to the normalized sheet. Parity with nfcore-rnaseq-wrapper.
        if sample not in sample_names:
            sample_names.append(sample)
        fastq_paths.extend(resolved_fastqs.values())

    return normalized_rows, sample_names, fastq_paths


def _normalize_samplesheet_row(
    row: dict[str, str],
    *,
    line_number: int,
    input_path: Path,
    output_columns: list[str],
    expected_cells_override: int | None,
    preset: str | None,
) -> tuple[dict[str, str], str, dict[str, Path]]:
    sample = _validate_sample_name(row, line_number)
    sample_type = _validate_sample_type(row, line_number, preset=preset)
    feature_type = _validate_feature_type(row, line_number, preset=preset)
    resolved_fastqs = _resolve_fastq_columns(
        row, line_number=line_number, input_path=input_path, preset=preset
    )
    _validate_preset_fastq_naming(
        resolved_fastqs, preset=preset, line_number=line_number
    )
    expected_cells = _validate_expected_cells(row, line_number, expected_cells_override)
    normalized = {column: str(row.get(column, "")).strip() for column in output_columns}
    normalized.update(
        {
            "sample": sample,
            "fastq_1": _resolved_output(resolved_fastqs["fastq_1"]),
            "fastq_2": _resolved_output(resolved_fastqs["fastq_2"]),
            "expected_cells": expected_cells,
            "seq_center": str(row.get("seq_center", "")).strip(),
        }
    )
    if sample_type:
        normalized["sample_type"] = sample_type
    if feature_type:
        normalized["feature_type"] = feature_type
    for optional_fastq in _optional_fastq_columns_for_preset(preset):
        if optional_fastq in resolved_fastqs:
            normalized[optional_fastq] = _resolved_output(resolved_fastqs[optional_fastq])
    return normalized, sample, resolved_fastqs


def _validate_sample_name(row: dict[str, str], line_number: int) -> str:
    raw_sample = str(row.get("sample", "")).strip()
    sample = _normalize_sample_name(raw_sample)
    if not sample or not _SAMPLE_NAME_RE.match(sample):
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="Sample names must be non-empty and cannot contain whitespace after normalization.",
            fix="Provide a non-empty sample value; whitespace is normalized to underscores before validation.",
            details={
                "line": line_number,
                "sample": raw_sample,
                "normalized_sample": sample,
            },
        )
    return sample


def _normalize_sample_name(raw_sample: str) -> str:
    return _SAMPLE_WHITESPACE_RE.sub("_", raw_sample.strip())


def _validate_sample_type(
    row: dict[str, str], line_number: int, *, preset: str | None
) -> str:
    sample_type = str(row.get("sample_type", "")).strip().lower()
    if preset == "cellrangerarc" and not sample_type:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="cellrangerarc rows require sample_type to identify ATAC versus GEX libraries.",
            fix="Set sample_type to either 'atac' or 'gex' for every cellrangerarc row.",
            details={"line": line_number, "column": "sample_type"},
        )
    if not sample_type:
        return sample_type
    # Enum is validated whenever the column carries a value, regardless of preset:
    # assets/schema_input.json declares sample_type as a property-level enum, so
    # nf-schema would reject an invalid value at runtime for ANY aligner. Fail fast
    # in preflight instead of deferring to a late Nextflow error (audit H-08). Only
    # the *presence* requirement above is preset-specific (cellrangerarc).
    if sample_type not in SUPPORTED_SAMPLE_TYPE_VALUES:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="sample_type must be one of the values supported by nf-core/scrnaseq.",
            fix=f"Use one of: {', '.join(sorted(SUPPORTED_SAMPLE_TYPE_VALUES))}.",
            details={"line": line_number, "sample_type": sample_type},
        )
    return sample_type


def _validate_feature_type(
    row: dict[str, str], line_number: int, *, preset: str | None
) -> str:
    feature_type = str(row.get("feature_type", "")).strip().lower()
    if preset == "cellrangermulti" and not feature_type:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="cellrangermulti rows require feature_type to identify each library type.",
            fix=f"Set feature_type to one of: {', '.join(sorted(SUPPORTED_FEATURE_TYPE_VALUES))}.",
            details={"line": line_number, "column": "feature_type"},
        )
    if not feature_type:
        return feature_type
    # Enum is validated whenever the column carries a value, regardless of preset:
    # assets/schema_input.json declares feature_type as a property-level enum, so
    # nf-schema would reject an invalid value at runtime for ANY aligner. Fail fast
    # in preflight instead of deferring to a late Nextflow error (audit H-08). Only
    # the *presence* requirement above is preset-specific (cellrangermulti).
    if feature_type not in SUPPORTED_FEATURE_TYPE_VALUES:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="feature_type must be one of the values supported by nf-core/scrnaseq.",
            fix=f"Use one of: {', '.join(sorted(SUPPORTED_FEATURE_TYPE_VALUES))}.",
            details={"line": line_number, "feature_type": feature_type},
        )
    return feature_type


def _resolve_fastq_columns(
    row: dict[str, str],
    *,
    line_number: int,
    input_path: Path,
    preset: str | None,
) -> dict[str, _ResolvedPath]:
    resolved = {
        column: _resolve_required_fastq(
            row, column, line_number=line_number, input_path=input_path
        )
        for column in _REQUIRED_FASTQ_COLUMNS
    }
    for column in _optional_fastq_columns_for_preset(preset):
        raw_value = str(row.get(column, "")).strip()
        if _row_requires_fastq_barcode(row, column=column, preset=preset):
            resolved[column] = _resolve_required_fastq(
                row, column, line_number=line_number, input_path=input_path
            )
            continue
        if _looks_like_fastq_path(raw_value):
            resolved[column] = _resolve_existing_fastq(
                raw_value,
                column,
                line_number=line_number,
                input_path=input_path,
            )
    return resolved


def _resolve_required_fastq(
    row: dict[str, str], column: str, *, line_number: int, input_path: Path
) -> _ResolvedPath:
    raw_value = str(row.get(column, "")).strip()
    if not raw_value:
        _raise_missing_fastq_column(column, line_number)
    return _resolve_existing_fastq(
        raw_value, column, line_number=line_number, input_path=input_path
    )


def _raise_missing_fastq_column(column: str, line_number: int) -> None:
    if column == "fastq_barcode":
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="CellRanger ARC ATAC rows require fastq_barcode to point to a barcode FASTQ file.",
            fix="Provide fastq_barcode for each cellrangerarc row with sample_type=atac.",
            details={"line": line_number, "column": column},
        )
    raise SkillError(
        stage="validation",
        error_code=ErrorCode.INVALID_SAMPLESHEET,
        message=(
            "Both fastq_1 and fastq_2 are required for every row: nf-core/scrnaseq 4.1.0 "
            "lists 'sample', 'fastq_1' and 'fastq_2' in the samplesheet schema's required "
            "columns (assets/schema_input.json)."
        ),
        fix="Provide both FASTQ mates for each row.",
        details={"line": line_number, "column": column},
    )


def _optional_fastq_columns_for_preset(preset: str | None) -> tuple[str, ...]:
    return (
        _OPTIONAL_FASTQ_COLUMNS
        if preset in {"cellrangerarc", "cellrangermulti"}
        else ()
    )


def _row_requires_fastq_barcode(
    row: dict[str, str], *, column: str, preset: str | None
) -> bool:
    sample_type = str(row.get("sample_type", "")).strip().lower()
    return (
        preset == "cellrangerarc"
        and column == "fastq_barcode"
        and sample_type == "atac"
    )


def _looks_like_fastq_path(value: str) -> bool:
    if not value:
        return False
    return value.endswith(_FASTQ_SUFFIXES)


def _resolve_existing_fastq(
    raw_path: str, column: str, *, line_number: int, input_path: Path
) -> _ResolvedPath:
    # Remote URIs (s3://, gs://, https://, …) are passed through verbatim: nf-core
    # supports remote inputs and Nextflow stages them in the execution context, so
    # the wrapper validates only the basename and defers existence to Nextflow.
    if _is_remote_uri(raw_path):
        _validate_fastq_path(raw_path, raw_path, column, line_number)
        return raw_path
    fastq_path = Path(raw_path).expanduser()
    fastq_path = (
        (input_path.parent / fastq_path).resolve()
        if not fastq_path.is_absolute()
        else fastq_path.resolve()
    )
    _validate_fastq_path(fastq_path, raw_path, column, line_number)
    return fastq_path


def _validate_fastq_path(
    value: _ResolvedPath, raw_path: str, column: str, line_number: int
) -> None:
    # Remote URI: validate the basename only; existence/readability are deferred to
    # Nextflow staging (it reads in the true execution context).
    if _is_remote_uri(raw_path):
        if not _FASTQ_BASENAME_RE.match(_resolved_name(value)):
            raise SkillError(
                stage="validation",
                error_code=ErrorCode.INVALID_FASTQ,
                message="FASTQ filenames must match the nf-core/scrnaseq schema.",
                fix="Use a basename without whitespace and with lowercase .fastq.gz or .fq.gz extension.",
                details={
                    "line": line_number,
                    "column": column,
                    "path": raw_path,
                    "filename": _resolved_name(value),
                },
            )
        return
    fastq_path = value  # local Path
    if not fastq_path.exists():
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.MISSING_FASTQ,
            message="A FASTQ file listed in the samplesheet does not exist.",
            fix="Correct the path in the samplesheet before re-running.",
            details={"line": line_number, "column": column, "path": raw_path},
        )
    if not _FASTQ_BASENAME_RE.match(fastq_path.name):
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_FASTQ,
            message="FASTQ filenames must match the nf-core/scrnaseq schema.",
            fix="Use a basename without whitespace and with lowercase .fastq.gz or .fq.gz extension.",
            details={
                "line": line_number,
                "column": column,
                "path": raw_path,
                "filename": fastq_path.name,
            },
        )
    if not fastq_path.is_file():
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.MISSING_FASTQ,
            message="A FASTQ path exists but does not point to a regular file.",
            fix="Ensure the FASTQ path refers to a regular file, not a directory.",
            details={"line": line_number, "column": column, "path": raw_path},
        )
    # NOTE: readability is intentionally NOT pre-checked here. Nextflow reads the FASTQ
    # data in the execution context (often a root container under the default Docker
    # profile), so an os.access(R_OK) pre-check by the launching user would false-block
    # valid runs. Existence + regular-file type are validated; readability is deferred to
    # Nextflow's staging. (Mirrors nfcore-rnaseq-wrapper's documented policy.)


def _validate_preset_fastq_naming(
    resolved_fastqs: dict[str, _ResolvedPath],
    *,
    preset: str | None,
    line_number: int,
) -> None:
    # FASTQ-naming enforcement is deliberately asymmetric across the Cell Ranger
    # family (audit F-5):
    #   * cellranger    → lenient R1/R2 pair-key check.
    #   * cellrangerarc → strict 10x Illumina naming (ATAC needs the barcode read).
    #   * cellrangermulti → NOT validated here, on purpose. Cell Ranger Multi maps
    #     libraries through its own multi samplesheet (--cellranger-multi-barcodes)
    #     and a [libraries] config, so the per-file 10x convention is resolved by
    #     the multi config rather than the wrapper's main samplesheet. Imposing the
    #     strict regex here would reject valid Multi inputs; the naming is left to
    #     Cell Ranger, which fails clearly if a name is wrong.
    if preset == "cellranger":
        _validate_cellranger_fastq_pair(resolved_fastqs, line_number=line_number)
    if preset == "cellrangerarc":
        _validate_cellrangerarc_fastq_names(resolved_fastqs, line_number=line_number)


def _validate_cellranger_fastq_pair(
    resolved_fastqs: dict[str, _ResolvedPath], *, line_number: int
) -> None:
    r1_name = _resolved_name(resolved_fastqs["fastq_1"])
    r2_name = _resolved_name(resolved_fastqs["fastq_2"])
    r1_key = _cellranger_pair_key(r1_name, expected_marker="_R1")
    r2_key = _cellranger_pair_key(r2_name, expected_marker="_R2")
    if r1_key and r1_key == r2_key:
        return
    raise SkillError(
        stage="validation",
        error_code=ErrorCode.INVALID_SAMPLESHEET,
        message="Cell Ranger FASTQ pairs must differ only by R1/R2.",
        fix="Use matched FASTQ basenames such as sample_S1_L001_R1_001.fastq.gz and sample_S1_L001_R2_001.fastq.gz.",
        details={
            "line": line_number,
            "preset": "cellranger",
            "fastq_1": r1_name,
            "fastq_2": r2_name,
        },
    )


def _cellranger_pair_key(filename: str, *, expected_marker: str) -> str:
    match = _CELLRANGER_READ_MARKER_RE.match(filename)
    if not match or match.group("marker") != expected_marker:
        return ""
    return f"{match.group('prefix')}<READ>{match.group('suffix')}"


def _validate_cellrangerarc_fastq_names(
    resolved_fastqs: dict[str, _ResolvedPath], *, line_number: int
) -> None:
    for column, fastq_path in resolved_fastqs.items():
        if _TENX_FASTQ_RE.match(_resolved_name(fastq_path)):
            continue
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_FASTQ,
            message="Cell Ranger ARC FASTQ filenames must follow the 10x naming convention.",
            fix=(
                "Use names like sample_S1_L001_R1_001.fastq.gz, sample_S1_L001_R2_001.fastq.gz, "
                "and sample_S1_L001_I2_001.fastq.gz."
            ),
            details={
                "line": line_number,
                "preset": "cellrangerarc",
                "column": column,
                "filename": _resolved_name(fastq_path),
            },
        )


def _validate_expected_cells(
    row: dict[str, str],
    line_number: int,
    expected_cells_override: int | None,
) -> str:
    expected_cells = (
        expected_cells_override
        if expected_cells_override is not None
        else row.get("expected_cells", "")
    )
    expected_cells_str = str(expected_cells).strip()
    if not expected_cells_str:
        return expected_cells_str
    try:
        expected_cells_int = int(expected_cells_str)
    except ValueError as exc:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="expected_cells must be an integer when provided.",
            fix="Use an integer value for expected_cells or leave it blank.",
            details={"line": line_number, "value": expected_cells_str},
        ) from exc
    if expected_cells_int < 1:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="expected_cells must be a positive integer (>= 1).",
            fix="Set expected_cells to a positive integer or leave the column blank.",
            details={"line": line_number, "value": expected_cells_int},
        )
    return expected_cells_str


def _reject_sample_name_collision(
    raw_samples_by_normalized_sample: dict[str, str],
    row: dict[str, str],
    line_number: int,
    normalized_sample: str,
) -> None:
    raw_sample = str(row.get("sample", "")).strip()
    normalized_raw_sample = _normalize_sample_name(raw_sample)
    if normalized_raw_sample != normalized_sample:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="Sample name normalization was inconsistent.",
            fix="Report this as a bug with the offending samplesheet row.",
            details={
                "line": line_number,
                "sample": raw_sample,
                "normalized_sample": normalized_sample,
                "renormalized_sample": normalized_raw_sample,
            },
        )
    previous_raw_sample = raw_samples_by_normalized_sample.setdefault(
        normalized_sample, raw_sample
    )
    if previous_raw_sample == raw_sample:
        return
    raise SkillError(
        stage="validation",
        error_code=ErrorCode.INVALID_SAMPLESHEET,
        message="Distinct sample names collapse to the same normalized nf-core sample identifier.",
        fix=(
            "Rename samples so they remain unique after whitespace is converted to underscores "
            "(for example, avoid using both 'sample A' and 'sample_A')."
        ),
        details={
            "line": line_number,
            "sample": raw_sample,
            "previous_sample": previous_raw_sample,
            "normalized_sample": normalized_sample,
        },
    )


def _reject_duplicate_fastq_row(
    seen_rows: set[tuple[str, str, str]],
    line_number: int,
    sample: str,
    resolved_fastqs: dict[str, _ResolvedPath],
) -> None:
    row_key = (
        sample,
        _resolved_output(resolved_fastqs["fastq_1"]),
        _resolved_output(resolved_fastqs["fastq_2"]),
    )
    if row_key in seen_rows:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="Samplesheet contains a duplicate FASTQ row.",
            fix="Remove exact duplicate rows; repeated sample names are allowed only for distinct FASTQ pairs.",
            details={"line": line_number, "sample": sample},
        )
    seen_rows.add(row_key)


def _reject_inconsistent_repeated_sample_metadata(
    metadata_by_sample: dict[str, dict[str, str]],
    line_number: int,
    sample: str,
    normalized_row: dict[str, str],
) -> None:
    current_metadata = {
        "expected_cells": str(normalized_row.get("expected_cells", "")).strip(),
        "seq_center": str(normalized_row.get("seq_center", "")).strip(),
    }
    previous_metadata = metadata_by_sample.setdefault(sample, current_metadata)
    for column, current_value in current_metadata.items():
        previous_value = previous_metadata.get(column, "")
        if current_value == previous_value:
            continue
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message=f"Repeated sample rows must use the same {column} value.",
            fix=(
                f"Use the same {column} for every row of sample {sample!r}, "
                "or split runs from different metadata groups into distinct sample names."
            ),
            details={
                "line": line_number,
                "sample": sample,
                "column": column,
                "previous_value": previous_value,
                "value": current_value,
            },
        )


def _write_normalized_samplesheet(
    output_path: Path,
    output_columns: list[str],
    normalized_rows: list[dict[str, str]],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    # newline="" + lineterminator="\n": csv.writer defaults to CRLF, which would
    # make the bundle samplesheet differ from every other (LF) artifact and vary
    # its checksum. Force LF on every OS for a byte-stable bundle.
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(normalized_rows)
