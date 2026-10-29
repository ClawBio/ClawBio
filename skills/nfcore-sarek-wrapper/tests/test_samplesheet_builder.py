"""Regression tests for Sarek-compatible samplesheet handling."""
from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from errors import ErrorCode, SkillError
from samplesheet_builder import validate_and_normalize_samplesheet


def _touch(path: Path, content: str = "x") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def _row(tmp_path: Path, **overrides: str) -> dict[str, str]:
    fastq_1 = _touch(tmp_path / "S1_R1.fastq.gz")
    fastq_2 = _touch(tmp_path / "S1_R2.fastq.gz")
    row = {
        "patient": "P1",
        "sample": "S1",
        "lane": "L001",
        "fastq_1": fastq_1.as_posix(),
        "fastq_2": fastq_2.as_posix(),
        "sex": "NA",
        "status": "0",
    }
    row.update(overrides)
    return row


def _write_delimited(path: Path, rows: list[dict[str, str]], delimiter: str) -> Path:
    fieldnames = list(dict.fromkeys(k for row in rows for k in row))
    lines = [delimiter.join(fieldnames)]
    for row in rows:
        lines.append(delimiter.join(str(row.get(k, "")) for k in fieldnames))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


@pytest.mark.parametrize("suffix,writer", [
    (".csv", lambda path, rows: _write_delimited(path, rows, ",")),
    (".tsv", lambda path, rows: _write_delimited(path, rows, "\t")),
    (".json", lambda path, rows: path.write_text(json.dumps(rows), encoding="utf-8") or path),
    (".yaml", lambda path, rows: path.write_text(yaml.safe_dump(rows), encoding="utf-8") or path),
])
def test_samplesheet_accepts_official_structured_formats(tmp_path, suffix, writer):
    input_path = tmp_path / f"samplesheet{suffix}"
    writer(input_path, [_row(tmp_path)])

    report = validate_and_normalize_samplesheet(
        input_path,
        tmp_path / "out.csv",
        step="mapping",
    )

    assert report["sample_names"] == ["S1"]


def test_samplesheet_allows_url_inputs_without_local_exists(tmp_path):
    input_path = tmp_path / "remote.csv"
    _write_delimited(input_path, [_row(tmp_path, fastq_1="https://example.org/S1_R1.fastq.gz")], ",")

    report = validate_and_normalize_samplesheet(input_path, tmp_path / "out.csv", step="mapping")

    assert "https://example.org/S1_R1.fastq.gz" in report["fastq_paths"]


def test_samplesheet_rejects_single_fastq_mapping_input(tmp_path):
    input_path = _write_delimited(tmp_path / "single_fastq.csv", [_row(tmp_path, fastq_2="")], ",")

    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(input_path, tmp_path / "out.csv", step="mapping")

    assert exc.value.error_code == ErrorCode.INVALID_SAMPLESHEET


def test_samplesheet_rejects_lane_with_whitespace(tmp_path):
    # schema_input.json constrains lane to `^\S+$` (no internal whitespace).
    input_path = _write_delimited(tmp_path / "ws_lane.csv", [_row(tmp_path, lane="lane 1")], ",")

    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(input_path, tmp_path / "out.csv", step="mapping")

    assert exc.value.error_code == ErrorCode.INVALID_SAMPLESHEET


def test_samplesheet_rejects_spring2_without_spring1(tmp_path):
    # schema_input.json declares spring_2 -> spring_1.
    spring2 = _touch(tmp_path / "S1.fastq.gz.spring")
    row = {"patient": "P1", "sample": "S1", "lane": "L001",
           "spring_2": spring2.as_posix(), "sex": "NA", "status": "0"}
    input_path = _write_delimited(tmp_path / "spring2only.csv", [row], ",")
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(input_path, tmp_path / "out.csv", step="mapping")
    assert exc.value.error_code == ErrorCode.INVALID_SAMPLESHEET


def test_mapping_bam_captures_optional_bai(tmp_path):
    # Sarek propagates [meta, bam, bai] for a mapping-step BAM input.
    bam = _touch(tmp_path / "S1.bam")
    bai = _touch(tmp_path / "S1.bai")
    row = {"patient": "P1", "sample": "S1", "lane": "L001",
           "bam": bam.as_posix(), "bai": bai.as_posix(), "sex": "NA", "status": "0"}
    input_path = _write_delimited(tmp_path / "ubam.csv", [row], ",")
    report = validate_and_normalize_samplesheet(input_path, tmp_path / "out.csv", step="mapping")
    assert report["sample_names"] == ["S1"]


def test_samplesheet_report_surfaces_rows_by_patient_for_preflight(tmp_path):
    input_path = tmp_path / "remote.csv"
    _write_delimited(
        input_path,
        [_row(tmp_path, sex="NA", status="1", fastq_1="https://example.org/T1_R1.fastq.gz")],
        ",",
    )

    report = validate_and_normalize_samplesheet(input_path, tmp_path / "out.csv", step="mapping")

    assert report["rows_by_patient"]["P1"][0]["sample"] == "S1"
    assert report["rows_by_patient"]["P1"][0]["sex"] == "NA"


def test_samplesheet_allows_single_spring_archive(tmp_path):
    spring = _touch(tmp_path / "S1.fastq.gz.spring")
    row = _row(tmp_path, fastq_1="", fastq_2="", spring_1=spring.as_posix(), spring_2="")
    input_path = _write_delimited(tmp_path / "spring.csv", [row], ",")

    report = validate_and_normalize_samplesheet(input_path, tmp_path / "out.csv", step="mapping")

    assert report["spring_paths"] == [spring]


def test_samplesheet_rejects_sample_whitespace_without_normalizing(tmp_path):
    input_path = _write_delimited(tmp_path / "bad.csv", [_row(tmp_path, sample="Sample One")], ",")

    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(input_path, tmp_path / "out.csv", step="mapping")

    assert exc.value.error_code == ErrorCode.INVALID_SAMPLESHEET


def test_samplesheet_rejects_duplicate_patient_sample_status_lane(tmp_path):
    rows = [_row(tmp_path), _row(tmp_path)]
    input_path = _write_delimited(tmp_path / "dupe.csv", rows, ",")

    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(input_path, tmp_path / "out.csv", step="mapping")

    assert exc.value.error_code == ErrorCode.INVALID_SAMPLESHEET


def test_samplesheet_rejects_same_sample_across_multiple_patients(tmp_path):
    rows = [_row(tmp_path), _row(tmp_path, patient="P2", lane="L002")]
    input_path = _write_delimited(tmp_path / "multi_patient.csv", rows, ",")

    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(input_path, tmp_path / "out.csv", step="mapping")

    assert exc.value.error_code == ErrorCode.INVALID_SAMPLESHEET


def test_samplesheet_allows_multiple_normals_unless_exactly_one_tumor(tmp_path):
    # Sarek rejects multiple normals only with exactly one tumor. One normal
    # plus multiple tumors is a valid somatic_paired configuration.
    rows = [
        _row(tmp_path, sample="N1", lane="L001", status="0", fastq_1=_touch(tmp_path / "N1.fastq.gz").as_posix()),
        _row(tmp_path, sample="T1", lane="L001", status="1", fastq_1=_touch(tmp_path / "T1.fastq.gz").as_posix()),
        _row(tmp_path, sample="T2", lane="L001", status="1", fastq_1=_touch(tmp_path / "T2.fastq.gz").as_posix()),
    ]
    input_path = _write_delimited(tmp_path / "multi_normals.csv", rows, ",")

    report = validate_and_normalize_samplesheet(input_path, tmp_path / "out.csv", step="mapping")

    assert report["analysis_mode"] == "somatic_paired"


def test_samplesheet_rejects_multiple_normals_with_exactly_one_tumor(tmp_path):
    rows = [
        _row(tmp_path, sample="N1", lane="L001", status="0", fastq_1=_touch(tmp_path / "N1.fastq.gz").as_posix()),
        _row(tmp_path, sample="N2", lane="L001", status="0", fastq_1=_touch(tmp_path / "N2.fastq.gz").as_posix()),
        _row(tmp_path, sample="T1", lane="L001", status="1", fastq_1=_touch(tmp_path / "T1.fastq.gz").as_posix()),
    ]
    input_path = _write_delimited(tmp_path / "bad_multi_normals.csv", rows, ",")

    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(input_path, tmp_path / "out.csv", step="mapping")

    assert exc.value.error_code == ErrorCode.INVALID_PAIRING


def test_every_reported_path_group_is_remote_gated(tmp_path):
    """The remote-input gate reads the builder's path groups from one list; a new
    group the builder reports must be added there or remote URIs in it go ungated."""
    from preflight import SAMPLESHEET_PATH_GROUPS

    report = validate_and_normalize_samplesheet(
        _write_delimited(tmp_path / "samples.csv", [_row(tmp_path)], ","), tmp_path / "out.csv"
    )
    reported = {key for key in report if key.endswith("_paths") or key == "tables"}
    assert reported == set(SAMPLESHEET_PATH_GROUPS)
