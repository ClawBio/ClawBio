"""Regression tests for Sarek-compatible samplesheet handling."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import pytest
import yaml
from hypothesis import HealthCheck, example, given, settings
from hypothesis import strategies as st

from errors import ErrorCode, SkillError
from clawbio.common.tests.nfcore_test_imports import load_tested_module

_samplesheet_builder = load_tested_module(Path(__file__).resolve().parent.parent, "samplesheet_builder")
validate_and_normalize_samplesheet = _samplesheet_builder.validate_and_normalize_samplesheet


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



# Rules nf-core/sarek 3.8.1 enforces itself when Nextflow starts, before any task runs:
# nf-schema checks assets/schema_input.json (required patient/sample, ^\S+$ patterns,
# the sex/status enums, contamination's type, file patterns and existence,
# dependentRequired, uniqueEntries), and subworkflows/local/samplesheet_to_channel
# refuses a row whose inputs do not fit --step, duplicate patient-sample-status-lane
# rows, a sample under two patients, and two normals with one tumour (and, in
# preflight's domain, varlociraptor tumours without contamination). The wrapper
# keeps no copy: such a sheet is written
# through (paths resolved, since Nextflow runs from --output) for the pipeline to name.
def _fq(tmp_path: Path, name: str) -> str:
    return _touch(tmp_path / name).as_posix()


_NF_CORE_CASES = {
    "missing patient column": dict(rows=[{"sample": "S1", "lane": "L1"}], fastq=True, drop="patient"),
    "sample with whitespace": dict(rows=[{"sample": "Sample One"}], fastq=True),
    "patient with surrounding whitespace": dict(rows=[{"patient": " P1"}], fastq=True),
    "sex outside the enum": dict(rows=[{"sex": "male"}], fastq=True),
    "status outside the enum": dict(rows=[{"status": "2"}], fastq=True),
    "status that is not a number": dict(rows=[{"status": "tumor"}], fastq=True),
    "lane with whitespace": dict(rows=[{"lane": "lane 1"}], fastq=True),
    "contamination that is not a number": dict(rows=[{"status": "1", "contamination": "high"}], fastq=True),
    "FASTQ that is not gzipped": dict(rows=[{"fastq_1": "R1.fastq", "fastq_2": "R2.fastq"}], files=("R1.fastq", "R2.fastq")),
    "FASTQ that does not exist": dict(rows=[{"fastq_1": "missing_R1.fastq.gz", "fastq_2": "missing_R2.fastq.gz"}]),
    "mapping row with only fastq_1": dict(rows=[{"fastq_1": "R1.fastq.gz"}], files=("R1.fastq.gz",)),
    "mapping row without a lane": dict(rows=[{"lane": ""}], fastq=True),
    "mapping row without inputs": dict(rows=[{}]),
    "spring_2 without spring_1": dict(rows=[{"spring_2": "S1.fastq.gz.spring"}], files=("S1.fastq.gz.spring",)),
    "uBAM with the wrong extension": dict(rows=[{"bam": "S1.sam"}], files=("S1.sam",)),
    "bai with the wrong extension": dict(rows=[{"bam": "S1.bam", "bai": "S1.idx"}], files=("S1.bam", "S1.idx")),
    "variant_calling row without BAM or CRAM": dict(rows=[{"lane": ""}], step="variant_calling"),
    "variant_calling row with a bai but no bam": dict(rows=[{"lane": "", "bai": "S1.bai"}], files=("S1.bai",), step="variant_calling"),
    "variant_calling CRAM that does not exist": dict(rows=[{"lane": "", "cram": "S1.cram", "crai": "S1.crai"}], files=("S1.crai",), step="variant_calling"),
    "recalibration table with the wrong extension": dict(rows=[{"lane": "", "cram": "S1.cram", "crai": "S1.crai", "table": "S1.txt"}], files=("S1.cram", "S1.crai", "S1.txt"), step="recalibrate"),
    "annotate row without a VCF": dict(rows=[{"lane": ""}], step="annotate"),
    "annotate VCF that does not exist": dict(rows=[{"lane": "", "vcf": "S1.vcf.gz"}], step="annotate"),
    "duplicate patient-sample-status-lane": dict(rows=[{}, {}], fastq=True),
    "one sample under two patients": dict(rows=[{}, {"patient": "P2", "lane": "L2"}], fastq=True),
    "two normals and one tumour": dict(rows=[{"sample": "N1"}, {"sample": "N2"}, {"sample": "T1", "status": "1"}], fastq=True),
}


@pytest.mark.parametrize("case", sorted(_NF_CORE_CASES))
def test_rules_nf_core_enforces_reach_it_unchanged(tmp_path, case):
    spec = _NF_CORE_CASES[case]
    for name in spec.get("files", ()):
        _touch(tmp_path / name)
    base = {"patient": "P1", "sample": "S1", "lane": "L1", "sex": "XX", "status": "0"}
    if spec.get("fastq"):
        base.update(fastq_1=_fq(tmp_path, "R1.fastq.gz"), fastq_2=_fq(tmp_path, "R2.fastq.gz"))
    rows = [{**base, **row} for row in spec["rows"]]
    rows = [{k: v for k, v in row.items() if k != spec.get("drop")} for row in rows]
    sheet = _write_delimited(tmp_path / "samplesheet.csv", rows, ",")
    out = tmp_path / "out.csv"

    validate_and_normalize_samplesheet(sheet, out, step=spec.get("step", "mapping"))

    written = list(csv.DictReader(out.read_text(encoding="utf-8").splitlines()))
    assert len(written) == len(rows)
    for given_row, row in zip(rows, written):
        for column, value in given_row.items():
            if column in _PATH_COLUMNS and value:
                assert row[column] == (tmp_path / value).resolve().as_posix()
            elif column not in ("sex", "status"):
                assert row[column] == value, column


_PATH_COLUMNS = ("fastq_1", "fastq_2", "spring_1", "spring_2", "bam", "bai", "cram", "crai", "table", "vcf")


# What nf-core/sarek 3.8.1 does not refuse, the wrapper still does: the pipeline
# silently picks one input mode for a row that mixes them, passes a null index for a
# BAM/CRAM row without one at non-mapping steps (refused only when a lane is set),
# routes recalibrate rows without a table as if they had one, and the schema puts no
# bounds on contamination.
@pytest.mark.parametrize("step,row,needle", [
    ("mapping", {"bam": "S1.bam"}, "mix FASTQ, Spring and BAM"),
    ("variant_calling", {"lane": "", "bam": "S1.bam", "bai": "S1.bai", "cram": "S1.cram", "crai": "S1.crai"}, "mix BAM and CRAM"),
    ("variant_calling", {"lane": "", "bam": "S1.bam"}, "`bai` is required"),
    ("markduplicates", {"lane": "", "cram": "S1.cram"}, "`crai` is required"),
    ("recalibrate", {"lane": "", "cram": "S1.cram", "crai": "S1.crai"}, "`table`"),
])
def test_rules_nf_core_does_not_check_are_still_refused(tmp_path, step, row, needle):
    for name in ("S1.bam", "S1.bai", "S1.cram", "S1.crai"):
        _touch(tmp_path / name)
    base = {"patient": "P1", "sample": "S1", "lane": "L1", "status": "0"}
    if step == "mapping":
        base.update(fastq_1=_fq(tmp_path, "R1.fastq.gz"), fastq_2=_fq(tmp_path, "R2.fastq.gz"))
    sheet = _write_delimited(tmp_path / "samplesheet.csv", [{**base, **row}], ",")
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(sheet, tmp_path / "out.csv", step=step)
    assert exc.value.error_code == ErrorCode.INVALID_SAMPLESHEET
    assert needle in exc.value.message


_STATUS = st.one_of(st.sampled_from(["", "0", "1", "2", "01", "-1", "tumor", "nan", "1.0"]), st.text(max_size=5))


@settings(suppress_health_check=[HealthCheck.function_scoped_fixture], max_examples=60, deadline=None)
@given(
    patient=st.text(alphabet="PQ1 _\t", min_size=0, max_size=4),
    sample=st.text(alphabet="ST2 _", min_size=0, max_size=4),
    sex=st.sampled_from(["", "XX", "XY", "NA", "male", " XX"]),
    status=_STATUS,
    lane=st.text(alphabet="L3 ", max_size=3),
)
@example(patient="", sample="", sex="", status="tumor", lane="")
@example(patient="", sample="", sex="", status="0\x1e0", lane="")
def test_lenient_parsing_never_refuses_what_nf_schema_judges(tmp_path, patient, sample, sex, status, lane):
    """Any value of a column nf-schema validates is written through, and the
    pairing summary the wrapper derives from it never refuses or crashes."""
    if any(ch in value for value in (patient, sample, sex, status, lane) for ch in "\r\n,\"\x00"):
        return
    row = {"patient": patient, "sample": sample, "sex": sex, "status": status, "lane": lane,
           "fastq_1": _fq(tmp_path, "R1.fastq.gz"), "fastq_2": _fq(tmp_path, "R2.fastq.gz")}
    sheet = _write_delimited(tmp_path / "samplesheet.csv", [row], ",")
    report = validate_and_normalize_samplesheet(sheet, tmp_path / "out.csv", step="mapping")
    # Read as CSV, not with str.splitlines(), which also breaks lines at \x1c-\x1e.
    with (tmp_path / "out.csv").open(encoding="utf-8", newline="") as handle:
        written = next(csv.DictReader(handle))
    assert (written["patient"], written["sample"], written["lane"]) == (patient, sample, lane)
    assert report["analysis_mode"] in {"germline", "tumor_only", "somatic_paired", "mixed"}


@settings(suppress_health_check=[HealthCheck.function_scoped_fixture], max_examples=60, deadline=None)
@given(st.one_of(st.floats(allow_nan=True, allow_infinity=True), st.floats(0, 1)))
@example(float("nan"))
@example(0.0)
@example(1.0)
def test_contamination_is_accepted_only_as_a_fraction(tmp_path, value):
    """nf-schema checks only that contamination is a number; NaN, infinity and
    values outside [0, 1] must never reach Varlociraptor as a confident purity."""
    row = {"patient": "P1", "sample": "T1", "lane": "L1", "status": "1", "contamination": repr(value),
           "fastq_1": _fq(tmp_path, "R1.fastq.gz"), "fastq_2": _fq(tmp_path, "R2.fastq.gz")}
    sheet = _write_delimited(tmp_path / "samplesheet.csv", [row], ",")
    if math.isfinite(value) and 0.0 <= value <= 1.0:
        validate_and_normalize_samplesheet(sheet, tmp_path / "out.csv", step="mapping")
    else:
        with pytest.raises(SkillError) as exc:
            validate_and_normalize_samplesheet(sheet, tmp_path / "out.csv", step="mapping")
        assert exc.value.error_code == ErrorCode.INVALID_SAMPLESHEET
