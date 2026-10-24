from __future__ import annotations

import csv
from pathlib import Path
import sys

import pytest

_SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_SKILL_DIR))


def _purge_foreign_modules(*names: str) -> None:
    for name in names:
        module = sys.modules.get(name)
        module_file = Path(getattr(module, "__file__", "") or "")
        if module is not None and _SKILL_DIR not in module_file.parents and module_file != _SKILL_DIR / f"{name}.py":
            sys.modules.pop(name, None)


def _purge_local_modules(*names: str) -> None:
    for name in names:
        module = sys.modules.get(name)
        module_file = Path(getattr(module, "__file__", "") or "")
        if module is not None and (module_file == _SKILL_DIR / f"{name}.py" or _SKILL_DIR in module_file.parents):
            sys.modules.pop(name, None)


def _remove_skill_dir_from_sys_path() -> None:
    while str(_SKILL_DIR) in sys.path:
        sys.path.remove(str(_SKILL_DIR))


_purge_foreign_modules("errors", "schemas", "samplesheet_builder")

from errors import ErrorCode, SkillError
from samplesheet_builder import validate_and_normalize_samplesheet

_purge_local_modules("errors", "schemas", "samplesheet_builder")
_remove_skill_dir_from_sys_path()


def _touch(path: Path, text: str = "x") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _write_sheet(path: Path, rows: list[dict[str, str]], fieldnames: list[str] | None = None) -> Path:
    if fieldnames is None:
        seen: list[str] = []
        for row in rows:
            for key in row:
                if key not in seen:
                    seen.append(key)
        fieldnames = seen
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return path


def _valid_row(tmp_path: Path, **overrides: str) -> dict[str, str]:
    r1 = _touch(tmp_path / "reads" / "sample_R1.fastq.gz")
    row = {"sample": "sampleA", "fastq_1": str(r1), "strandedness": "auto"}
    row.update(overrides)
    return row


def test_accepts_single_end_fastq_and_writes_optional_fastq_2_blank(tmp_path):
    src = _write_sheet(tmp_path / "samplesheet.csv", [_valid_row(tmp_path)])
    out = tmp_path / "normalized.csv"

    result = validate_and_normalize_samplesheet(src, out)

    row = next(csv.DictReader(out.open(encoding="utf-8")))
    assert result["sample_count"] == 1
    assert row["fastq_2"] == ""
    assert Path(row["fastq_1"]).is_absolute()


def test_resolves_relative_fastq_against_samplesheet_directory(tmp_path):
    sheet_dir = tmp_path / "inputs"
    r1 = _touch(sheet_dir / "reads" / "sample_R1.fastq.gz")
    src = _write_sheet(sheet_dir / "samplesheet.csv", [{"sample": "sampleA", "fastq_1": "reads/sample_R1.fastq.gz", "strandedness": "auto"}])
    out = tmp_path / "normalized.csv"

    validate_and_normalize_samplesheet(src, out)

    row = next(csv.DictReader(out.open(encoding="utf-8")))
    assert row["fastq_1"] == r1.resolve().as_posix()


def test_preserves_remote_fastq_uri_without_local_existence_check(tmp_path):
    uri = "s3://bucket.example.org/reads/sample_R1.fastq.gz"
    src = _write_sheet(
        tmp_path / "samplesheet.csv",
        [{"sample": "sampleA", "fastq_1": uri, "strandedness": "auto"}],
    )
    out = tmp_path / "normalized.csv"

    result = validate_and_normalize_samplesheet(src, out)

    row = next(csv.DictReader(out.open(encoding="utf-8")))
    assert row["fastq_1"] == uri
    assert result["fastq_paths"] == [uri]


def test_normalizes_sample_whitespace_to_underscores(tmp_path):
    src = _write_sheet(tmp_path / "samplesheet.csv", [_valid_row(tmp_path, sample="sample  A")])

    result = validate_and_normalize_samplesheet(src, tmp_path / "normalized.csv")

    assert result["sample_names"] == ["sample_A"]


def test_rejects_sample_collision_after_whitespace_normalization(tmp_path):
    r1a = _touch(tmp_path / "a_R1.fastq.gz")
    r1b = _touch(tmp_path / "b_R1.fastq.gz")
    src = _write_sheet(
        tmp_path / "samplesheet.csv",
        [
            {"sample": "sample A", "fastq_1": str(r1a), "strandedness": "auto"},
            {"sample": "sample_A", "fastq_1": str(r1b), "strandedness": "auto"},
        ],
    )

    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "normalized.csv")

    assert exc.value.error_code == ErrorCode.INVALID_SAMPLESHEET


def test_rejects_duplicate_fastq_rows(tmp_path):
    r1 = _touch(tmp_path / "sample_R1.fastq.gz")
    src = _write_sheet(
        tmp_path / "samplesheet.csv",
        [
            {"sample": "sampleA", "fastq_1": str(r1), "strandedness": "auto"},
            {"sample": "sampleA", "fastq_1": str(r1), "strandedness": "auto"},
        ],
    )

    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "normalized.csv")

    assert exc.value.error_code == ErrorCode.INVALID_SAMPLESHEET


def test_rejects_row_mixing_fastq_and_bam_columns(tmp_path):
    gbam = _touch(tmp_path / "sample.genome.bam")
    src = _write_sheet(
        tmp_path / "samplesheet.csv",
        [_valid_row(tmp_path, genome_bam=str(gbam))],  # fastq_1 set by default
    )

    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "normalized.csv")

    assert exc.value.error_code == ErrorCode.INVALID_SAMPLESHEET
    assert "FASTQ and BAM" in exc.value.message


def test_reports_unknown_columns_and_preserves_them(tmp_path):
    src = _write_sheet(tmp_path / "samplesheet.csv", [_valid_row(tmp_path, condition="treated")])
    out = tmp_path / "normalized.csv"

    result = validate_and_normalize_samplesheet(src, out)

    row = next(csv.DictReader(out.open(encoding="utf-8")))
    assert result["unknown_columns"] == ["condition"]
    assert row["condition"] == "treated"


def test_rejects_empty_samplesheet(tmp_path):
    (tmp_path / "samplesheet.csv").write_text("sample,fastq_1,strandedness\n", encoding="utf-8")

    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(tmp_path / "samplesheet.csv", tmp_path / "normalized.csv")

    assert exc.value.error_code == ErrorCode.INVALID_SAMPLESHEET


def test_rejects_missing_samplesheet_file(tmp_path):
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(tmp_path / "missing.csv", tmp_path / "normalized.csv")

    assert exc.value.error_code == ErrorCode.MISSING_INPUT


def test_normalized_output_has_stable_rnaseq_columns(tmp_path):
    src = _write_sheet(tmp_path / "samplesheet.csv", [_valid_row(tmp_path, condition="treated")])
    out = tmp_path / "normalized.csv"

    validate_and_normalize_samplesheet(src, out)

    reader = csv.DictReader(out.open(encoding="utf-8"))
    assert reader.fieldnames[:6] == [
        "sample",
        "fastq_1",
        "fastq_2",
        "strandedness",
        "genome_bam",
        "transcriptome_bam",
    ]


# ── mixed PE/SE / BAM coverage / metadata pass-through ──────────────────────


def test_accepts_mixed_paired_and_single_end_rows_in_same_sheet(tmp_path):
    pe_r1 = _touch(tmp_path / "reads" / "pe_R1.fastq.gz")
    pe_r2 = _touch(tmp_path / "reads" / "pe_R2.fastq.gz")
    se_r1 = _touch(tmp_path / "reads" / "se_R1.fastq.gz")
    rows = [
        {"sample": "PE", "fastq_1": str(pe_r1), "fastq_2": str(pe_r2), "strandedness": "auto"},
        {"sample": "SE", "fastq_1": str(se_r1), "fastq_2": "", "strandedness": "forward"},
    ]
    src = _write_sheet(tmp_path / "samplesheet.csv", rows)
    out = tmp_path / "normalized.csv"

    result = validate_and_normalize_samplesheet(src, out)

    assert result["sample_count"] == 2
    assert result["sample_names"] == ["PE", "SE"]
    rows_out = list(csv.DictReader(out.open(encoding="utf-8")))
    assert rows_out[0]["fastq_2"] == str(pe_r2.resolve())
    assert rows_out[1]["fastq_2"] == ""


def test_bam_reprocessing_preserves_official_samplesheet_fastq_columns(tmp_path):
    gbam = _touch(tmp_path / "sample.genome.bam")
    r1 = _touch(tmp_path / "reads" / "sample_R1.fastq.gz")
    src = _write_sheet(
        tmp_path / "samplesheet.csv",
        [{"sample": "sampleA", "fastq_1": str(r1), "strandedness": "auto", "genome_bam": str(gbam)}],
    )
    out = tmp_path / "normalized.csv"
    result = validate_and_normalize_samplesheet(src, out, skip_alignment=True)
    import csv as _csv
    row = next(_csv.DictReader(out.open(encoding="utf-8")))
    assert row["fastq_1"] == str(r1.resolve())
    assert result["fastq_paths"] == [r1.resolve()]
    assert result["bam_paths"] == [gbam.resolve()]


def test_preserves_remote_bam_uri_without_local_existence_check(tmp_path):
    bam_uri = "s3://bucket.example.org/bams/sample.genome.bam"
    r1 = _touch(tmp_path / "reads" / "sample_R1.fastq.gz")
    src = _write_sheet(
        tmp_path / "samplesheet.csv",
        [{"sample": "sampleA", "fastq_1": str(r1), "strandedness": "auto", "genome_bam": bam_uri}],
    )
    out = tmp_path / "normalized.csv"

    result = validate_and_normalize_samplesheet(src, out, skip_alignment=True)

    row = next(csv.DictReader(out.open(encoding="utf-8")))
    assert row["genome_bam"] == bam_uri
    assert result["bam_paths"] == [bam_uri]


def test_non_csv_extension_rejected(tmp_path):
    """Samplesheet with .tsv extension must be rejected — schema: input must be a .csv file."""
    tsv = tmp_path / "samplesheet.tsv"
    tsv.write_text("sample,fastq_1,strandedness\nsampleA,/reads/A_R1.fastq.gz,auto\n", encoding="utf-8")
    out = tmp_path / "normalized.csv"
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(tsv, out)
    assert exc.value.error_code == ErrorCode.INVALID_SAMPLESHEET


# Rules nf-core/rnaseq 3.26.0 enforces itself when Nextflow starts: nf-schema checks
# assets/schema_input.json (required columns, patterns, the strandedness enum, file
# existence, percent_mapped bounds) and checkSamplesAfterGrouping refuses runs of one
# sample with different strandedness. Checked live: nf-schema refused
# `strandedness: sideways`, `seq_center: My Center`, a missing FASTQ and
# `percent_mapped: 150`, naming each entry.
_NF_SCHEMA_CASES = {
    "missing strandedness column": dict(drop_column="strandedness"),
    "strandedness outside the enum": dict(row={"strandedness": "sideways"}),
    "whitespace in seq_center": dict(row={"seq_center": "My Center"}),
    "FASTQ with a .txt suffix": dict(fastq="sample_R1.txt"),
    "FASTQ basename with a space": dict(fastq="sample R1.fastq.gz"),
    "FASTQ that does not exist": dict(fastq="missing.fastq.gz", create=False),
    "runs of one sample with different strandedness": dict(second_run={"strandedness": "reverse"}),
    "BAM that does not exist": dict(row={"genome_bam": "missing.bam"}, skip_alignment=True),
    "percent_mapped above 100": dict(row={"genome_bam": "a.bam", "percent_mapped": "150"}, skip_alignment=True, bam=True),
    "empty strandedness": dict(row={"strandedness": ""}),
    "strandedness in the wrong case": dict(row={"strandedness": "AUTO"}),
}


@pytest.mark.parametrize("case", sorted(_NF_SCHEMA_CASES))
def test_rules_nf_schema_enforces_reach_it_unchanged(tmp_path, case):
    """The wrapper keeps no copy of upstream's samplesheet rules: such a row is
    normalized (paths resolved) and left for nf-schema, which names the entry."""
    spec = _NF_SCHEMA_CASES[case]
    fastq = tmp_path / spec.get("fastq", "sample_R1.fastq.gz")
    if spec.get("create", True):
        _touch(fastq)
    if spec.get("bam"):
        _touch(tmp_path / "a.bam")
    row = {"sample": "sampleA", "fastq_1": fastq.name, "strandedness": "forward", **spec.get("row", {})}
    rows = [row]
    if "second_run" in spec:
        other = _touch(tmp_path / "sample_L002_R1.fastq.gz")
        rows.append({**row, "fastq_1": other.name, **spec["second_run"]})
    fieldnames = [name for name in dict.fromkeys(key for r in rows for key in r) if name != spec.get("drop_column")]
    src = _write_sheet(tmp_path / "samplesheet.csv", [{k: r[k] for k in fieldnames} for r in rows], fieldnames=fieldnames)
    result = validate_and_normalize_samplesheet(
        src, tmp_path / "normalized.csv", skip_alignment=spec.get("skip_alignment", False)
    )
    written = list(csv.DictReader((tmp_path / "normalized.csv").read_text().splitlines()))
    assert result["sample_count"] == 1 and len(written) == len(rows)
    assert written[0]["fastq_1"] == (tmp_path / fastq.name).resolve().as_posix()
    for key, value in spec.get("row", {}).items():
        expected = (tmp_path / value).resolve().as_posix() if key == "genome_bam" else value
        assert written[0][key] == expected


@pytest.mark.parametrize("fastq_1", ["", "   "])
def test_a_row_without_fastq_1_is_left_for_nf_schema(tmp_path, fastq_1):
    """fastq_1 is required by schema_input.json, BAM reprocessing rows included;
    nf-schema names the entry that lacks it."""
    bam = _touch(tmp_path / "a.bam")
    src = _write_sheet(tmp_path / "samplesheet.csv",
                       [{"sample": "sampleA", "fastq_1": fastq_1, "strandedness": "auto", "genome_bam": bam.name}])
    validate_and_normalize_samplesheet(src, tmp_path / "normalized.csv", skip_alignment=True)
    written = list(csv.DictReader((tmp_path / "normalized.csv").read_text().splitlines()))
    assert written[0]["fastq_1"] == "" and written[0]["genome_bam"] == bam.resolve().as_posix()
