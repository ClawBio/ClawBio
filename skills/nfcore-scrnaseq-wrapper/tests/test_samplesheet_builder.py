from __future__ import annotations

import csv
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from errors import ErrorCode, SkillError
from samplesheet_builder import (
    _looks_like_fastq_path,
    validate_and_normalize_samplesheet,
)


def _touch_fastq(path: Path) -> None:
    path.write_text("x", encoding="utf-8")


def test_validate_and_normalize_samplesheet(tmp_path):
    r1 = tmp_path / "a_R1.fastq.gz"
    r2 = tmp_path / "a_R2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2,expected_cells\nsampleA,{r1},{r2},1000\n",
        encoding="utf-8",
    )
    out = tmp_path / "normalized.csv"
    result = validate_and_normalize_samplesheet(src, out)
    assert result["sample_count"] == 1
    assert out.exists()


def test_validate_and_normalize_samplesheet_rejects_missing_fastq(tmp_path):
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        "sample,fastq_1,fastq_2\nsampleA,missing_R1.fastq.gz,missing_R2.fastq.gz\n",
        encoding="utf-8",
    )
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "normalized.csv")
    assert exc.value.error_code == "MISSING_FASTQ"


def test_validate_rejects_fastq_basename_with_space(tmp_path):
    r1 = tmp_path / "sample A_R1.fastq.gz"
    r2 = tmp_path / "sampleA_R2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2\nsampleA,{r1},{r2}\n",
        encoding="utf-8",
    )
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "normalized.csv")
    assert exc.value.error_code == "INVALID_FASTQ"
    assert exc.value.details["column"] == "fastq_1"


def test_validate_accepts_upstream_sample_name_with_non_whitespace_symbols(tmp_path):
    r1 = tmp_path / "sample_R1.fastq.gz"
    r2 = tmp_path / "sample_R2.fastq.gz"
    _touch_fastq(r1)
    _touch_fastq(r2)
    src = tmp_path / "samplesheet.csv"
    out = tmp_path / "normalized.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2\nsample+alpha:1,{r1},{r2}\n",
        encoding="utf-8",
    )

    summary = validate_and_normalize_samplesheet(src, out, preset="star")

    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert summary["sample_names"] == ["sample+alpha:1"]
    assert rows[0]["sample"] == "sample+alpha:1"


def test_validate_rejects_fastq_extension_case_not_matching_nfcore_schema(tmp_path):
    r1 = tmp_path / "sample_R1.FASTQ.GZ"
    r2 = tmp_path / "sample_R2.fastq.gz"
    _touch_fastq(r1)
    _touch_fastq(r2)
    src = tmp_path / "samplesheet.csv"
    src.write_text(f"sample,fastq_1,fastq_2\nsampleA,{r1},{r2}\n", encoding="utf-8")

    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "normalized.csv", preset="star")

    assert exc.value.error_code == ErrorCode.INVALID_FASTQ
    assert exc.value.details["filename"] == "sample_R1.FASTQ.GZ"


def test_validate_allows_parent_directory_with_space(tmp_path):
    reads_dir = tmp_path / "reads with space"
    reads_dir.mkdir()
    r1 = reads_dir / "sampleA_R1.fastq.gz"
    r2 = reads_dir / "sampleA_R2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2\nsampleA,{r1},{r2}\n",
        encoding="utf-8",
    )
    out = tmp_path / "normalized.csv"
    result = validate_and_normalize_samplesheet(src, out)
    assert result["sample_count"] == 1


def test_validate_rejects_fastq_barcode_basename_with_space_for_arc(tmp_path):
    r1 = tmp_path / "sampleA_R1.fastq.gz"
    r2 = tmp_path / "sampleA_R2.fastq.gz"
    barcode = tmp_path / "sample A_I2.fastq.gz"
    for path in (r1, r2, barcode):
        path.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2,sample_type,fastq_barcode\nsampleA,{r1},{r2},atac,{barcode}\n",
        encoding="utf-8",
    )
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "normalized.csv", preset="cellrangerarc")
    assert exc.value.error_code == "INVALID_FASTQ"
    assert exc.value.details["column"] == "fastq_barcode"


def test_validate_rejects_required_columns_out_of_order(tmp_path):
    r1 = tmp_path / "a_R1.fastq.gz"
    r2 = tmp_path / "a_R2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"fastq_1,sample,fastq_2\n{r1},sampleA,{r2}\n",
        encoding="utf-8",
    )
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "normalized.csv")
    assert exc.value.error_code == "INVALID_SAMPLESHEET"
    assert exc.value.details["expected_first_columns"] == ["sample", "fastq_1", "fastq_2"]


def test_validate_normalizes_sample_name_spaces_to_underscores(tmp_path):
    r1 = tmp_path / "a_R1.fastq.gz"
    r2 = tmp_path / "a_R2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2\nsample A,{r1},{r2}\n",
        encoding="utf-8",
    )
    out = tmp_path / "normalized.csv"
    result = validate_and_normalize_samplesheet(src, out)

    import csv as _csv
    row = next(_csv.DictReader(out.open(encoding="utf-8")))
    assert row["sample"] == "sample_A"
    assert result["sample_names"] == ["sample_A"]


def test_validate_normalizes_sample_name_whitespace_runs_to_single_underscore(tmp_path):
    r1 = tmp_path / "a_R1.fastq.gz"
    r2 = tmp_path / "a_R2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2\nsample\t  A,{r1},{r2}\n",
        encoding="utf-8",
    )
    out = tmp_path / "normalized.csv"
    result = validate_and_normalize_samplesheet(src, out)

    import csv as _csv
    row = next(_csv.DictReader(out.open(encoding="utf-8")))
    assert row["sample"] == "sample_A"
    assert result["sample_names"] == ["sample_A"]


def test_validate_rejects_sample_names_that_collide_after_space_normalization(tmp_path):
    r1a = tmp_path / "a_R1.fastq.gz"
    r2a = tmp_path / "a_R2.fastq.gz"
    r1b = tmp_path / "b_R1.fastq.gz"
    r2b = tmp_path / "b_R2.fastq.gz"
    for path in (r1a, r2a, r1b, r2b):
        path.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2\nsample A,{r1a},{r2a}\nsample_A,{r1b},{r2b}\n",
        encoding="utf-8",
    )
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "normalized.csv")
    assert exc.value.error_code == "INVALID_SAMPLESHEET"
    assert exc.value.details["normalized_sample"] == "sample_A"


def test_validate_rejects_sample_names_that_collide_after_whitespace_normalization(tmp_path):
    r1a = tmp_path / "a_R1.fastq.gz"
    r2a = tmp_path / "a_R2.fastq.gz"
    r1b = tmp_path / "b_R1.fastq.gz"
    r2b = tmp_path / "b_R2.fastq.gz"
    for path in (r1a, r2a, r1b, r2b):
        path.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2\nsample\tA,{r1a},{r2a}\nsample_A,{r1b},{r2b}\n",
        encoding="utf-8",
    )
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "normalized.csv")
    assert exc.value.error_code == "INVALID_SAMPLESHEET"
    assert exc.value.details["normalized_sample"] == "sample_A"


def test_validate_accepts_upstream_sample_name_with_special_chars(tmp_path):
    r1 = tmp_path / "a_R1.fastq.gz"
    r2 = tmp_path / "a_R2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2\nsample@bad!,{r1},{r2}\n",
        encoding="utf-8",
    )
    result = validate_and_normalize_samplesheet(src, tmp_path / "normalized.csv")
    assert result["sample_names"] == ["sample@bad!"]


def test_validate_accepts_valid_sample_name_formats(tmp_path):
    r1 = tmp_path / "a_R1.fastq.gz"
    r2 = tmp_path / "a_R2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    for valid_name in ("sample1", "Sample_A", "sample.1", "sample-1"):
        src = tmp_path / "samplesheet.csv"
        src.write_text(
            f"sample,fastq_1,fastq_2\n{valid_name},{r1},{r2}\n",
            encoding="utf-8",
        )
        result = validate_and_normalize_samplesheet(src, tmp_path / "normalized.csv")
        assert result["sample_count"] == 1


def test_normalized_csv_writes_absolute_posix_paths(tmp_path):
    """FASTQ paths in the normalized CSV must be absolute and use forward slashes."""
    r1 = tmp_path / "reads" / "s_R1.fastq.gz"
    r2 = tmp_path / "reads" / "s_R2.fastq.gz"
    r1.parent.mkdir()
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(f"sample,fastq_1,fastq_2\nsampleA,{r1},{r2}\n", encoding="utf-8")
    out = tmp_path / "norm.csv"
    validate_and_normalize_samplesheet(src, out)
    import csv as _csv
    rows = list(_csv.DictReader(out.open(encoding="utf-8")))
    assert Path(rows[0]["fastq_1"]).is_absolute()
    assert "\\" not in rows[0]["fastq_1"]
    assert "\\" not in rows[0]["fastq_2"]


def test_relative_fastq_path_resolved_against_csv_directory(tmp_path):
    """A relative FASTQ path should be resolved against the CSV's directory, not CWD."""
    reads_dir = tmp_path / "data"
    reads_dir.mkdir()
    r1 = reads_dir / "s_R1.fastq.gz"
    r2 = reads_dir / "s_R2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    # Write relative paths (relative to the CSV's location = tmp_path)
    src = tmp_path / "samplesheet.csv"
    src.write_text("sample,fastq_1,fastq_2\nsampleA,data/s_R1.fastq.gz,data/s_R2.fastq.gz\n", encoding="utf-8")
    out = tmp_path / "norm.csv"
    result = validate_and_normalize_samplesheet(src, out)
    assert result["sample_count"] == 1
    import csv as _csv
    rows = list(_csv.DictReader(out.open(encoding="utf-8")))
    assert rows[0]["fastq_1"] == r1.resolve().as_posix()


def test_exact_duplicate_fastq_rows_rejected(tmp_path):
    r1 = tmp_path / "r1.fastq.gz"
    r2 = tmp_path / "r2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2\nsampleA,{r1},{r2}\nsampleA,{r1},{r2}\n",
        encoding="utf-8",
    )
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "norm.csv")
    assert exc.value.error_code == "INVALID_SAMPLESHEET"
    assert "duplicate" in exc.value.message.lower() or "duplicate" in str(exc.value.details).lower()


def test_repeated_sample_names_allowed_for_distinct_fastqs(tmp_path):
    r1a = tmp_path / "lane1_R1.fastq.gz"
    r2a = tmp_path / "lane1_R2.fastq.gz"
    r1b = tmp_path / "lane2_R1.fastq.gz"
    r2b = tmp_path / "lane2_R2.fastq.gz"
    for path in (r1a, r2a, r1b, r2b):
        path.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2\nsampleA,{r1a},{r2a}\nsampleA,{r1b},{r2b}\n",
        encoding="utf-8",
    )
    result = validate_and_normalize_samplesheet(src, tmp_path / "norm.csv")
    assert result["sample_count"] == 2


def test_repeated_sample_requires_consistent_expected_cells(tmp_path):
    r1a = tmp_path / "lane1_R1.fastq.gz"
    r2a = tmp_path / "lane1_R2.fastq.gz"
    r1b = tmp_path / "lane2_R1.fastq.gz"
    r2b = tmp_path / "lane2_R2.fastq.gz"
    for path in (r1a, r2a, r1b, r2b):
        path.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        "sample,fastq_1,fastq_2,expected_cells\n"
        f"sampleA,{r1a},{r2a},1000\n"
        f"sampleA,{r1b},{r2b},2000\n",
        encoding="utf-8",
    )

    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "norm.csv")

    assert exc.value.error_code == "INVALID_SAMPLESHEET"
    assert exc.value.details["column"] == "expected_cells"
    assert exc.value.details["sample"] == "sampleA"


def test_repeated_sample_requires_consistent_seq_center(tmp_path):
    r1a = tmp_path / "lane1_R1.fastq.gz"
    r2a = tmp_path / "lane1_R2.fastq.gz"
    r1b = tmp_path / "lane2_R1.fastq.gz"
    r2b = tmp_path / "lane2_R2.fastq.gz"
    for path in (r1a, r2a, r1b, r2b):
        path.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        "sample,fastq_1,fastq_2,seq_center\n"
        f"sampleA,{r1a},{r2a},CoreA\n"
        f"sampleA,{r1b},{r2b},CoreB\n",
        encoding="utf-8",
    )

    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "norm.csv")

    assert exc.value.error_code == "INVALID_SAMPLESHEET"
    assert exc.value.details["column"] == "seq_center"
    assert exc.value.details["sample"] == "sampleA"


def test_repeated_sample_allows_matching_metadata(tmp_path):
    r1a = tmp_path / "lane1_R1.fastq.gz"
    r2a = tmp_path / "lane1_R2.fastq.gz"
    r1b = tmp_path / "lane2_R1.fastq.gz"
    r2b = tmp_path / "lane2_R2.fastq.gz"
    for path in (r1a, r2a, r1b, r2b):
        path.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        "sample,fastq_1,fastq_2,expected_cells,seq_center\n"
        f"sampleA,{r1a},{r2a},1000,CoreA\n"
        f"sampleA,{r1b},{r2b},1000,CoreA\n",
        encoding="utf-8",
    )

    result = validate_and_normalize_samplesheet(src, tmp_path / "norm.csv")

    assert result["sample_count"] == 2


def test_expected_cells_zero_rejected(tmp_path):
    r1 = tmp_path / "r1.fastq.gz"
    r2 = tmp_path / "r2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(f"sample,fastq_1,fastq_2,expected_cells\nsampleA,{r1},{r2},0\n", encoding="utf-8")
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "norm.csv")
    assert exc.value.error_code == "INVALID_SAMPLESHEET"


def test_expected_cells_negative_rejected(tmp_path):
    r1 = tmp_path / "r1.fastq.gz"
    r2 = tmp_path / "r2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(f"sample,fastq_1,fastq_2,expected_cells\nsampleA,{r1},{r2},-500\n", encoding="utf-8")
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "norm.csv")
    assert exc.value.error_code == "INVALID_SAMPLESHEET"


def test_missing_fastq_2_rejected(tmp_path):
    """A row with fastq_1 but empty fastq_2 must raise INVALID_SAMPLESHEET."""
    r1 = tmp_path / "r1.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(f"sample,fastq_1,fastq_2\nsampleA,{r1},\n", encoding="utf-8")
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "norm.csv")
    assert exc.value.error_code == "INVALID_SAMPLESHEET"


def test_bom_prefixed_samplesheet_accepted(tmp_path):
    """UTF-8 BOM (Excel export) must not cause a MISSING_COLUMNS error."""
    r1 = tmp_path / "r1.fastq.gz"
    r2 = tmp_path / "r2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    # Write BOM + CSV header explicitly
    src.write_bytes(
        "﻿sample,fastq_1,fastq_2\n".encode("utf-8")
        + f"sampleA,{r1},{r2}\n".encode("utf-8")
    )
    result = validate_and_normalize_samplesheet(src, tmp_path / "norm.csv")
    assert result["sample_count"] == 1


def test_unknown_columns_reported(tmp_path):
    """Extra columns in the samplesheet must be surfaced in the result dict."""
    r1 = tmp_path / "r1.fastq.gz"
    r2 = tmp_path / "r2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(f"sample,fastq_1,fastq_2,library_type\nsampleA,{r1},{r2},GEX\n", encoding="utf-8")
    result = validate_and_normalize_samplesheet(src, tmp_path / "norm.csv")
    assert "library_type" in result["unknown_columns"]


def test_unknown_columns_preserved_in_normalized_samplesheet(tmp_path):
    r1 = tmp_path / "sampleA_gex_S1_L001_R1_001.fastq.gz"
    r2 = tmp_path / "sampleA_gex_S1_L001_R2_001.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2,sample_type,fastq_barcode,feature_type\n"
        f"sampleA,{r1},{r2},gex,SI-GA-A1,gex\n",
        encoding="utf-8",
    )
    out = tmp_path / "norm.csv"
    validate_and_normalize_samplesheet(src, out, preset="cellrangerarc")
    import csv as _csv
    row = next(_csv.DictReader(out.open(encoding="utf-8")))
    assert row["sample_type"] == "gex"
    assert row["fastq_barcode"] == "SI-GA-A1"
    assert row["feature_type"] == "gex"


def test_cellrangerarc_requires_sample_type_and_fastq_barcode_columns(tmp_path):
    r1 = tmp_path / "r1.fastq.gz"
    r2 = tmp_path / "r2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(f"sample,fastq_1,fastq_2\nsampleA,{r1},{r2}\n", encoding="utf-8")
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "norm.csv", preset="cellrangerarc")
    assert exc.value.error_code == "INVALID_SAMPLESHEET"
    assert exc.value.details["missing_columns"] == ["sample_type", "fastq_barcode"]


def test_cellrangerarc_rejects_invalid_sample_type(tmp_path):
    r1 = tmp_path / "r1.fastq.gz"
    r2 = tmp_path / "r2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2,sample_type,fastq_barcode\nsampleA,{r1},{r2},rna,\n",
        encoding="utf-8",
    )
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "norm.csv", preset="cellrangerarc")
    assert exc.value.error_code == "INVALID_SAMPLESHEET"
    assert exc.value.details["sample_type"] == "rna"


def test_cellrangermulti_requires_feature_type_column(tmp_path):
    r1 = tmp_path / "r1.fastq.gz"
    r2 = tmp_path / "r2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(f"sample,fastq_1,fastq_2\nsampleA,{r1},{r2}\n", encoding="utf-8")
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "norm.csv", preset="cellrangermulti")
    assert exc.value.error_code == "INVALID_SAMPLESHEET"
    assert exc.value.details["missing_columns"] == ["feature_type"]


def test_cellrangermulti_rejects_invalid_feature_type(tmp_path):
    r1 = tmp_path / "r1.fastq.gz"
    r2 = tmp_path / "r2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2,feature_type\nsampleA,{r1},{r2},beam\n",
        encoding="utf-8",
    )
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "norm.csv", preset="cellrangermulti")
    assert exc.value.error_code == "INVALID_SAMPLESHEET"
    assert exc.value.details["feature_type"] == "beam"


def test_fastq_barcode_relative_path_is_normalized_for_arc(tmp_path):
    reads = tmp_path / "reads"
    reads.mkdir()
    r1 = reads / "sampleA_gex_S1_L001_R1_001.fastq.gz"
    r2 = reads / "sampleA_gex_S1_L001_R2_001.fastq.gz"
    barcode = reads / "sampleA_atac_S1_L001_I2_001.fastq.gz"
    for path in (r1, r2, barcode):
        path.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        "sample,fastq_1,fastq_2,sample_type,fastq_barcode\n"
        "sampleA,reads/sampleA_gex_S1_L001_R1_001.fastq.gz,"
        "reads/sampleA_gex_S1_L001_R2_001.fastq.gz,gex,reads/sampleA_atac_S1_L001_I2_001.fastq.gz\n",
        encoding="utf-8",
    )
    out = tmp_path / "norm.csv"
    validate_and_normalize_samplesheet(src, out, preset="cellrangerarc")
    import csv as _csv
    row = next(_csv.DictReader(out.open(encoding="utf-8")))
    assert row["fastq_barcode"] == barcode.resolve().as_posix()


def test_arc_atac_row_requires_fastq_barcode(tmp_path):
    r1 = tmp_path / "gex_R1.fastq.gz"
    r2 = tmp_path / "gex_R2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2,sample_type,fastq_barcode\nsampleA,{r1},{r2},atac,\n",
        encoding="utf-8",
    )
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "norm.csv", preset="cellrangerarc")
    assert exc.value.error_code == "INVALID_SAMPLESHEET"
    assert exc.value.details["column"] == "fastq_barcode"


def test_arc_atac_row_rejects_symbolic_fastq_barcode(tmp_path):
    r1 = tmp_path / "gex_R1.fastq.gz"
    r2 = tmp_path / "gex_R2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2,sample_type,fastq_barcode\nsampleA,{r1},{r2},atac,SI-GA-A1\n",
        encoding="utf-8",
    )
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "norm.csv", preset="cellrangerarc")
    assert exc.value.error_code == "MISSING_FASTQ"
    assert exc.value.details["column"] == "fastq_barcode"


def test_invalid_fastq_extension_rejected(tmp_path):
    r1 = tmp_path / "r1.txt"
    r2 = tmp_path / "r2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(f"sample,fastq_1,fastq_2\nsampleA,{r1},{r2}\n", encoding="utf-8")
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "norm.csv")
    assert exc.value.error_code == "INVALID_FASTQ"


def test_looks_like_fastq_rejects_path_with_slash_but_no_fastq_suffix():
    assert not _looks_like_fastq_path("/some/path/to/file.bam")
    assert not _looks_like_fastq_path("some/relative/path")
    assert not _looks_like_fastq_path("SI-GA-A1")


def test_looks_like_fastq_accepts_fastq_extensions():
    assert _looks_like_fastq_path("/data/sample_R1.fastq.gz")
    assert _looks_like_fastq_path("relative/sample.fq.gz")
    assert _looks_like_fastq_path("sample_1.fastq.gz")


def test_sample_name_accepts_leading_dot_allowed_by_nfcore_schema(tmp_path):
    r1 = tmp_path / "r1.fastq.gz"
    r2 = tmp_path / "r2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    ss = tmp_path / "samplesheet.csv"
    ss.write_text(f"sample,fastq_1,fastq_2\n.hidden,{r1},{r2}\n", encoding="utf-8")
    out = tmp_path / "out.csv"
    result = validate_and_normalize_samplesheet(ss, out)
    assert result["sample_names"] == [".hidden"]


def test_cellranger_rejects_fastq_pair_that_differs_by_more_than_read_marker(tmp_path):
    r1 = tmp_path / "sampleA_S1_L001_R1_001.fastq.gz"
    r2 = tmp_path / "other_S1_L001_R2_001.fastq.gz"
    _touch_fastq(r1)
    _touch_fastq(r2)
    src = tmp_path / "samplesheet.csv"
    src.write_text(f"sample,fastq_1,fastq_2\nsampleA,{r1},{r2}\n", encoding="utf-8")

    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "normalized.csv", preset="cellranger")

    assert exc.value.error_code == ErrorCode.INVALID_SAMPLESHEET
    assert exc.value.details["preset"] == "cellranger"
    assert exc.value.details["fastq_1"] == r1.name
    assert exc.value.details["fastq_2"] == r2.name


def test_cellranger_accepts_fastq_pair_that_only_differs_by_r1_r2(tmp_path):
    r1 = tmp_path / "sampleA_S1_L001_R1_001.fastq.gz"
    r2 = tmp_path / "sampleA_S1_L001_R2_001.fastq.gz"
    _touch_fastq(r1)
    _touch_fastq(r2)
    src = tmp_path / "samplesheet.csv"
    out = tmp_path / "normalized.csv"
    src.write_text(f"sample,fastq_1,fastq_2\nsampleA,{r1},{r2}\n", encoding="utf-8")

    validate_and_normalize_samplesheet(src, out, preset="cellranger")

    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert rows[0]["fastq_1"] == r1.as_posix()
    assert rows[0]["fastq_2"] == r2.as_posix()


def test_cellrangerarc_rejects_fastqs_not_following_10x_arc_naming(tmp_path):
    r1 = tmp_path / "arc_R1.fastq.gz"
    r2 = tmp_path / "arc_R2.fastq.gz"
    barcode = tmp_path / "arc_I2.fastq.gz"
    _touch_fastq(r1)
    _touch_fastq(r2)
    _touch_fastq(barcode)
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2,fastq_barcode,sample_type\n"
        f"sampleA,{r1},{r2},{barcode},atac\n",
        encoding="utf-8",
    )

    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "normalized.csv", preset="cellrangerarc")

    assert exc.value.error_code == ErrorCode.INVALID_FASTQ
    assert exc.value.details["preset"] == "cellrangerarc"


def test_cellrangerarc_accepts_10x_arc_naming_for_atac_and_gex_rows(tmp_path):
    atac_r1 = tmp_path / "sampleA_atac_S1_L001_R1_001.fastq.gz"
    atac_r2 = tmp_path / "sampleA_atac_S1_L001_R2_001.fastq.gz"
    atac_i2 = tmp_path / "sampleA_atac_S1_L001_I2_001.fastq.gz"
    gex_r1 = tmp_path / "sampleB_gex_S1_L001_R1_001.fastq.gz"
    gex_r2 = tmp_path / "sampleB_gex_S1_L001_R2_001.fastq.gz"
    for path in (atac_r1, atac_r2, atac_i2, gex_r1, gex_r2):
        _touch_fastq(path)
    src = tmp_path / "samplesheet.csv"
    out = tmp_path / "normalized.csv"
    src.write_text(
        "sample,fastq_1,fastq_2,fastq_barcode,sample_type\n"
        f"sampleA,{atac_r1},{atac_r2},{atac_i2},atac\n"
        f"sampleB,{gex_r1},{gex_r2},,gex\n",
        encoding="utf-8",
    )

    validate_and_normalize_samplesheet(src, out, preset="cellrangerarc")

    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert rows[0]["fastq_barcode"] == atac_i2.as_posix()
    assert rows[1]["sample_type"] == "gex"
