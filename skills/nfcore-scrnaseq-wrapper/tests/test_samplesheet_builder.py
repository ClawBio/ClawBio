from __future__ import annotations

import csv
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from hypothesis import HealthCheck, example, given, settings
from hypothesis import strategies as st

from errors import ErrorCode, SkillError
from samplesheet_builder import validate_and_normalize_samplesheet


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


def test_validate_accepts_remote_fastq_urls_verbatim(tmp_path):
    """Remote FASTQ URIs are passed through verbatim (parity with nf-core +
    nfcore-sarek/rnaseq): existence is deferred to Nextflow staging, only the
    basename is validated, and the URI is written unchanged to the normalized
    samplesheet (no Path() mangling of the `//` scheme separator)."""
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        "sample,fastq_1,fastq_2\n"
        "sampleA,https://example.org/S1_R1.fastq.gz,s3://bucket/S1_R2.fastq.gz\n",
        encoding="utf-8",
    )
    out = tmp_path / "normalized.csv"

    summary = validate_and_normalize_samplesheet(src, out)

    assert "https://example.org/S1_R1.fastq.gz" in summary["fastq_paths"]
    assert "s3://bucket/S1_R2.fastq.gz" in summary["fastq_paths"]
    written = out.read_text(encoding="utf-8")
    assert "https://example.org/S1_R1.fastq.gz" in written
    assert "s3://bucket/S1_R2.fastq.gz" in written


def test_validate_accepts_required_columns_in_any_order_and_normalizes(tmp_path):
    """nf-core validates the samplesheet by column NAME, not position. A header
    with the required columns present but in a different order is valid and must
    be accepted; the normalized output always leads with sample,fastq_1,fastq_2."""
    r1 = tmp_path / "a_R1.fastq.gz"
    r2 = tmp_path / "a_R2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"fastq_1,sample,fastq_2\n{r1},sampleA,{r2}\n",
        encoding="utf-8",
    )
    out = tmp_path / "normalized.csv"
    result = validate_and_normalize_samplesheet(src, out)

    assert result["sample_count"] == 1
    assert result["sample_names"] == ["sampleA"]
    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert list(rows[0].keys())[:3] == ["sample", "fastq_1", "fastq_2"]
    assert rows[0]["sample"] == "sampleA"


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
    assert (
        "duplicate" in exc.value.message.lower()
        or "duplicate" in str(exc.value.details).lower()
    )


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


def test_expected_cells_zero_rejected(tmp_path):
    r1 = tmp_path / "r1.fastq.gz"
    r2 = tmp_path / "r2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2,expected_cells\nsampleA,{r1},{r2},0\n",
        encoding="utf-8",
    )
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "norm.csv")
    assert exc.value.error_code == "INVALID_SAMPLESHEET"


def test_expected_cells_override_rejected_for_multiple_samples(tmp_path):
    r1a = tmp_path / "a_R1.fastq.gz"
    r2a = tmp_path / "a_R2.fastq.gz"
    r1b = tmp_path / "b_R1.fastq.gz"
    r2b = tmp_path / "b_R2.fastq.gz"
    for path in (r1a, r2a, r1b, r2b):
        path.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2\nsampleA,{r1a},{r2a}\nsampleB,{r1b},{r2b}\n",
        encoding="utf-8",
    )

    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(
            src, tmp_path / "norm.csv", expected_cells_override=1000
        )

    assert exc.value.error_code == "INVALID_SAMPLESHEET"
    assert exc.value.details["field"] == "expected_cells_override"
    assert exc.value.details["sample_count"] == 2


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
        validate_and_normalize_samplesheet(
            src, tmp_path / "normalized.csv", preset="cellrangerarc"
        )

    assert exc.value.error_code == ErrorCode.INVALID_FASTQ
    assert exc.value.details["preset"] == "cellrangerarc"


def test_cellrangerarc_accepts_10x_naming_with_fq_gz_extension(tmp_path):
    """The upstream samplesheet schema allows both .fastq.gz and .fq.gz; the ARC
    10x-naming check must accept .fq.gz too rather than rejecting a valid file by
    extension alone (audit F-4)."""
    r1 = tmp_path / "sampleA_S1_L001_R1_001.fq.gz"
    r2 = tmp_path / "sampleA_S1_L001_R2_001.fq.gz"
    barcode = tmp_path / "sampleA_S1_L001_I2_001.fq.gz"
    for fastq in (r1, r2, barcode):
        _touch_fastq(fastq)
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2,fastq_barcode,sample_type\n"
        f"sampleA,{r1},{r2},{barcode},atac\n",
        encoding="utf-8",
    )

    result = validate_and_normalize_samplesheet(
        src, tmp_path / "normalized.csv", preset="cellrangerarc"
    )

    assert result["sample_count"] == 1


def test_cellranger_family_warns_when_expected_cells_supplied(tmp_path, capsys):
    """nf-core docs advise against expected_cells since cellranger v7; the wrapper
    surfaces that advisory (stdout) for cellranger-family presets, without failing."""
    r1 = tmp_path / "r1.fastq.gz"
    r2 = tmp_path / "r2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2,feature_type,expected_cells\nsampleA,{r1},{r2},gex,5000\n",
        encoding="utf-8",
    )
    validate_and_normalize_samplesheet(
        src, tmp_path / "norm.csv", preset="cellrangermulti"
    )
    out = capsys.readouterr().out
    assert "expected_cells" in out and "Cell Ranger v7" in out


def test_non_cellranger_preset_does_not_warn_about_expected_cells(tmp_path, capsys):
    r1 = tmp_path / "r1.fastq.gz"
    r2 = tmp_path / "r2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2,expected_cells\nsampleA,{r1},{r2},5000\n",
        encoding="utf-8",
    )
    validate_and_normalize_samplesheet(src, tmp_path / "norm.csv", preset="standard")
    assert "Cell Ranger v7" not in capsys.readouterr().out


def test_protocol_column_is_reported_and_omitted_from_upstream_sheet(tmp_path):
    """Protocol is a global 4.1.0 parameter, not an input-schema column."""
    r1 = tmp_path / "a_R1.fastq.gz"
    r2 = tmp_path / "a_R2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    src = tmp_path / "samplesheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2,protocol,expected_cells\nsampleA,{r1},{r2},10XV2,1000\n",
        encoding="utf-8",
    )
    normalized = tmp_path / "normalized.csv"
    result = validate_and_normalize_samplesheet(src, normalized)
    assert "protocol" in result["unknown_columns"]
    assert "protocol" not in normalized.read_text(encoding="utf-8").splitlines()[0]


def test_unknown_metadata_columns_are_omitted_from_upstream_sheet(tmp_path):
    """The reproducibility sheet is an upstream input, not a metadata archive."""
    r1 = tmp_path / "S_R1.fastq.gz"
    r2 = tmp_path / "S_R2.fastq.gz"
    r1.write_bytes(b"")
    r2.write_bytes(b"")
    sheet = tmp_path / "samples.csv"
    sheet.write_text(
        "sample,fastq_1,fastq_2,cohort\n"
        f"S,{r1.name},{r2.name},discovery\n",
        encoding="utf-8",
    )
    normalized = tmp_path / "normalized.csv"

    result = validate_and_normalize_samplesheet(sheet, normalized)

    assert result["unknown_columns"] == ["cohort"]
    assert "cohort" not in normalized.read_text(encoding="utf-8").splitlines()[0]


def test_multi_lane_sample_counted_once(tmp_path):
    """A sample split across lanes is ONE sample, not one per row.

    nf-core/scrnaseq documents re-sequencing the same sample across lanes by repeating
    the `sample` identifier on multiple rows: "The sample identifiers have to be the same
    when you have re-sequenced the same sample more than once", and the pipeline
    "concatenate[s] the raw reads before performing any downstream analysis".

    The summary counted raw CSV rows (`len(normalized_rows)`) and appended to
    `sample_names` with no duplicate guard, so a 2-sample sheet with one sample split
    over two lanes reported 3 samples and listed the split sample twice. That count is
    not merely cosmetic: it is written into the provenance bundle (inputs.json
    `sample_count`) and is preferred by report.md over the samples actually detected in
    the outputs — so the audit trail claimed more samples than the run had.

    nfcore-rnaseq-wrapper already deduplicates (`if sample not in sample_names`) and
    counts `len(sample_names)`; this pins the same contract here.
    """
    lane1_r1 = tmp_path / "s1_L001_R1.fastq.gz"
    lane1_r2 = tmp_path / "s1_L001_R2.fastq.gz"
    lane2_r1 = tmp_path / "s1_L002_R1.fastq.gz"
    lane2_r2 = tmp_path / "s1_L002_R2.fastq.gz"
    other_r1 = tmp_path / "s2_R1.fastq.gz"
    other_r2 = tmp_path / "s2_R2.fastq.gz"
    for f in (lane1_r1, lane1_r2, lane2_r1, lane2_r2, other_r1, other_r2):
        _touch_fastq(f)

    src = tmp_path / "samplesheet.csv"
    src.write_text(
        "sample,fastq_1,fastq_2\n"
        f"sampleA,{lane1_r1},{lane1_r2}\n"
        f"sampleA,{lane2_r1},{lane2_r2}\n"
        f"sampleB,{other_r1},{other_r2}\n",
        encoding="utf-8",
    )
    out = tmp_path / "normalized.csv"
    result = validate_and_normalize_samplesheet(src, out)

    assert result["sample_count"] == 2, "a lane-split sample must count once, not per row"
    assert result["sample_names"] == ["sampleA", "sampleB"], "sample_names must not repeat"

    # All three rows must still reach the pipeline: nf-core concatenates them.
    rows = list(csv.DictReader(out.read_text(encoding="utf-8").splitlines()))
    assert len(rows) == 3, "every lane row must be preserved for the pipeline to concatenate"
    # And every FASTQ is still staged.
    assert len(result["fastq_paths"]) == 6



# Rules nf-core/scrnaseq 4.1.0 enforces itself when Nextflow starts, before any task
# runs: nf-schema checks assets/schema_input.json (required sample, fastq_1 and fastq_2,
# the FASTQ pattern and existence, file-path refusing a directory, the integer type of
# expected_cells, the sample_type and feature_type enums), and PIPELINE_INITIALISATION
# refuses a cellrangerarc row whose sample_type is not gex/atac or an ATAC row without
# fastq_barcode. The wrapper keeps no copy: the row is written through (paths resolved,
# sample names normalized) for the pipeline to name.
_NF_CORE_CASES = {
    "missing fastq_2 column": dict(drop="fastq_2"),
    "empty sample": dict(row={"sample": ""}),
    "FASTQ that does not exist": dict(row={"fastq_1": "missing_R1.fastq.gz"}),
    "FASTQ that is a directory": dict(row={"fastq_1": "reads_R1.fastq.gz"}, directory="reads_R1.fastq.gz"),
    "FASTQ basename with a space": dict(row={"fastq_1": "sample A_R1.fastq.gz"}, files=("sample A_R1.fastq.gz",)),
    "FASTQ extension in upper case": dict(row={"fastq_1": "s_R1.FASTQ.GZ"}, files=("s_R1.FASTQ.GZ",)),
    "remote FASTQ that is not a FASTQ": dict(row={"fastq_2": "https://example.org/notes.txt"}),
    "expected_cells that is not an integer": dict(row={"expected_cells": "many"}),
    "sample_type outside the enum": dict(row={"sample_type": "beam"}),
    "feature_type outside the enum": dict(row={"feature_type": "beam"}),
    "cellrangerarc row without sample_type": dict(preset="cellrangerarc", arc=True),
    "cellrangerarc ATAC row without fastq_barcode": dict(preset="cellrangerarc", arc=True, row={"sample_type": "atac"}),
}


@pytest.mark.parametrize("case", sorted(_NF_CORE_CASES))
def test_rules_nf_core_enforces_reach_it_unchanged(tmp_path, case):
    spec = _NF_CORE_CASES[case]
    names = ("S1_S1_L001_R1_001.fastq.gz", "S1_S1_L001_R2_001.fastq.gz") if spec.get("arc") else ("s_R1.fastq.gz", "s_R2.fastq.gz")
    for name in (*names, *spec.get("files", ())):
        _touch_fastq(tmp_path / name)
    if "directory" in spec:
        (tmp_path / spec["directory"]).mkdir()
    row = {"sample": "sampleA", "fastq_1": names[0], "fastq_2": names[1], **spec.get("row", {})}
    row = {k: v for k, v in row.items() if k != spec.get("drop")}
    src = tmp_path / "samplesheet.csv"
    src.write_text(",".join(row) + "\n" + ",".join(row.values()) + "\n", encoding="utf-8")
    out = tmp_path / "normalized.csv"

    validate_and_normalize_samplesheet(src, out, preset=spec.get("preset"))

    written = list(csv.DictReader(out.read_text(encoding="utf-8").splitlines()))
    assert len(written) == 1
    for column, value in row.items():
        if column.startswith("fastq") and value and "://" not in value:
            assert written[0][column] == (tmp_path / value).resolve().as_posix()
        elif column in written[0]:
            assert written[0][column] == value, column


def test_cellrangermulti_rows_still_need_a_feature_type(tmp_path):
    """align_cellrangermulti routes rows by feature_type with no default branch, so
    nf-core/scrnaseq 4.1.0 silently drops a row without one; the wrapper refuses it."""
    for name in ("r1.fastq.gz", "r2.fastq.gz"):
        _touch_fastq(tmp_path / name)
    src = tmp_path / "samplesheet.csv"
    src.write_text("sample,fastq_1,fastq_2\nsampleA,r1.fastq.gz,r2.fastq.gz\n", encoding="utf-8")
    with pytest.raises(SkillError) as exc:
        validate_and_normalize_samplesheet(src, tmp_path / "norm.csv", preset="cellrangermulti")
    assert exc.value.error_code == ErrorCode.INVALID_SAMPLESHEET
    assert exc.value.details["column"] == "feature_type"


@settings(suppress_health_check=[HealthCheck.function_scoped_fixture], max_examples=60, deadline=None)
@given(st.one_of(st.integers(-5, 5), st.text(alphabet="0123456789.-e x", max_size=6)))
@example("0")
@example("1.5")
@example(" 7 ")
def test_expected_cells_is_refused_only_as_an_integer_below_one(tmp_path, value):
    """The schema checks only that expected_cells is an integer; the wrapper refuses an
    integer below 1 (an aligner cannot expect no cells) and passes anything else on."""
    text = str(value)
    if "," in text:
        return
    for name in ("r1.fastq.gz", "r2.fastq.gz"):
        _touch_fastq(tmp_path / name)
    src = tmp_path / "samplesheet.csv"
    src.write_text(f"sample,fastq_1,fastq_2,expected_cells\nsampleA,r1.fastq.gz,r2.fastq.gz,{text}\n", encoding="utf-8")
    out = tmp_path / "norm.csv"
    try:
        number = int(text.strip())
    except ValueError:
        number = None
    if number is not None and number < 1:
        with pytest.raises(SkillError) as exc:
            validate_and_normalize_samplesheet(src, out)
        assert exc.value.error_code == ErrorCode.INVALID_SAMPLESHEET
    else:
        validate_and_normalize_samplesheet(src, out)
        written = next(csv.DictReader(out.read_text(encoding="utf-8").splitlines()))
        assert written["expected_cells"] == text.strip()
