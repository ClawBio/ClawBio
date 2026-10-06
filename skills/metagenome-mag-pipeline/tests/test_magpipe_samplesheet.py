"""Samplesheet loading, discovery and writing.

Every rule here is upstream's own, from `read_samplesheet()` in
`workflow/rules/helpers.smk` at the pinned commit. The wrapper repeats them
because upstream only checks them once snakemake has started, and a cohort should
not be lost to a typo in column three.

The one deliberate difference: FASTQ paths resolve against the samplesheet's own
directory and are made absolute, because upstream's `resolve()` joins relative
paths to the pipeline root, not to the sheet.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

from magpipe_errors import ErrorCode, SkillError
from magpipe_samplesheet import (
    Sample,
    build_samplesheet,
    discover_samples,
    load_samplesheet,
    write_samplesheet,
)
from magpipe_schemas import SAMPLESHEET_COLUMNS

COLUMNS = list(SAMPLESHEET_COLUMNS)


def make_fastq(path: Path, records: int = 2) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        for index in range(records):
            handle.write(f"@read{index}\nACGTACGTAC\n+\nIIIIIIIIII\n")
    return path


def write_sheet(path: Path, rows: list[dict], columns: list[str] | None = None) -> Path:
    columns = columns or COLUMNS
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join(str(row.get(column, "") or "") for column in columns))
    path.write_text("\n".join(lines) + "\n")
    return path


def paired_row(tmp_path: Path, sample: str = "S1", lanes: int = 1) -> dict:
    first, second = [], []
    for lane in range(1, lanes + 1):
        suffix = f"_lane{lane}" if lanes > 1 else ""
        first.append(str(make_fastq(tmp_path / f"{sample}{suffix}_R1.fastq.gz")))
        second.append(str(make_fastq(tmp_path / f"{sample}{suffix}_R2.fastq.gz")))
    return {
        "sample": sample,
        "fastq_1": ";".join(first),
        "fastq_2": ";".join(second),
    }


@pytest.fixture
def paired_sheet(tmp_path: Path) -> Path:
    return write_sheet(tmp_path / "samplesheet.tsv", [paired_row(tmp_path, "S1", lanes=2)])


class TestLoadValid:
    def test_valid_paired(self, tmp_path: Path) -> None:
        sheet = write_sheet(tmp_path / "s.tsv", [paired_row(tmp_path)])
        (sample,) = load_samplesheet(sheet)
        assert sample.sample == "S1"
        assert sample.paired is True
        assert len(sample.fastq_1) == len(sample.fastq_2) == 1
        assert sample.group == "S1"
        assert sample.sra_run is None

    def test_valid_single_end(self, tmp_path: Path) -> None:
        row = {
            "sample": "S1",
            "fastq_1": str(make_fastq(tmp_path / "S1_R1.fastq.gz")),
            "layout": "single",
        }
        (sample,) = load_samplesheet(write_sheet(tmp_path / "s.tsv", [row]))
        assert sample.paired is False
        assert sample.fastq_2 == ()

    def test_paired_inferred_when_layout_absent(self, tmp_path: Path) -> None:
        (sample,) = load_samplesheet(write_sheet(tmp_path / "s.tsv", [paired_row(tmp_path)]))
        assert sample.paired is True

    def test_multi_lane_joined_in_order(self, tmp_path: Path) -> None:
        sheet = write_sheet(tmp_path / "s.tsv", [paired_row(tmp_path, lanes=3)])
        (sample,) = load_samplesheet(sheet)
        assert [p.name for p in sample.fastq_1] == [
            "S1_lane1_R1.fastq.gz",
            "S1_lane2_R1.fastq.gz",
            "S1_lane3_R1.fastq.gz",
        ]
        assert len(sample.fastq_2) == 3

    def test_group_taken_from_the_column(self, tmp_path: Path) -> None:
        row = paired_row(tmp_path)
        row["group"] = "GROUP_A"
        (sample,) = load_samplesheet(write_sheet(tmp_path / "s.tsv", [row]))
        assert sample.group == "GROUP_A"

    def test_group_defaults_to_the_sample(self, tmp_path: Path) -> None:
        row = paired_row(tmp_path, sample="S9")
        assert "group" not in row
        (sample,) = load_samplesheet(write_sheet(tmp_path / "s.tsv", [row]))
        assert sample.group == "S9"

    def test_upstream_demo_columns_are_accepted(self, tmp_path: Path) -> None:
        """The upstream generator writes sample/fastq_1/fastq_2/group/assembly."""
        columns = ["sample", "fastq_1", "fastq_2", "group", "assembly"]
        row = paired_row(tmp_path)
        row["assembly"] = "individual"
        (sample,) = load_samplesheet(write_sheet(tmp_path / "s.tsv", [row], columns))
        assert sample.paired is True

    def test_sra_run_row(self, tmp_path: Path) -> None:
        row = {"sample": "SRR1761673", "sra_run": "SRR1761673", "layout": "paired"}
        (sample,) = load_samplesheet(write_sheet(tmp_path / "sra.tsv", [row]))
        assert sample.sra_run == "SRR1761673"
        assert sample.fastq_1 == ()
        assert sample.paired is True

    def test_sra_run_single_layout(self, tmp_path: Path) -> None:
        row = {"sample": "SRR1761673", "sra_run": "SRR1761673", "layout": "single"}
        (sample,) = load_samplesheet(write_sheet(tmp_path / "sra.tsv", [row]))
        assert sample.paired is False


class TestLoadRejects:
    @pytest.mark.parametrize(
        ("mutate", "fragment"),
        [
            (lambda row: row.update(sample="bad sample"), "invalid sample identifier"),
            (lambda row: row.update(sample=""), "empty sample identifier"),
            (lambda row: row.update(group="bad group"), "invalid assembly group"),
            (lambda row: row.update(layout="PE"), "layout must be"),
            (lambda row: row.update(layout="paired", fastq_2=""), "layout=paired"),
            (lambda row: row.update(layout="single"), "layout=single"),
        ],
    )
    def test_field_rules(self, tmp_path: Path, mutate, fragment: str) -> None:
        row = paired_row(tmp_path)
        mutate(row)
        with pytest.raises(SkillError) as excinfo:
            load_samplesheet(write_sheet(tmp_path / "s.tsv", [row]))
        assert excinfo.value.error_code == ErrorCode.INVALID_SAMPLESHEET
        assert fragment in excinfo.value.message

    def test_duplicate_sample(self, tmp_path: Path) -> None:
        rows = [paired_row(tmp_path, "S1"), paired_row(tmp_path, "S2")]
        rows[1]["sample"] = "S1"
        with pytest.raises(SkillError) as excinfo:
            load_samplesheet(write_sheet(tmp_path / "s.tsv", rows))
        assert excinfo.value.error_code == ErrorCode.INVALID_SAMPLESHEET
        assert "duplicate" in excinfo.value.message

    def test_unequal_lane_counts(self, tmp_path: Path) -> None:
        row = paired_row(tmp_path, lanes=2)
        row["fastq_2"] = row["fastq_2"].split(";")[0]
        with pytest.raises(SkillError) as excinfo:
            load_samplesheet(write_sheet(tmp_path / "s.tsv", [row]))
        assert "unequal numbers of R1/R2 lanes" in excinfo.value.message

    def test_repeated_path(self, tmp_path: Path) -> None:
        path = str(make_fastq(tmp_path / "S1_R1.fastq.gz"))
        row = {"sample": "S1", "fastq_1": f"{path};{path}", "layout": "single"}
        with pytest.raises(SkillError) as excinfo:
            load_samplesheet(write_sheet(tmp_path / "s.tsv", [row]))
        assert "repeated FASTQ path" in excinfo.value.message

    def test_sra_run_with_paths(self, tmp_path: Path) -> None:
        row = paired_row(tmp_path)
        row["sra_run"] = "SRR1761673"
        row["layout"] = "paired"
        with pytest.raises(SkillError) as excinfo:
            load_samplesheet(write_sheet(tmp_path / "s.tsv", [row]))
        assert "mutually exclusive" in excinfo.value.message

    def test_sra_run_without_layout(self, tmp_path: Path) -> None:
        row = {"sample": "SRR1761673", "sra_run": "SRR1761673"}
        with pytest.raises(SkillError) as excinfo:
            load_samplesheet(write_sheet(tmp_path / "sra.tsv", [row]))
        assert "require an explicit layout" in excinfo.value.message

    def test_malformed_accession(self, tmp_path: Path) -> None:
        row = {"sample": "RUN1", "sra_run": "not-an-accession", "layout": "paired"}
        with pytest.raises(SkillError) as excinfo:
            load_samplesheet(write_sheet(tmp_path / "sra.tsv", [row]))
        assert "is not a run accession" in excinfo.value.message

    def test_missing_required_column(self, tmp_path: Path) -> None:
        columns = ["sample", "fastq_2", "group"]
        row = {"sample": "S1", "group": "S1"}
        with pytest.raises(SkillError) as excinfo:
            load_samplesheet(write_sheet(tmp_path / "s.tsv", [row], columns))
        assert "missing required column" in excinfo.value.message
        assert "fastq_1" in excinfo.value.message

    def test_empty_sheet(self, tmp_path: Path) -> None:
        with pytest.raises(SkillError) as excinfo:
            load_samplesheet(write_sheet(tmp_path / "s.tsv", []))
        assert excinfo.value.error_code == ErrorCode.INVALID_SAMPLESHEET
        assert "no samples" in excinfo.value.message

    def test_missing_fastq_file(self, tmp_path: Path) -> None:
        row = {"sample": "S1", "fastq_1": str(tmp_path / "gone_R1.fastq.gz")}
        with pytest.raises(SkillError) as excinfo:
            load_samplesheet(write_sheet(tmp_path / "s.tsv", [row]))
        assert excinfo.value.error_code == ErrorCode.MISSING_FASTQ
        assert "S1" in excinfo.value.message

    def test_missing_sheet_file(self, tmp_path: Path) -> None:
        with pytest.raises(SkillError) as excinfo:
            load_samplesheet(tmp_path / "nope.tsv")
        assert excinfo.value.error_code == ErrorCode.MISSING_INPUT


class TestPathResolution:
    def test_relative_paths_become_absolute(self, tmp_path: Path) -> None:
        make_fastq(tmp_path / "reads" / "S1_R1.fastq.gz")
        make_fastq(tmp_path / "reads" / "S1_R2.fastq.gz")
        row = {"sample": "S1", "fastq_1": "reads/S1_R1.fastq.gz", "fastq_2": "reads/S1_R2.fastq.gz"}
        (sample,) = load_samplesheet(write_sheet(tmp_path / "s.tsv", [row]))
        assert sample.fastq_1[0].is_absolute()
        assert sample.fastq_1[0] == (tmp_path / "reads" / "S1_R1.fastq.gz").resolve()

    def test_relative_paths_resolve_against_the_sheet_not_the_cwd(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        (tmp_path / "reads").mkdir()
        make_fastq(tmp_path / "reads" / "S1_R1.fastq.gz")
        sheet = write_sheet(
            tmp_path / "reads" / "s.tsv", [{"sample": "S1", "fastq_1": "S1_R1.fastq.gz"}]
        )
        elsewhere = tmp_path / "elsewhere"
        elsewhere.mkdir()
        monkeypatch.chdir(elsewhere)
        (sample,) = load_samplesheet(sheet)
        assert sample.fastq_1[0].parent == (tmp_path / "reads").resolve()


class TestDiscovery:
    def test_flat_r1_r2(self, tmp_path: Path) -> None:
        make_fastq(tmp_path / "S1_R1.fastq.gz")
        make_fastq(tmp_path / "S1_R2.fastq.gz")
        (sample,) = discover_samples(tmp_path)
        assert sample.sample == "S1" and sample.paired is True

    def test_flat_underscore_mates(self, tmp_path: Path) -> None:
        make_fastq(tmp_path / "TEST_001_1.fastq.gz")
        make_fastq(tmp_path / "TEST_001_2.fastq.gz")
        (sample,) = discover_samples(tmp_path)
        assert sample.sample == "TEST_001"
        assert [p.name for p in sample.fastq_1] == ["TEST_001_1.fastq.gz"]

    def test_flat_multi_lane(self, tmp_path: Path) -> None:
        for lane in (1, 2):
            make_fastq(tmp_path / f"S1_lane{lane}_R1.fastq.gz")
            make_fastq(tmp_path / f"S1_lane{lane}_R2.fastq.gz")
        (sample,) = discover_samples(tmp_path)
        assert sample.sample == "S1"
        assert [p.name for p in sample.fastq_1] == ["S1_lane1_R1.fastq.gz", "S1_lane2_R1.fastq.gz"]

    def test_one_folder_per_sample(self, tmp_path: Path) -> None:
        for name in ("S1", "S2"):
            make_fastq(tmp_path / name / f"{name}_R1.fastq.gz")
            make_fastq(tmp_path / name / f"{name}_R2.fastq.gz")
        samples = discover_samples(tmp_path)
        assert [s.sample for s in samples] == ["S1", "S2"]
        assert all(s.paired for s in samples)

    def test_flat_single_end(self, tmp_path: Path) -> None:
        make_fastq(tmp_path / "S1_R1.fastq.gz")
        (sample,) = discover_samples(tmp_path)
        assert sample.paired is False
        assert sample.fastq_2 == ()

    def test_unpaired_filename_without_mate_token(self, tmp_path: Path) -> None:
        make_fastq(tmp_path / "solo.fastq.gz")
        (sample,) = discover_samples(tmp_path)
        assert sample.sample == "solo" and sample.paired is False

    def test_no_reads(self, tmp_path: Path) -> None:
        (tmp_path / "empty").mkdir()
        with pytest.raises(SkillError) as excinfo:
            discover_samples(tmp_path / "empty")
        assert excinfo.value.error_code == ErrorCode.MISSING_INPUT
        assert "No FASTQ files" in excinfo.value.message

    def test_missing_directory(self, tmp_path: Path) -> None:
        with pytest.raises(SkillError) as excinfo:
            discover_samples(tmp_path / "nope")
        assert excinfo.value.error_code == ErrorCode.MISSING_INPUT

    def test_unequal_flat_mates(self, tmp_path: Path) -> None:
        make_fastq(tmp_path / "S1_lane1_R1.fastq.gz")
        make_fastq(tmp_path / "S1_lane2_R1.fastq.gz")
        make_fastq(tmp_path / "S1_lane1_R2.fastq.gz")
        with pytest.raises(SkillError) as excinfo:
            discover_samples(tmp_path)
        assert excinfo.value.error_code == ErrorCode.INVALID_SAMPLESHEET
        assert "R2 files" in excinfo.value.message

    def test_r2_without_r1(self, tmp_path: Path) -> None:
        make_fastq(tmp_path / "S1_R2.fastq.gz")
        with pytest.raises(SkillError) as excinfo:
            discover_samples(tmp_path)
        assert "without any R1" in excinfo.value.message

    def test_mixed_layouts_refused(self, tmp_path: Path) -> None:
        make_fastq(tmp_path / "loose_R1.fastq.gz")
        make_fastq(tmp_path / "S1" / "S1_R1.fastq.gz")
        with pytest.raises(SkillError) as excinfo:
            discover_samples(tmp_path)
        assert "mixes per-sample subdirectories" in excinfo.value.message

    def test_non_fastq_files_ignored(self, tmp_path: Path) -> None:
        make_fastq(tmp_path / "S1_R1.fastq.gz")
        (tmp_path / "checksums.md5").write_text("x\n")
        (sample,) = discover_samples(tmp_path)
        assert sample.sample == "S1"


class TestWriteAndRoundTrip:
    def test_writes_the_six_columns(self, tmp_path: Path) -> None:
        sample = Sample(
            sample="S1",
            fastq_1=(Path("/abs/S1_R1.fastq.gz"),),
            fastq_2=(Path("/abs/S1_R2.fastq.gz"),),
            group="G1",
            sra_run=None,
            paired=True,
        )
        dest = write_samplesheet([sample], tmp_path / "out" / "samplesheet.tsv")
        header, row = dest.read_text().splitlines()[:2]
        assert header.split("\t") == COLUMNS
        cells = row.split("\t")
        assert cells[0] == "S1" and cells[3] == "G1" and cells[5] == "paired"

    def test_lanes_joined_with_semicolons(self, tmp_path: Path) -> None:
        sample = Sample(
            sample="S1",
            fastq_1=(Path("/a/1_R1.fq.gz"), Path("/a/2_R1.fq.gz")),
            fastq_2=(Path("/a/1_R2.fq.gz"), Path("/a/2_R2.fq.gz")),
            group="S1",
            sra_run=None,
            paired=True,
        )
        dest = write_samplesheet([sample], tmp_path / "samplesheet.tsv")
        assert dest.read_text().splitlines()[1].split("\t")[1] == "/a/1_R1.fq.gz;/a/2_R1.fq.gz"

    def test_round_trip(self, tmp_path: Path) -> None:
        sheet = write_sheet(
            tmp_path / "in.tsv",
            [paired_row(tmp_path, "S1", lanes=2), paired_row(tmp_path, "S2", lanes=1)],
        )
        original = load_samplesheet(sheet)
        dest = write_samplesheet(original, tmp_path / "out.tsv")
        assert load_samplesheet(dest) == original

    def test_round_trip_preserves_groups(self, tmp_path: Path) -> None:
        first, second = paired_row(tmp_path, "S1"), paired_row(tmp_path, "S2")
        first["group"] = second["group"] = "GROUP_A"
        original = load_samplesheet(write_sheet(tmp_path / "in.tsv", [first, second]))
        dest = write_samplesheet(original, tmp_path / "out.tsv")
        assert [s.group for s in load_samplesheet(dest)] == ["GROUP_A", "GROUP_A"]


class TestBuildSamplesheet:
    def test_dispatches_to_a_file(self, tmp_path: Path) -> None:
        sheet = write_sheet(tmp_path / "in.tsv", [paired_row(tmp_path)])
        dest = tmp_path / "out" / "samplesheet.tsv"
        samples = build_samplesheet(sheet, dest)
        assert [s.sample for s in samples] == ["S1"]
        assert dest.is_file()

    def test_dispatches_to_a_directory(self, tmp_path: Path) -> None:
        reads = tmp_path / "reads"
        make_fastq(reads / "S1_R1.fastq.gz")
        make_fastq(reads / "S1_R2.fastq.gz")
        dest = tmp_path / "out" / "samplesheet.tsv"
        samples = build_samplesheet(reads, dest)
        assert [s.sample for s in samples] == ["S1"]
        assert samples[0].fastq_1[0].is_absolute()

    def test_missing_input(self, tmp_path: Path) -> None:
        with pytest.raises(SkillError) as excinfo:
            build_samplesheet(tmp_path / "nope", tmp_path / "out.tsv")
        assert excinfo.value.error_code == ErrorCode.MISSING_INPUT
