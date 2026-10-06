"""Samplesheet loading, discovery and writing.

Every rule enforced here is upstream's own, from `read_samplesheet()` in
`workflow/rules/helpers.smk` at the pinned commit. The wrapper repeats them
because upstream only checks them after snakemake has started, and a cohort of
runs should not be lost to a typo in column three.

The one deliberate difference: FASTQ paths are resolved against the samplesheet's
own directory and made absolute, because upstream's `resolve()` joins relative
paths to the pipeline root, not to the sheet.
"""

from __future__ import annotations

import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

from magpipe_errors import ErrorCode, SkillError
from magpipe_schemas import (
    FASTQ_SUFFIXES,
    SAMPLE_ID_PATTERN,
    SAMPLESHEET_COLUMNS,
    SRA_ACCESSION_PATTERN,
)

STAGE = "samplesheet"

_SAMPLE_ID = re.compile(SAMPLE_ID_PATTERN)
_SRA_ACCESSION = re.compile(SRA_ACCESSION_PATTERN)
# The mate token sits at the end of the stem: `SAMPLE_01_R1`, `SAMPLE_01_1`,
# `SRR1761673.2`. A stem ending in a bare digit with no separator before it
# (`TEST_001`) is a sample name, not a mate, so the separator is required.
_MATE_TAIL = re.compile(r"[._-](?:(?P<r>[Rr])(?P<mr>[12])|(?P<m>[12]))$")
# `lane1`, `lane_01`, `Lane2`: a lane label inside a flat sample name.
_LANE_TOKEN = re.compile(r"lane[0-9_]*\Z", re.IGNORECASE)


@dataclass(frozen=True)
class Sample:
    sample: str
    fastq_1: tuple[Path, ...]
    fastq_2: tuple[Path, ...]
    group: str
    sra_run: str | None
    paired: bool


def _fail(code: str, message: str, fix: str, **details: object) -> SkillError:
    return SkillError(
        stage=STAGE, error_code=code, message=message, fix=fix, details=dict(details)
    )


def _is_fastq(path: Path) -> bool:
    return path.is_file() and path.name.lower().endswith(FASTQ_SUFFIXES)


def _stem_of(path: Path) -> str:
    """`SAMPLE_01_lane1_R1.fastq.gz` -> `SAMPLE_01_lane1_R1`."""
    name = path.name
    for suffix in FASTQ_SUFFIXES:
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def _split_lanes(value: str) -> list[str]:
    """Semicolon-separated lanes, in matching mate order."""
    return [part.strip() for part in (value or "").split(";") if part.strip()]


def _absolute(path_text: str, base: Path, sample: str) -> Path:
    candidate = Path(path_text).expanduser()
    resolved = candidate if candidate.is_absolute() else (base / candidate)
    if not resolved.exists():
        raise _fail(
            ErrorCode.MISSING_FASTQ,
            f"{sample}: FASTQ not found: {path_text}",
            "Check the path in the samplesheet. Relative paths resolve against the "
            "samplesheet's own directory.",
            sample=sample,
            path=str(resolved),
        )
    return resolved.resolve()


def load_samplesheet(path: Path) -> list[Sample]:
    """Read and validate a tab-separated samplesheet.

    Raises `SkillError` with `INVALID_SAMPLESHEET` naming the sample and the rule
    it broke, or `MISSING_FASTQ` when a named file does not exist.
    """
    path = Path(path)
    if not path.is_file():
        raise _fail(
            ErrorCode.MISSING_INPUT,
            f"Samplesheet not found: {path}",
            "Pass --input <samplesheet.tsv> or a directory of FASTQ files.",
            path=str(path),
        )

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = {"sample", "fastq_1"} - set(reader.fieldnames or [])
        if missing:
            raise _fail(
                ErrorCode.INVALID_SAMPLESHEET,
                f"{path} is missing required column(s): {sorted(missing)}",
                "The sheet must have at least `sample` and `fastq_1`; see "
                "examples/samplesheet.example.tsv.",
                missing=sorted(missing),
            )
        rows = list(reader)

    if not rows:
        raise _fail(
            ErrorCode.INVALID_SAMPLESHEET,
            f"{path} contains no samples",
            "Add at least one row, or run with --demo to generate synthetic reads.",
            path=str(path),
        )

    base = path.parent.resolve()
    samples: list[Sample] = []
    seen: set[str] = set()
    for row in rows:
        name = (row.get("sample") or "").strip()
        if not name:
            raise _fail(
                ErrorCode.INVALID_SAMPLESHEET,
                f"{path} has a row with an empty sample identifier",
                "Every row needs a unique `sample`.",
                path=str(path),
            )
        if not _SAMPLE_ID.fullmatch(name):
            raise _fail(
                ErrorCode.INVALID_SAMPLESHEET,
                f"invalid sample identifier {name!r}",
                "Sample identifiers may hold letters, digits, underscore, dot and "
                "hyphen only.",
                sample=name,
            )
        if name in seen:
            raise _fail(
                ErrorCode.INVALID_SAMPLESHEET,
                f"duplicate sample {name!r} in {path}",
                "Sample identifiers must be unique; merge lanes with `;` instead.",
                sample=name,
            )
        seen.add(name)

        layout = (row.get("layout") or "").strip().lower()
        if layout not in ("", "paired", "single"):
            raise _fail(
                ErrorCode.INVALID_SAMPLESHEET,
                f"{name}: layout must be 'paired' or 'single', got {layout!r}",
                "Leave layout empty to infer it from the mates present.",
                sample=name,
                layout=layout,
            )

        accession = (row.get("sra_run") or "").strip()
        first_text = (row.get("fastq_1") or "").strip()
        second_text = (row.get("fastq_2") or "").strip()

        if accession:
            if first_text or second_text:
                raise _fail(
                    ErrorCode.INVALID_SAMPLESHEET,
                    f"{name}: sra_run and fastq_1/fastq_2 are mutually exclusive",
                    "A row either names a public run or holds paths, never both.",
                    sample=name,
                )
            if not _SRA_ACCESSION.fullmatch(accession):
                raise _fail(
                    ErrorCode.INVALID_SAMPLESHEET,
                    f"{name}: {accession!r} is not a run accession",
                    "Use an SRA, ENA or DRA run accession such as SRR1761673.",
                    sample=name,
                    sra_run=accession,
                )
            if not layout:
                raise _fail(
                    ErrorCode.INVALID_SAMPLESHEET,
                    f"{name}: sra_run rows require an explicit layout",
                    "sra-tools only discovers the layout while downloading, which is "
                    "too late for the workflow, so state `paired` or `single`.",
                    sample=name,
                )
            first: tuple[Path, ...] = ()
            second = ()
            paired = layout == "paired"
        else:
            first_lanes = _split_lanes(first_text)
            second_lanes = _split_lanes(second_text)
            if not first_lanes or (second_lanes and len(first_lanes) != len(second_lanes)):
                raise _fail(
                    ErrorCode.INVALID_SAMPLESHEET,
                    f"{name}: missing R1 or unequal numbers of R1/R2 lanes",
                    "One sample may hold several lanes, but both mates must list the "
                    "same number of them in matching order.",
                    sample=name,
                    fastq_1_lanes=len(first_lanes),
                    fastq_2_lanes=len(second_lanes),
                )
            if len(set(first_lanes + second_lanes)) != len(first_lanes + second_lanes):
                raise _fail(
                    ErrorCode.INVALID_SAMPLESHEET,
                    f"{name}: repeated FASTQ path",
                    "The same file may not fill two lanes of one sample.",
                    sample=name,
                )
            if layout == "paired" and not second_lanes:
                raise _fail(
                    ErrorCode.INVALID_SAMPLESHEET,
                    f"{name}: layout=paired but fastq_2 is empty",
                    "Either supply the R2 lanes or declare layout=single.",
                    sample=name,
                )
            if layout == "single" and second_lanes:
                raise _fail(
                    ErrorCode.INVALID_SAMPLESHEET,
                    f"{name}: layout=single but fastq_2 is set",
                    "A single-end sample may not carry R2 lanes.",
                    sample=name,
                )
            first = tuple(_absolute(lane, base, name) for lane in first_lanes)
            second = tuple(_absolute(lane, base, name) for lane in second_lanes)
            paired = layout == "paired" if layout else bool(second)

        group = (row.get("group") or name).strip()
        if not _SAMPLE_ID.fullmatch(group):
            raise _fail(
                ErrorCode.INVALID_SAMPLESHEET,
                f"invalid assembly group {group!r}",
                "Assembly groups may hold letters, digits, underscore, dot and "
                "hyphen only.",
                sample=name,
                group=group,
            )

        samples.append(
            Sample(
                sample=name,
                fastq_1=first,
                fastq_2=second,
                group=group,
                sra_run=accession or None,
                paired=paired,
            )
        )
    return samples


def _name_and_lane(path: Path) -> tuple[str, str, int | None]:
    """(sample, lane label, mate) for one FASTQ filename.

    The mate token must be at the end of the stem and preceded by a separator, so
    `TEST_001_1.fastq.gz` is mate 1 of `TEST_001` while `TEST_001.fastq.gz` is an
    unpaired read file for a sample called `TEST_001`.
    """
    stem = _stem_of(path)
    match = _MATE_TAIL.search(stem)
    if match is None:
        return stem, "", None
    base = stem[: match.start()]
    mate = int(match.group("mr") or match.group("m"))
    lane = ""
    if "." in base or "_" in base or "-" in base:
        head, _, tail = base.replace("-", "_").replace(".", "_").rpartition("_")
        if head and _LANE_TOKEN.fullmatch(tail):
            base, lane = head, tail
    return base, lane, mate


def _samples_in(directory: Path, forced_name: str | None) -> list[Sample]:
    """Build samples from the FASTQ files directly inside `directory`."""
    entries: dict[str, dict[str, list[tuple[str, Path]]]] = {}
    for path in sorted(directory.iterdir()):
        if not _is_fastq(path):
            continue
        name, lane, mate = _name_and_lane(path)
        if forced_name is not None:
            # One folder per sample: the folder names the sample, and the lane
            # label inside the filename still orders that sample's lanes.
            name = forced_name
        # A filename carrying no mate token is a single-end read file, i.e. R1.
        bucket = entries.setdefault(name, {"1": [], "2": []})
        bucket["1" if mate is None else str(mate)].append((lane, path))

    samples: list[Sample] = []
    for name in sorted(entries):
        if not _SAMPLE_ID.fullmatch(name):
            raise _fail(
                ErrorCode.INVALID_SAMPLESHEET,
                f"invalid sample identifier {name!r} discovered in {directory}",
                "Rename the file or folder, or write a samplesheet by hand.",
                sample=name,
                directory=str(directory),
            )
        bucket = entries.setdefault(name, {"1": [], "2": []})
        first_lanes = sorted(bucket["1"], key=lambda item: (item[0], item[1].name))
        second_lanes = sorted(bucket["2"], key=lambda item: (item[0], item[1].name))
        if not first_lanes:
            raise _fail(
                ErrorCode.INVALID_SAMPLESHEET,
                f"{name}: R2 reads found without any R1 in {directory}",
                "Every paired-end library needs both mates.",
                sample=name,
                directory=str(directory),
            )
        if second_lanes and len(first_lanes) != len(second_lanes):
            raise _fail(
                ErrorCode.INVALID_SAMPLESHEET,
                f"{name}: {len(first_lanes)} R1 files but {len(second_lanes)} R2 files "
                f"in {directory}",
                "Lanes of one sample must come in matching pairs, named in the same "
                "lane order.",
                sample=name,
                directory=str(directory),
            )
        samples.append(
            Sample(
                sample=name,
                fastq_1=tuple(path.resolve() for _, path in first_lanes),
                fastq_2=tuple(path.resolve() for _, path in second_lanes),
                group=name,
                sra_run=None,
                paired=bool(second_lanes),
            )
        )
    return samples


def discover_samples(reads_dir: Path) -> list[Sample]:
    """Build samples from a directory of FASTQ files.

    Two layouts are accepted, and mixing them is refused rather than guessed at:
    flat files whose names carry the mate (`<id>_R1`/`<id>_R2`, `<id>_1`/`<id>_2`),
    or one subdirectory per sample with the same naming inside it. Several lanes
    of one sample are sorted by name and joined into one sample.
    """
    reads_dir = Path(reads_dir)
    if not reads_dir.is_dir():
        raise _fail(
            ErrorCode.MISSING_INPUT,
            f"Reads directory not found: {reads_dir}",
            "Pass a directory of FASTQ files, or a samplesheet TSV.",
            path=str(reads_dir),
        )

    subdirs = [entry for entry in sorted(reads_dir.iterdir()) if entry.is_dir()]
    with_reads = [entry for entry in subdirs if any(_is_fastq(f) for f in entry.iterdir())]
    loose = [entry for entry in sorted(reads_dir.iterdir()) if _is_fastq(entry)]

    if with_reads and loose:
        raise _fail(
            ErrorCode.INVALID_SAMPLESHEET,
            f"{reads_dir} mixes per-sample subdirectories with loose FASTQ files",
            "Point --input at a directory in one layout only, or write a samplesheet "
            "listing every sample explicitly.",
            directory=str(reads_dir),
        )

    samples: list[Sample] = []
    if with_reads:
        for entry in with_reads:
            samples.extend(_samples_in(entry, forced_name=entry.name))
    elif loose:
        samples = _samples_in(reads_dir, forced_name=None)
    else:
        raise _fail(
            ErrorCode.MISSING_INPUT,
            f"No FASTQ files found in {reads_dir}",
            "Expected FASTQ files (*.fastq.gz, *.fq.gz, *.fastq, *.fq), either flat "
            "or one directory per sample.",
            directory=str(reads_dir),
        )

    if not samples:
        raise _fail(
            ErrorCode.MISSING_INPUT,
            f"No samples could be derived from {reads_dir}",
            "Check the file naming, or write a samplesheet TSV by hand.",
            directory=str(reads_dir),
        )
    return sorted(samples, key=lambda sample: sample.sample)


def write_samplesheet(samples: list[Sample], dest: Path) -> Path:
    """Write the six upstream columns, tab-separated, with absolute paths."""
    if not samples:
        raise _fail(
            ErrorCode.INVALID_SAMPLESHEET,
            "Refusing to write an empty samplesheet",
            "Provide at least one sample, or run with --demo.",
        )
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(SAMPLESHEET_COLUMNS)
        for sample in samples:
            writer.writerow(
                [
                    sample.sample,
                    ";".join(str(path) for path in sample.fastq_1),
                    ";".join(str(path) for path in sample.fastq_2),
                    sample.group,
                    sample.sra_run or "",
                    "paired" if sample.paired else "single",
                ]
            )
    return dest


def build_samplesheet(input_path: Path, dest: Path) -> list[Sample]:
    """Load `--input` (a file or a directory) and write the sheet to `dest`."""
    path = Path(input_path)
    if path.is_dir():
        samples = discover_samples(path)
    elif path.is_file():
        samples = load_samplesheet(path)
    else:
        raise _fail(
            ErrorCode.MISSING_INPUT,
            f"Input not found: {path}",
            "Pass a samplesheet TSV, a directory of FASTQ files, or --demo to "
            "generate synthetic reads.",
            path=str(path),
        )
    write_samplesheet(samples, dest)
    return samples
