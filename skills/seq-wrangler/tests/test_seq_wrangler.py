import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import seq_wrangler as sw
from pathlib import Path
import json
import types
import pytest
from inline_snapshot import snapshot

from clawbio.common import demo_runner
import seq_wrangler as sw


def test_load_samplesheet_mixed_se_pe(tmp_path: Path):
    p = tmp_path / "samples.csv"
    p.write_text(
        "sample,fastq_1,fastq_2\n"
        "CTRL,a_R1.fastq.gz,a_R2.fastq.gz\n"
        "TREAT,b.fastq.gz,\n",
        encoding="utf-8"
    )
    rows = sw.load_samplesheet(p)
    assert len(rows) == 2
    assert rows[0]["sample"] == "CTRL"
    assert rows[0]["fastq_2"] == "a_R2.fastq.gz"
    assert rows[1]["sample"] == "TREAT"
    assert rows[1]["fastq_2"] == ""


def test_load_samplesheet_missing_required_columns(tmp_path: Path):
    p = tmp_path / "bad.csv"
    p.write_text("sample,foo\nS1,x.fastq.gz\n", encoding="utf-8")
    with pytest.raises(ValueError):
        sw.load_samplesheet(p)


def test_detect_threads_positive():
    assert sw.detect_threads() >= 1


def test_check_tool_false_for_nonexistent():
    assert sw.check_tool("definitely_not_a_real_binary_123") is False


def test_sha256_file(tmp_path: Path):
    p = tmp_path / "x.txt"
    p.write_text("abc", encoding="utf-8")
    h = sw.sha256_file(p)
    assert len(h) == 64
    assert h == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"


def test_run_command_non_critical_failure():
    result = sw.run_command(
        [sys.executable, "-c", "import sys; sys.exit(3)"],
        "fail-non-critical",
        critical=False,
    )
    assert result.returncode == 3


def test_validate_args_demo_ok():
    args = types.SimpleNamespace(
        demo=True, index=None, samplesheet=None, r1=None, r2=None
    )
    sw.validate_args(args)  # no exception


def test_validate_args_requires_index_when_not_demo():
    args = types.SimpleNamespace(
        demo=False, index=None, samplesheet=None, r1="a.fastq.gz", r2=None
    )
    with pytest.raises(ValueError):
        sw.validate_args(args)


def test_validate_args_samplesheet_xor_single():
    args = types.SimpleNamespace(
        demo=False, index="/ref/hg38", samplesheet="samples.csv", r1="a.fastq.gz", r2=None
    )
    with pytest.raises(ValueError):
        sw.validate_args(args)


def test_generate_report_and_demo_outputs(tmp_path: Path):
    args = types.SimpleNamespace(
        output=str(tmp_path),
        aligner="bwa",
        mapq=20,
        remove_duplicates=False,
        samplesheet=None,
        r1=None,
        r2=None,
        index=None,
        threads=2,
        trim=False,
        run_fastqc=False,
        run_multiqc=False,
        demo=True,
        genome_build="GRCh38",
    )

    sw.run_demo(args)

    report = tmp_path / "report.md"
    summary = tmp_path / "summary.json"
    repro_cmd = tmp_path / "reproducibility" / "commands.sh"
    repro_env = tmp_path / "reproducibility" / "environment.yml"

    assert report.exists()
    assert summary.exists()
    assert repro_cmd.exists()
    assert repro_env.exists()

    data = json.loads(summary.read_text(encoding="utf-8"))
    assert isinstance(data, list)
    assert len(data) == 2
    assert data[0]["sample"] in {"CTRL_REP1", "TREAT_REP1"}

# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/seq-wrangler/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("seq-wrangler", tmp_path) == snapshot({'report.md': """\
# Seq Wrangler Report

**Date**: <TIMESTAMP> UTC  \n\
**Tool**: ClawBio Seq Wrangler v0.1.0  \n\
**Mode**: Demo (synthetic data)

---

## Run parameters

| Parameter | Value |
|-----------|-------|
| Aligner | bwa |
| MAPQ threshold | ≥20 |
| Duplicate handling | Mark only |
| Samples processed | 2 |

---

## Results per sample

### CTRL_REP1 (Paired-end)

**Alignment statistics**

| Metric | Value |
|--------|-------|
| Total reads | 10,000,000 |
| Mapped reads | 9,750,000 (97.5%) |
| Duplicates | 850,000 |
| Properly paired | 9,500,000 |

**Coverage (top chromosomes)**

| Chromosome | Mean depth | Breadth |
|------------|------------|---------|
| chr1 | 28.4x | 97.0% |
| chr2 | 27.9x | 96.0% |
| chrX | 14.2x | 91.0% |
| chrM | 1842.0x | 100.0% |

**Insert size**: mean 312.4 bp | median 305.0 bp | std 48.2 bp

**BAM**: `<OUTPUT_DIR>/bam/CTRL_REP1_sorted.bam`

### TREAT_REP1 (Single-end)

**Alignment statistics**

| Metric | Value |
|--------|-------|
| Total reads | 8,500,000 |
| Mapped reads | 8,200,000 (96.5%) |
| Duplicates | 600,000 |

**Coverage (top chromosomes)**

| Chromosome | Mean depth | Breadth |
|------------|------------|---------|
| chr1 | 24.1x | 95.0% |
| chr2 | 23.8x | 94.0% |
| chrX | 23.5x | 94.0% |
| chrM | 1560.0x | 100.0% |

**BAM**: `<OUTPUT_DIR>/bam/TREAT_REP1_sorted.bam`


---

## Methods

- **QC**: FastQC per-sample quality assessment
- **Trimming**: fastp (if --trim flag used)
- **Alignment**: bwa with MAPQ ≥ 20 filter
- **BAM processing**: samtools sort -n → fixmate → coordinate-sort → markdup (mark only) → index
- **Coverage**: samtools coverage (per-chromosome mean depth and breadth)
- **Insert size**: samtools stats IS histogram (paired-end only)
- **MultiQC**: aggregated QC report across samples (if --run-multiqc)

## Input checksums

- `demo_ctrl_R1.fastq.gz`: demo data (no file on disk)- `demo_treat_R1.fastq.gz`: demo data (no file on disk)

## Disclaimer

> ClawBio Seq Wrangler is a research and educational tool. Results must be validated before use in clinical or production settings.
""", 'summary.json': """\
[
  {
    "sample": "CTRL_REP1",
    "is_pe": true,
    "bam": "<OUTPUT_DIR>/bam/CTRL_REP1_sorted.bam",
    "flagstat": {
      "total": 10000000,
      "mapped": 9750000,
      "mapped_pct": 97.5,
      "duplicates": 850000,
      "dup_pct": 8.7,
      "paired": true,
      "properly_paired": 9500000
    },
    "coverage": [
      {
        "chr": "chr1",
        "mean_depth": 28.4,
        "breadth": 0.97
      },
      {
        "chr": "chr2",
        "mean_depth": 27.9,
        "breadth": 0.96
      },
      {
        "chr": "chrX",
        "mean_depth": 14.2,
        "breadth": 0.91
      },
      {
        "chr": "chrM",
        "mean_depth": 1842.0,
        "breadth": 1.0
      }
    ],
    "insert_size": {
      "mean": 312.4,
      "std": 48.2,
      "median": 305.0
    }
  },
  {
    "sample": "TREAT_REP1",
    "is_pe": false,
    "bam": "<OUTPUT_DIR>/bam/TREAT_REP1_sorted.bam",
    "flagstat": {
      "total": 8500000,
      "mapped": 8200000,
      "mapped_pct": 96.5,
      "duplicates": 600000,
      "dup_pct": 7.3,
      "paired": false,
      "properly_paired": 0
    },
    "coverage": [
      {
        "chr": "chr1",
        "mean_depth": 24.1,
        "breadth": 0.95
      },
      {
        "chr": "chr2",
        "mean_depth": 23.8,
        "breadth": 0.94
      },
      {
        "chr": "chrX",
        "mean_depth": 23.5,
        "breadth": 0.94
      },
      {
        "chr": "chrM",
        "mean_depth": 1560.0,
        "breadth": 1.0
      }
    ],
    "insert_size": null
  }
]\
"""})
