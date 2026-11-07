"""Shared fixtures for the metagenome-mag-pipeline tests.

Nothing here touches the network or runs real snakemake. The `fake_pipeline`
fixture is a minimal tree satisfying `PIPELINE_REQUIRED_FILES`, the `stub_snakemake`
and `stub_conda` fixtures put fake executables on PATH, and `fake_results` is a
hand-written `results/` tree for the `full` preset.

The upstream `workflow/config/config.yaml` and `config_test.yaml` are copied
verbatim into `tests/fixtures/` at the pinned commit; they are the schema this
wrapper renders against, so paraphrasing them would make the tests meaningless.
"""

from __future__ import annotations

import importlib
import json
import os
import shutil
import stat
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"

if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))


def load(name: str):
    """Import `magpipe_<name>` from the skill directory."""
    if str(SKILL_DIR) not in sys.path:
        sys.path.insert(0, str(SKILL_DIR))
    return importlib.import_module(f"magpipe_{name}")


def _make_executable(path: Path) -> None:
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def write_rows(
    path: Path,
    columns: list[str],
    rows: list[list[str]],
    delimiter: str = "\t",
) -> Path:
    """Write a delimited table with a header, creating parent directories."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [delimiter.join(columns)]
    lines += [
        delimiter.join("" if value is None else str(value) for value in row) for row in rows
    ]
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def fake_pipeline(tmp_path: Path) -> Path:
    """A pipeline checkout that satisfies PIPELINE_REQUIRED_FILES.

    The two config files are verbatim copies of the upstream files at the pinned
    commit. `test/make_test_dataset.py` is a stand-in with upstream's CLI, writing
    two synthetic samples of gzipped FASTQ and a samplesheet with the same columns
    upstream writes.
    """
    schemas = load("schemas")

    root = tmp_path / "metagenomics-workflow"
    (root / "workflow" / "config").mkdir(parents=True)
    (root / "workflow" / "rules").mkdir(parents=True)
    (root / "workflow" / "resources").mkdir(parents=True)
    (root / "test").mkdir(parents=True)

    (root / "workflow" / "Snakefile").write_text(
        '# Minimal stand-in for the upstream Snakefile.\nmin_version("9.0")\n'
    )
    shutil.copyfile(
        FIXTURES_DIR / "upstream_config.yaml", root / "workflow" / "config" / "config.yaml"
    )
    shutil.copyfile(
        FIXTURES_DIR / "upstream_config_test.yaml",
        root / "workflow" / "config" / "config_test.yaml",
    )
    (root / "workflow" / "rules" / "helpers.smk").write_text("# stand-in\n")
    (root / "workflow" / "resources" / "illumina_adapters.fa").write_text(">adapter\nACGT\n")
    shutil.copyfile(
        FIXTURES_DIR / "stub_make_test_dataset.py", root / "test" / "make_test_dataset.py"
    )

    missing = [name for name in schemas.PIPELINE_REQUIRED_FILES if not (root / name).exists()]
    assert not missing, f"fake_pipeline is missing {missing}"
    return root


@pytest.fixture
def stub_snakemake(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Put a fake `snakemake` first on PATH and return its binary directory.

    The stub answers `--version`, records its argv to $STUB_SNAKEMAKE_ARGV, exits
    with $STUB_SNAKEMAKE_EXIT when set, and otherwise writes the pipeline's final
    tables under the results_dir named in `--configfile` — unless `--dry-run`.
    """
    bindir = tmp_path / "stub_bin"
    bindir.mkdir()
    stub = bindir / "snakemake"
    shutil.copyfile(FIXTURES_DIR / "stub_snakemake.py", stub)
    _make_executable(stub)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("STUB_SNAKEMAKE_ARGV", str(tmp_path / "snakemake_argv.json"))
    return bindir


@pytest.fixture
def stub_conda(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Put a fake `conda` first on PATH and return its binary directory."""
    bindir = tmp_path / "stub_conda_bin"
    bindir.mkdir()
    stub = bindir / "conda"
    shutil.copyfile(FIXTURES_DIR / "stub_conda.py", stub)
    _make_executable(stub)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")
    return bindir


@pytest.fixture
def no_pipeline_tools(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Take snakemake, conda and mamba off PATH for the duration of a test."""
    empty = tmp_path / "empty_path"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))


@pytest.fixture
def recorded_argv(stub_snakemake: Path) -> callable:
    """Return a reader for the argv the stub snakemake recorded."""

    def read() -> list[str]:
        path = Path(os.environ["STUB_SNAKEMAKE_ARGV"])
        return json.loads(path.read_text()) if path.exists() else []

    return read


@pytest.fixture
def fake_results(tmp_path: Path) -> Path:
    """A populated `results/` tree for the `full` preset.

    Hand-written tables with upstream's column names, including a blank cell in
    `bottleneck.tsv` and in `strain_heterogeneity.tsv` so the "blank is not zero"
    rule is exercised by a real fixture.
    """
    results = tmp_path / "results"

    write_rows(
        results / "final" / "qc" / "read_counts.tsv",
        ["sample", "layout", "retained_reads", "retained_bases", "retained_pairs"],
        [
            ["S1", "paired", "4000", "1200000", "2000"],
            ["S2", "single", "2500", "750000", ""],
        ],
    )
    write_rows(
        results / "final" / "qc" / "sourmash_similarities.csv",
        ["sample", "S1", "S2"],
        [["S1", "1.0", "0.31"], ["S2", "0.31", "1.0"]],
        delimiter=",",
    )
    # MetaPhlAn format: version line, header line, then percentage relative abundance.
    (results / "final" / "profile").mkdir(parents=True, exist_ok=True)
    (results / "final" / "profile" / "profile.tsv").write_text(
        "mpa_vJan26_CHOCOPhlAnSGB_202605\n"
        "clade_name\tS1\tS2\n"
        "k__Bacteria\t98.5000\t97.9000\n"
        "k__Bacteria|p__Firmicutes\t61.2500\t12.0000\n"
        "k__Bacteria|p__Proteobacteria\t20.0000\t70.5000\n"
    )
    write_rows(
        results / "final" / "assembly" / "quast.tsv",
        ["Assembly", "S1"],
        [
            ["# contigs (>= 0 bp)", "412"],
            ["# contigs", "331"],
            ["Largest contig", "41200"],
            ["Total length", "1284000"],
            ["GC (%)", "51.20"],
        ],
    )
    (results / "reports" / "assembly").mkdir(parents=True, exist_ok=True)
    (results / "reports" / "assembly" / "quast.html").write_text("<html>quast</html>\n")
    write_rows(
        results / "final" / "binning" / "checkm2" / "quality_summary.tsv",
        ["genome", "completeness", "contamination", "passes_qa"],
        [
            ["S1.bin.1.fasta", "94.53", "2.10", "True"],
            ["S1.bin.2.fasta", "71.20", "8.40", "True"],
            ["S2.bin.1.fasta", "42.00", "12.30", "False"],
        ],
    )
    write_rows(
        results / "final" / "binning" / "drep" / "GenomeInfo.csv",
        ["genome", "completeness", "contamination", "secondary_cluster"],
        [["S1.bin.1.fasta", "94.53", "2.10", "NA"]],
        delimiter=",",
    )
    write_rows(
        results / "final" / "binning" / "gtdb" / "gtdb_summary.tsv",
        [
            "user_genome",
            "classification",
            "closest_genome_reference",
            "closest_genome_reference_score",
            "closest_genome_taxonomy",
        ],
        [
            [
                "S1.bin.1",
                "d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus_x",
                "GCA_000195955.1",
                "98.50",
                "d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus",
            ],
            ["S1.bin.2", "", "", "", ""],
        ],
    )
    write_rows(
        results / "final" / "binning" / "bins.tsv",
        [
            "genome",
            "completeness",
            "contamination",
            "passes_qa",
            "classification",
            "closest_genome_reference",
            "closest_genome_reference_score",
            "closest_genome_taxonomy",
            "sample",
        ],
        [
            [
                "S1.bin.1",
                "94.53",
                "2.10",
                "True",
                "d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus",
                "GCA_000195955.1",
                "98.50",
                "d__Bacteria;p__Bacillota;c__Bacilli",
                "S1",
            ],
            ["S1.bin.2", "71.20", "8.40", "True", "", "", "", "", "S1"],
            ["S2.bin.1", "42.00", "12.30", "False", "", "", "", "", "S2"],
        ],
    )
    write_rows(
        results / "intermediate" / "binning" / "dastool" / "S1" / "contig2bin.tsv",
        ["contig", "bin", "bin_type"],
        [["contig_1", "bin.1", "dastool"], ["contig_2", "bin.2", "dastool"]],
    )
    write_rows(
        results / "intermediate" / "binning" / "dastool" / "S2" / "contig2bin.tsv",
        ["contig", "bin", "bin_type"],
        [["contig_1", "bin.1", "dastool"]],
    )
    # assembly_capture_pct is blank for the single-end sample: not measured, never zero.
    write_rows(
        results / "final" / "diagnostic" / "bottleneck.tsv",
        [
            "sample",
            "layout",
            "reads",
            "assembly_unique",
            "assembly_capture_pct",
            "catalogue_unique",
            "catalogue_capture_pct",
        ],
        [
            ["S1", "paired", "4000", "3100", "77.50", "2600", "65.00"],
            ["S2", "single", "2500", "", "", "1500", "60.00"],
        ],
    )
    write_rows(
        results / "final" / "taxonomy" / "contig_taxonomy.tsv",
        ["contig", "taxonomy", "assigned_genes", "total_genes", "disparity"],
        [
            ["contig_1", "k__Bacteria;p__Bacillota", "18", "19", "0.05"],
            ["contig_2", "unclassified", "3", "14", "0.79"],
        ],
    )
    write_rows(
        results / "final" / "annotation" / "function_abundance.tsv",
        ["bin", "source", "function", "genes"],
        [
            ["unassigned", "kegg", "ko:K00001", "31"],
            ["S1.bin.1", "cog", "COG0001", "12"],
            ["S1.bin.1", "pfam", "PF00001", "4"],
        ],
    )
    # heterogeneity_pct is blank for S2.bin.1: no coding position passed the depth filter.
    write_rows(
        results / "final" / "heterogeneity" / "strain_heterogeneity.tsv",
        [
            "sample",
            "bin",
            "genes",
            "covered_positions",
            "non_synonymous_codons",
            "synonymous_codons",
            "ambiguous_codons",
            "heterogeneity_pct",
        ],
        [
            ["S1", "S1.bin.1", "412", "88400", "132", "640", "0", "0.1493"],
            ["S1", "S1.bin.2", "118", "23100", "88", "210", "0", "0.3810"],
            ["S2", "S2.bin.1", "97", "0", "0", "0", "0", ""],
        ],
    )
    write_rows(
        results / "final" / "benchmarks" / "runtime.tsv",
        ["rule", "h:m:s", "s", "max_rss", "max_vmem", "cpu_time", "benchmark_file", "label_source"],
        [
            ["bbtools", "0:02:13", "133.0", "1200000", "2000000", "133.0", "b.tsv", "run"],
            ["spades", "1:40:02", "6002.0", "9000000", "12000000", "2100.0", "s.tsv", "run"],
        ],
    )
    write_rows(
        results / "provenance" / "software_versions.tsv",
        ["tool", "version", "evidence"],
        [
            ["snakemake", "9.11.2", "importlib.metadata"],
            ["bbtools", "39.01", "workflow/envs/bbtools.yaml"],
            ["metaphlan", "4.2.5", "workflow/envs/metaphlan4.yaml"],
        ],
    )
    write_rows(
        results / "provenance" / "pipeline_revision.tsv",
        ["field", "value"],
        [
            ["commit", "7351a702d29857801963cae7aeff6600d59ebe77"],
            ["commit_date", "2026-06-14T09:12:44+00:00"],
            ["branch", "main"],
            ["remote", "https://github.com/mobashirrahman/metagenomics-workflow.git"],
            ["dirty", "false"],
        ],
    )
    (results / "provenance" / "reference_db.md5").write_text(
        "resource\tmd5\n"
        "host_grch38.fa\t6f1a1f0c5d1c2b3a4d5e6f708192a3b4\n"
        "checkm2_uniref100.dmnd\t0a1b2c3d4e5f60718293a4b5c6d7e8f9\n"
    )
    (results / "provenance" / "run_manifest.json").write_text(
        '{"effective_config": {"results_dir": "results"},'
        ' "source_sha256": {"workflow/Snakefile": "abc"},'
        ' "samplesheet_sha256": "def"}\n'
    )
    (results / "reports" / "multiqc").mkdir(parents=True, exist_ok=True)
    (results / "reports" / "multiqc" / "multiqc_report.html").write_text(
        "<html>multiqc</html>\n"
    )
    return results
