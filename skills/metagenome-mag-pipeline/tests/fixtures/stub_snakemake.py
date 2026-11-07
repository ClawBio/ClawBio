#!/usr/bin/env python3
"""A stand-in for the `snakemake` executable, for tests only.

Behaviour:
  * `--version`           prints 9.11.2
  * always                records its own argv to $STUB_SNAKEMAKE_ARGV as JSON
  * $STUB_SNAKEMAKE_EXIT  exits with that status without writing anything
  * `--dry-run`           prints a dry-run notice and writes nothing
  * otherwise             writes the pipeline's `final/`, `reports/` and
                          `provenance/` artifacts under the `results_dir` named in
                          `--configfile`, keyed on that config's enabled stages

The tables it writes use the upstream column names (see the pinned commit's
`workflow/scripts/`), so a parser that guesses a column name fails here.
"""

import json
import os
import sys
from pathlib import Path

import yaml

SAMPLE_COLUMNS = ["sample", "layout", "retained_reads", "retained_bases", "retained_pairs"]


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def table(columns: list[str], rows: list[list[object]]) -> str:
    lines = ["\t".join(columns)]
    lines += ["\t".join("" if value is None else str(value) for value in row) for row in rows]
    return "\n".join(lines) + "\n"


def read_samples(samplesheet: Path) -> list[tuple[str, bool]]:
    """(sample, paired) from a tab-separated samplesheet."""
    rows = []
    lines = samplesheet.read_text().splitlines()
    if not lines:
        return rows
    header = lines[0].split("\t")
    for line in lines[1:]:
        if not line.strip():
            continue
        cells = dict(zip(header, line.split("\t")))
        rows.append((cells.get("sample", ""), bool(cells.get("fastq_2"))))
    return rows


def main() -> int:
    argv = sys.argv[1:]

    record = os.environ.get("STUB_SNAKEMAKE_ARGV")
    if record:
        Path(record).write_text(json.dumps(argv))

    if "--version" in argv:
        print("9.11.2")
        return 0

    forced = os.environ.get("STUB_SNAKEMAKE_EXIT")
    if forced:
        print("stub snakemake: forced failure", file=sys.stderr)
        return int(forced)

    if "--dry-run" in argv:
        print("Job stats:\ndry run: 4 jobs, 0 total")
        print("Nothing to execute; this is a dry run.")
        return 0

    config_path = Path(argv[argv.index("--configfile") + 1]) if "--configfile" in argv else None
    if config_path is None or not config_path.is_file():
        print("stub snakemake: no readable --configfile", file=sys.stderr)
        return 1
    config = yaml.safe_load(config_path.read_text()) or {}

    results = Path(config["results_dir"])
    samples = read_samples(Path(config["samplesheet"])) or [("TEST_001", True)]

    # Always present.
    write(
        results / "final" / "qc" / "read_counts.tsv",
        table(
            SAMPLE_COLUMNS,
            [
                [name, "paired" if paired else "single", 4000, 1200000, 2000 if paired else None]
                for name, paired in samples
            ],
        ),
    )
    write(
        results / "final" / "diagnostic" / "bottleneck.tsv",
        table(
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
                [
                    name,
                    "paired" if paired else "single",
                    4000,
                    3100 if paired else None,
                    "77.50" if paired else None,
                    None,
                    None,
                ]
                for name, paired in samples
            ],
        ),
    )
    write(
        results / "final" / "benchmarks" / "runtime.tsv",
        table(
            ["rule", "h:m:s", "s", "max_rss", "max_vmem", "cpu_time", "benchmark_file", "label_source"],
            [["spades", "0:04:11", "251.0", "9000000", "12000000", "251.0", "b.tsv", "run"]],
        ),
    )
    write(
        results / "provenance" / "software_versions.tsv",
        table(
            ["tool", "version", "evidence"],
            [
                ["snakemake", "9.11.2", "importlib.metadata"],
                ["bbtools", "39.01", "workflow/envs/bbtools.yaml"],
            ],
        ),
    )
    write(
        results / "provenance" / "pipeline_revision.tsv",
        table(
            ["field", "value"],
            [
                ["commit", "7351a702d29857801963cae7aeff6600d59ebe77"],
                ["branch", "main"],
                ["dirty", "false"],
            ],
        ),
    )
    write(
        results / "provenance" / "reference_db.md5",
        "# resource\tmd5\nhost_grch38.fa\t6f1a1f0c5d1c2b3a4d5e6f708192a3b4\n",
    )
    write(
        results / "provenance" / "run_manifest.json",
        json.dumps(
            {
                "effective_config": {"results_dir": str(results)},
                "source_sha256": {"workflow/Snakefile": "0" * 64},
                "samplesheet_sha256": "1" * 64,
            }
        )
        + "\n",
    )

    if config.get("qc", {}).get("sourmash", {}).get("enabled"):
        names = [name for name, _ in samples]
        write(
            results / "final" / "qc" / "sourmash_similarities.csv",
            ",".join(["sample"] + names)
            + "\n"
            + "\n".join(
                ",".join([name] + ["1.0" if other == name else "0.42" for other in names])
                for name in names
            )
            + "\n",
        )

    if config.get("profile", {}).get("profiler", "none") != "none":
        names = [name for name, _ in samples]
        write(
            results / "final" / "profile" / "profile.tsv",
            "mpa_vJan26_CHOCOPhlAnSGB_202605\n"
            + "\t".join(["clade_name"] + names)
            + "\n"
            + "k__Bacteria\t"
            + "\t".join(["98.5000"] * len(names))
            + "\n"
            + "k__Bacteria|p__Bacillota\t"
            + "\t".join(["61.2500", "12.0000"][: len(names)])
            + "\n",
        )

    assembly = bool(config.get("assembly", {}).get("enabled", True))
    if assembly:
        write(
            results / "final" / "assembly" / "quast.tsv",
            table(
                ["Assembly"] + [name for name, _ in samples],
                [
                    ["# contigs", "331"],
                    ["Largest contig", "41200"],
                    ["Total length", "1284000"],
                    ["GC (%)", "51.20"],
                ],
            ),
        )
        write(results / "reports" / "assembly" / "quast.html", "<html>quast</html>\n")

    if config.get("binning", {}).get("enabled", True):
        for name, _ in samples:
            write(
                results / "intermediate" / "binning" / "dastool" / name / "contig2bin.tsv",
                table(
                    ["contig", "bin", "bin_type"],
                    [["contig_1", "bin.1", "dastool"], ["contig_2", "bin.2", "dastool"]],
                ),
            )

    if config.get("mag", {}).get("enabled", True):
        write(
            results / "final" / "binning" / "checkm2" / "quality_summary.tsv",
            table(
                ["genome", "completeness", "contamination", "passes_qa"],
                [["TEST_001.bin.1.fasta", "94.53", "2.10", "True"]],
            ),
        )
        write(
            results / "final" / "binning" / "drep" / "GenomeInfo.csv",
            "genome,completeness,contamination,secondary_cluster\n"
            "TEST_001.bin.1.fasta,94.53,2.10,NA\n",
        )
        write(
            results / "final" / "binning" / "gtdb" / "gtdb_summary.tsv",
            table(
                [
                    "user_genome",
                    "classification",
                    "closest_genome_reference",
                    "closest_genome_reference_score",
                    "closest_genome_taxonomy",
                ],
                [["TEST_001.bin.1", "d__Bacteria", "GCA_000195955.1", "98.50", "d__Bacteria"]],
            ),
        )
        write(
            results / "final" / "binning" / "bins.tsv",
            table(
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
                [["TEST_001.bin.1", "94.53", "2.10", "True", "d__Bacteria", "", "", "", "TEST_001"]],
            ),
        )

    if config.get("taxonomy", {}).get("enabled", True) and assembly:
        write(
            results / "final" / "taxonomy" / "contig_taxonomy.tsv",
            table(
                ["contig", "taxonomy", "assigned_genes", "total_genes", "disparity"],
                [
                    ["contig_1", "k__Bacteria;p__Bacillota", "18", "19", "0.05"],
                    ["contig_2", "unclassified", "3", "14", "0.79"],
                ],
            ),
        )

    if config.get("annotation", {}).get("enabled", True) and assembly:
        write(
            results / "final" / "annotation" / "function_abundance.tsv",
            table(
                ["bin", "source", "function", "genes"],
                [["unassigned", "kegg", "ko:K00001", "31"]],
            ),
        )
        if config.get("annotation", {}).get("minpath", {}).get("data_dir"):
            for name, _ in samples:
                write(
                    results / "intermediate" / "annotation" / "minpath" / f"{name}.minpath",
                    "pathway\tgenes\nmap00010\t3\n",
                )

    if config.get("heterogeneity", {}).get("enabled", True) and config.get("binning", {}).get(
        "enabled", True
    ):
        write(
            results / "final" / "heterogeneity" / "strain_heterogeneity.tsv",
            table(
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
                [["TEST_001", "TEST_001.bin.1", "412", "88400", "132", "640", "0", "0.1493"]],
            ),
        )

    if config.get("report", {}).get("multiqc", True):
        write(
            results / "reports" / "multiqc" / "multiqc_report.html", "<html>multiqc</html>\n"
        )

    print("stub snakemake: wrote results to " + str(results))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
