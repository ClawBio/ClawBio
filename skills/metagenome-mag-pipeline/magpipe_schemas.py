"""Constants shared by every module of the metagenome-mag-pipeline wrapper.

Nothing here computes anything. These are the facts that must not drift: where
the upstream pipeline lives, the commit this skill is pinned to, the identifiers
upstream's `read_samplesheet()` accepts, and the tables it writes. Every value
here is read from the upstream checkout at `PINNED_REF`, not from memory.
"""

from __future__ import annotations

SKILL_NAME = "metagenome-mag-pipeline"
SKILL_VERSION = "0.1.0"
CLI_ALIAS = "mag-pipeline"

UPSTREAM_REPO = "https://github.com/mobashirrahman/metagenomics-workflow.git"
UPSTREAM_LICENSE = "MIT"

# Upstream release tag cut on the pinned commit (2026-10-06). The commit stays the
# operational pin — tags can move, commits cannot — and the tag is recorded for
# provenance and human-readable references.
PINNED_REF = "7351a702d29857801963cae7aeff6600d59ebe77"
PINNED_TAG = "v0.2.0"

# The upstream Snakefile declares min_version("9.0") and pins snakemake=9.11.2
# in its environment.yml.
MIN_SNAKEMAKE: tuple[int, int] = (9, 0)
SNAKEMAKE_PIN = "9.11.2"

# Files whose presence marks a directory as a usable checkout of the upstream
# workflow. Checked before anything is launched.
PIPELINE_REQUIRED_FILES: tuple[str, ...] = (
    "workflow/Snakefile",
    "workflow/config/config.yaml",
    "workflow/config/config_test.yaml",
    "test/make_test_dataset.py",
)

# Presets and the stage switches they apply. The mapping itself lives in
# magpipe_config.PRESET_SWITCHES; this tuple is the allowlist argparse validates
# --preset against.
PRESETS: tuple[str, ...] = ("qc", "profile", "assembly", "full")

# binning.binners entries upstream accepts. helpers.smk rejects anything else at
# parse time and needs at least two for a consensus to mean anything.
BINNERS: tuple[str, ...] = ("metabat", "maxbin", "vamb")

# taxonomy.sensitivity must name one of DIAMOND's own settings; upstream raises
# at config-parse time rather than letting the tool reject it hours into a cohort.
DIAMOND_SENSITIVITIES: tuple[str, ...] = (
    "fast",
    "mid-sensitive",
    "sensitive",
    "more-sensitive",
    "very-sensitive",
    "ultra-sensitive",
)

# The columns write_samplesheet() emits, in upstream's order.
SAMPLESHEET_COLUMNS: tuple[str, ...] = (
    "sample",
    "fastq_1",
    "fastq_2",
    "group",
    "sra_run",
    "layout",
)

# Upstream's own patterns, from workflow/rules/helpers.smk.
SAMPLE_ID_PATTERN = r"[A-Za-z0-9_.-]+"
SRA_ACCESSION_PATTERN = r"(?:SR[RXZP]|EN[RXZP]|DR[RXZP]|ER[RXZP])\d{6,}"

FASTQ_SUFFIXES: tuple[str, ...] = (".fastq.gz", ".fq.gz", ".fastq", ".fq")

# Logical name -> path relative to the pipeline's results_dir. Only the paths
# upstream writes under final/, reports/ and provenance/ belong here; the
# intermediate/ trees are working files, not report material.
FINAL_OUTPUTS: dict[str, str] = {
    "read_counts": "final/qc/read_counts.tsv",
    "sourmash_similarities": "final/qc/sourmash_similarities.csv",
    "profile": "final/profile/profile.tsv",
    "quast": "final/assembly/quast.tsv",
    "quast_html": "reports/assembly/quast.html",
    "bins": "final/binning/bins.tsv",
    "checkm2_quality": "final/binning/checkm2/quality_summary.tsv",
    "drep_genome_info": "final/binning/drep/GenomeInfo.csv",
    "gtdb_summary": "final/binning/gtdb/gtdb_summary.tsv",
    "bottleneck": "final/diagnostic/bottleneck.tsv",
    "contig_taxonomy": "final/taxonomy/contig_taxonomy.tsv",
    "function_abundance": "final/annotation/function_abundance.tsv",
    "strain_heterogeneity": "final/heterogeneity/strain_heterogeneity.tsv",
    "runtime": "final/benchmarks/runtime.tsv",
    "multiqc": "reports/multiqc/multiqc_report.html",
    "software_versions": "provenance/software_versions.tsv",
    "pipeline_revision": "provenance/pipeline_revision.tsv",
    "reference_db_md5": "provenance/reference_db.md5",
    "run_manifest": "provenance/run_manifest.json",
}

# The four provenance files plus the timing table the upstream `all` rule always
# requests, whatever stages are enabled.
ALWAYS_EXPECTED_OUTPUTS: tuple[str, ...] = (
    FINAL_OUTPUTS["runtime"],
    FINAL_OUTPUTS["software_versions"],
    FINAL_OUTPUTS["pipeline_revision"],
    FINAL_OUTPUTS["reference_db_md5"],
    FINAL_OUTPUTS["run_manifest"],
)

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device and "
    "does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)

DEFAULT_TIMEOUT_HOURS = 12.0
