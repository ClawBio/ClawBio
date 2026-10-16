from __future__ import annotations

import re
from pathlib import Path



# resourceLimits of the pinned release's own -profile test (the `test` profile in nextflow.config). --demo keeps
# them: the launcher may lower them on a smaller machine but never raises them.
DEMO_RESOURCE_LIMITS = {"cpus": 4, "memory_gb": 15}

SKILL_NAME = "nfcore-rnaseq-wrapper"
SKILL_ALIAS = "rnaseq-pipeline"
SKILL_VERSION = "0.1.0"
SKILL_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SKILL_DIR.parent.parent
REPO_PARENT = PROJECT_ROOT.parent
DEFAULT_LOCAL_PIPELINE_DIR = REPO_PARENT / "rnaseq"
DEFAULT_REMOTE_PIPELINE = "nf-core/rnaseq"
DEFAULT_PIPELINE_VERSION = "3.26.0"
DEFAULT_PROFILE = "docker"
# Default wall-clock ceiling for a local Nextflow run. Overridable per-run via
# --timeout-hours so large cohorts are not killed mid-flight (audit F-02).
DEFAULT_TIMEOUT_HOURS = 12
DEFAULT_TIMEOUT_SECONDS = DEFAULT_TIMEOUT_HOURS * 60 * 60
SUPPORTED_PROFILES = {
    # Execution backends
    "docker",
    "conda",
    "mamba",
    "singularity",
    "apptainer",
    "podman",
    "shifter",
    "charliecloud",
    "wave",
    # Scientific modifier profiles
    "prokaryotic",
    "rapid_quant",
    # Architecture modifiers
    "arm64",
    "emulate_amd64",        # Run amd64 containers on ARM hardware via QEMU emulation
    # Self-contained nf-core test profiles (ship with params.input + refs)
    "test",
    "test_prokaryotic",
    "test_full",
    "test_full_aws",
    "test_full_gcp",
    "test_full_azure",
    "test_gpu",
    # Debug / development
    "debug",
}
JAVA_MIN_VERSION = 17
NEXTFLOW_MIN_VERSION = (25, 4, 3)

# Environment Nextflow needs for the pinned release. nf-core/rnaseq 3.26.0 parses
# under both the v1 and the strict v2 syntax parser, so it needs nothing.
NEXTFLOW_ENV: dict[str, str] = {}
PIPELINE_REQUIRED_FILES = ("main.nf", "nextflow.config", "assets/schema_input.json")

SUPPORTED_ALIGNERS = {"star_salmon", "star_rsem", "hisat2", "bowtie2_salmon"}
SUPPORTED_PSEUDO_ALIGNERS = {"salmon", "kallisto"}
SUPPORTED_TRIMMERS = {"trimgalore", "fastp"}
# Single source of truth for the default trimmer (nf-core/rnaseq 3.26.0 default).
# Consumed by params_builder and the wrapper's argument parser so the default
# lives in exactly one place.
DEFAULT_TRIMMER = "trimgalore"
SUPPORTED_RIBO_TOOLS = {"sortmerna", "ribodetector", "bowtie2"}
# nf-core/rnaseq 3.26.0 nextflow_schema.json contaminant_screening enum. Single
# source shared by the argument parser and mirrored by both reproducibility JSONs.
SUPPORTED_CONTAMINANT_SCREENING = {"kraken2", "kraken2_bracken", "sylph"}
# RSeQC modules accepted by nf-core/rnaseq 3.26.0 (nextflow_schema.json rseqc_modules).
SUPPORTED_RSEQC_MODULES = {
    "bam_stat",
    "inner_distance",
    "infer_experiment",
    "junction_annotation",
    "junction_saturation",
    "read_distribution",
    "read_duplication",
    "tin",
}
SUPPORTED_UMI_TOOLS = {"umitools", "umicollapse"}
SUPPORTED_STRANDEDNESS = {"forward", "reverse", "unstranded", "auto"}
# nf-core/rnaseq 3.26.0 nextflow_schema.json enum values
SUPPORTED_UMI_EXTRACT_METHODS = {"string", "regex"}
SUPPORTED_UMI_GROUPING_METHODS = {"unique", "percentile", "cluster", "adjacency", "directional"}
SUPPORTED_SALMON_LIB_TYPES = {
    "A", "IS", "ISF", "ISR", "IU", "MS", "MSF", "MSR",
    "MU", "OS", "OSF", "OSR", "OU", "SF", "SR", "U",
}
SUPPORTED_PUBLISH_DIR_MODES = {"symlink", "rellink", "link", "copy", "copyNoFollow", "move"}
# Bracken taxonomic precision levels, ordered domain→species as upstream lists them
# (nf-core/rnaseq 3.26.0 nextflow_schema.json bracken_precision enum). Tuple, not set,
# so --help preserves the meaningful taxonomic order.
SUPPORTED_BRACKEN_PRECISION = ("D", "P", "C", "O", "F", "G", "S")

# Shared regexes — defined once here so every module enforces the identical rule.
# FASTQ_BASENAME_RE: basename ends in .fq/.fastq(.gz) and contains no whitespace or
# path separator (nf-core schema fastq pattern). WHITESPACE_RE: any whitespace char,
# used to reject paths the nf-core ^\S+ schema would later refuse.
FASTQ_BASENAME_RE = re.compile(r"^[^\s/]+\.f(?:ast)?q(?:\.gz)?$")
WHITESPACE_RE = re.compile(r"\s")
SUPPORTED_SAMPLE_COLUMNS = {
    "sample",
    "fastq_1",
    "fastq_2",
    "strandedness",
    "genome_bam",
    "transcriptome_bam",
    "seq_platform",
    "seq_center",
    "percent_mapped",
}
REQUIRED_SAMPLE_COLUMNS = ("sample", "fastq_1", "strandedness")

# Top-level genome keys verified from nf-core/rnaseq 3.26.0
# conf/igenomes.config.
SUPPORTED_IGENOMES_NAMES = frozenset(
    {
        "AGPv3",
        "BDGP6",
        "CHIMP2.1.4",
        "CHM13",
        "CanFam3.1",
        "EB1",
        "EB2",
        "EF2",
        "EquCab2",
        "Galgal4",
        "Gm01",
        "GRCh37",
        "GRCh38",
        "GRCm38",
        "GRCz10",
        "IRGSP-1.0",
        "Mmul_1",
        "R64-1-1",
        "Rnor_5.0",
        "Rnor_6.0",
        "Sbi1",
        "Sscrofa10.2",
        "TAIR10",
        "UMD3.1",
        "WBcel235",
        "bosTau8",
        "canFam3",
        "ce10",
        "danRer10",
        "dm6",
        "equCab2",
        "galGal4",
        "hg19",
        "hg38",
        "mm10",
        "panTro4",
        "rn6",
        "sacCer3",
        "susScr3",
    }
)
