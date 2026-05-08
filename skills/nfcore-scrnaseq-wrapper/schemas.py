from __future__ import annotations

from pathlib import Path


SKILL_NAME = "nfcore-scrnaseq-wrapper"
SKILL_ALIAS = "scrnaseq-pipeline"
SKILL_VERSION = "0.1.0"
SKILL_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SKILL_DIR.parent.parent
REPO_PARENT = PROJECT_ROOT.parent
DEFAULT_LOCAL_PIPELINE_DIR = REPO_PARENT / "scrnaseq"
DEFAULT_REMOTE_PIPELINE = "nf-core/scrnaseq"
DEFAULT_PIPELINE_VERSION = "4.1.0"
DEFAULT_PROFILE = "docker"
DEFAULT_PRESET = "standard"
DEFAULT_TIMEOUT_SECONDS = 60 * 60 * 12
SUPPORTED_PRESETS = {"standard", "star", "kallisto", "cellranger", "cellrangerarc", "cellrangermulti"}
SUPPORTED_PROFILES = {
    "docker", "conda", "mamba",
    "singularity", "apptainer",
    "podman", "shifter", "charliecloud",
}
JAVA_MIN_VERSION = 17
NEXTFLOW_MIN_VERSION = (25, 4, 0)
PIPELINE_REQUIRED_FILES = ("main.nf", "nextflow.config", "assets/schema_input.json")
SUPPORTED_SAMPLE_COLUMNS = {
    "sample",
    "fastq_1",
    "fastq_2",
    "expected_cells",
    "seq_center",
    "sample_type",
    "fastq_barcode",
    "feature_type",
}
REQUIRED_SAMPLE_COLUMNS = ("sample", "fastq_1", "fastq_2")
SUPPORTED_SAMPLE_TYPE_VALUES = {"atac", "gex"}
SUPPORTED_FEATURE_TYPE_VALUES = {"gex", "vdj", "ab", "crispr", "cmo"}
PREFERRED_H5AD_ORDER = (
    "combined_cellbender_filter_matrix.h5ad",
    "combined_filtered_matrix.h5ad",
    "combined_raw_matrix.h5ad",
)

COMMON_GENOMES = (
    "GRCh38",
    "GRCh37",
    "hg38",
    "hg19",
    "GRCm38",
    "GRCm39",
    "mm10",
    "mm9",
    "BDGP6",
    "WBcel235",
    "R64-1-1",
    "CanFam3.1",
    "UMD3.1",
    "Mmul_8.0.1",
    "Sscrofa11.1",
)

PRESET_ALIGNERS = {
    "standard": "simpleaf",
    "star": "star",
    "kallisto": "kallisto",
    "cellranger": "cellranger",
    "cellrangerarc": "cellrangerarc",
    "cellrangermulti": "cellrangermulti",
}

PRESET_REQUIREMENTS = {
    "standard": {
        "requires_any": [("genome",), ("simpleaf_index",), ("fasta", "gtf"), ("transcript_fasta", "txp2gene")],
    },
    "star": {
        "requires_any": [("genome",), ("star_index",), ("fasta", "gtf")],
    },
    "kallisto": {
        "requires_any": [("genome",), ("kallisto_index",), ("fasta", "gtf")],
    },
    "cellranger": {
        "requires_any": [("genome",), ("cellranger_index",), ("fasta", "gtf")],
    },
    "cellrangerarc": {
        "requires_any": [("genome",), ("cellranger_index",), ("fasta", "gtf")],
    },
    "cellrangermulti": {
        "requires_any": [("genome",), ("cellranger_index",), ("fasta", "gtf")],
    },
}
