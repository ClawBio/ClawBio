from __future__ import annotations

from pathlib import Path



# resourceLimits of the pinned release's own -profile test (conf/test.config). --demo keeps
# them: the launcher may lower them on a smaller machine but never raises them.
DEMO_RESOURCE_LIMITS = {"cpus": 4, "memory_gb": 15}

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
SUPPORTED_PRESETS = {
    "standard",
    "star",
    "kallisto",
    "cellranger",
    "cellrangerarc",
    "cellrangermulti",
}
SUPPORTED_PROFILES = {
    "docker",
    "conda",
    "mamba",
    "singularity",
    "apptainer",
    "podman",
    "shifter",
    "charliecloud",
    # Wave containers / GPU execution — both defined in nf-core/scrnaseq's
    # nextflow.config profiles block (kept in sync with pinned_versions.json).
    "wave",
    "gpu",
    # nf-core/scrnaseq 4.1.0 Nextflow-native modifier/test profiles.
    "debug",
    "arm64",
    "emulate_amd64",
    "test",
    "test_full",
    "test_cellrangermulti",
    "test_multiome",
}
JAVA_MIN_VERSION = 17
# Comparison tuple (numeric) and the canonical human-readable string. Nextflow
# releases use a zero-padded YY.MM.patch format ("25.04.0"), which is the ONLY
# correct way to show the version to users — the tuple renders to "25.4.0",
# which is not a real Nextflow version. Keep both in sync (guarded by a test).
NEXTFLOW_MIN_VERSION = (25, 4, 0)
NEXTFLOW_MIN_VERSION_DISPLAY = "25.04.0"

# Environment Nextflow needs for the pinned release, applied to the run, shown in
# the recorded command and exported by commands.sh. nf-core/scrnaseq 4.1.0 predates
# Nextflow's strict syntax, which Nextflow 26.04 applies by default: it rejects the
# release's include of conf/test_multiome.config, a file the tag does not ship. The
# v1 parser runs the release unmodified.
NEXTFLOW_ENV: dict[str, str] = {"NXF_SYNTAX_PARSER": "v1"}
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

PRESET_ALIGNERS = {
    "standard": "simpleaf",
    "star": "star",
    "kallisto": "kallisto",
    "cellranger": "cellranger",
    "cellrangerarc": "cellrangerarc",
    "cellrangermulti": "cellrangermulti",
}

# ── Reference field categories (single source of truth) ───────────────────────
# Both params_builder and preflight import these so the rules for "what counts as
# a genome reference" cannot drift between modules. Presence of a GENOME reference
# means iGenomes is unused (→ igenomes_ignore) and conflicts with --genome.
# Auxiliary files (barcode whitelists, CMO/probe/feature references, primers,
# multiplexed-sample samplesheets) are local inputs that are NOT genome
# references: they must never suppress iGenomes nor conflict with --genome.
GENOME_REFERENCE_FIELDS = (
    "fasta",
    "gtf",
    "transcript_fasta",
    "txp2gene",
    "simpleaf_index",
    "kallisto_index",
    "star_index",
    "cellranger_index",
)
AUXILIARY_PATH_FIELDS = (
    "barcode_whitelist",
    "kb_t1c",
    "kb_t2c",
    "motifs",
    "cellrangerarc_config",
    "cellranger_vdj_index",
    "gex_frna_probe_set",
    "gex_target_panel",
    "gex_cmo_set",
    "fb_reference",
    "vdj_inner_enrichment_primers",
    "gex_barcode_sample_assignment",
    "cellranger_multi_barcodes",
)
# Symbolic identifiers (not local paths; never existence-checked).
SYMBOLIC_REFERENCE_FIELDS = ("genome", "cellrangerarc_reference")
# Every local path written to params.yaml / existence-checked in preflight.
ALL_REFERENCE_PATH_FIELDS = GENOME_REFERENCE_FIELDS + AUXILIARY_PATH_FIELDS


def profile_components(profile: str) -> list[str]:
    return [part.strip() for part in str(profile).split(",") if part.strip()]


def profile_includes(profile: str, component: str) -> bool:
    return component in profile_components(profile)


PRESET_REQUIREMENTS = {
    "standard": {
        "requires_any": [
            ("genome",),
            ("simpleaf_index",),
            ("fasta", "gtf"),
            ("transcript_fasta", "txp2gene"),
        ],
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
    # ARC/Multi reference building is documented from a prebuilt cellranger_index
    # or from fasta+gtf (ARC: auto-build + optional motif/config; Multi: GEX/VDJ
    # built on the fly from fasta+gtf). --genome (iGenomes) is a valid path because
    # it resolves to fasta+gtf for the pipeline. The genuine gap — a GEX-only
    # cellranger_index that cannot build a VDJ reference — is enforced separately by
    # _check_cellrangermulti_vdj_reference_policy (audit finding #3).
    "cellrangerarc": {
        "requires_any": [("genome",), ("cellranger_index",), ("fasta", "gtf")],
    },
    "cellrangermulti": {
        "requires_any": [("genome",), ("cellranger_index",), ("fasta", "gtf")],
    },
}
