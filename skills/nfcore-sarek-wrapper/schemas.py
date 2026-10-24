from __future__ import annotations

from pathlib import Path



# resourceLimits of the pinned release's own -profile test (conf/test.config). --demo keeps
# them: the launcher may lower them on a smaller machine but never raises them.
DEMO_RESOURCE_LIMITS = {"cpus": 4, "memory_gb": 15}

SKILL_NAME = "nfcore-sarek-wrapper"
SKILL_ALIAS = "sarek-pipeline"
SKILL_VERSION = "0.1.0"
SKILL_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SKILL_DIR.parent.parent
REPO_PARENT = PROJECT_ROOT.parent
DEFAULT_LOCAL_PIPELINE_DIR = REPO_PARENT / "sarek"
DEFAULT_REMOTE_PIPELINE = "nf-core/sarek"
DEFAULT_PIPELINE_VERSION = "3.8.1"
DEFAULT_PROFILE = "docker"
DEFAULT_TIMEOUT_SECONDS = 60 * 60 * 24   # 24h — sarek WGS+variant calling can be long
DEFAULT_STEP = "mapping"
DEFAULT_ALIGNER = "bwa-mem"
DEFAULT_GENOME = "GATK.GRCh38"
# nf-core/sarek 3.8.1 nextflow.config: igenomes_base = 's3://ngi-igenomes/igenomes/'.
# When this default base is in effect, --genome must be one of the keys defined in
# conf/igenomes.config (SUPPORTED_IGENOMES_NAMES); a custom base may define others.
DEFAULT_IGENOMES_BASE = "s3://ngi-igenomes/igenomes/"

JAVA_MIN_VERSION = 17
# From nf-core/sarek 3.8.1 nextflow.config: manifest.nextflowVersion = '!>=25.10.2'.
NEXTFLOW_MIN_VERSION = (25, 10, 2)

# Environment Nextflow needs for the pinned release, applied to the run, shown in
# the recorded command and exported by commands.sh. nf-core/sarek 3.8.1 predates
# Nextflow's strict syntax, which Nextflow 26.04 applies by default and which
# rejects its configs ("Unexpected input: ':'" in conf/modules); the v1 parser
# runs the release unmodified.
NEXTFLOW_ENV: dict[str, str] = {"NXF_SYNTAX_PARSER": "v1"}
PIPELINE_REQUIRED_FILES = ("main.nf", "nextflow.config", "assets/schema_input.json")


# 25 profiles from nf-core/sarek 3.8.1 nextflow.config (verified against
# https://github.com/nf-core/sarek/blob/3.8.1/nextflow.config — profiles block).
SUPPORTED_PROFILES = {
    # Container runtimes (9)
    "docker", "podman", "singularity", "apptainer",
    "shifter", "charliecloud", "wave", "conda", "mamba",
    # Architecture modifiers (2)
    "arm64", "emulate_amd64",
    # Runtime modifiers: gpu, spark. `mutect` is a TEST-DATA profile
    # (conf/test_mutect2.config hardcodes Mutect2 `--normal-sample normal`); the
    # token stays allowed for the test profile but is gated to --demo in
    # _validate_wrapper_flags so it never corrupts a real paired run.
    "gpu", "spark", "mutect",
    # Test profiles (10)
    "test", "test_aws", "test_azure",
    "test_full", "test_full_aws", "test_full_azure",
    "test_full_germline", "test_full_germline_aws", "test_full_germline_azure",
    "test_full_germline_ncbench_agilent",
    # Development (1)
    "debug",
}

SUPPORTED_ALIGNERS = {"bwa-mem", "bwa-mem2", "dragmap", "sentieon-bwamem", "parabricks"}

SUPPORTED_STEPS = {
    "mapping",
    "markduplicates",
    "prepare_recalibration",
    "recalibrate",
    "variant_calling",
    "annotate",
}

# Values recognised by --tools (Sarek 3.8.1 nextflow_schema.json `tools` pattern,
# 28 tokens). Used by preflight to reject unknown tools, so this MUST be complete.
SUPPORTED_TOOLS = {
    # Germline callers
    "haplotypecaller", "sentieon_haplotyper", "sentieon_dnascope",
    "deepvariant", "freebayes", "mpileup",
    # Somatic / tumor-only callers
    "mutect2", "muse", "lofreq", "strelka", "varlociraptor",
    "sentieon_tnscope",
    # Structural variants
    "manta", "tiddit",
    # CNV / ploidy
    "ascat", "cnvkit", "controlfreec", "indexcov",
    # MSI
    "msisensor2", "msisensorpro",
    # Preprocessing / dedup
    "sentieon_dedup",
    # Annotators
    "snpeff", "vep", "snpsift", "bcfann", "merge",
    # Misc QC / filtering
    "ngscheckmate", "bbsplit",
}

# Values recognised by --skip_tools (Sarek 3.8.1 `skip_tools` pattern, 15 tokens).
SUPPORTED_SKIP_TOOLS = {
    "baserecalibrator", "baserecalibrator_report",
    "bcftools", "dnascope_filter", "documentation",
    "fastqc", "haplotypecaller_filter", "haplotyper_filter",
    "markduplicates", "markduplicates_report",
    "mosdepth", "multiqc", "samtools", "vcftools", "versions",
}

# fgbio GroupReadsByUmi strategies (schema enum: Identity, Edit, Adjacency, Paired).
SUPPORTED_GROUP_BY_UMI = {"Identity", "Edit", "Adjacency", "Paired"}
SUPPORTED_UMI_LOCATIONS = {"read1", "read2", "per_read", "index1", "index2", "per_index"}

# GATK PCR indel model values. `gatk_pcr_indel_model` has no schema enum (free
# string) and is left unvalidated; `sentieon_dnascope_pcr_indel_model` is gated
# against this set in preflight (single source of truth — no inline duplication).
SUPPORTED_GATK_PCR_INDEL_MODEL = {"NONE", "HOSTILE", "AGGRESSIVE", "CONSERVATIVE"}
# Sentieon Haplotyper/DNAscope emit modes. The 3.8.1 schema declares these as
# free strings (no enum), but only the Sentieon-documented base modes and their
# `gvcf` pairings are valid; preflight rejects anything else. Single source of
# truth for both emit-mode params (no inline regex/list duplication).
_SENTIEON_EMIT_BASE_MODES = ("variant", "confident", "all", "gvcf")
SUPPORTED_SENTIEON_EMIT_MODE = frozenset(
    _SENTIEON_EMIT_BASE_MODES
    + tuple(f"gvcf,{m}" for m in ("all", "confident", "variant"))
    + tuple(f"{m},gvcf" for m in ("all", "confident", "variant"))
)
SUPPORTED_USE_GATK_SPARK = {"baserecalibrator", "markduplicates"}

SUPPORTED_VEP_OUT_FORMAT = {"json", "tab", "vcf"}
SUPPORTED_ASCAT_GENOME = {"hg19", "hg38"}

SUPPORTED_PUBLISH_DIR_MODES = {"symlink", "rellink", "link", "copy", "copyNoFollow", "move"}

# The 42 iGenomes keys accepted by --genome, mirrored from nf-core/sarek 3.8.1 conf/igenomes.config
# https://raw.githubusercontent.com/nf-core/sarek/3.8.1/conf/igenomes.config
# Bare "GRCh37"/"GRCh38"/"GRCm39"/"mm39" are NOT iGenomes
# keys — use the prefixed variants (GATK.GRCh38, Ensembl.GRCh37, NCBI.GRCh38).
SUPPORTED_IGENOMES_NAMES = frozenset({
    "AGPv3",
    "BDGP6",
    "CHIMP2.1.4",
    "CHM13",
    "CanFam3.1",
    "EB1",
    "EB2",
    "EF2",
    "Ensembl.GRCh37",
    "EquCab2",
    "GATK.GRCh37",
    "GATK.GRCh38",
    "GRCm38",
    "GRCz10",
    "Galgal4",
    "Gm01",
    "IRGSP-1.0",
    "Mmul_1",
    "NCBI.GRCh38",
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
    "testdata.nf-core.sarek",
    # Some upstream entries are duplicated under different cases (e.g. galGal4
    # vs Galgal4, equCab2 vs EquCab2, canFam3 vs CanFam3.1) — both are kept
    # because both resolve in conf/igenomes.config.
})

# Samplesheet columns recognised across all --step values
# (step-aware required-column logic lives in samplesheet_builder.py).
SUPPORTED_SAMPLE_COLUMNS = {
    "patient", "sample", "lane", "sex", "status",
    "fastq_1", "fastq_2", "spring_1", "spring_2",
    "bam", "bai", "cram", "crai", "table",
    "vcf", "variantcaller",
    "contamination",   # required by varlociraptor for tumor-only/somatic
}
