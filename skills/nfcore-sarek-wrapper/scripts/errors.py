from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class SkillError(Exception):
    stage: str
    error_code: str
    message: str
    fix: str
    details: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        super().__init__(self.message)

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": "error",
            "ok": False,
            "stage": self.stage,
            "error_code": self.error_code,
            "message": self.message,
            "fix": self.fix,
            "details": self.details,
        }


class ErrorCode:
    MISSING_INPUT = "MISSING_INPUT"
    MISSING_FASTQ = "MISSING_FASTQ"
    INVALID_FASTQ = "INVALID_FASTQ"
    # NOTE: input readability is intentionally NOT pre-checked. The wrapper does not
    # read the FASTQ/BAM/CRAM data — Nextflow does, and under the default Docker profile
    # the container often runs as root and can read files the launching user cannot, so an
    # os.access(R_OK) pre-check by the launcher would false-block valid runs. Output
    # *writability* IS checked (OUTPUT_DIR_NOT_WRITABLE) because the wrapper itself writes
    # there. Existence of inputs is validated (MISSING_FASTQ); readability is deferred to
    # Nextflow's staging, which reads in the true execution context.
    INVALID_SAMPLESHEET = "INVALID_SAMPLESHEET"
    # Remote input/reference URIs (s3://, gs://, https://, ftp://, …) are rejected by
    # default (local-first); they are only allowed with --allow-remote-inputs, which
    # also emits a runtime warning that data is being fetched over the network.
    REMOTE_INPUT_NOT_ALLOWED = "REMOTE_INPUT_NOT_ALLOWED"
    # --demo composes nf-core's upstream `-profile test`, whose inputs/references are
    # remote GitHub URLs by design. Under NXF_OFFLINE the demo cannot fetch them, so
    # preflight fails fast with this code instead of a cryptic Nextflow abort.
    DEMO_REQUIRES_NETWORK = "DEMO_REQUIRES_NETWORK"
    OUTPUT_DIR_NOT_EMPTY = "OUTPUT_DIR_NOT_EMPTY"
    OUTPUT_DIR_NOT_WRITABLE = "OUTPUT_DIR_NOT_WRITABLE"
    MISSING_REFERENCE = "MISSING_REFERENCE"
    CONFLICTING_REFERENCES = "CONFLICTING_REFERENCES"
    INVALID_GENOME = "INVALID_GENOME"
    INVALID_PRESET_CONFIGURATION = "INVALID_PRESET_CONFIGURATION"
    UNSUPPORTED_MODE = "UNSUPPORTED_MODE"
    INVALID_PROFILE = "INVALID_PROFILE"
    MISSING_DOCKER = "MISSING_DOCKER"
    DOCKER_NOT_RUNNING = "DOCKER_NOT_RUNNING"
    MISSING_CONDA = "MISSING_CONDA"
    MISSING_SINGULARITY = "MISSING_SINGULARITY"
    MISSING_PODMAN = "MISSING_PODMAN"
    PODMAN_NOT_RUNNING = "PODMAN_NOT_RUNNING"
    MISSING_HPC_RUNTIME = "MISSING_HPC_RUNTIME"
    MISSING_JAVA = "MISSING_JAVA"
    MISSING_NEXTFLOW = "MISSING_NEXTFLOW"
    NEXTFLOW_VERSION_TOO_OLD = "NEXTFLOW_VERSION_TOO_OLD"
    PIPELINE_SOURCE_INVALID = "PIPELINE_SOURCE_INVALID"
    INVALID_RESUME_STATE = "INVALID_RESUME_STATE"
    EXECUTION_FAILED = "EXECUTION_FAILED"
    EXPECTED_OUTPUTS_NOT_FOUND = "EXPECTED_OUTPUTS_NOT_FOUND"
    UNEXPECTED_ERROR = "UNEXPECTED_ERROR"
    INVALID_ALIGNER = "INVALID_ALIGNER"
    CONFLICTING_FASTA_ARGS = "CONFLICTING_FASTA_ARGS"
    REFERENCE_PATH_NOT_FOUND = "REFERENCE_PATH_NOT_FOUND"
    INVALID_BAM = "INVALID_BAM"
    INVALID_STEP = "INVALID_STEP"
    INVALID_TOOLS = "INVALID_TOOLS"
    INVALID_PAIRING = "INVALID_PAIRING"
    CONFLICTING_STEPS = "CONFLICTING_STEPS"
    MISSING_PON = "MISSING_PON"
    MISSING_GERMLINE_RESOURCE = "MISSING_GERMLINE_RESOURCE"
    MISSING_VEP_CACHE = "MISSING_VEP_CACHE"
    MISSING_SNPEFF_DB = "MISSING_SNPEFF_DB"
    INVALID_ASCAT_RESOURCES = "INVALID_ASCAT_RESOURCES"
    UMI_SPARK_INCOMPATIBLE = "UMI_SPARK_INCOMPATIBLE"
    VARLOCIRAPTOR_MISSING_CONTAMINATION = "VARLOCIRAPTOR_MISSING_CONTAMINATION"
    INVALID_SEX = "INVALID_SEX"
    INVALID_STATUS = "INVALID_STATUS"
    INVALID_INTERVALS = "INVALID_INTERVALS"
    INVALID_EMAIL = "INVALID_EMAIL"
    OUTPUT_DIR_INSIDE_REPO = "OUTPUT_DIR_INSIDE_REPO"
    NEXTFLOW_NOT_FOUND = "NEXTFLOW_NOT_FOUND"
    JAVA_VERSION_TOO_OLD = "JAVA_VERSION_TOO_OLD"
    BACKEND_UNAVAILABLE = "BACKEND_UNAVAILABLE"
    MISSING_ANNOTATION_CACHE = "MISSING_ANNOTATION_CACHE"
    RESUME_DRIFT = "RESUME_DRIFT"
    INVALID_FLAG_COMBINATION = "INVALID_FLAG_COMBINATION"
    INVALID_CONSENSUS_CALLING = "INVALID_CONSENSUS_CALLING"
    MISSING_INTERVALS = "MISSING_INTERVALS"
    MISSING_DBNSFP = "MISSING_DBNSFP"
    INVALID_BCFTOOLS_FILTER = "INVALID_BCFTOOLS_FILTER"
    MISSING_SNPSIFT_DATABASES = "MISSING_SNPSIFT_DATABASES"
    BUILD_ONLY_INDEX_CONFLICT = "BUILD_ONLY_INDEX_CONFLICT"
    GENOME_REQUIRED = "GENOME_REQUIRED"
    PARABRICKS_CONDA_CONFLICT = "PARABRICKS_CONDA_CONFLICT"
    MISSING_SENTIEON_LICENSE = "MISSING_SENTIEON_LICENSE"
    JOINT_GERMLINE_REQUIREMENT = "JOINT_GERMLINE_REQUIREMENT"
    JOINT_MUTECT2_REQUIREMENT = "JOINT_MUTECT2_REQUIREMENT"
    SPARK_OUTPUT_INCOMPATIBLE = "SPARK_OUTPUT_INCOMPATIBLE"
