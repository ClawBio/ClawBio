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
    # read the FASTQ data — Nextflow does, and under the default Docker profile the
    # container often runs as root and can read files the launching user cannot, so an
    # os.access(R_OK) pre-check by the launcher would false-block valid runs. Output
    # *writability* IS checked (OUTPUT_DIR_NOT_WRITABLE) because the wrapper itself writes
    # there. Existence of inputs is validated (MISSING_FASTQ); readability is deferred to
    # Nextflow's staging, which reads in the true execution context.
    INVALID_SAMPLESHEET = "INVALID_SAMPLESHEET"
    # Remote input/reference URIs (s3://, gs://, https://, ftp://, …) are rejected by
    # default (local-first); they are only allowed with --allow-remote-inputs, which
    # also emits a runtime warning that data is being fetched over the network.
    REMOTE_INPUT_NOT_ALLOWED = "REMOTE_INPUT_NOT_ALLOWED"
    # --demo delegates to nf-core's upstream `-profile test`, whose inputs/references
    # are remote GitHub URLs by design. Under NXF_OFFLINE the demo cannot fetch them,
    # so preflight fails fast with this code instead of a cryptic Nextflow abort.
    DEMO_REQUIRES_NETWORK = "DEMO_REQUIRES_NETWORK"
    OUTPUT_DIR_NOT_EMPTY = "OUTPUT_DIR_NOT_EMPTY"
    OUTPUT_DIR_NOT_WRITABLE = "OUTPUT_DIR_NOT_WRITABLE"
    # Writing pipeline outputs inside the ClawBio source tree would pollute the
    # repository, so it is refused with a dedicated, self-describing code (parity
    # with nfcore-sarek / nfcore-rnaseq) rather than the generic _NOT_WRITABLE.
    OUTPUT_DIR_INSIDE_REPO = "OUTPUT_DIR_INSIDE_REPO"
    MISSING_REFERENCE = "MISSING_REFERENCE"
    CONFLICTING_REFERENCES = "CONFLICTING_REFERENCES"
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
    JAVA_VERSION_TOO_OLD = "JAVA_VERSION_TOO_OLD"
    MISSING_NEXTFLOW = "MISSING_NEXTFLOW"
    NEXTFLOW_VERSION_TOO_OLD = "NEXTFLOW_VERSION_TOO_OLD"
    PIPELINE_SOURCE_INVALID = "PIPELINE_SOURCE_INVALID"
    INVALID_RESUME_STATE = "INVALID_RESUME_STATE"
    EXECUTION_FAILED = "EXECUTION_FAILED"
    EXPECTED_OUTPUTS_NOT_FOUND = "EXPECTED_OUTPUTS_NOT_FOUND"
    UNEXPECTED_ERROR = "UNEXPECTED_ERROR"
