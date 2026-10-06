"""Structured failures for the metagenome-mag-pipeline wrapper.

Every failure the wrapper raises on purpose is a `SkillError`: it names the
stage that failed, a stable error code, what went wrong, and the command or flag
that fixes it. `main()` prints `to_dict()` as JSON on stderr and exits 1, so an
agent reads the fix out of the error rather than guessing.
"""

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
    INVALID_SAMPLESHEET = "INVALID_SAMPLESHEET"
    MISSING_FASTQ = "MISSING_FASTQ"
    # A samplesheet row naming a public SRA/ENA/DRA run downloads data. That is
    # network egress under a local-first skill, so it needs --allow-remote-inputs.
    REMOTE_INPUT_NOT_ALLOWED = "REMOTE_INPUT_NOT_ALLOWED"
    INVALID_PRESET = "INVALID_PRESET"
    INVALID_CONFIG = "INVALID_CONFIG"
    # taxonomy, annotation and the MAG stage need external databases that upstream
    # does not bundle. A partially supplied set is an error, never a silent skip.
    MISSING_DATABASE = "MISSING_DATABASE"
    # metaSPAdes needs both mates; upstream refuses a single-end library outright.
    ASSEMBLY_REQUIRES_PAIRED = "ASSEMBLY_REQUIRES_PAIRED"
    MISSING_SNAKEMAKE = "MISSING_SNAKEMAKE"
    SNAKEMAKE_VERSION_TOO_OLD = "SNAKEMAKE_VERSION_TOO_OLD"
    MISSING_CONDA = "MISSING_CONDA"
    PIPELINE_SOURCE_INVALID = "PIPELINE_SOURCE_INVALID"
    PIPELINE_FETCH_FAILED = "PIPELINE_FETCH_FAILED"
    DEMO_REQUIRES_NETWORK = "DEMO_REQUIRES_NETWORK"
    OUTPUT_DIR_NOT_EMPTY = "OUTPUT_DIR_NOT_EMPTY"
    OUTPUT_DIR_NOT_WRITABLE = "OUTPUT_DIR_NOT_WRITABLE"
    EXECUTION_FAILED = "EXECUTION_FAILED"
    EXECUTION_TIMEOUT = "EXECUTION_TIMEOUT"
    EXPECTED_OUTPUTS_NOT_FOUND = "EXPECTED_OUTPUTS_NOT_FOUND"
    UNEXPECTED_ERROR = "UNEXPECTED_ERROR"
