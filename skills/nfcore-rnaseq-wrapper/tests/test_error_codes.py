from __future__ import annotations

from pathlib import Path
import inspect
import json
import sys

_SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_SKILL_DIR))


def _purge_foreign_modules(*names: str) -> None:
    for name in names:
        module = sys.modules.get(name)
        module_file = Path(getattr(module, "__file__", "") or "")
        if module is not None and _SKILL_DIR not in module_file.parents and module_file != _SKILL_DIR / f"{name}.py":
            sys.modules.pop(name, None)


def _purge_local_modules(*names: str) -> None:
    for name in names:
        module = sys.modules.get(name)
        module_file = Path(getattr(module, "__file__", "") or "")
        if module is not None and (module_file == _SKILL_DIR / f"{name}.py" or _SKILL_DIR in module_file.parents):
            sys.modules.pop(name, None)


def _remove_skill_dir_from_sys_path() -> None:
    while str(_SKILL_DIR) in sys.path:
        sys.path.remove(str(_SKILL_DIR))


_purge_foreign_modules("errors", "schemas")

from errors import ErrorCode, SkillError
from schemas import (
    DEFAULT_LOCAL_PIPELINE_DIR,
    DEFAULT_PIPELINE_VERSION,
    DEFAULT_PROFILE,
    DEFAULT_REMOTE_PIPELINE,
    JAVA_MIN_VERSION,
    NEXTFLOW_MIN_VERSION,
    SKILL_ALIAS,
    SKILL_NAME,
    SUPPORTED_ALIGNERS,
    SUPPORTED_CONTAMINANT_SCREENING,
    SUPPORTED_PROFILES,
    SUPPORTED_PSEUDO_ALIGNERS,
    SUPPORTED_RIBO_TOOLS,
    SUPPORTED_TRIMMERS,
    SUPPORTED_UMI_TOOLS,
)

_purge_local_modules("errors", "schemas")
_remove_skill_dir_from_sys_path()


def test_all_error_codes_are_string_constants():
    attrs = {k: v for k, v in inspect.getmembers(ErrorCode) if k.isupper()}
    assert attrs, "ErrorCode must define at least one constant"
    assert all(isinstance(v, str) for v in attrs.values())


def test_rnaseq_specific_error_codes_are_present():
    expected = {
        "INVALID_ALIGNER",
        "INVALID_PSEUDO_ALIGNER",
        "INVALID_TRIMMER",
        "INVALID_RIBO_TOOL",
        "INVALID_UMI_TOOL",
        "INCOMPATIBLE_ALIGNER_ARGS",
        "CONFLICTING_FASTA_ARGS",
        "BAM_REPROCESSING_INCOMPLETE",
        "HISAT2_NO_QUANTIFICATION",
        "REFERENCE_PATH_NOT_FOUND",
    }
    assert expected <= {v for k, v in inspect.getmembers(ErrorCode) if k.isupper()}


def test_skill_error_to_dict_shape():
    err = SkillError(
        stage="validation",
        error_code=ErrorCode.INVALID_ALIGNER,
        message="bad aligner",
        fix="choose a supported aligner",
        details={"aligner": "bad"},
    )
    assert err.to_dict() == {
        # Shared cross-wrapper failure contract: status + ok discriminators.
        "status": "error",
        "ok": False,
        "stage": "validation",
        "error_code": "INVALID_ALIGNER",
        "message": "bad aligner",
        "fix": "choose a supported aligner",
        "details": {"aligner": "bad"},
    }


def test_pinned_versions_match_runtime_constants():
    """reproducibility/pinned_versions.json must stay in sync with schemas.py
    constants. Drift here means provenance auditors and the wrapper itself
    disagree about what's supported — a documentation-vs-implementation bug
    that's invisible until someone reads the JSON expecting authoritative info.
    """
    path = Path(__file__).resolve().parent.parent / "reproducibility" / "pinned_versions.json"
    pinned = json.loads(path.read_text(encoding="utf-8"))
    assert pinned["pipeline"]["default_version"] == DEFAULT_PIPELINE_VERSION
    assert int(pinned["minimum_versions"]["java"]) == JAVA_MIN_VERSION
    assert tuple(map(int, pinned["minimum_versions"]["nextflow"].split("."))) == NEXTFLOW_MIN_VERSION
    # Profiles are grouped by purpose inside pinned_versions; flatten before comparing.
    all_profiles_in_json: set[str] = set()
    for group in pinned["supported_profiles"].values():
        all_profiles_in_json.update(group)
    assert all_profiles_in_json == SUPPORTED_PROFILES
    assert set(pinned["supported_aligners"]) == SUPPORTED_ALIGNERS
    assert set(pinned["supported_pseudo_aligners"]) == SUPPORTED_PSEUDO_ALIGNERS
    assert set(pinned["supported_trimmers"]) == SUPPORTED_TRIMMERS
    assert set(pinned["supported_ribo_tools"]) == SUPPORTED_RIBO_TOOLS
    assert set(pinned["supported_umi_tools"]) == SUPPORTED_UMI_TOOLS
    assert set(pinned["supported_contaminant_screening"]) == SUPPORTED_CONTAMINANT_SCREENING


def test_contaminant_screening_enum_value_is_canonical():
    assert SUPPORTED_CONTAMINANT_SCREENING == {"kraken2", "kraken2_bracken", "sylph"}


def test_rnaseq_runtime_constants_match_phase_1_contract():
    assert SKILL_NAME == "nfcore-rnaseq-wrapper"
    assert SKILL_ALIAS == "rnaseq-pipeline"
    assert DEFAULT_REMOTE_PIPELINE == "nf-core/rnaseq"
    assert DEFAULT_PIPELINE_VERSION == "3.26.0"
    assert DEFAULT_LOCAL_PIPELINE_DIR.name == "rnaseq"
    assert DEFAULT_PROFILE == "docker"
    assert JAVA_MIN_VERSION == 17
    assert NEXTFLOW_MIN_VERSION == (25, 4, 3)
    assert SUPPORTED_ALIGNERS == {"star_salmon", "star_rsem", "hisat2", "bowtie2_salmon"}
    assert SUPPORTED_PSEUDO_ALIGNERS == {"salmon", "kallisto"}
    assert SUPPORTED_TRIMMERS == {"trimgalore", "fastp"}
    assert SUPPORTED_PROFILES == {
        # Execution backends
        "docker", "conda", "mamba", "singularity", "apptainer",
        "podman", "shifter", "charliecloud", "wave",
        # Scientific modifiers
        "prokaryotic", "rapid_quant",
        # Architecture modifier
        "arm64",
        # Self-contained nf-core test profiles
        "test", "test_prokaryotic", "test_full",
        "test_full_aws", "test_full_gcp", "test_full_azure", "test_gpu",
        # Debug / development
        "debug", "emulate_amd64",
    }
