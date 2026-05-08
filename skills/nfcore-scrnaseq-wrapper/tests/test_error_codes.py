from __future__ import annotations

from pathlib import Path
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from errors import ErrorCode, SkillError
from schemas import (
    DEFAULT_PIPELINE_VERSION,
    JAVA_MIN_VERSION,
    NEXTFLOW_MIN_VERSION,
    SUPPORTED_PRESETS,
    SUPPORTED_PROFILES,
)


def test_all_error_codes_are_string_constants():
    import inspect
    attrs = {k: v for k, v in inspect.getmembers(ErrorCode) if not k.startswith("_")}
    assert attrs, "ErrorCode must define at least one constant"
    assert all(isinstance(v, str) for v in attrs.values()), "All ErrorCode values must be strings"


def test_skill_error_accepts_error_code_constant():
    err = SkillError(
        stage="preflight",
        error_code=ErrorCode.MISSING_JAVA,
        message="java missing",
        fix="install java",
    )
    assert err.error_code == "MISSING_JAVA"


def test_error_response_shape():
    payload = SkillError(
        stage="preflight",
        error_code="MISSING_JAVA",
        message="java missing",
        fix="install java",
    ).to_dict()
    assert payload["ok"] is False
    assert payload["error_code"] == "MISSING_JAVA"


def test_skill_error_to_dict():
    err = SkillError("validation", "INVALID_SAMPLESHEET", "bad", "fix")
    payload = err.to_dict()
    assert payload["stage"] == "validation"


def test_pinned_versions_match_runtime_constants():
    path = Path(__file__).resolve().parent.parent / "reproducibility" / "pinned_versions.json"
    pinned = json.loads(path.read_text(encoding="utf-8"))
    assert pinned["pipeline"]["default_version"] == DEFAULT_PIPELINE_VERSION
    assert int(pinned["minimum_versions"]["java"]) == JAVA_MIN_VERSION
    assert tuple(map(int, pinned["minimum_versions"]["nextflow"].split("."))) == NEXTFLOW_MIN_VERSION
    assert set(pinned["supported_profiles"]) == SUPPORTED_PROFILES
    assert set(pinned["supported_presets"]) == SUPPORTED_PRESETS


def test_compatibility_policy_matches_runtime_execution_defaults():
    path = Path(__file__).resolve().parent.parent / "reproducibility" / "compatibility_policy.json"
    policy = json.loads(path.read_text(encoding="utf-8"))
    assert policy["resume_policy"]["requires_same_params_checksum"] is True
    assert policy["execution_policy"]["work_dir"] == "output/upstream/work"
