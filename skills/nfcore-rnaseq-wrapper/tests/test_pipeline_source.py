from __future__ import annotations

from pathlib import Path
import shutil
import subprocess
import sys

import pytest

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


_purge_foreign_modules("errors", "schemas", "pipeline_source")

from errors import SkillError
from pipeline_source import resolve_pipeline_source

_purge_local_modules("errors", "schemas", "pipeline_source")
_remove_skill_dir_from_sys_path()

_requires_git = pytest.mark.skipif(shutil.which("git") is None, reason="git not on PATH")


def _make_valid_local_checkout(path: Path, *, manifest_version: str | None = None) -> None:
    path.mkdir(parents=True, exist_ok=True)
    (path / "main.nf").write_text("// main", encoding="utf-8")
    if manifest_version is None:
        config_text = "// config"
    else:
        config_text = (
            "manifest {\n"
            "    name = 'nf-core/rnaseq'\n"
            f"    version = '{manifest_version}'\n"
            "    nextflowVersion = '!>=25.04.3'\n"
            "}\n"
        )
    (path / "nextflow.config").write_text(config_text, encoding="utf-8")
    assets = path / "assets"
    assets.mkdir()
    (assets / "schema_input.json").write_text("{}", encoding="utf-8")


def test_local_checkout_parses_manifest_version(tmp_path):
    local = tmp_path / "rnaseq"
    _make_valid_local_checkout(local, manifest_version="3.26.0")
    result = resolve_pipeline_source(requested_version="3.26.0", local_pipeline_dir=local)
    assert result["manifest_version"] == "3.26.0"


def test_local_checkout_manifest_version_ignores_nextflow_version(tmp_path):
    local = tmp_path / "rnaseq"
    _make_valid_local_checkout(local, manifest_version="3.20.0")
    result = resolve_pipeline_source(requested_version="3.26.0", local_pipeline_dir=local)
    # Must read manifest.version, not manifest.nextflowVersion.
    assert result["manifest_version"] == "3.20.0"


def test_local_checkout_manifest_version_ignores_custom_config_version(tmp_path):
    """A real nf-core config declares `custom_config_version` above the manifest block.

    The synthetic fixture used by the tests above does not, which is why the previous
    regex could read `master` from every real checkout while these tests stayed green.
    """
    local = tmp_path / "rnaseq"
    local.mkdir(parents=True)
    (local / "main.nf").write_text("// main", encoding="utf-8")
    (local / "nextflow.config").write_text(
        "params {\n"
        "    version                    = false\n"
        "    custom_config_version      = 'master'\n"
        "    custom_config_base         = \"https://example.invalid/${params.custom_config_version}\"\n"
        "}\n"
        "manifest {\n"
        "    name            = 'nf-core/rnaseq'\n"
        "    version         = '3.26.0'\n"
        "    nextflowVersion = '!>=25.04.3'\n"
        "}\n",
        encoding="utf-8",
    )
    assets = local / "assets"
    assets.mkdir()
    (assets / "schema_input.json").write_text("{}", encoding="utf-8")
    result = resolve_pipeline_source(requested_version="3.26.0", local_pipeline_dir=local)
    assert result["manifest_version"] == "3.26.0"


def test_local_checkout_manifest_version_empty_when_absent(tmp_path):
    local = tmp_path / "rnaseq"
    _make_valid_local_checkout(local)  # config has no manifest block
    result = resolve_pipeline_source(requested_version="3.26.0", local_pipeline_dir=local)
    assert result["manifest_version"] == ""


def test_remote_manifest_version_mirrors_requested(tmp_path):
    result = resolve_pipeline_source(requested_version="3.26.0", local_pipeline_dir=tmp_path / "missing")
    assert result["manifest_version"] == "3.26.0"


def test_resolves_to_remote_when_no_local_dir(tmp_path):
    absent = tmp_path / "nonexistent_rnaseq"
    result = resolve_pipeline_source(
        requested_version="3.26.0",
        local_pipeline_dir=absent,
    )
    assert result["source_kind"] == "remote_repo"
    assert result["source_ref"] == "nf-core/rnaseq"
    assert result["resolved_version"] == "3.26.0"
    assert result["dirty"] is False
    assert result["branch"] == ""


def test_local_checkout_with_whitespace_path_records_rejection_diagnostic(tmp_path):
    # A sibling checkout whose path contains whitespace is rejected (Docker on macOS
    # cannot reliably run scripts from spaced paths) and the wrapper falls back to the
    # remote pipeline. The rejection must be recorded so provenance/upstream.json can
    # explain why the local checkout was not used (provenance.build_upstream_payload
    # reads these keys).
    local = tmp_path / "has space" / "rnaseq"
    local.mkdir(parents=True)
    result = resolve_pipeline_source(requested_version="3.26.0", local_pipeline_dir=local)
    assert result["source_kind"] == "remote_repo"
    assert result["local_attempted"] == str(local.resolve())
    assert "whitespace" in str(result["local_rejected_reason"]).lower()


def test_remote_without_local_rejection_has_no_diagnostic(tmp_path):
    result = resolve_pipeline_source(
        requested_version="3.26.0",
        local_pipeline_dir=tmp_path / "absent_rnaseq",
    )
    assert "local_attempted" not in result


def test_resolves_to_local_when_valid_checkout(tmp_path):
    local = tmp_path / "rnaseq"
    _make_valid_local_checkout(local)
    result = resolve_pipeline_source(
        requested_version="3.26.0",
        local_pipeline_dir=local,
    )
    assert result["source_kind"] == "local_checkout"
    assert result["source_ref"] == str(local.resolve())
    assert result["resolved_version"] == "3.26.0"
    assert result["branch"] == ""
    assert result["dirty"] is False


def test_raises_when_local_dir_missing_required_files(tmp_path):
    local = tmp_path / "rnaseq"
    local.mkdir()
    (local / "main.nf").write_text("// main", encoding="utf-8")
    with pytest.raises(SkillError) as exc:
        resolve_pipeline_source(
            requested_version="3.26.0",
            local_pipeline_dir=local,
        )
    assert exc.value.error_code == "PIPELINE_SOURCE_INVALID"
    assert exc.value.details["missing_files"] == ["nextflow.config", "assets/schema_input.json"]
    assert exc.value.details["path"] == str(local.resolve())


def test_local_checkout_result_contains_expected_keys(tmp_path):
    local = tmp_path / "rnaseq"
    _make_valid_local_checkout(local)
    result = resolve_pipeline_source(
        requested_version="3.26.0",
        local_pipeline_dir=local,
    )
    assert set(result) == {"source_kind", "source_ref", "resolved_version", "manifest_version", "branch", "dirty"}


def test_remote_result_contains_expected_keys(tmp_path):
    result = resolve_pipeline_source(
        requested_version="3.26.0",
        local_pipeline_dir=tmp_path / "missing",
    )
    assert set(result) == {"source_kind", "source_ref", "resolved_version", "manifest_version", "branch", "dirty"}


@_requires_git
def test_local_checkout_uses_git_commit_when_available(tmp_path):
    local = tmp_path / "rnaseq"
    _make_valid_local_checkout(local)
    subprocess.run(["git", "init"], cwd=local, check=True, capture_output=True, text=True)
    subprocess.run(["git", "config", "user.email", "test@example.com"], cwd=local, check=True)
    subprocess.run(["git", "config", "user.name", "Test User"], cwd=local, check=True)
    subprocess.run(["git", "add", "."], cwd=local, check=True)
    subprocess.run(["git", "commit", "-m", "init"], cwd=local, check=True, capture_output=True, text=True)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=local, check=True, capture_output=True, text=True).stdout.strip()
    result = resolve_pipeline_source(
        requested_version="3.26.0",
        local_pipeline_dir=local,
    )
    assert result["resolved_version"] == commit
    assert result["dirty"] is False
