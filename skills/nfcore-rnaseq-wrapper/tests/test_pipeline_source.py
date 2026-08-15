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


# Mirrors the shape of a real nf-core nextflow.config: the standard template sets
# `custom_config_version` (and pipelines such as sarek add tool-version keys like
# `vep_version`) in the params block, ~300 lines ABOVE the manifest block. A
# file-wide version scan therefore reads one of those instead of manifest.version
# (ClawBio#333). Keeping this preamble in the shared fixture means every test in
# this module exercises the real config shape rather than a thin idealised one.
_NFCORE_PARAMS_PREAMBLE = (
    "params {\n"
    "    version                 = false\n"
    "    vep_version             = \"111.0-0\"\n"
    "    custom_config_version   = 'master'\n"
    "    custom_config_base      = \"https://raw.githubusercontent.com/nf-core/configs/${params.custom_config_version}\"\n"
    "}\n"
)


def _make_valid_local_checkout(path: Path, *, manifest_version: str | None = None) -> None:
    path.mkdir(parents=True, exist_ok=True)
    (path / "main.nf").write_text("// main", encoding="utf-8")
    config_text = _NFCORE_PARAMS_PREAMBLE
    if manifest_version is not None:
        config_text += (
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
    """ClawBio#333: `custom_config_version = 'master'` precedes the manifest block.

    A file-wide scan returns 'master' and the pinned-version gate then rejects a
    genuinely correct checkout, pushing users towards
    --allow-pipeline-version-override, which would also mask a real mismatch.
    """
    local = tmp_path / "rnaseq"
    _make_valid_local_checkout(local, manifest_version="3.26.0")
    config_text = (local / "nextflow.config").read_text(encoding="utf-8")
    assert "custom_config_version   = 'master'" in config_text, "fixture must carry the real trap"
    assert config_text.index("custom_config_version") < config_text.index("manifest {")

    result = resolve_pipeline_source(requested_version="3.26.0", local_pipeline_dir=local)
    assert result["manifest_version"] == "3.26.0"


def test_local_checkout_manifest_version_ignores_tool_version_keys(tmp_path):
    """nf-core/sarek carries `vep_version = "111.0-0"` above the manifest block.

    Worse than the 'master' case: '111.0-0' reads as a plausible pipeline version,
    so the resulting mismatch error misdirects debugging rather than looking absurd.
    """
    local = tmp_path / "rnaseq"
    _make_valid_local_checkout(local, manifest_version="3.26.0")
    result = resolve_pipeline_source(requested_version="3.26.0", local_pipeline_dir=local)
    assert result["manifest_version"] != "111.0-0"
    assert result["manifest_version"] == "3.26.0"


def test_local_checkout_manifest_version_ignores_dotted_assignment(tmp_path):
    """`params.version = '9.9.9'` is not the manifest version.

    A leading-boundary lookbehind alone does not exclude this, because the
    preceding character is '.', which is why the fix scopes to the manifest block.
    """
    local = tmp_path / "rnaseq"
    _make_valid_local_checkout(local, manifest_version="3.26.0")
    config = local / "nextflow.config"
    config.write_text("params.version = '9.9.9'\n" + config.read_text(encoding="utf-8"), encoding="utf-8")
    result = resolve_pipeline_source(requested_version="3.26.0", local_pipeline_dir=local)
    assert result["manifest_version"] == "3.26.0"


def test_local_checkout_manifest_version_ignores_commented_out_version(tmp_path):
    """A commented-out version above the manifest block must not win."""
    local = tmp_path / "rnaseq"
    _make_valid_local_checkout(local, manifest_version="3.26.0")
    config = local / "nextflow.config"
    config.write_text("// version = '0.0.1-dev'\n" + config.read_text(encoding="utf-8"), encoding="utf-8")
    result = resolve_pipeline_source(requested_version="3.26.0", local_pipeline_dir=local)
    assert result["manifest_version"] == "3.26.0"


def test_local_checkout_manifest_version_ignores_comment_inside_manifest(tmp_path):
    """A commented-out version *inside* the manifest block must not win either."""
    local = tmp_path / "rnaseq"
    _make_valid_local_checkout(local, manifest_version="3.26.0")
    config = local / "nextflow.config"
    config.write_text(
        config.read_text(encoding="utf-8").replace(
            "manifest {\n", "manifest {\n    // version = '0.0.1-dev'\n"
        ),
        encoding="utf-8",
    )
    result = resolve_pipeline_source(requested_version="3.26.0", local_pipeline_dir=local)
    assert result["manifest_version"] == "3.26.0"


def test_local_checkout_manifest_version_handles_nested_braces(tmp_path):
    """The manifest scan must be brace-balanced, not stop at the first '}'."""
    local = tmp_path / "rnaseq"
    _make_valid_local_checkout(local, manifest_version="3.26.0")
    config = local / "nextflow.config"
    config.write_text(
        config.read_text(encoding="utf-8").replace(
            "manifest {\n    name = 'nf-core/rnaseq'\n",
            "manifest {\n    name = 'nf-core/rnaseq'\n    defaultBranch = \"${x ? 'a' : 'b'}\"\n",
        ),
        encoding="utf-8",
    )
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
