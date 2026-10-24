"""Tests for clawbio.common.report."""

import json
import tempfile
from pathlib import Path

import pytest
from hypothesis import given, strategies as st

from clawbio.common.report import (
    generate_report_header,
    skill_version,
    version_from_skill_md,
    write_audit_log,
    write_result_json,
)


def test_write_audit_log_appends_jsonl_record(tmp_path):
    log_path = tmp_path / "audit.log"
    write_audit_log(
        skill="pharmgx",
        version="0.2.0",
        input_checksum="abc123",
        output_dir=str(tmp_path / "out"),
        log_path=log_path,
    )
    records = [json.loads(l) for l in log_path.read_text().splitlines()]
    assert len(records) == 1
    assert records[0]["skill"] == "pharmgx"
    assert records[0]["version"] == "0.2.0"
    assert records[0]["input_checksum"] == "abc123"
    assert "timestamp" in records[0]


def test_write_audit_log_appends_multiple_runs(tmp_path):
    log_path = tmp_path / "audit.log"
    write_audit_log(skill="pharmgx", version="0.2.0", input_checksum="aaa", output_dir="/tmp/a", log_path=log_path)
    write_audit_log(skill="gwas-prs", version="0.1.0", input_checksum="bbb", output_dir="/tmp/b", log_path=log_path)
    records = [json.loads(l) for l in log_path.read_text().splitlines()]
    assert len(records) == 2
    assert records[1]["skill"] == "gwas-prs"


def test_write_audit_log_creates_parent_dir(tmp_path):
    log_path = tmp_path / "nested" / "audit.log"
    write_audit_log(skill="pharmgx", version="0.2.0", input_checksum="x", output_dir="/tmp", log_path=log_path)
    assert log_path.exists()


def test_write_result_json_includes_datasets(tmp_path):
    write_result_json(
        output_dir=tmp_path,
        skill="pharmgx",
        version="0.2.0",
        summary={},
        data={},
        datasets={"pgx_panel": "2024-01", "cpic_guidelines": "v1.19"},
    )
    result = json.loads((tmp_path / "result.json").read_text())
    assert result["datasets"] == {"pgx_panel": "2024-01", "cpic_guidelines": "v1.19"}


def test_write_result_json_datasets_defaults_to_empty(tmp_path):
    write_result_json(
        output_dir=tmp_path,
        skill="pharmgx",
        version="0.2.0",
        summary={},
        data={},
    )
    result = json.loads((tmp_path / "result.json").read_text())
    assert result["datasets"] == {}


def test_generate_report_header_includes_skill_version():
    result = generate_report_header(
        title="Test Report",
        skill_name="pharmgx",
        skill_version="1.2.3",
    )
    assert "**Version**: 1.2.3" in result


def test_generate_report_header_allows_an_omitted_skill_version():
    """`skill_version` is optional: 1a91430 gave it a default so existing
    callers kept working. This test asserted a TypeError instead, and CI ran
    no test under clawbio/common/tests/, so it was red on main unnoticed.
    Omitting the version still renders the header, with the Version line
    blank, so callers that can pass a version should."""
    result = generate_report_header(title="Test Report", skill_name="pharmgx")
    assert "# Test Report" in result
    assert "**Version**: \n" in result


def test_generate_report_header_includes_skill_name():
    result = generate_report_header(
        title="Test Report",
        skill_name="pharmgx",
        skill_version="0.1.0",
    )
    assert "**Skill**: pharmgx" in result
    assert "**Version**: 0.1.0" in result


SEMVER = st.tuples(*[st.integers(0, 999)] * 3).map(lambda t: ".".join(map(str, t)))


def _skill_md(version_line: str, body: str = "") -> str:
    return f"---\nname: demo-skill\nmetadata:\n  author: x\n{version_line}\n---\n\n# Demo\n{body}"


@given(version=SEMVER, decoy=SEMVER)
def test_version_from_skill_md_reads_frontmatter_never_the_body(version, decoy):
    text = _skill_md(f"  version: {version}", body=f"\nmetadata:\n  version: {decoy}\n")
    assert version_from_skill_md(text) == version


def test_version_from_skill_md_without_a_version_raises():
    with pytest.raises(ValueError, match="metadata.version"):
        version_from_skill_md(_skill_md("  tags: [a]"))


def test_skill_version_reads_the_skill_md_next_to_the_script(tmp_path):
    (tmp_path / "SKILL.md").write_text(_skill_md("  version: 1.4.2"))
    assert skill_version(tmp_path / "my_skill.py") == "1.4.2"


def test_skill_version_walks_up_from_a_subfolder_module(tmp_path):
    (tmp_path / "SKILL.md").write_text(_skill_md("  version: 0.3.0"))
    (tmp_path / "core").mkdir()
    assert skill_version(tmp_path / "core" / "report.py") == "0.3.0"


def test_skill_version_without_a_skill_md_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="SKILL.md"):
        skill_version(tmp_path / "orphan.py")
