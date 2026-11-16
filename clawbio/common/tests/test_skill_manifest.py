"""Tests for clawbio.common.skill_manifest: the one reader of SKILL.md frontmatter."""

import sys
from pathlib import Path

import pytest
from hypothesis import given, strategies as st

from clawbio.common.skill_manifest import SkillManifestError, read_skill_manifest

SKILLS_DIR = Path(__file__).resolve().parents[3] / "skills"

SEMVER = st.tuples(*[st.integers(0, 999)] * 3).map(lambda t: ".".join(map(str, t)))
LICENCE = st.sampled_from(["MIT", "Apache-2.0", "GPL-3.0", "cc-by-nc-sa-4.0"])


def _skill(tmp_path: Path, frontmatter: str, body: str = "") -> Path:
    (tmp_path / "SKILL.md").write_text(f"---\n{frontmatter}\n---\n\n# Demo\n{body}", encoding="utf-8")
    return tmp_path


@given(version=SEMVER, decoy=SEMVER, licence=LICENCE, decoy_licence=LICENCE)
def test_reads_frontmatter_never_the_body(tmp_path_factory, version, decoy, licence, decoy_licence):
    folder = _skill(
        tmp_path_factory.mktemp("s"),
        f"name: demo\nlicense: {licence}\nmetadata:\n  version: {version}",
        body=f"\n---\nlicense: {decoy_licence}\nmetadata:\n  version: {decoy}\n---\n",
    )
    manifest = read_skill_manifest(folder)
    assert manifest.version == version
    assert manifest.license == licence


def test_reads_named_fields_from_metadata(tmp_path):
    folder = _skill(tmp_path, "\n".join([
        "name: demo",
        "description: A demo skill",
        "license: MIT",
        "metadata:",
        "  version: 1.2.3",
        "  tags: [a, b]",
        "  model_license: cc-by-nc-sa-4.0",
        "  data_license: CC0-1.0",
    ]))
    m = read_skill_manifest(folder)
    assert (m.name, m.description, m.license, m.version) == ("demo", "A demo skill", "MIT", "1.2.3")
    assert m.tags == ("a", "b")
    assert (m.model_license, m.data_license) == ("cc-by-nc-sa-4.0", "CC0-1.0")


def test_absent_artefact_licences_are_empty_not_permissive(tmp_path):
    """No third-party artefact declared is not a claim that one is permissive."""
    m = read_skill_manifest(_skill(tmp_path, "name: demo\nlicense: MIT\nmetadata:\n  version: 0.1.0"))
    assert (m.model_license, m.data_license, m.tags) == ("", "", ())


def test_top_level_keys_other_than_the_spec_four_are_not_read(tmp_path):
    """Legacy top-level placement is gone: every skill nests under `metadata:`."""
    m = read_skill_manifest(_skill(tmp_path, "name: demo\nversion: 9.9.9\nmodel_license: x\nmetadata:\n  version: 1.0.0"))
    assert m.version == "1.0.0"
    assert m.model_license == ""


def test_missing_version_raises_only_when_version_is_read(tmp_path):
    m = read_skill_manifest(_skill(tmp_path, "name: demo\nlicense: MIT\nmetadata:\n  tags: [x]"))
    assert (m.license, m.tags) == ("MIT", ("x",))
    with pytest.raises(SkillManifestError, match="metadata.version"):
        m.version


@pytest.mark.parametrize("text", [
    "# no frontmatter\n",
    "---\n- just\n- a list\n---\n",
    "---\nname: s\ndescription: foo: bar\n---\n",  # invalid YAML: callers catch one error type
])
def test_missing_or_non_mapping_frontmatter_raises_on_read(tmp_path, text):
    (tmp_path / "SKILL.md").write_text(text, encoding="utf-8")
    with pytest.raises(SkillManifestError):
        read_skill_manifest(tmp_path)


def test_walks_up_from_a_file_in_a_subfolder(tmp_path):
    _skill(tmp_path, "name: demo\nmetadata:\n  version: 0.3.0")
    (tmp_path / "core").mkdir()
    assert read_skill_manifest(tmp_path / "core" / "report.py").version == "0.3.0"


def test_without_a_skill_md_raises(tmp_path):
    with pytest.raises(SkillManifestError, match="SKILL.md"):
        read_skill_manifest(tmp_path / "orphan.py")


def test_has_no_regex_fallback_without_pyyaml(tmp_path, monkeypatch):
    """Without PyYAML it fails loudly rather than falling back to regexes that
    backtrack catastrophically on crafted frontmatter (ReDoS)."""
    folder = _skill(tmp_path, "name: x\ndependencies:\n" + "  a: 1\n" * 5000)
    monkeypatch.setitem(sys.modules, "yaml", None)
    with pytest.raises(ImportError):
        read_skill_manifest(folder)


def test_every_shipped_skill_has_a_readable_manifest_with_a_version():
    for skill_md in sorted(SKILLS_DIR.glob("*/SKILL.md")):
        assert read_skill_manifest(skill_md.parent).version, skill_md.parent.name


def test_folds_a_block_scalar_description(tmp_path):
    """Moved from bigquery-public's tests, which pinned the old catalog parser."""
    m = read_skill_manifest(_skill(tmp_path, "name: demo\ndescription: >-\n  First line\n  second line\nmetadata:\n  tags: [one, two]"))
    assert m.description == "First line second line"
    assert m.tags == ("one", "two")
