"""Tests for scaffold_skill.py.

A scaffolded skill must produce the reproducibility bundle its own SKILL.md
promises. Before this suite, the scaffolder emitted SKILL.md prose describing
`reproducibility/` and generated no code to write it, so every new skill
started life with documented-but-absent outputs (see #372).
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import scaffold_skill  # noqa: E402


@pytest.fixture
def scaffolded(tmp_path, monkeypatch):
    """Scaffold a throwaway skill into tmp_path and return its directory."""
    skills_dir = tmp_path / "skills"
    skills_dir.mkdir()
    monkeypatch.setattr(scaffold_skill, "SKILLS_DIR", skills_dir)
    scaffold_skill.scaffold(
        "demo-repro-skill",
        "Synthetic skill used to test the scaffolder.",
        bench_dir=tmp_path / "bench",
    )
    return skills_dir / "demo-repro-skill"


def _run_demo(skill_dir: Path, out_dir: Path):
    script = skill_dir / "demo_repro_skill.py"
    # A real skill resolves the repo root as parents[2]; a skill scaffolded into
    # tmp_path cannot, so put it on PYTHONPATH for the subprocess instead.
    env = dict(os.environ, PYTHONPATH=str(REPO_ROOT))
    return subprocess.run(
        [sys.executable, str(script), "--demo", "--output", str(out_dir)],
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
        env=env,
    )


class TestReproducibilityBundle:
    def test_demo_run_writes_bundle(self, scaffolded, tmp_path):
        out = tmp_path / "out"
        result = _run_demo(scaffolded, out)
        assert result.returncode == 0, f"demo run failed: {result.stderr}"
        repro = out / "reproducibility"
        assert repro.is_dir(), "scaffolded skill wrote no reproducibility/ directory"
        for artefact in ("commands.sh", "environment.yml", "checksums.sha256"):
            assert (repro / artefact).exists(), f"missing {artefact}"

    def test_bundle_uses_shared_layer(self, scaffolded):
        src = (scaffolded / "demo_repro_skill.py").read_text()
        assert "clawbio.common.reproducibility" in src, (
            "generated skill must import the shared layer, not hand-roll writers"
        )
        for helper in ("write_commands_sh", "write_environment_yml", "write_checksums"):
            assert helper in src, f"generated skill does not call {helper}"

    def test_checksums_resolve_from_output_dir(self, scaffolded, tmp_path):
        """Labels must be relative to output_dir, per docs/reproducibility.md."""
        out = tmp_path / "out"
        assert _run_demo(scaffolded, out).returncode == 0
        lines = (out / "reproducibility" / "checksums.sha256").read_text().strip().splitlines()
        assert lines, "checksum manifest is empty"
        for line in lines:
            _, label = line.split("  ", 1)
            assert (out / label).exists(), f"label {label!r} does not resolve from output_dir"

    def test_skill_md_promises_the_bundle(self, scaffolded):
        """The Output Structure tree must list what the skill actually writes."""
        skill_md = (scaffolded / "SKILL.md").read_text()
        tree = skill_md.split("## Output Structure", 1)[1].split("```")[1]
        for artefact in ("commands.sh", "environment.yml", "checksums.sha256"):
            assert artefact in tree, f"Output Structure tree omits {artefact}"

    def test_output_contract_test_is_generated(self, scaffolded):
        tests = (scaffolded / "tests" / "test_demo_repro_skill.py").read_text()
        assert "class TestOutputContract" in tests, (
            "scaffolded tests must include the output-contract guard"
        )


class TestConformanceChecklist:
    def test_checklist_includes_reproducibility(self, scaffolded, capsys, tmp_path, monkeypatch):
        """The printed checklist mirrors CLAUDE.md and must cover the bundle."""
        monkeypatch.setattr(scaffold_skill, "SKILLS_DIR", scaffolded.parent)
        scaffold_skill.scaffold(
            "demo-repro-skill",
            "Synthetic skill used to test the scaffolder.",
            force=True,
            bench_dir=tmp_path / "bench2",
        )
        out = capsys.readouterr().out
        assert "reproducibility" in out.lower(), (
            "conformance checklist does not mention the reproducibility bundle"
        )
