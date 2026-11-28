"""`clawbio run` launches a skill that has no SKILLS entry, from its own folder.

The skill gets only the standard contract (--input, --output, --demo); every
extra flag is dropped, so the INT-001 allowlist still holds. Skills that cannot
honour the contract are refused with a reason, never run half-way. Nothing here
reads skills/catalog.json, so the catalog can stop being committed.
"""

import sys
import types

import pytest

from clawbio import cli

FULL = (
    "import argparse\np = argparse.ArgumentParser()\n"
    "p.add_argument('--input')\np.add_argument('--output')\np.add_argument('--demo', action='store_true')\n"
)


@pytest.fixture
def captured_run(monkeypatch):
    seen = {}

    def fake_run(cmd, *args, **kwargs):
        seen["cmd"] = [str(part) for part in cmd]
        seen["cwd"] = kwargs.get("cwd")
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(cli.subprocess, "run", fake_run)
    return seen


@pytest.fixture
def skills_dir(tmp_path, monkeypatch):
    """A skills/ folder a test fills, with no catalog.json in it."""
    root = tmp_path / "skills"
    root.mkdir()
    monkeypatch.setattr(cli, "SKILLS_DIR", root)
    return root


def _skill(root, name, files, version="1.2.3"):
    d = root / name
    d.mkdir()
    (d / "SKILL.md").write_text(f"---\nname: {name}\ndescription: test skill\nmetadata:\n  version: {version}\n---\n")
    for rel, body in files.items():
        (d / rel).write_text(body)
    return d


def test_unregistered_skill_runs_its_demo(captured_run, tmp_path):
    assert "vcf-annotator" not in cli.SKILLS
    result = cli.run_skill("vcf-annotator", demo=True, output_dir=str(tmp_path))

    assert result["success"] is True
    cmd = captured_run["cmd"]
    assert cmd[1].endswith("skills/vcf-annotator/vcf_annotator.py")
    assert cmd[2:] == ["--demo", "--output", str(tmp_path.resolve())]


def test_runs_with_no_catalog_file(skills_dir, captured_run, tmp_path):
    _skill(skills_dir, "my-skill", {"my_skill.py": FULL})
    result = cli.run_skill("my-skill", demo=True, output_dir=str(tmp_path / "out"))

    assert not (skills_dir / "catalog.json").exists()
    assert result["success"] is True
    assert captured_run["cmd"][1] == str((skills_dir / "my-skill" / "my_skill.py").resolve())


def test_unregistered_skill_takes_input_as_an_absolute_path(captured_run, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "x.vcf").write_text("##fileformat=VCFv4.2\n")
    result = cli.run_skill("vcf-annotator", input_path="x.vcf", output_dir=str(tmp_path / "out"))

    assert result["success"] is True
    assert captured_run["cmd"][2:4] == ["--input", str((tmp_path / "x.vcf").resolve())]


def test_extra_flags_are_dropped_for_an_unregistered_skill(captured_run, tmp_path):
    cli.run_skill(
        "vcf-annotator", demo=True, output_dir=str(tmp_path),
        extra_args=["--upload", "https://example.org", "--threads", "8"],
    )
    cmd = captured_run["cmd"]
    assert "--upload" not in cmd and "https://example.org" not in cmd and "--threads" not in cmd


def test_demo_keeps_the_skills_extra_demo_arguments(captured_run, tmp_path):
    # ancestry-risk-profiler's demo needs the synthetic patient's ancestry.
    cli.run_skill("ancestry-risk-profiler", demo=True, output_dir=str(tmp_path))
    assert captured_run["cmd"][2:5] == ["--demo", "--ancestry", "SAS"]


def test_runs_from_the_script_folder(captured_run, tmp_path):
    cli.run_skill("vcf-annotator", demo=True, output_dir=str(tmp_path))
    assert captured_run["cwd"].endswith("skills/vcf-annotator")


NO_OUT = "import argparse\np = argparse.ArgumentParser()\np.add_argument('--input')\np.add_argument('--demo')\n"
NO_IN = "import argparse\np = argparse.ArgumentParser()\np.add_argument('--output')\np.add_argument('--demo')\n"
NO_DEMO = "import argparse\np = argparse.ArgumentParser()\np.add_argument('--output')\np.add_argument('--input')\n"


@pytest.mark.parametrize(
    "files, kwargs, reason",
    [
        ({"r_skill.R": "x\n"}, {"demo": True}, "Python"),
        ({"s.py": NO_OUT}, {"demo": True}, "--output"),
        ({"s.py": NO_IN}, {"input_path": "x.vcf"}, "--input"),
        ({"s.py": NO_IN}, {"profile_path": "PT001.json"}, "--input"),  # the profile's file goes in as --input
        ({"s.py": NO_DEMO}, {"demo": True}, "--demo"),
    ],
)
def test_a_skill_that_breaks_the_contract_is_refused(skills_dir, captured_run, tmp_path, files, kwargs, reason):
    _skill(skills_dir, "bad", files)
    result = cli.run_skill("bad", output_dir=str(tmp_path / "out"), **kwargs)

    assert result["success"] is False
    assert reason in result["stderr"]
    assert "SKILL.md" in result["stderr"]
    assert "cmd" not in captured_run  # never started


def test_a_skill_listed_as_having_no_demo_is_refused(skills_dir, captured_run, tmp_path, monkeypatch):
    _skill(skills_dir, "nodemo", {"nodemo.py": FULL})
    monkeypatch.setattr(cli, "NO_DEMO", {"nodemo"})
    result = cli.run_skill("nodemo", demo=True, output_dir=str(tmp_path / "out"))

    assert result["success"] is False
    assert "demo" in result["stderr"]
    assert "cmd" not in captured_run


def test_a_script_linked_from_outside_the_skill_folder_is_refused(skills_dir, captured_run, tmp_path):
    outside = tmp_path / "elsewhere.py"
    outside.write_text(FULL)
    d = _skill(skills_dir, "escape", {})
    (d / "escape.py").symlink_to(outside)
    result = cli.run_skill("escape", demo=True, output_dir=str(tmp_path / "out"))

    assert result["success"] is False
    assert "cmd" not in captured_run


def test_a_skill_with_no_script_still_points_at_its_skill_md(captured_run, tmp_path):
    result = cli.run_skill("bgpt-mcp", demo=True, output_dir=str(tmp_path))

    assert result["success"] is False
    assert "SKILL.md" in result["stderr"]
    assert "cmd" not in captured_run


def test_the_audit_record_names_the_skill_md_version(skills_dir):
    _skill(skills_dir, "v-skill", {"v_skill.py": FULL}, version="4.5.6")
    assert cli._skill_version("v-skill") == "4.5.6"


def test_cli_run_launches_an_unregistered_skill(captured_run, monkeypatch, tmp_path):
    monkeypatch.setattr(sys, "argv", ["clawbio", "run", "vcf-annotator", "--demo", "--output", str(tmp_path)])
    with pytest.raises(SystemExit) as exit_info:
        cli.main()
    assert exit_info.value.code in (0, None)
    assert captured_run["cmd"][1].endswith("vcf_annotator.py")


def _list_line(capsys, name):
    cli.list_skills()
    lines = capsys.readouterr().out.splitlines()
    i = next(i for i, line in enumerate(lines) if f" {name} " in f"{line} ")
    return lines[i + 1]


def test_list_says_an_unregistered_skill_runs_through_clawbio_run(capsys):
    assert "standard flags only" in _list_line(capsys, "vcf-annotator")


def test_list_says_a_refused_skill_does_not(capsys):
    # turingdb-graph's script defines no --output.
    assert "not via clawbio run" in _list_line(capsys, "turingdb-graph")
