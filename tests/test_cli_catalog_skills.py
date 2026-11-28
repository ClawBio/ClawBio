"""`clawbio run` launches a catalog skill that has no SKILLS entry.

The skill gets only the standard contract (--input, --output, --demo); every
extra flag is dropped, so the INT-001 allowlist still holds. Skills that cannot
honour the contract are refused with a reason, never run half-way.
"""

import sys
import types

import pytest

from clawbio import cli


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
def fake_catalog(monkeypatch):
    """Replace the catalog with entries a test defines."""
    entries = {}
    monkeypatch.setattr(cli, "_catalog_entry", lambda name: entries.get(name))
    return entries


def _entry(name, script, flags, demo_command=None, version="1.2.3"):
    return {
        "name": name,
        "script": script,
        "run_flags": flags,
        "demo_command": demo_command,
        "description": "test skill",
        "version": version,
    }


def test_unregistered_skill_runs_its_demo(captured_run, tmp_path):
    assert "vcf-annotator" not in cli.SKILLS
    result = cli.run_skill("vcf-annotator", demo=True, output_dir=str(tmp_path))

    assert result["success"] is True
    cmd = captured_run["cmd"]
    assert cmd[1].endswith("skills/vcf-annotator/vcf_annotator.py")
    assert cmd[2:] == ["--demo", "--output", str(tmp_path.resolve())]


def test_unregistered_skill_takes_input_as_an_absolute_path(captured_run, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "x.vcf").write_text("##fileformat=VCFv4.2\n")
    out = tmp_path / "out"
    result = cli.run_skill("vcf-annotator", input_path="x.vcf", output_dir=str(out))

    assert result["success"] is True
    cmd = captured_run["cmd"]
    assert cmd[2:4] == ["--input", str((tmp_path / "x.vcf").resolve())]


def test_extra_flags_are_dropped_for_an_unregistered_skill(captured_run, tmp_path):
    cli.run_skill(
        "vcf-annotator", demo=True, output_dir=str(tmp_path),
        extra_args=["--upload", "https://example.org", "--threads", "8"],
    )
    cmd = captured_run["cmd"]
    assert "--upload" not in cmd and "https://example.org" not in cmd and "--threads" not in cmd


def test_demo_keeps_the_catalog_demo_arguments(captured_run, tmp_path):
    # The catalog demo for ancestry-risk-profiler supplies the synthetic patient's ancestry.
    cli.run_skill("ancestry-risk-profiler", demo=True, output_dir=str(tmp_path))
    cmd = captured_run["cmd"]
    assert cmd[2:5] == ["--demo", "--ancestry", "SAS"]


def test_runs_from_the_script_folder(captured_run, tmp_path):
    cli.run_skill("vcf-annotator", demo=True, output_dir=str(tmp_path))
    assert captured_run["cwd"].endswith("skills/vcf-annotator")


@pytest.mark.parametrize(
    "entry, kwargs, reason",
    [
        (_entry("r-skill", "skills/r-skill/r_skill.R", []), {"demo": True}, "Python"),
        (_entry("no-out", "skills/no-out/no_out.py", ["--demo", "--input"], "python skills/no-out/no_out.py --demo"),
         {"demo": True}, "--output"),
        (_entry("no-in", "skills/no-in/no_in.py", ["--demo", "--output"], "python skills/no-in/no_in.py --demo"),
         {"input_path": "x.vcf"}, "--input"),
        (_entry("no-demo", "skills/no-demo/no_demo.py", ["--input", "--output"]), {"demo": True}, "--demo"),
    ],
)
def test_a_skill_that_breaks_the_contract_is_refused(fake_catalog, captured_run, tmp_path, entry, kwargs, reason):
    fake_catalog[entry["name"]] = entry
    result = cli.run_skill(entry["name"], output_dir=str(tmp_path / "out"), **kwargs)

    assert result["success"] is False
    assert reason in result["stderr"]
    assert "SKILL.md" in result["stderr"]
    assert "cmd" not in captured_run  # never started


def test_a_script_outside_the_skills_folder_is_refused(fake_catalog, captured_run, tmp_path):
    fake_catalog["escape"] = _entry("escape", "skills/../clawbio.py", ["--demo", "--input", "--output"],
                                    "python skills/../clawbio.py --demo")
    result = cli.run_skill("escape", demo=True, output_dir=str(tmp_path))

    assert result["success"] is False
    assert "cmd" not in captured_run


def test_a_skill_with_no_script_still_points_at_its_skill_md(captured_run, tmp_path):
    result = cli.run_skill("bgpt-mcp", demo=True, output_dir=str(tmp_path))

    assert result["success"] is False
    assert "SKILL.md" in result["stderr"]
    assert "cmd" not in captured_run


def test_the_audit_record_names_the_catalog_version(fake_catalog):
    fake_catalog["v-skill"] = _entry("v-skill", "skills/v-skill/v.py", ["--demo", "--input", "--output"],
                                     version="4.5.6")
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


def test_list_says_a_catalog_skill_runs_through_clawbio_run(capsys):
    assert "standard flags only" in _list_line(capsys, "vcf-annotator")


def test_list_says_a_refused_skill_does_not(capsys):
    # turingdb-graph's script defines no --output.
    assert "not via clawbio run" in _list_line(capsys, "turingdb-graph")


def test_a_profile_input_is_refused_when_the_script_takes_no_input(fake_catalog, captured_run, tmp_path):
    # --profile resolves to the profile's stored input file, which is passed as --input.
    fake_catalog["no-in"] = _entry("no-in", "skills/no-in/no_in.py", ["--demo", "--output"],
                                   "python skills/no-in/no_in.py --demo")
    result = cli.run_skill("no-in", profile_path=str(tmp_path / "PT001.json"), output_dir=str(tmp_path / "out"))

    assert result["success"] is False
    assert "--input" in result["stderr"]
    assert "cmd" not in captured_run


def test_a_catalog_script_missing_from_this_install_says_so(fake_catalog, captured_run, tmp_path):
    # A pip install leaves out licence-excluded skills, but the shipped catalog still lists them.
    fake_catalog["gone"] = _entry("gone", "skills/gone/gone.py", ["--demo", "--input", "--output"],
                                  "python skills/gone/gone.py --demo")
    result = cli.run_skill("gone", demo=True, output_dir=str(tmp_path))

    assert result["success"] is False
    assert "repository checkout" in result["stderr"]
    assert "cmd" not in captured_run
