from __future__ import annotations

import json
from importlib import import_module
import os
from pathlib import Path

import pytest
from hypothesis import HealthCheck, example, given, settings, strategies as st

# The runner engine lives in clawbio.cli; import it directly so monkeypatching
# its subprocess attribute affects run_skill's behaviour.
from clawbio import cli as clawbio_runner

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def test_run_skill_passes_absolute_output_to_subprocess(monkeypatch, tmp_path):
    captured: dict[str, object] = {}

    class Proc:
        returncode = 0
        stdout = ""
        stderr = ""

    def fake_run(cmd, capture_output, text, timeout, cwd):
        captured["cmd"] = cmd
        captured["cwd"] = cwd
        return Proc()

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(clawbio_runner.subprocess, "run", fake_run)

    result = clawbio_runner.run_skill(
        skill_name="gwas",
        output_dir="relative_gwas_output",
        extra_args=["--rsid", "rs861539", "--no-cache"],
    )

    expected_output = tmp_path / "relative_gwas_output"
    assert result["success"] is True
    assert result["output_dir"] == str(expected_output)
    cmd = captured["cmd"]
    assert cmd[cmd.index("--output") + 1] == str(expected_output)
    assert captured["cwd"] == str(PROJECT_ROOT / "skills" / "gwas-lookup")


def test_run_skill_passes_absolute_relative_input_to_subprocess(monkeypatch, tmp_path):
    captured: dict[str, object] = {}

    class Proc:
        returncode = 0
        stdout = ""
        stderr = ""

    def fake_run(cmd, capture_output=True, text=True, timeout=None, cwd=None):
        captured["cmd"] = cmd
        captured["cwd"] = cwd
        return Proc()

    samplesheet = tmp_path / "samplesheet.csv"
    samplesheet.write_text("patient,sample,lane,fastq_1,fastq_2\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    launcher = import_module("clawbio.common.nfcore.launcher")

    monkeypatch.setattr(clawbio_runner.subprocess, "run", fake_run)
    monkeypatch.setattr(launcher, "run_wrapper", fake_run)
    monkeypatch.setattr(
        clawbio_runner.subprocess, "Popen",
        lambda *args, **kwargs: pytest.fail("A routing test must not start a real skill"),
    )

    result = clawbio_runner.run_skill(
        skill_name="sarek-pipeline",
        input_path="samplesheet.csv",
        output_dir=str(tmp_path / "sarek_out"),
        extra_args=["--check", "--no-banner"],
    )

    assert result["success"] is True
    cmd = captured["cmd"]
    assert cmd[cmd.index("--input") + 1] == str(samplesheet)
    assert captured["cwd"] == str(PROJECT_ROOT / "skills" / "nfcore-sarek-wrapper")


def _run_gwas(monkeypatch, tmp_path, returncode):
    class Proc:
        stdout = ""
        stderr = "boom"

    Proc.returncode = returncode
    monkeypatch.setattr(clawbio_runner.subprocess, "run", lambda *a, **k: Proc())
    clawbio_runner.run_skill(skill_name="gwas", output_dir=str(tmp_path / "out"), extra_args=["--rsid", "rs861539"])
    log = Path(os.environ["CLAWBIO_AUDIT_LOG"])
    return [json.loads(line) for line in log.read_text().splitlines()]


def test_run_skill_writes_an_audit_record(monkeypatch, tmp_path):
    (record,) = _run_gwas(monkeypatch, tmp_path, 0)
    assert record["event"] == "skill_run"
    assert record["gen_ai.agent.id"] == "gwas"
    assert record["status"] == "OK"
    assert record["output.value"] == str(tmp_path / "out")


def test_run_skill_audit_record_marks_a_failed_run(monkeypatch, tmp_path):
    (record,) = _run_gwas(monkeypatch, tmp_path, 2)
    assert record["status"] == "ERROR"
    assert "2" in record["error"]


def test_run_skill_audit_record_carries_the_skill_md_version(monkeypatch, tmp_path):
    (record,) = _run_gwas(monkeypatch, tmp_path, 0)
    frontmatter = (PROJECT_ROOT / "skills" / "gwas-lookup" / "SKILL.md").read_text().split("---")[1]
    version = next(l.split(":", 1)[1].strip() for l in frontmatter.splitlines() if l.strip().startswith("version:"))
    assert record["gen_ai.agent.version"] == version


_NFCORE_ALIASES = ("rnaseq-pipeline", "sarek-pipeline", "scrnaseq-pipeline")


def _runner_timeout_for(monkeypatch, tmp_path, *argv):
    """The timeout `clawbio run` gives run_skill for the given command line."""
    import sys

    seen: dict[str, object] = {}

    def fake_run_skill(**kwargs):
        seen["timeout"] = kwargs["timeout"]
        return {"skill": kwargs["skill_name"], "success": True, "exit_code": 0,
                "output_dir": str(tmp_path), "files": [], "stdout": "", "stderr": "",
                "duration_seconds": 0}

    monkeypatch.setattr(clawbio_runner, "run_skill", fake_run_skill)
    monkeypatch.setattr(sys, "argv", ["clawbio", "run", *argv, "--output", str(tmp_path)])
    with pytest.raises(SystemExit):
        clawbio_runner.main()
    return seen["timeout"]



@pytest.mark.parametrize("alias", _NFCORE_ALIASES)
def test_runner_outlasts_the_wrappers_own_timeout(alias, monkeypatch, tmp_path):
    """Let the wrapper reach its own deadline before the fallback runner cap.

    The historical fixed 12h10m cap stopped --timeout-hours 48 or 0 too early;
    the runner deadline must still follow the requested wrapper deadline.
    """
    default = clawbio_runner.SKILLS[alias]["default_timeout_seconds"]
    assert _runner_timeout_for(monkeypatch, tmp_path, alias, "--demo") == default
    assert _runner_timeout_for(monkeypatch, tmp_path, alias, "--demo", "--timeout-hours", "48") == 48 * 3600 + 600
    assert _runner_timeout_for(monkeypatch, tmp_path, alias, "--demo", "--timeout-hours=0.5") == 1800 + 600
    assert _runner_timeout_for(monkeypatch, tmp_path, alias, "--demo", "--timeout-hours", "0") is None
    # An explicit runner --timeout still wins; a bad value is the wrapper's usage error.
    assert _runner_timeout_for(monkeypatch, tmp_path, alias, "--demo", "--timeout-hours", "48", "--timeout", "999") == 999
    assert _runner_timeout_for(monkeypatch, tmp_path, alias, "--demo", "--timeout-hours", "nan") == default



@pytest.mark.parametrize("alias", _NFCORE_ALIASES)
@settings(deadline=None, max_examples=10,
          suppress_health_check=[HealthCheck.function_scoped_fixture])
@example(seconds=300)
@given(seconds=st.integers(min_value=0, max_value=100000))
def test_explicit_runner_timeout_always_wins(alias, monkeypatch, tmp_path, seconds):
    assert _runner_timeout_for(
        monkeypatch, tmp_path, alias, "--demo", "--timeout-hours", "48",
        "--timeout", str(seconds),
    ) == seconds



@pytest.mark.parametrize("alias", _NFCORE_ALIASES)
@settings(deadline=None, max_examples=16,
          suppress_health_check=[HealthCheck.function_scoped_fixture])
@example(hours="0", flag="--timeout_hours", joined=False)
@example(hours="48", flag="--timeout_hours", joined=True)
@given(hours=st.sampled_from(("0", "0.5", "12", "48")),
       flag=st.sampled_from(("--timeout-hours", "--timeout_hours")),
       joined=st.booleans())
def test_wrapper_timeout_spellings_share_runner_deadline(
    alias, monkeypatch, tmp_path, hours, flag, joined,
):
    expected = {"0": None, "0.5": 2400, "12": 43800, "48": 173400}
    arguments = [f"{flag}={hours}"] if joined else [flag, hours]
    assert _runner_timeout_for(
        monkeypatch, tmp_path, alias, "--demo", *arguments,
    ) == expected[hours]


def _run_main_with(monkeypatch, tmp_path, capsys, *, success, stderr):
    import sys

    def fake_run_skill(**kwargs):
        return {"skill": kwargs["skill_name"], "success": success, "exit_code": 0 if success else 1,
                "output_dir": str(tmp_path), "files": [], "stdout": "", "stderr": stderr,
                "duration_seconds": 0}

    monkeypatch.setattr(clawbio_runner, "run_skill", fake_run_skill)
    monkeypatch.setattr(sys, "argv", ["clawbio", "run", "rnaseq-pipeline", "--demo", "--output", str(tmp_path)])
    with pytest.raises(SystemExit):
        clawbio_runner.main()
    return capsys.readouterr().out


def test_successful_run_relays_the_skills_warnings(monkeypatch, tmp_path, capsys):
    """A skill's stderr was shown only on failure, so `clawbio run rnaseq-pipeline
    --demo --skip-qc` succeeded without saying --skip-qc had no effect."""
    stderr = (
        "progress line\n"
        "WARNING: --demo ignores tuning flags (--skip-qc); the upstream test profile ...\n"
        "WARNING: --demo forces --aligner star_salmon (requested: 'hisat2').\n"
        "WARNING: --demo forces --aligner star_salmon (requested: 'hisat2').\n"
    )
    out = _run_main_with(monkeypatch, tmp_path, capsys, success=True, stderr=stderr)
    assert out.count("WARNING: --demo ignores tuning flags (--skip-qc)") == 1
    assert out.count("WARNING: --demo forces --aligner star_salmon") == 1
    assert "progress line" not in out


def test_successful_run_without_warnings_prints_none(monkeypatch, tmp_path, capsys):
    out = _run_main_with(monkeypatch, tmp_path, capsys, success=True, stderr="progress line\n")
    assert "WARNING" not in out
