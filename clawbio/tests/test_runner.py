from __future__ import annotations

import json
import os
from pathlib import Path

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

    def fake_run(cmd, capture_output, text, timeout, cwd):
        captured["cmd"] = cmd
        captured["cwd"] = cwd
        return Proc()

    samplesheet = tmp_path / "samplesheet.csv"
    samplesheet.write_text("patient,sample,lane,fastq_1,fastq_2\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(clawbio_runner.subprocess, "run", fake_run)

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
