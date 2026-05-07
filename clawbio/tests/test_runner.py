from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
_RUNNER_SPEC = importlib.util.spec_from_file_location("clawbio_runner", PROJECT_ROOT / "clawbio.py")
clawbio_runner = importlib.util.module_from_spec(_RUNNER_SPEC)
sys.modules["clawbio_runner"] = clawbio_runner
assert _RUNNER_SPEC.loader is not None
_RUNNER_SPEC.loader.exec_module(clawbio_runner)


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
