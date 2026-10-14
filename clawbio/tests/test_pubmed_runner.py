"""Runner forwards the PubMed summary options without allowing arbitrary flags."""
from unittest.mock import Mock
import io
import sys

import pytest

from clawbio import cli


def test_pubmed_runner_forwards_summary_options(tmp_path, monkeypatch):
    process = Mock(return_value=Mock(returncode=0, stdout="", stderr=""))
    monkeypatch.setattr(cli.subprocess, "run", process)
    options = ["--query", "BRCA1", "--summary-method", "llm", "--provider", "ollama",
               "--model", "test-model", "--model-params", '{"reasoning_effort":"none"}',
               "--llm-timeout", "120", "--summary-max-tokens", "200",
               "--base-url", "http://localhost:11434/v1", "--max-results", "2"]
    result = cli.run_skill("pubmed-summariser", output_dir=str(tmp_path),
                          extra_args=[*options, "--arbitrary", "blocked"])
    assert result["success"]
    command = process.call_args.args[0]
    assert command[-len(options):] == options
    assert "--arbitrary" not in command


def test_pubmed_report_preview_with_legacy_windows_output(tmp_path, monkeypatch):
    (tmp_path / "report.md").write_text("# Study of α\n", encoding="utf-8")
    monkeypatch.setattr(cli, "run_skill", Mock(return_value={
        "success": True, "exit_code": 0, "duration_seconds": 1,
        "output_dir": str(tmp_path), "files": ["report.md"], "stdout": "", "stderr": "",
    }))
    monkeypatch.setattr(sys, "argv", ["clawbio", "run", "pubmed-summariser", "--demo"])
    buffer = io.BytesIO()
    stream = io.TextIOWrapper(buffer, encoding="cp1252", errors="strict")
    monkeypatch.setattr(sys, "stdout", stream)
    try:
        with pytest.raises(SystemExit) as exc:
            cli.main()
        assert exc.value.code == 0
        stream.flush()
        assert b"\\u03b1" in buffer.getvalue()
    finally:
        stream.close()
