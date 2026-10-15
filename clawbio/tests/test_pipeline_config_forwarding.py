"""Regression tests for nf-core pipeline Nextflow-config forwarding via clawbio.cli.

Repository bug being guarded:

``clawbio.cli`` historically exposed two *equivalent* argparse options for the
Nextflow ``-c`` config file — ``-c/--config`` (stored in ``args.extra_config``)
and ``--nextflow-config`` (stored in ``args.nextflow_config``). Each pipeline
forwarded only ONE of those two destinations:

* ``scrnaseq-pipeline`` forwarded ``extra_config`` and therefore silently dropped
  any ``--nextflow-config`` the user supplied.
* ``rnaseq-pipeline`` and ``sarek-pipeline`` forwarded ``nextflow_config`` and
  therefore silently dropped any ``-c/--config`` the user supplied.

All three wrappers accept ``-c``/``--config``/``--nextflow-config`` as aliases of
a single repeatable option, so the launcher must accept and forward them the same
way. After the fix every spelling is normalised and forwarded to every pipeline
wrapper as ``--nextflow-config <path>``.
"""
from __future__ import annotations

import pytest

from clawbio import cli


PIPELINES = ["scrnaseq-pipeline", "rnaseq-pipeline", "sarek-pipeline"]
SPELLINGS = ["-c", "--config", "--nextflow-config"]


def _drive_main(monkeypatch, argv):
    """Run ``cli.main()`` with ``run_skill`` stubbed, returning its ``extra_args``."""
    captured: dict[str, object] = {}

    def fake_run_skill(**kwargs):
        captured.update(kwargs)
        return {
            "output_dir": kwargs.get("output_dir"),
            "success": True,
            "exit_code": 0,
            "stdout": "",
            "stderr": "",
            "duration_seconds": 0,
            "files": [],
        }

    monkeypatch.setattr(cli, "run_skill", fake_run_skill)
    monkeypatch.setattr(cli.sys, "argv", argv)
    with pytest.raises(SystemExit):
        cli.main()
    return list(captured.get("extra_args") or [])


@pytest.mark.parametrize("skill", PIPELINES)
@pytest.mark.parametrize("spelling", SPELLINGS)
def test_config_flag_forwarded_for_every_pipeline(monkeypatch, tmp_path, skill, spelling):
    cfg = tmp_path / "hpc.config"
    cfg.write_text("// test config\n", encoding="utf-8")
    out = tmp_path / "out"
    argv = ["clawbio", "run", skill, "--demo", spelling, str(cfg), "--output", str(out)]

    extra = _drive_main(monkeypatch, argv)

    assert "--nextflow-config" in extra, f"{skill} via {spelling} dropped the config: {extra!r}"
    idx = extra.index("--nextflow-config")
    assert extra[idx + 1] == str(cfg), f"{skill} via {spelling} forwarded wrong path: {extra!r}"


@pytest.mark.parametrize("skill", PIPELINES)
def test_repeated_and_mixed_config_spellings_all_forwarded(monkeypatch, tmp_path, skill):
    c1 = tmp_path / "a.config"
    c1.write_text("// a\n", encoding="utf-8")
    c2 = tmp_path / "b.config"
    c2.write_text("// b\n", encoding="utf-8")
    out = tmp_path / "out"
    argv = [
        "clawbio", "run", skill, "--demo",
        "-c", str(c1),
        "--nextflow-config", str(c2),
        "--output", str(out),
    ]

    extra = _drive_main(monkeypatch, argv)
    joined = " ".join(extra)

    assert f"--nextflow-config {c1}" in joined, f"{skill}: first config dropped: {joined!r}"
    assert f"--nextflow-config {c2}" in joined, f"{skill}: second config dropped: {joined!r}"


@pytest.mark.parametrize("skill", PIPELINES)
@pytest.mark.parametrize("work_dir", ["", "work O'Brien $literal", "s3://analysis-bucket/work"])
def test_work_directory_reaches_every_pipeline_launcher(monkeypatch, tmp_path, skill, work_dir):
    extra = _drive_main(monkeypatch, [
        "clawbio", "run", skill, "--demo", "--resume", "--check",
        "--work-dir", work_dir, "--output", str(tmp_path / "out"),
    ])

    assert extra.count("--work-dir") == 1, extra
    assert extra[extra.index("--work-dir") + 1] == work_dir
    assert "--resume" in extra and "--check" in extra


def test_work_directory_is_not_forwarded_to_unrelated_skills(monkeypatch, tmp_path):
    extra = _drive_main(monkeypatch, [
        "clawbio", "run", "pharmgx", "--demo", "--work-dir", "work",
        "--output", str(tmp_path / "out"),
    ])
    assert "--work-dir" not in extra


@pytest.mark.parametrize("skill", ["sarek-pipeline", "scrnaseq-pipeline"])
def test_trace_report_suffix_reaches_supported_pipelines(monkeypatch, tmp_path, skill):
    suffix = "validation-run"
    extra = _drive_main(monkeypatch, [
        "clawbio", "run", skill, "--demo", "--trace-report-suffix", suffix,
        "--output", str(tmp_path / "out"),
    ])
    assert extra.count("--trace-report-suffix") == 1, extra
    assert extra[extra.index("--trace-report-suffix") + 1] == suffix
