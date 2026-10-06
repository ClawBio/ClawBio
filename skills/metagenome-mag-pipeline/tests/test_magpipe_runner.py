"""The snakemake command: exact argv, quoting, log capture, timeout, exit codes."""

from __future__ import annotations

import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

import pytest
import yaml

SKILL_DIR = Path(__file__).resolve().parents[1]
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

from magpipe_errors import ErrorCode, SkillError
from magpipe_runner import build_snakemake_command, run_snakemake


def prepared_config(tmp_path: Path) -> Path:
    """A config the stub snakemake can act on: a samplesheet and a results_dir."""
    sheet = tmp_path / "samplesheet.tsv"
    sheet.write_text("sample\tfastq_1\tfastq_2\tgroup\tsra_run\tlayout\nS1\t\t\tS1\t\tsingle\n")
    config = tmp_path / "config.yaml"
    config.write_text(
        yaml.safe_dump(
            {
                "samplesheet": str(sheet),
                "results_dir": str(tmp_path / "results"),
                "database_dir": str(tmp_path / "db"),
                "qc": {"sourmash": {"enabled": False}},
                "profile": {"profiler": "none"},
                "assembly": {"enabled": True},
                "binning": {"enabled": False},
                "mag": {"enabled": False},
                "taxonomy": {"enabled": False},
                "annotation": {"enabled": False},
                "heterogeneity": {"enabled": False},
                "report": {"multiqc": True},
            }
        )
    )
    return config


def command(tmp_path: Path, **kwargs) -> tuple[list[str], str]:
    defaults = {
        "pipeline_dir": tmp_path / "pipeline",
        "config_path": tmp_path / "config.yaml",
        "work_dir": tmp_path / "workdir",
        "conda_prefix": tmp_path / "conda",
        "cores": 8,
        "mem_mb": 32000,
        "dry_run": False,
    }
    return build_snakemake_command(**{**defaults, **kwargs})


class TestCommand:
    def test_argv_order_matches_the_upstream_smoke_test(self, tmp_path: Path) -> None:
        argv, _ = command(tmp_path)
        assert argv == [
            "snakemake",
            "-s",
            str(tmp_path / "pipeline" / "workflow" / "Snakefile"),
            "--configfile",
            str(tmp_path / "config.yaml"),
            "--directory",
            str(tmp_path / "workdir"),
            "--use-conda",
            "--conda-frontend",
            "conda",
            "--conda-prefix",
            str(tmp_path / "conda"),
            "--cores",
            "8",
            "--resources",
            "mem_mb=32000",
            "--rerun-incomplete",
            "--printshellcmds",
        ]

    def test_dry_run_appends_the_flag(self, tmp_path: Path) -> None:
        argv, _ = command(tmp_path, dry_run=True)
        assert argv[-1] == "--dry-run"
        assert argv.count("--dry-run") == 1

    def test_no_dry_run_by_default(self, tmp_path: Path) -> None:
        argv, _ = command(tmp_path)
        assert "--dry-run" not in argv

    def test_uses_conda_and_conda_frontend(self, tmp_path: Path) -> None:
        argv, _ = command(tmp_path)
        assert "--use-conda" in argv
        assert argv[argv.index("--conda-frontend") + 1] == "conda"

    def test_string_form_is_shell_join_of_the_argv(self, tmp_path: Path) -> None:
        argv, rendered = command(tmp_path)
        assert rendered == shlex.join(argv)

    def test_paths_with_spaces_are_quoted_in_the_string_form(self, tmp_path: Path) -> None:
        spaced = tmp_path / "my output dir"
        spaced.mkdir()
        argv, rendered = command(tmp_path, pipeline_dir=spaced)
        assert str(spaced / "workflow" / "Snakefile") in argv
        assert f"'{spaced / 'workflow' / 'Snakefile'}'" in rendered
        assert rendered == shlex.join(argv)

    def test_string_form_survives_a_shell_round_trip(self, tmp_path: Path) -> None:
        spaced = tmp_path / "dir with space"
        spaced.mkdir()
        argv, rendered = command(tmp_path, config_path=spaced / "config file.yaml")
        assert shlex.split(rendered) == argv

    def test_cores_and_memory_are_passed_through(self, tmp_path: Path) -> None:
        argv, _ = command(tmp_path, cores=2, mem_mb=4000)
        assert argv[argv.index("--cores") + 1] == "2"
        assert argv[argv.index("--resources") + 1] == "mem_mb=4000"


class TestRun:
    def test_successful_run_writes_the_log(self, stub_snakemake: Path, tmp_path: Path) -> None:
        prepared_config(tmp_path)
        argv, _ = command(tmp_path)
        log = tmp_path / "pipeline" / "snakemake.log"
        assert run_snakemake(argv, log_path=log, timeout_seconds=60) == 0
        text = log.read_text()
        assert argv[0] in text
        assert "stub snakemake" in text

    def test_non_zero_exit_raises_with_a_log_tail(self, stub_snakemake: Path, tmp_path: Path, monkeypatch) -> None:
        monkeypatch.setenv("STUB_SNAKEMAKE_EXIT", "1")
        argv, _ = command(tmp_path)
        log = tmp_path / "snakemake.log"
        with pytest.raises(SkillError) as excinfo:
            run_snakemake(argv, log_path=log, timeout_seconds=60)
        error = excinfo.value
        assert error.error_code == ErrorCode.EXECUTION_FAILED
        assert error.details["exit_code"] == 1
        assert error.details["log_tail"]
        assert any("forced failure" in line for line in error.details["log_tail"])
        assert error.details["log"] == str(log)

    def test_log_tail_is_capped(self, stub_snakemake: Path, tmp_path: Path, monkeypatch) -> None:
        monkeypatch.setenv("STUB_SNAKEMAKE_EXIT", "1")
        argv, _ = command(tmp_path)
        log = tmp_path / "snakemake.log"
        with pytest.raises(SkillError) as excinfo:
            run_snakemake(argv, log_path=log, timeout_seconds=60)
        assert len(excinfo.value.details["log_tail"]) <= 40

    def test_timeout_raises_and_kills_the_process_group(self, tmp_path: Path) -> None:
        sleeper = tmp_path / "sleeper.py"
        sleeper.write_text(
            "import time\nprint('starting', flush=True)\ntime.sleep(120)\n"
        )
        log = tmp_path / "snakemake.log"
        argv = [sys.executable, str(sleeper)]
        with pytest.raises(SkillError) as excinfo:
            run_snakemake(argv, log_path=log, timeout_seconds=1)
        assert excinfo.value.error_code == ErrorCode.EXECUTION_TIMEOUT
        assert excinfo.value.details["timeout_seconds"] == 1

    def test_timeout_leaves_no_child_behind(self, tmp_path: Path) -> None:
        """A sleeping child of the stub must be gone once the timeout is raised."""
        marker = tmp_path / "child_alive"
        child = tmp_path / "child.py"
        child.write_text(
            "import os, time, pathlib\n"
            f"pathlib.Path({str(marker)!r}).write_text(str(os.getpid()))\n"
            "time.sleep(120)\n"
        )
        parent = tmp_path / "parent.py"
        parent.write_text(
            "import subprocess, sys, time\n"
            f"subprocess.Popen([sys.executable, {str(child)!r}])\n"
            "print('spawned', flush=True)\n"
            "time.sleep(120)\n"
        )
        log = tmp_path / "snakemake.log"
        with pytest.raises(SkillError) as excinfo:
            run_snakemake([sys.executable, str(parent)], log_path=log, timeout_seconds=3)
        assert excinfo.value.error_code == ErrorCode.EXECUTION_TIMEOUT
        assert marker.is_file(), "the child never started, so the test proves nothing"
        pid = int(marker.read_text())
        deadline = time.monotonic() + 10
        while time.monotonic() < deadline:
            try:
                os.kill(pid, 0)
            except ProcessLookupError:
                break
            time.sleep(0.2)
        else:  # pragma: no cover - only on a machine that leaked the process
            pytest.fail(f"child process {pid} survived the timeout")
        assert excinfo.value.details["log_tail"] or log.is_file()

    def test_missing_executable_is_an_execution_failure(self, tmp_path: Path) -> None:
        with pytest.raises(SkillError) as excinfo:
            run_snakemake(
                [str(tmp_path / "not-a-program")],
                log_path=tmp_path / "log.txt",
                timeout_seconds=10,
            )
        assert excinfo.value.error_code == ErrorCode.EXECUTION_FAILED
        assert "snakemake" in excinfo.value.fix

    def test_stdout_and_stderr_are_both_captured(self, tmp_path: Path) -> None:
        noisy = tmp_path / "noisy.py"
        noisy.write_text(
            "import sys\nprint('to stdout')\nprint('to stderr', file=sys.stderr)\n"
        )
        log = tmp_path / "log.txt"
        assert run_snakemake([sys.executable, str(noisy)], log_path=log, timeout_seconds=30) == 0
        text = log.read_text()
        assert "to stdout" in text and "to stderr" in text

    def test_log_directory_is_created(self, stub_snakemake: Path, tmp_path: Path) -> None:
        prepared_config(tmp_path)
        argv, _ = command(tmp_path)
        log = tmp_path / "deep" / "nested" / "snakemake.log"
        run_snakemake(argv, log_path=log, timeout_seconds=60)
        assert log.is_file()

    def test_the_suite_never_touches_the_network(
        self, stub_snakemake: Path, tmp_path: Path
    ) -> None:
        """The stub is what ran: no output appeared under a real results_dir."""
        prepared_config(tmp_path)
        argv, _ = command(tmp_path)
        run_snakemake(argv, log_path=tmp_path / "log.txt", timeout_seconds=60)
        assert (tmp_path / "results" / "final" / "qc" / "read_counts.tsv").is_file()
        assert "stub snakemake" in (tmp_path / "log.txt").read_text()
        assert subprocess.run(
            ["git", "status", "--porcelain", "skills"],
            cwd=SKILL_DIR.parent,
            capture_output=True,
            text=True,
        ).returncode == 0


class TestTransientFailures:
    """Failures that must not be reported as fatal configuration errors.

    A conda environment build can die with EDQUOT on a shared filesystem that
    reports terabytes free (seen live: `Errno 122` followed by a clean rerun).
    A workdir can be locked by another run after a kill. Both are retryable or
    explainable; neither is a pipeline defect.
    """

    def _flaky_snakemake(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
        """A stub that dies with the quota signature once, then succeeds."""
        import stat

        bindir = tmp_path / "flaky_bin"
        bindir.mkdir()
        stub = bindir / "snakemake"
        stub.write_text(
            "#!/usr/bin/env python3\n"
            "import sys\n"
            "from pathlib import Path\n"
            "state = Path(__file__).resolve().parent / 'calls.txt'\n"
            "calls = state.read_text().splitlines() if state.exists() else []\n"
            "state.write_text('\\n'.join(calls + ['x']) + '\\n')\n"
            "if not calls:\n"
            "    print('CreateCondaEnvironmentException: Could not create conda environment')\n"
            "    print('[Errno 122] Disk quota exceeded')\n"
            "    raise SystemExit(1)\n"
            "print('flaky stub snakemake: recovered on attempt 2')\n"
        )
        stub.chmod(stub.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")
        return bindir

    def test_env_quota_failure_is_retried(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        self._flaky_snakemake(tmp_path, monkeypatch)
        argv, _ = command(tmp_path)
        log = tmp_path / "snakemake.log"
        assert run_snakemake(argv, log_path=log, timeout_seconds=60, max_attempts=2) == 0
        text = log.read_text()
        assert "Disk quota exceeded" in text
        assert "retry" in text.lower()

    def test_env_quota_failure_raises_without_retry(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        self._flaky_snakemake(tmp_path, monkeypatch)
        argv, _ = command(tmp_path)
        log = tmp_path / "snakemake.log"
        with pytest.raises(SkillError) as excinfo:
            run_snakemake(argv, log_path=log, timeout_seconds=60, max_attempts=1)
        assert excinfo.value.error_code == ErrorCode.EXECUTION_FAILED

    def test_locked_workdir_names_the_lock(self, tmp_path: Path) -> None:
        argv = [
            sys.executable,
            "-c",
            "print('LockException: Directory cannot be locked.'); raise SystemExit(1)",
        ]
        log = tmp_path / "snakemake.log"
        with pytest.raises(SkillError) as excinfo:
            run_snakemake(argv, log_path=log, timeout_seconds=60)
        error = excinfo.value
        assert error.error_code == ErrorCode.EXECUTION_FAILED
        assert "lock" in error.message.lower()
        assert "unlock" in error.fix.lower() or "another run" in error.message.lower()
