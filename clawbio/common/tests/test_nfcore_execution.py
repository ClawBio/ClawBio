"""Real subprocess execution and injectable failure-policy contracts."""

import importlib
import sys
from pathlib import Path

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st


def _run(api, command, output, **overrides):
    return api.execute_nextflow(
        command,
        cwd=output,
        output_dir=output,
        timeout_seconds=overrides.pop("timeout_seconds", 10),
        timeout_fix="site-specific timeout advice",
        **overrides,
    )


def test_execution_captures_real_stdout_stderr(tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.execution")
    result = _run(
        api,
        [
            sys.executable,
            "-c",
            'import sys; print("hello"); print("oops",file=sys.stderr)',
        ],
        tmp_path,
        timeout_seconds=None,
    )
    assert result == {
        "exit_code": 0,
        "stdout_path": str(tmp_path / "logs/stdout.txt"),
        "stderr_path": str(tmp_path / "logs/stderr.txt"),
    }
    assert (tmp_path / "logs/stdout.txt").read_text() == "hello\n"
    assert (tmp_path / "logs/stderr.txt").read_text() == "oops\n"


def test_launch_failure_removes_empty_logs(tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.execution")
    errors = importlib.import_module("clawbio.common.nfcore.errors")
    with pytest.raises(errors.SkillError) as failure:
        _run(api, [str(tmp_path / "missing-executable")], tmp_path)
    assert failure.value.stage == "preflight"
    assert failure.value.error_code == "EXECUTION_FAILED"
    assert not list((tmp_path / "logs").iterdir())


def test_nonzero_execution_reads_nextflow_log_hint(tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.execution")
    errors = importlib.import_module("clawbio.common.nfcore.errors")
    (tmp_path / ".nextflow.log").write_text("java.net.UnknownHostException")
    with pytest.raises(errors.SkillError) as failure:
        _run(api, [sys.executable, "-c", "raise SystemExit(7)"], tmp_path)
    assert failure.value.details["exit_code"] == 7
    assert "IPv6" in failure.value.fix


def test_timeout_uses_supplied_advice_and_cleans_up(tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.execution")
    errors = importlib.import_module("clawbio.common.nfcore.errors")
    with pytest.raises(errors.SkillError) as failure:
        _run(
            api,
            [sys.executable, "-c", "import time; time.sleep(60)"],
            tmp_path,
            timeout_seconds=0.02,
        )
    assert failure.value.fix == "site-specific timeout advice"
    assert failure.value.details == {"timeout_seconds": 0.02}


def test_pipeline_hint_follows_the_shared_environment_hints(tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.execution")
    errors = importlib.import_module("clawbio.common.nfcore.errors")
    (tmp_path / ".nextflow.log").write_text("Network is unreachable")
    with pytest.raises(errors.SkillError) as failure:
        _run(
            api,
            [sys.executable, "-c", "raise SystemExit(1)"],
            tmp_path,
            pipeline_failure_hint=lambda stdout, stderr: " PIPELINE-HINT",
        )
    fix = failure.value.fix
    assert fix.startswith("Inspect logs/stdout.txt")
    assert fix.index("preferIPv6Addresses") < fix.index(" PIPELINE-HINT")
    assert fix.endswith(" PIPELINE-HINT")


@pytest.mark.parametrize(
    ("platform", "under_tmp", "expected"),
    [("darwin", True, True), ("darwin", False, False), ("linux", True, False)],
)
def test_macos_tmp_hint_only_for_macos_output_under_tmp(
    platform, under_tmp, expected, monkeypatch
):
    api = importlib.import_module("clawbio.common.nfcore.execution")
    monkeypatch.setattr(sys, "platform", platform)
    monkeypatch.setattr(api, "is_under_tmp", lambda _path: under_tmp)
    hint = api.macos_tmp_failure_hint(Path("/tmp/run"))
    assert bool(hint) is expected
    if expected:
        assert "/tmp" in hint and "HOME" in hint


def test_pipeline_environment_overrides_the_callers(tmp_path, monkeypatch):
    """A pipeline that needs a specific Nextflow setting (e.g. the v1 syntax parser
    for releases predating strict syntax) must get it even if the user's shell
    exports the opposite; everything else is inherited."""
    api = importlib.import_module("clawbio.common.nfcore.execution")
    monkeypatch.setenv("NXF_SYNTAX_PARSER", "v2")
    monkeypatch.setenv("CLAWBIO_INHERITED", "yes")
    _run(
        api,
        ["sh", "-c", 'printf "%s %s" "$NXF_SYNTAX_PARSER" "$CLAWBIO_INHERITED"'],
        tmp_path,
        env={"NXF_SYNTAX_PARSER": "v1"},
    )
    assert (tmp_path / "logs/stdout.txt").read_text() == "v1 yes"


def test_log_tail_reads_only_the_tail_of_a_large_log(tmp_path, monkeypatch):
    """Failure hints scan the end of Nextflow's logs, which can run to gigabytes on
    a long run; reading them whole just to keep the last 64 KiB wastes memory."""
    api = importlib.import_module("clawbio.common.nfcore.execution")
    log = tmp_path / "stdout.txt"
    log.write_bytes(b"x" * 1_000_000 + "é tail marker".encode())
    read_sizes = []
    real_open = Path.open

    def spy_open(self, *args, **kwargs):
        handle = real_open(self, *args, **kwargs)
        real_read = handle.read

        def read(size=-1):
            read_sizes.append(size)
            return real_read(size)

        handle.read = read
        return handle

    monkeypatch.setattr(Path, "open", spy_open)
    tail = api.read_log_tail(log, limit=1024)
    assert tail.endswith("é tail marker")
    assert len(tail.encode()) <= 1024
    assert read_sizes and all(0 <= size <= 1024 for size in read_sizes)


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX signals and process groups")
def test_sigterm_stops_nextflow_before_the_wrapper_exits(tmp_path):
    """`kill <wrapper>` (or an agent or scheduler stopping it) sends SIGTERM, which
    by default ends Python at once and leaves Nextflow running in its own session.
    The runner must stop Nextflow first and then exit 128+SIGTERM."""
    import os
    import signal
    import subprocess
    import textwrap
    import time

    pid_file = tmp_path / "child.pid"
    driver = textwrap.dedent(
        f"""
        from pathlib import Path
        from clawbio.common.nfcore.execution import execute_nextflow
        execute_nextflow(
            ["sh", "-c", "echo $$ > {pid_file}; exec sleep 60"],
            cwd=Path({str(tmp_path)!r}), output_dir=Path({str(tmp_path)!r}),
            timeout_seconds=None, timeout_fix="",
        )
        """
    )
    root = Path(__file__).resolve().parents[3]
    wrapper = subprocess.Popen(
        [sys.executable, "-c", driver], env={**os.environ, "PYTHONPATH": str(root)}
    )
    deadline = time.monotonic() + 30
    while not pid_file.exists() and time.monotonic() < deadline:
        time.sleep(0.05)
    assert pid_file.exists(), "the stand-in for Nextflow never started"
    time.sleep(0.2)
    wrapper.send_signal(signal.SIGTERM)
    assert wrapper.wait(timeout=30) == 128 + signal.SIGTERM
    with pytest.raises(ProcessLookupError):
        os.kill(int(pid_file.read_text()), 0)


# nf-schema's refusal, as Nextflow 26.04.3 printed it for nf-core/sarek 3.8.1.
_NF_SCHEMA_REFUSAL = """ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --clip_r1 (-1): -1 is less than 0
* --wes (false): Value is [string] but should be [boolean]

 -- Check script 'subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68
"""


def test_nf_schema_refusal_names_every_invalid_value(tmp_path):
    """nf-schema is the only parameter validator, so its refusal must reach the
    user as clearly as the wrapper's own errors do."""
    api = importlib.import_module("clawbio.common.nfcore.execution")
    errors = importlib.import_module("clawbio.common.nfcore.errors")
    (tmp_path / ".nextflow.log").write_text(_NF_SCHEMA_REFUSAL)
    with pytest.raises(errors.SkillError) as failure:
        _run(api, [sys.executable, "-c", "raise SystemExit(1)"], tmp_path)
    assert "nf-schema" in failure.value.fix
    assert "--clip_r1 (-1): -1 is less than 0" in failure.value.fix
    assert "--wes (false): Value is [string] but should be [boolean]" in failure.value.fix


@settings(suppress_health_check=[HealthCheck.function_scoped_fixture], max_examples=25)
@given(st.lists(
    st.tuples(st.from_regex(r"[a-z][a-z_]{0,12}", fullmatch=True), st.text(alphabet="abc019 -.[]", min_size=1, max_size=20)),
    min_size=1, max_size=4, unique_by=lambda item: item[0],
))
def test_each_refused_value_is_listed_once_across_logs(tmp_path, refused):
    api = importlib.import_module("clawbio.common.nfcore.execution")
    listed = [f"--{name} ({value}): {value} is invalid" for name, value in refused]
    text = "The following invalid input values have been detected:\n\n" + "".join(f"* {line}\n" for line in listed)
    logs = [tmp_path / "stdout.txt", tmp_path / "stderr.txt", tmp_path / ".nextflow.log"]
    for log in logs:  # Nextflow repeats the refusal in its console output and its log
        log.write_text(text + "\nunrelated line\n* --not_a_refusal (x): listed elsewhere\n")
    hint = api.environment_failure_hints(*logs)
    assert [hint.count(line) for line in listed] == [1] * len(listed)
    assert "not_a_refusal" not in hint


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX signals and process groups")
def test_tasks_are_stopped_even_after_nextflow_itself_has_exited(tmp_path):
    """The group was looked up with getpgid(Nextflow's pid), which fails once
    Nextflow has exited and been reaped, so the tasks it left in its session kept
    running. The session leader's pid is the group id."""
    import contextlib
    import os
    import subprocess
    import time

    api = importlib.import_module("clawbio.common.nfcore.execution")
    pid_file = tmp_path / "task.pid"
    leader = subprocess.Popen(
        ["sh", "-c", f"sleep 60 & echo $! > {pid_file}; exit 0"], start_new_session=True
    )
    leader.wait()
    task = int(pid_file.read_text())
    try:
        api.terminate_process_tree(leader, grace_seconds=2)
        deadline = time.monotonic() + 10
        while time.monotonic() < deadline:
            try:
                os.kill(task, 0)
            except ProcessLookupError:
                break
            time.sleep(0.1)
        else:
            pytest.fail("the task process outlived terminate_process_tree")
    finally:
        # Never leave the stand-in task behind, even when the assertion failed.
        with contextlib.suppress(ProcessLookupError):
            os.kill(task, 9)
