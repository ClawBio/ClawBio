"""Tests for the headless Claude Code trace parser used by the agent evals.

These run in normal CI: they need no model, only recorded stream-json lines.
"""
import importlib.util
import json
import sys
from pathlib import Path

import pytest
from hypothesis import example, given, strategies as st


def _load():
    path = Path(__file__).resolve().parent / "agent_trace.py"
    spec = importlib.util.spec_from_file_location("clawbio_agent_trace", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module  # @dataclass reads sys.modules[cls.__module__]
    spec.loader.exec_module(module)
    return module


MOD = _load()


def _assistant(*blocks):
    return json.dumps({"type": "assistant", "message": {"content": list(blocks)}})


def _tool(name, **inp):
    return {"type": "tool_use", "id": f"t_{name}", "name": name, "input": inp}


RESULT = json.dumps({"type": "result", "subtype": "success", "result": "Done.",
                     "total_cost_usd": 0.12, "num_turns": 3})


def test_parses_tool_calls_and_final_answer():
    lines = [
        json.dumps({"type": "system", "subtype": "init"}),
        _assistant({"type": "text", "text": "Running it."}, _tool("Bash", command="python clawbio.py run pharmgx --demo")),
        json.dumps({"type": "system", "subtype": "permission_denied", "tool_name": "Bash"}),
        RESULT,
    ]
    trace = MOD.parse_stream(lines)
    assert [c.command for c in trace.calls] == ["python clawbio.py run pharmgx --demo"]
    assert trace.final == "Done."


def _bash(tool_id, command):
    return _assistant({"type": "tool_use", "id": tool_id, "name": "Bash", "input": {"command": command}})


def _denied(tool_id):
    # Shape Claude Code 2.1.302 emits in dontAsk mode.
    return json.dumps({"type": "system", "subtype": "permission_denied", "tool_name": "Bash",
                       "tool_use_id": tool_id, "decision_reason_type": "mode"})


def test_a_denied_bash_call_is_attempted_but_not_run():
    """A denied skill call must not count as the agent running the skill."""
    lines = [
        _bash("t1", "python clawbio.py run pharmgx --demo --output /tmp/x 2>&1 | tail; ls /tmp/x"),
        _denied("t1"),
        _bash("t2", "python clawbio.py run pharmgx --demo --output output/x"),
        _result("t2", "done"),
        RESULT,
    ]
    trace = MOD.parse_stream(lines)
    assert [c.outcome for c in trace.calls] == ["denied", "ok"]
    assert not trace.ran_ok("/tmp/x")


def test_paths_touched_covers_read_grep_and_glob():
    lines = [
        _assistant(_tool("Read", file_path="/x/demo_report/report.md")),
        _assistant(_tool("Grep", pattern="def main", path="skills/pharmgx-reporter")),
        _assistant(_tool("Glob", pattern="skills/*/SKILL.md")),
        RESULT,
    ]
    assert MOD.parse_stream(lines).paths_touched == [
        "/x/demo_report/report.md", "skills/pharmgx-reporter", "skills/*/SKILL.md",
    ]


def test_missing_result_event_means_no_final_answer():
    """A run that died mid-stream must not read as an empty-but-valid answer."""
    trace = MOD.parse_stream([_assistant({"type": "text", "text": "partial"})])
    assert trace.final is None


def test_blank_and_non_json_lines_are_skipped():
    trace = MOD.parse_stream(["", "Warning: no stdin data received", RESULT])
    assert trace.final == "Done."


noise = st.sampled_from([
    json.dumps({"type": "system", "subtype": "thinking_tokens"}),
    json.dumps({"type": "rate_limit_event"}),
    _assistant({"type": "text", "text": "thinking aloud"}),
])


@given(commands=st.lists(st.text(min_size=1, max_size=20), max_size=6),
       gaps=st.lists(st.lists(noise, max_size=3), min_size=7, max_size=7))
def test_bash_commands_keep_order_whatever_events_sit_between(commands, gaps):
    lines = []
    for gap, cmd in zip(gaps, commands):
        lines += gap + [_assistant(_tool("Bash", command=cmd))]
    lines += gaps[-1] + [RESULT]
    assert [c.command for c in MOD.parse_stream(lines).calls] == commands


def _result(tool_id, content, is_error=False):
    block = {"type": "tool_result", "tool_use_id": tool_id, "content": content, "is_error": is_error}
    return json.dumps({"type": "user", "message": {"content": [block]}})


# Outcome of each Bash call, as Claude Code 2.1.303 reports it: a denial is a
# permission_denied event plus an is_error result; a non-zero exit is an is_error result.
OUTCOMES = {
    "ok": lambda i: [_result(i, "done")],
    "failed": lambda i: [_result(i, "Exit code 1\nModuleNotFoundError: No module named 'sushie'", True)],
    "denied": lambda i: [_denied(i), _result(i, "Permission to use Bash has been denied.", True)],
}


@given(calls=st.lists(st.tuples(st.text(min_size=1, max_size=20), st.sampled_from(sorted(OUTCOMES))), max_size=6))
def test_ran_and_failed_follow_each_calls_outcome(calls):
    """A crashed skill ran but must not read as a working run; a denied one never ran."""
    lines = []
    for n, (cmd, outcome) in enumerate(calls):
        lines += [_bash(f"t{n}", cmd), *OUTCOMES[outcome](f"t{n}")]
    trace = MOD.parse_stream(lines + [RESULT])
    assert [(c.command, c.outcome) for c in trace.calls] == calls
    assert [c.output for c in trace.calls if c.outcome == "ok"] == ["done" for _, o in calls if o == "ok"]


VIOLATION = ("HTTP/1.1 403 Forbidden\n<sandbox_violations>\n"
             "deny network-outbound example.com:443 (user denied)\n"
             "deny network-outbound api.ncbi.nlm.nih.gov:443 (user denied)\n</sandbox_violations>")


@pytest.mark.parametrize("content", [VIOLATION, [{"type": "text", "text": VIOLATION}]])
def test_sandbox_network_violations_are_recorded(content):
    """A skill reaching an unlisted host is a local-first finding, even when the command exits 0."""
    trace = MOD.parse_stream([_bash("t1", "python clawbio.py run x --demo"), _result("t1", content), RESULT])
    assert trace.network_blocked == ["example.com:443", "api.ncbi.nlm.nih.gov:443"]
    assert [c.outcome for c in trace.calls] == ["ok"]


def test_each_run_starts_from_a_fresh_clone(monkeypatch, tmp_path):
    """Leftover files, git history and auto-memory from earlier runs must not reach the agent."""
    seen = {}
    real_run = MOD.subprocess.run

    def fake_run(cmd, **kw):
        if "status" in cmd:
            return MOD.subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")
        if cmd[0] != "claude":
            return real_run(cmd, **kw)
        seen["cwd"] = Path(kw["cwd"])
        seen["has_agents_md"] = (seen["cwd"] / "AGENTS.md").is_file()
        seen["env"] = kw["env"]
        return MOD.subprocess.CompletedProcess(cmd, 0, stdout=RESULT + "\n", stderr="")

    monkeypatch.setattr(MOD.subprocess, "run", fake_run)
    trace = MOD.run_agent("hi", tmp_path / "log.jsonl")
    assert seen["cwd"] != MOD.REPO and seen["has_agents_md"]
    assert not seen["cwd"].exists()  # removed after the run
    # `uv run` in the clone reuses the repo's venv instead of building one per run.
    assert seen["env"]["UV_NO_SYNC"] == "1"
    assert trace.final == "Done."


def test_run_trials_returns_one_trace_per_trial_with_its_own_log(monkeypatch, tmp_path):
    logs = []
    monkeypatch.setattr(MOD, "run_agent", lambda prompt, log_path, timeout=600: logs.append(log_path) or MOD.Trace())
    traces = MOD.run_trials("hi", tmp_path, n=3)
    assert len(traces) == 3 and len(set(logs)) == 3


def test_refuses_to_run_with_uncommitted_changes(monkeypatch, tmp_path):
    """The clone only has committed files, so an uncommitted AGENTS.md edit would go untested."""
    def fake_run(cmd, **kw):
        assert cmd[0] == "git" and "status" in cmd, f"ran {cmd} on a dirty checkout"
        return MOD.subprocess.CompletedProcess(cmd, 0, stdout=" M AGENTS.md\n", stderr="")

    monkeypatch.setattr(MOD.subprocess, "run", fake_run)
    with pytest.raises(RuntimeError, match="AGENTS.md"):
        MOD.run_agent("hi", tmp_path / "log.jsonl")


def test_refuses_to_write_a_transcript_inside_the_repo(monkeypatch):
    """Raw logs can hold local paths and the user's ~/.claude/CLAUDE.md; they must not get committed."""
    monkeypatch.setattr(MOD.subprocess, "run", lambda cmd, **kw: pytest.fail(f"ran {cmd}"))
    with pytest.raises(ValueError, match="outside the repo"):
        MOD.run_agent("hi", MOD.REPO / "tests" / "agent_evals" / "log.jsonl")


# A pipe reports its last command, so `skill | tail` exits 0 when the skill crashed; the output still shows it.
CRASHED = {
    "clean": "Status: OK\nreport written",
    "traceback": "Traceback (most recent call last):\n  ...\nKeyError: 'gene'",
    "clawbio_failed": "Status:   FAILED\nexit code 1",
}


@given(calls=st.lists(st.tuples(st.sampled_from(["python clawbio.py run pharmgx", "python clawbio.py run nutrigx",
                                                  "cat skills/pharmgx-reporter/SKILL.md"]),
                                st.sampled_from(sorted(OUTCOMES)), st.sampled_from(sorted(CRASHED))), max_size=5))
@example(calls=[("python clawbio.py run pharmgx", "failed", "clean"),
                ("python clawbio.py run pharmgx", "ok", "clean")])  # a retry after a failure
@example(calls=[("cat skills/pharmgx-reporter/SKILL.md", "ok", "clean")])  # reading is not running
def test_ran_ok_means_some_matching_call_exited_0_with_clean_output(calls):
    lines = []
    for n, (cmd, outcome, out) in enumerate(calls):
        events = OUTCOMES[outcome](f"t{n}")
        if outcome == "ok":
            events = [_result(f"t{n}", CRASHED[out])]
        lines += [_bash(f"t{n}", cmd), *events]
    trace = MOD.parse_stream(lines + [RESULT])
    want = any(c == "python clawbio.py run pharmgx" and o == "ok" and out == "clean" for c, o, out in calls)
    assert trace.ran_ok("pharmgx") == want


def test_only_bash_output_can_report_a_blocked_host():
    """Reading a file that quotes a sandbox message (like this one) is not a network leak."""
    read = _assistant({"type": "tool_use", "id": "r1", "name": "Read", "input": {"file_path": "x.py"}})
    assert MOD.parse_stream([read, _result("r1", VIOLATION), RESULT]).network_blocked == []


def test_a_grep_pattern_is_not_a_path():
    """Searching for the word demo_report is not opening the stale report."""
    trace = MOD.parse_stream([_assistant(_tool("Grep", pattern="demo_report")), RESULT])
    assert trace.paths_touched == []


def _fake_claude(monkeypatch, result):
    def fake_run(cmd, **kw):
        if "status" in cmd:
            return MOD.subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")
        if cmd[0] != "claude":
            return MOD.subprocess.CompletedProcess(cmd, 0)  # git clone
        if isinstance(result, Exception):
            raise result
        return result
    monkeypatch.setattr(MOD.subprocess, "run", fake_run)


def test_a_claude_error_lands_in_the_log(monkeypatch, tmp_path):
    """An expired sign-in must say so in the log, not just leave it empty."""
    _fake_claude(monkeypatch, MOD.subprocess.CompletedProcess(["claude"], 1, stdout="", stderr="Not logged in"))
    MOD.run_agent("hi", tmp_path / "log.jsonl")
    assert "Not logged in" in (tmp_path / "log.jsonl").read_text()


def test_a_timed_out_run_keeps_its_partial_log(monkeypatch, tmp_path):
    partial = _assistant({"type": "text", "text": "still going"})
    _fake_claude(monkeypatch, MOD.subprocess.TimeoutExpired(["claude"], 600, output=partial.encode()))
    with pytest.raises(MOD.subprocess.TimeoutExpired):
        MOD.run_agent("hi", tmp_path / "log.jsonl")
    assert "still going" in (tmp_path / "log.jsonl").read_text()
