"""turingdb 3.x host handling and the shared reproducibility bundle."""

from __future__ import annotations

import argparse
import importlib.util
import shlex
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL_DIR))

import turingdb_graph as tg  # noqa: E402

from clawbio.common.checksums import sha256_file  # noqa: E402


@pytest.fixture(autouse=True)
def _cwd_in_tmp(tmp_path: Path, monkeypatch) -> None:
    """Run each test, and any subprocess it starts, from its own tmp_path (AGENTS.md)."""
    monkeypatch.chdir(tmp_path)


def _checksums(out: Path) -> dict[str, str]:
    lines = (out / "reproducibility" / "checksums.sha256").read_text().splitlines()
    return {label: digest for digest, label in (line.split("  ", 1) for line in lines)}


@pytest.mark.parametrize(
    ("host", "expected"),
    [
        ("localhost:6666", ("localhost", 6666)),
        ("localhost", ("localhost", 6666)),
        ("http://localhost:6666", ("localhost", 6666)),
        ("db.example.org:7000", ("db.example.org", 7000)),
    ],
)
def test_split_host(host: str, expected: tuple[str, int]) -> None:
    assert tg.split_host(host) == expected


def test_output_and_out_both_set_the_output_dir(tmp_path: Path, monkeypatch) -> None:
    seen = []
    monkeypatch.setattr(tg.Connection, "connect", lambda self, auto_start=True: None)
    monkeypatch.setattr(tg, "_run_subcommand", lambda parser, args, conn: seen.append(args.out))
    monkeypatch.setattr(tg, "write_reproducibility", lambda args, out: None)
    tg.main(["--demo", "cohort", "--output", str(tmp_path / "a")])
    tg.main(["--demo", "cohort", "--out", str(tmp_path / "b")])
    assert seen == [tmp_path / "a", tmp_path / "b"]


def test_bundle_checksums_every_output(tmp_path: Path) -> None:
    out = tmp_path / "out"
    (out / "analysis").mkdir(parents=True)
    (out / "report.md").write_text("# report\n")
    (out / "analysis" / "summary.json").write_text("{}\n")
    args = argparse.Namespace(
        host=tg.DEFAULT_HOST, no_auto_start=False, demo=None, build=False, query=True,
        analyse_cohort=False, graph="g1", cypher="MATCH (n) WHERE n.x = 'a' RETURN n",
        result_format="md",
    )
    tg.write_reproducibility(args, out)

    entries = _checksums(out)
    assert set(entries) == {
        "report.md",
        "analysis/summary.json",
        "reproducibility/commands.sh",
        "reproducibility/environment.yml",
    }
    for label, digest in entries.items():
        assert digest == sha256_file(out / label), label

    commands = (out / "reproducibility" / "commands.sh").read_text()
    assert '--output \\\n  "$OUTPUT_DIR"' in commands
    # The Cypher string is shell-quoted, so the shell hands it back unchanged.
    cypher_line = next(line for line in commands.splitlines() if "MATCH" in line)
    assert shlex.split(cypher_line.rstrip(" \\")) == [args.cypher]
    assert "turingdb>=3.1,<4" in (out / "reproducibility" / "environment.yml").read_text()


@pytest.mark.skipif(
    importlib.util.find_spec("turingdb") is None,
    reason="needs the turingdb-graph extra (uv sync --extra turingdb-graph)",
)
def test_demo_runs_in_process_and_writes_bundle(tmp_path: Path) -> None:
    # A subprocess, as the sweep runs it: turingdb allows one embedded engine per process.
    out = tmp_path / "out"
    subprocess.run(
        [sys.executable, str(SKILL_DIR / "turingdb_graph.py"), "--demo", "cohort", "--output", str(out)],
        check=True, capture_output=True, text=True, cwd=tmp_path,
    )
    entries = _checksums(out)
    assert "demo-cohort/analysis/report.md" in entries
    for label, digest in entries.items():
        assert digest == sha256_file(out / label), label
    assert sorted(p.name for p in tmp_path.iterdir()) == ["out"], "demo wrote outside --output"
    commands = (out / "reproducibility" / "commands.sh").read_text()
    assert "--host" not in commands and "--no-auto-start" not in commands


# ---------------------------------------------------------------------------
# Properties
# ---------------------------------------------------------------------------

from hypothesis import HealthCheck, example, given, settings  # noqa: E402
from hypothesis import strategies as st  # noqa: E402

_label = st.from_regex(r"[a-z0-9]([a-z0-9-]{0,14}[a-z0-9])?", fullmatch=True)
_hostname = st.one_of(
    st.lists(_label, min_size=1, max_size=4).map(".".join),
    st.tuples(*[st.integers(0, 255)] * 4).map(lambda t: ".".join(map(str, t))),
)


@given(hostname=_hostname, port=st.integers(1, 65535), scheme=st.sampled_from(["", "http://"]))
def test_split_host_recovers_host_and_port(hostname: str, port: int, scheme: str) -> None:
    assert tg.split_host(f"{scheme}{hostname}:{port}") == (hostname, port)
    assert tg.split_host(f"{scheme}{hostname}") == (hostname, tg.DEFAULT_PORT)


# Any text a user can pass on the command line: no NUL, no lone surrogates.
_cli_text = st.text(st.characters(blacklist_categories=("Cs",), blacklist_characters="\x00"), min_size=1)


@settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(graph=_cli_text, cypher=_cli_text, host=_cli_text)
# The LF-only writer behind commands.sh turned a quoted "\r" into "\n" on replay.
@example(graph="\r", cypher="0", host="0")
def test_commands_sh_hands_free_text_back_unchanged(tmp_path_factory, graph, cypher, host) -> None:
    """bash, running commands.sh, passes the original strings to the skill."""
    root = tmp_path_factory.mktemp("prop")
    out = root / "out"
    out.mkdir()
    args = argparse.Namespace(
        host=host, no_auto_start=False, demo=None, build=False, query=True,
        analyse_cohort=False, graph=graph, cypher=cypher, result_format="md",
    )
    tg.write_reproducibility(args, out)
    # Stand in for python: print the argv bash built, one NUL-terminated field each.
    fake_python = root / "fake_python"
    fake_python.write_text('#!/usr/bin/env bash\nprintf "%s\\0" "$@"\n')
    fake_python.chmod(0o755)
    before = set(root.iterdir())
    run = subprocess.run(
        ["bash", str(out / "reproducibility" / "commands.sh")],
        env={"PATH": "/usr/bin:/bin", "PYTHON": str(fake_python)},
        capture_output=True, check=True, cwd=root,
    )
    assert set(root.iterdir()) == before, "replay created files next to the output"
    argv = run.stdout.decode("utf-8", "surrogateescape").split("\0")[:-1]
    assert argv[argv.index("--host") + 1] == host
    assert argv[argv.index("--graph") + 1] == graph
    assert argv[argv.index("--cypher") + 1] == cypher
    assert argv[-2:] == ["--output", str(out)]

