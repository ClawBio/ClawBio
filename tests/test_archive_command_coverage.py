"""Every command an archive skill advertises must actually be reachable.

This is the test whose absence let `arrayexpress-fetch --command download-script`
ship: it was listed in `COMMANDS`, had no branch in `_to_upstream_argv`, and no
subparser in the vendored CLI, so it reached argparse as an unknown subcommand
and died with `SystemExit(2)`.

The check is deliberately cross-skill rather than per-skill. The defect is a
mismatch between two lists that live in different files, and a per-skill test
would have to be remembered five times.

Two ways a command can be legitimately implemented:

* the vendored CLI has a subparser for it, or
* the entry point handles it itself, before delegating -- these are declared in
  the module's `LOCAL_COMMANDS`, which exists so this test can tell "handled
  here" apart from "handled nowhere".
"""
from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
import socket
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SKILLS = ["arrayexpress-fetch", "biostudies-fetch", "ena-fetch",
          "geo-fetch", "pride-fetch"]


def _load_entry_point(skill: str):
    """Import a skill's entry point by path, as the skill's own tests do."""
    module_name = skill.replace("-", "_")
    path = ROOT / "skills" / skill / f"{module_name}.py"
    skill_dir = str(path.parent)
    if skill_dir not in sys.path:
        sys.path.insert(0, skill_dir)
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


class _NetworkBlocked(BaseException):
    """A BaseException, so the vendored `except Exception` retry loops cannot
    swallow it and sleep before retrying."""


class _Parsed(BaseException):
    """Raised once argparse has resolved the subcommand, so no command runs."""


@pytest.fixture
def parse_only(monkeypatch):
    """Stop the vendored `main()` straight after `parse_args` succeeds.

    Reaching that point is the whole claim under test. Letting the command body
    run as well made `ena-fetch fields`, which takes no arguments, query
    ebi.ac.uk for real.
    """
    real = argparse.ArgumentParser.parse_args

    def parse_then_stop(self, *args, **kwargs):
        real(self, *args, **kwargs)
        raise _Parsed

    monkeypatch.setattr(argparse.ArgumentParser, "parse_args", parse_then_stop)


@pytest.fixture
def no_network(monkeypatch):
    """Record, and refuse, every outbound connection the test attempts."""
    attempts = []

    def refuse(self, address, *args, **kwargs):
        attempts.append(address)
        raise _NetworkBlocked(address)

    monkeypatch.setattr(socket.socket, "connect", refuse)
    return attempts


@pytest.mark.parametrize("skill", SKILLS)
def test_every_advertised_command_is_implemented(skill, parse_only, no_network):
    """`--command X` must never reach the vendored CLI as an unknown subcommand.

    Probing with the bare subcommand is enough: argparse resolves the subparser
    before it validates that subparser's own arguments, so a *known* command
    either parses or fails with "the following arguments are required", and an
    *unknown* one fails with "invalid choice". `parse_only` stops a command that
    parses before it runs; `no_network` proves that nothing reached the network.
    """
    app = _load_entry_point(skill)
    local = set(getattr(app, "LOCAL_COMMANDS", ()))

    unreachable = []
    for cmd in app.COMMANDS:
        if cmd in local:
            continue
        err = io.StringIO()
        with contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(err), \
                contextlib.suppress(SystemExit, _Parsed):
            app.api.main([cmd])
        if "invalid choice" in err.getvalue():
            unreachable.append(cmd)

    assert not no_network, (
        f"{skill}: probing COMMANDS attempted network connections {no_network}; "
        "this suite must stay hermetic")

    assert not unreachable, (
        f"{skill} advertises {unreachable} in COMMANDS but the vendored CLI has "
        f"no subparser for them and they are not in LOCAL_COMMANDS. Either "
        f"implement them, handle them in the entry point, or stop advertising them.")


@pytest.mark.parametrize("skill", SKILLS)
def test_locally_handled_commands_are_advertised(skill):
    """The converse: LOCAL_COMMANDS must not name something COMMANDS omits,
    which would be a branch no user can reach."""
    app = _load_entry_point(skill)
    local = set(getattr(app, "LOCAL_COMMANDS", ()))
    orphaned = sorted(local - set(app.COMMANDS))
    assert not orphaned, f"{skill}: LOCAL_COMMANDS names {orphaned}, absent from COMMANDS"


@pytest.mark.parametrize("skill", SKILLS)
def test_search_still_defaults_to_twenty_hits(skill, tmp_path):
    """The shared `--limit` default is a sentinel so each command can pick its
    own. Search must be unaffected by that change.

    The sentinel exists because `--limit` meant two things: a search hit cap,
    and a row cap on `ena-fetch --command report`, whose vendored default is
    `0 = no limit`. A global default of 20 silently truncated a 95-run report.
    Supplying 20 at each search call site keeps the behaviour and makes it
    deliberate.
    """
    app = _load_entry_point(skill)
    if "search" not in app.COMMANDS:
        pytest.skip(f"{skill} has no search command")

    args = app._build_parser().parse_args(
        ["--command", "search", "--query", "microglia"])
    argv = app._to_upstream_argv(args, tmp_path)

    assert "--limit" in argv, f"{skill} search sends no --limit"
    assert argv[argv.index("--limit") + 1] == "20", (
        f"{skill} search default changed; the sentinel must not alter search")


@pytest.mark.parametrize("skill", SKILLS)
def test_every_runner_allowlisted_flag_exists_on_the_skill(skill):
    """`clawbio.py run` forwards only allowlisted flags (INT-001). A flag on the
    list that the skill's parser lacks passes the runner and then dies in
    argparse, so the allowlist would be advertising something that cannot work."""
    import clawbio.cli as cli

    app = _load_entry_point(skill)
    known = set(app._build_parser()._option_string_actions)
    entry = cli.SKILLS[skill]
    allowed = set(entry.get("allowed_extra_flags", ())) | set(
        entry.get("allowed_extra_flags_without_values", ()))
    missing = sorted(allowed - known)
    assert not missing, f"clawbio/cli.py allowlists {missing} for {skill}, which its parser lacks"


@pytest.mark.parametrize("skill", SKILLS)
def test_every_allowlisted_flag_survives_the_runner(skill, monkeypatch, tmp_path):
    """Pass each allowlisted flag through the real INT-001 filter in run_skill.

    `allowed_extra_flags_without_values` only says which allowed flags take no
    value. A flag listed there but not in `allowed_extra_flags` is dropped
    without a word, so `clawbio.py run geo-fetch --suppl` quietly downloaded
    the series matrix instead of the supplementary files.
    """
    import subprocess

    import clawbio.cli as cli

    entry = cli.SKILLS[skill]
    no_value = set(entry.get("allowed_extra_flags_without_values", ()))
    with_value = set(entry.get("allowed_extra_flags", ())) - no_value
    extra = []
    for flag in sorted(with_value):
        extra += [flag, "x"]
    extra += sorted(no_value)

    seen = []

    def fake_run(cmd, *args, **kwargs):
        seen.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(cli.subprocess, "run", fake_run)
    cli.run_skill(skill_name=skill, output_dir=str(tmp_path), extra_args=extra)
    forwarded = set(seen[0])
    dropped = sorted((with_value | no_value) - forwarded)
    assert not dropped, f"clawbio.py run silently drops {dropped} for {skill}"


@pytest.mark.parametrize("skill", SKILLS)
def test_boolean_flags_never_swallow_the_next_token(skill, monkeypatch, tmp_path):
    """A store_true flag registered as taking a value makes the runner forward
    the next token as its argument, so `--json stray` smuggles `stray` through
    the allowlist."""
    import subprocess

    import clawbio.cli as cli

    app = _load_entry_point(skill)
    parser = app._build_parser()
    boolean = sorted(
        opt for action in parser._actions if action.nargs == 0
        for opt in action.option_strings
        if opt in cli.SKILLS[skill].get("allowed_extra_flags", ()))
    assert boolean, f"{skill}: expected at least --json among the allowlisted booleans"

    seen = []
    monkeypatch.setattr(cli.subprocess, "run", lambda cmd, *a, **k: (
        seen.append(cmd), subprocess.CompletedProcess(cmd, 0, "", ""))[1])
    for flag in boolean:
        seen.clear()
        cli.run_skill(skill_name=skill, output_dir=str(tmp_path),
                      extra_args=[flag, "stray-token"])
        assert "stray-token" not in seen[0], f"{skill}: {flag} swallowed the next token"
