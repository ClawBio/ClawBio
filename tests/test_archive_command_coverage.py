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

import contextlib
import importlib.util
import io
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


@pytest.mark.parametrize("skill", SKILLS)
def test_every_advertised_command_is_implemented(skill):
    """`--command X` must never reach the vendored CLI as an unknown subcommand.

    Probing with the bare subcommand is enough: argparse resolves the subparser
    before it validates that subparser's own arguments, so a *known* command
    fails with "the following arguments are required" and an *unknown* one
    fails with "invalid choice". Neither reaches the network.
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
                contextlib.suppress(SystemExit, Exception):
            app.api.main([cmd])
        if "invalid choice" in err.getvalue():
            unreachable.append(cmd)

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
