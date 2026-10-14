"""Every flag in a skill's INT-001 allowlist must exist in that skill's argparse parser.

``run_skill`` forwards allowlisted extra flags to the skill script. A flag that is
allowlisted but not defined by the script's parser passes the filter and then
makes argparse exit with ``unrecognized arguments``, so the allowlist promises a
flag the skill cannot take.

Each script's parser is read by running the script's ``__main__`` in a subprocess
with ``ArgumentParser.parse_args`` patched to hand back the parser instead of
parsing, so nothing past argument parsing runs.
"""
from __future__ import annotations

import argparse
import json
import os
import runpy
import subprocess
import sys
from pathlib import Path

import pytest

from clawbio.cli import SKILLS

ROOT = Path(__file__).resolve().parents[1]

# Allowlist entries that do not match the script today. Strict xfail: once an
# entry is fixed its test XPASSes and fails, so remove it from here.
KNOWN_DRIFT = {
    # Registered alongside the IRB approval note (e1037c12) but never committed;
    # whether to keep this entry is a maintainer decision.
    "llm-bench": "registered script skills/llm-biobank-bench/ does not exist",
}


class _Captured(BaseException):  # BaseException: some scripts wrap main() in `except Exception`
    pass


def _probe(script: str) -> None:
    """Subprocess mode: print the JSON list of option strings the script's parser defines."""

    def grab(self, *args, **kwargs):
        raise _Captured(self)

    argparse.ArgumentParser.parse_args = grab
    argparse.ArgumentParser.parse_known_args = grab
    sys.argv = [script]
    sys.path.insert(0, str(Path(script).parent))
    try:
        runpy.run_path(script, run_name="__main__")
    except _Captured as captured:
        print(json.dumps({"flags": sorted(s for a in captured.args[0]._actions for s in a.option_strings)}))
    except ModuleNotFoundError as exc:
        print(json.dumps({"missing_module": exc.name}))
    except SystemExit as exc:  # a Click/Typer CLI exiting on its own validation
        print(json.dumps({"no_parser": repr(exc)}))


def _cases():
    for name, info in SKILLS.items():
        allowed = set(info.get("allowed_extra_flags") or ()) | set(
            info.get("allowed_extra_flags_without_values") or ()
        )
        if not allowed:
            continue
        marks = [pytest.mark.xfail(reason=KNOWN_DRIFT[name], strict=True)] if name in KNOWN_DRIFT else []
        yield pytest.param(name, Path(info["script"]), allowed, id=name, marks=marks)


@pytest.mark.parametrize("name, script, allowed", list(_cases()))
def test_allowlisted_flags_exist_in_parser(name, script, allowed):
    assert script.is_file(), f"{name}: registered script {script} does not exist"
    proc = subprocess.run(
        [sys.executable, __file__, str(script)],
        capture_output=True, text=True, timeout=120, cwd=ROOT,
        env={**os.environ, "PYTHONPATH": str(ROOT)},
    )
    lines = proc.stdout.strip().splitlines()
    assert lines and lines[-1].startswith("{"), f"{name}: probe produced no result\n{proc.stderr[-2000:]}"
    result = json.loads(lines[-1])
    if "missing_module" in result:
        pytest.skip(f"{name}: optional dependency {result['missing_module']!r} not installed")
    if "no_parser" in result:
        pytest.skip(f"{name}: no argparse parser reached: {result['no_parser']}")
    unknown = sorted(allowed - set(result["flags"]))
    assert not unknown, f"{name}: allowlisted flags not accepted by {script.name}: {unknown}"


if __name__ == "__main__":
    _probe(sys.argv[1])
