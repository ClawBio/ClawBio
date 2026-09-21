"""Shared machinery for the archive-fetch skills.

`biostudies-fetch`, `ena-fetch`, `geo-fetch`, `arrayexpress-fetch` and
`pride-fetch` all wrap a vendored, subcommand-driven upstream client the same
way: expose a runner-reachable `--command` form, capture what the upstream
prints, and write ClawBio's output contract around it.

This is a toolkit, not a framework. Each skill keeps its own `main()` and its
own argv translation, because that is where they genuinely differ; everything
below is what would otherwise be copy-pasted five times.

Note the split from the vendored modules: upstream code stays close to its
source so it can be re-synced, and *new* shared code lives here.
"""

from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
import shlex
import sys
from pathlib import Path

from clawbio.common.report import DISCLAIMER
from clawbio.common.reproducibility import (
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)

# Artifacts every archive skill promises in its Output Structure tree.
ARTIFACTS = (
    "report.md",
    "result.json",
    "tables/metadata.tsv",
    "reproducibility/commands.sh",
    "reproducibility/environment.yml",
    "reproducibility/checksums.sha256",
)


def load_sibling(skill_dir: Path, name: str):
    """Import a module from the skill directory without a package structure."""
    spec = importlib.util.spec_from_file_location(name, skill_dir / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def common_parser(prog: str, description: str, commands: tuple[str, ...]) -> argparse.ArgumentParser:
    """The flags every archive skill shares.

    `--command` exists because `clawbio/cli.py` filters extra arguments against
    a per-skill allowlist of *flags*; a bare positional subcommand is silently
    dropped, so without it a registered alias could only ever run `--demo`.
    """
    p = argparse.ArgumentParser(
        prog=prog, description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--demo", action="store_true",
                   help="run offline against the bundled fixture")
    p.add_argument("--output", metavar="DIR", help="output directory")
    p.add_argument("--command", choices=commands,
                   help="upstream subcommand to run (runner-reachable form)")
    p.add_argument("--accession", help="archive accession")
    p.add_argument("--query", help="search terms (--command search)")
    p.add_argument("--limit", type=int, default=20, help="search hits to return")
    p.add_argument("--json", action="store_true", help="emit raw JSON in the report")
    p.add_argument("--out", metavar="PATH",
                   help="upstream --out alias; resolved under --output, not cwd")
    return p


def expand_positional(argv: list[str], commands: tuple[str, ...],
                      query_commands: tuple[str, ...] = ("search",)) -> list[str]:
    """Translate the upstream positional form into the --command form.

    `metadata S-BSST2074` becomes `--command metadata --accession S-BSST2074`,
    so anyone following upstream's documentation is not stranded.
    """
    if not argv or argv[0] not in commands:
        return argv
    command = argv[0]
    rest = argv[1:]
    flag = "--query" if command in query_commands else "--accession"
    if rest and not rest[0].startswith("-"):
        return ["--command", command, flag, rest[0], *rest[1:]]
    return ["--command", command, *rest]


def resolve_out(user_out: str | None, output_dir: Path, default_rel: str) -> Path:
    """Anchor upstream's cwd-relative --out under --output.

    Upstream defaults land in the working directory (`metadata.tsv`,
    `./ena_out`, `--out .`). Under ClawBio nothing may be written outside the
    directory the caller named. An absolute --out is still honoured.
    """
    if user_out:
        candidate = Path(user_out).expanduser()
        resolved = candidate if candidate.is_absolute() else output_dir / candidate
    else:
        resolved = output_dir / default_rel
    resolved.parent.mkdir(parents=True, exist_ok=True)
    return resolved


def scrub(text: str, output_dir: Path) -> str:
    """Replace the absolute output path with a stable placeholder.

    Upstream prints the paths it wrote. Left alone those absolute paths land in
    report.md and result.json, so the same run in two directories yields two
    different reports and two different checksums. The report should describe
    what was produced, not where this invocation happened to put it.
    """
    return text.replace(str(output_dir), "<output>")


def run_upstream(api, argv: list[str], output_dir: Path) -> tuple[str, str]:
    """Run the vendored CLI, capturing what it prints."""
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        api.main(argv)
    return scrub(out.getvalue(), output_dir), scrub(err.getvalue(), output_dir)


def warn_if_overwriting(output_dir: Path, artifacts: tuple[str, ...] = ARTIFACTS) -> None:
    """CLAUDE.md Safety Rule 4: warn before overwriting existing reports."""
    if any((output_dir / name).exists() for name in artifacts):
        print(f"[warning] Existing output files will be overwritten in {output_dir}",
              file=sys.stderr)


def write_report(output_dir: Path, *, title: str, source_url: str,
                 sections: list[tuple[str, str]]) -> Path:
    """Write report.md: one fenced block per command, plus the disclaimer."""
    lines = [f"# {title}", "", f"Source: {source_url}", ""]
    for heading, body in sections:
        lines += [f"## {heading}", "", "```", body.rstrip() or "(no output)", "```", ""]
    lines += ["---", "", f"*{DISCLAIMER}*", ""]
    path = output_dir / "report.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def write_bundle(output_dir: Path, *, script: Path, env_name: str,
                 command: str, written: list[Path]) -> None:
    """Write reproducibility/ through the shared helpers.

    Never hand-roll these: the bundle layout, line endings and checksum format
    have to match every other ClawBio skill.
    """
    cmd = ["python", str(script), "--command", command, "--output", str(output_dir)]
    commands_sh = write_commands_sh(output_dir, shlex.join(cmd))
    environment_yml = write_environment_yml(
        output_dir, env_name,
        [],  # the vendored clients are standard library only
        python_version=f"{sys.version_info.major}.{sys.version_info.minor}",
    )
    write_checksums([*written, commands_sh, environment_yml], output_dir,
                    anchor=output_dir)
