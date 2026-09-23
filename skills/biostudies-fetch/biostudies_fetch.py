#!/usr/bin/env python3
"""biostudies-fetch — ClawBio entry point for EMBL-EBI BioStudies.

Ported from UKDRI/informatics_data_skills @ 7cc3e6e.
Copyright (c) 2026 UK Dementia Research Institute. Licensed MIT.

The archive logic lives in `biostudies_fetch_api.py`, vendored close to
upstream. The shared ClawBio machinery lives in `clawbio.common.archive_fetch`.
This module holds only what is specific to BioStudies: its subcommands, how
they map onto the upstream parser, and how `--demo` is served offline.

    python skills/biostudies-fetch/biostudies_fetch.py --demo --output /tmp/bs
    python skills/biostudies-fetch/biostudies_fetch.py \
        --command metadata --accession S-BSST2074 --output /tmp/bs
    python skills/biostudies-fetch/biostudies_fetch.py metadata S-BSST2074 --output /tmp/bs
"""

from __future__ import annotations

import sys
from pathlib import Path

_SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SKILL_DIR.parent.parent
for _p in (str(_PROJECT_ROOT), str(_SKILL_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from clawbio.common import archive_fetch as af  # noqa: E402
from clawbio.common.report import write_result_json  # noqa: E402

api = af.load_sibling(_SKILL_DIR, "biostudies_fetch_api")

SKILL = "biostudies-fetch"
VERSION = "0.1.0"
COMMANDS = ("metadata", "files", "download", "search", "metadata-table")
# Search hits shown in report.md. Explicit here rather than inherited from
# the shared parser, which has no default -- see archive_fetch.common_parser.
SEARCH_LIMIT = 20

DEMO_ACCESSION = "S-BSST2074"
DEMO_STUDY_FIXTURE = _SKILL_DIR / "examples" / f"demo_{DEMO_ACCESSION}.json"
DEMO_SEARCH_FIXTURE = _SKILL_DIR / "examples" / "demo_search.json"
DEMO_COMMANDS = ("metadata", "files", "metadata-table", "search")


def _build_parser():
    p = af.common_parser("biostudies_fetch.py", __doc__, COMMANDS)
    p.add_argument("--collection", help="restrict a search to one collection")
    p.add_argument("--match", help="only files whose path contains this substring")
    return p


def _to_upstream_argv(args, output_dir: Path) -> list[str]:
    """Translate the --command form into the vendored parser's argv."""
    cmd = args.command
    if cmd == "search":
        if not args.query:
            raise SystemExit("--command search needs --query")
        argv = ["search", args.query]
        if args.collection:
            argv += ["--collection", args.collection]
        argv += ["--limit", str(SEARCH_LIMIT if args.limit is None else args.limit)]
        if args.json:
            argv.append("--json")
        return argv

    if not args.accession:
        raise SystemExit(f"--command {cmd} needs --accession")
    argv = [cmd, args.accession]
    if cmd in ("metadata", "files") and args.json:
        argv.append("--json")
    elif cmd == "download":
        if args.match:
            argv += ["--match", args.match]
        argv += ["--out", str(af.resolve_out(args.out, output_dir, "downloads"))]
    elif cmd == "metadata-table":
        argv += ["--out", str(af.resolve_out(args.out, output_dir, "tables/metadata.tsv"))]
    return argv


def _install_demo_transport() -> None:
    """Serve the vendored module's HTTP layer from the committed fixtures.

    Patching `http_get` rather than the command functions keeps `--demo` on the
    same code path as a live run, so the demo exercises the real parsing.
    """
    study = DEMO_STUDY_FIXTURE.read_bytes()
    search = DEMO_SEARCH_FIXTURE.read_bytes()

    def _offline_get(url: str, retries: int = 3) -> bytes:
        if "/search" in url:
            return search
        if "/studies/" in url:
            return study
        if "/biosamples/" in url:
            return b"{}"
        raise SystemExit(f"demo mode has no fixture for {url}")

    api.http_get = _offline_get


def main(argv: list[str] | None = None) -> int:
    argv = af.expand_positional(
        list(sys.argv[1:] if argv is None else argv), COMMANDS)
    parser = _build_parser()
    args = parser.parse_args(argv)

    if not args.demo and not args.command:
        parser.error("give --demo, or --command with its arguments")
    if not args.output:
        parser.error("--output is required")

    output_dir = Path(args.output).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    af.warn_if_overwriting(output_dir)

    if args.demo:
        _install_demo_transport()
        args.accession = args.accession or DEMO_ACCESSION
        args.query = args.query or "N-masked reference genome"

    accession = args.accession or DEMO_ACCESSION
    # A demo run exercises the whole documented tree, not the cheapest command.
    commands = list(DEMO_COMMANDS) if args.demo else [args.command]

    sections: list[tuple[str, str]] = []
    for cmd in commands:
        args.command = cmd
        stdout, stderr = af.run_upstream(api, _to_upstream_argv(args, output_dir), output_dir)
        sections.append((cmd, stdout or stderr))

    report = af.write_report(
        output_dir,
        title=f"BioStudies report — {accession}",
        source_url=f"[EMBL-EBI BioStudies](https://www.ebi.ac.uk/biostudies/studies/{accession})",
        sections=sections,
    )
    written = [report]

    table = output_dir / "tables" / "metadata.tsv"
    if table.exists():
        written.append(table)

    written.append(write_result_json(
        output_dir,
        skill=SKILL,
        version=VERSION,
        summary={"accession": accession, "commands": commands, "demo": bool(args.demo)},
        data={"sections": dict(sections)},
        status="ok",
        ok=True,
    ))

    af.write_bundle(output_dir, script=Path(__file__).resolve(),
                    env_name="clawbio-biostudies-fetch",
                    command=commands[0], written=written)
    print(f"Wrote {len(written)} artifact(s) to {output_dir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
