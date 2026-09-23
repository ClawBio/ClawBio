#!/usr/bin/env python3
"""pride-fetch — ClawBio entry point for the PRIDE Archive.

Ported from UKDRI/informatics_data_skills @ 7cc3e6e.
Copyright (c) 2026 UK Dementia Research Institute. Licensed MIT.

The archive logic lives in `pride_fetch_api.py`, vendored close to upstream.
The shared ClawBio machinery lives in `clawbio.common.archive_fetch`. This
module holds only what is specific to PRIDE.

Note PRIDE keeps upstream's own `download-script`: it is driven by the project
file list and can unzip, unlike the samplesheet-driven emitter the nucleotide
archives share (`clawbio.common.download_script`).

    python skills/pride-fetch/pride_fetch.py --demo --output /tmp/pride
    python skills/pride-fetch/pride_fetch.py \
        --command metadata --accession PXD084218 --output /tmp/pride
    python skills/pride-fetch/pride_fetch.py metadata PXD084218 --output /tmp/pride
"""

from __future__ import annotations

import json
import sys
import urllib.parse
from pathlib import Path

_SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SKILL_DIR.parent.parent
for _p in (str(_PROJECT_ROOT), str(_SKILL_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from clawbio.common import archive_fetch as af  # noqa: E402
from clawbio.common.report import write_result_json  # noqa: E402

api = af.load_sibling(_SKILL_DIR, "pride_fetch_api")

SKILL = "pride-fetch"
VERSION = "0.1.0"
COMMANDS = ("metadata", "files", "download", "search", "samplesheet",
            "metadata-table", "download-script")
# Search hits shown in report.md. Explicit here rather than inherited from
# the shared parser, which has no default -- see archive_fetch.common_parser.
SEARCH_LIMIT = 20

DEMO_ACCESSION = "PXD084218"
_E = _SKILL_DIR / "examples"
DEMO_PROJECT = _E / f"demo_{DEMO_ACCESSION}_project.json"
DEMO_FILES = _E / f"demo_{DEMO_ACCESSION}_files.json"
DEMO_SDRF_LINKS = _E / f"demo_{DEMO_ACCESSION}_sdrf_links.json"
DEMO_SEARCH = _E / "demo_search.json"
DEMO_COMMANDS = ("metadata", "files", "metadata-table", "samplesheet", "download-script")


def _build_parser():
    p = af.common_parser("pride_fetch.py", __doc__, COMMANDS)
    p.add_argument("--ext", help="filter by extension, e.g. raw, mzid, mzML, mgf")
    p.add_argument("--from", dest="source", choices=["auto", "pride", "generate"],
                   default="auto", help="SDRF source for --command samplesheet")
    p.add_argument("--acquisition", choices=["dia", "dda"], default="dia")
    p.add_argument("--local-dir")
    p.add_argument("--tool", choices=["curl", "wget"], default="curl",
                   help="download-script transfer tool; curl is the default "
                        "because it also retries HTTP 403 and aborts a stalled "
                        "transfer. Both resume a partial file.")
    p.add_argument("--outdir", default="pride_data")
    p.add_argument("--unzip", action="store_true")
    p.add_argument("--no-slurm", action="store_true")
    p.add_argument("--job-name", default="pride_download")
    p.add_argument("--partition")
    p.add_argument("--account")
    p.add_argument("--cpus", default="1")
    p.add_argument("--mem", default="4G")
    p.add_argument("--time", default="24:00:00")
    p.add_argument("--email")
    return p


def _to_upstream_argv(args, output_dir: Path) -> list[str]:
    cmd = args.command
    if cmd == "search":
        if not args.query:
            raise SystemExit("--command search needs --query")
        return ["search", args.query, "--limit", str(SEARCH_LIMIT if args.limit is None else args.limit)] + (
            ["--json"] if args.json else [])

    if not args.accession:
        raise SystemExit(f"--command {cmd} needs --accession")

    if cmd == "metadata":
        return ["metadata", args.accession] + (["--json"] if args.json else [])
    if cmd == "files":
        argv = ["files", args.accession]
        if args.ext:
            argv += ["--ext", args.ext]
        return argv + (["--json"] if args.json else [])
    if cmd == "download":
        argv = ["download", args.accession,
                "--out", str(af.resolve_out(args.out, output_dir, "downloads"))]
        if args.ext:
            argv += ["--ext", args.ext]
        return argv
    if cmd == "metadata-table":
        return ["metadata-table", args.accession,
                "--out", str(af.resolve_out(args.out, output_dir, "tables/metadata.tsv"))]
    if cmd == "samplesheet":
        return ["samplesheet", args.accession, "--from", args.source,
                "--acquisition", args.acquisition,
                "--out", str(af.resolve_out(
                    args.out, output_dir, f"{args.accession}.sdrf.tsv"))] + (
            ["--local-dir", args.local_dir] if args.local_dir else [])
    if cmd == "download-script":
        argv = ["download-script", args.accession, "--tool", args.tool,
                "--outdir", args.outdir, "--job-name", args.job_name,
                "--cpus", args.cpus, "--mem", args.mem, "--time", args.time,
                "--out", str(af.resolve_out(args.out, output_dir, "download_pride.sh"))]
        if args.ext:
            argv += ["--ext", args.ext]
        for flag, value in (("--partition", args.partition),
                            ("--account", args.account), ("--email", args.email)):
            if value:
                argv += [flag, value]
        if args.unzip:
            argv.append("--unzip")
        if args.no_slurm:
            argv.append("--no-slurm")
        return argv
    raise SystemExit(f"unsupported command: {cmd}")


def _install_demo_transport() -> None:
    """Serve the vendored client's single HTTP entry point from fixtures."""
    project = DEMO_PROJECT.read_bytes()
    files = DEMO_FILES.read_bytes()
    sdrf_links = DEMO_SDRF_LINKS.read_bytes()
    search = DEMO_SEARCH.read_bytes()

    def _offline_get(url: str, retries: int = 3) -> bytes:
        path = urllib.parse.urlparse(url).path
        query = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)
        if "/search/projects" in path:
            return search
        if "/files/sdrf/" in path:
            # The demo project genuinely has no submitter SDRF, so `auto` falls
            # through to generating a minimal one from the data files. That
            # keeps the demo offline without pretending: the submitter SDRF,
            # when one exists, lives on ftp.pride.ebi.ac.uk.
            return sdrf_links
        if path.endswith("/files"):
            # Upstream pages until a short batch comes back.
            return files if int(query.get("page", ["0"])[0]) == 0 else b"[]"
        if "/projects/" in path:
            return project
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
        args.query = args.query or "Arabidopsis"

    accession = args.accession or DEMO_ACCESSION
    commands = list(DEMO_COMMANDS) if (args.demo and not args.command) else [args.command]

    sections: list[tuple[str, str]] = []
    for cmd in commands:
        args.command = cmd
        stdout, stderr = af.run_upstream(api, _to_upstream_argv(args, output_dir), output_dir)
        sections.append((cmd, stdout or stderr))

    report = af.write_report(
        output_dir,
        title=f"PRIDE report — {accession}",
        source_url=f"[PRIDE Archive](https://www.ebi.ac.uk/pride/archive/projects/{accession})",
        sections=sections,
    )
    written = [report]
    for rel in ("tables/metadata.tsv", f"{accession}.sdrf.tsv", "download_pride.sh"):
        if (output_dir / rel).exists():
            written.append(output_dir / rel)

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
                    env_name="clawbio-pride-fetch", command=commands[0], written=written)
    print(f"Wrote {len(written)} artifact(s) to {output_dir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
