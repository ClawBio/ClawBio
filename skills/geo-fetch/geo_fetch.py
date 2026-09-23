#!/usr/bin/env python3
"""geo-fetch — ClawBio entry point for the NCBI Gene Expression Omnibus.

Ported from UKDRI/informatics_data_skills @ 7cc3e6e.
Copyright (c) 2026 UK Dementia Research Institute. Licensed MIT.

The archive logic lives in `geo_fetch_api.py`, vendored close to upstream. The
shared ClawBio machinery lives in `clawbio.common.archive_fetch`. This module
holds only what is specific to GEO.

Credentials: `NCBI_EMAIL` and `NCBI_API_KEY` are read from the environment but
are **never sent** unless `--use-ncbi-credentials` is passed. The presence of
an environment variable is not consent; that rule is carried over from upstream
deliberately.

    python skills/geo-fetch/geo_fetch.py --demo --output /tmp/geo
    python skills/geo-fetch/geo_fetch.py \
        --command samplesheet --accession GSE30720 --assay bulk --output /tmp/geo
    python skills/geo-fetch/geo_fetch.py metadata GSE30720 --output /tmp/geo
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

_SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SKILL_DIR.parent.parent
for _p in (str(_PROJECT_ROOT), str(_SKILL_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from clawbio.common import archive_fetch as af  # noqa: E402
from clawbio.common.download_script import (  # noqa: E402
    SlurmOptions,
    urls_from_samplesheet,
    write_download_script,
)
from clawbio.common.report import write_result_json  # noqa: E402

api = af.load_sibling(_SKILL_DIR, "geo_fetch_api")

SKILL = "geo-fetch"
VERSION = "0.1.0"
COMMANDS = ("metadata", "samples", "files", "download", "search",
            "metadata-table", "runtable", "samplesheet", "download-script")
# Search hits shown in report.md. Explicit here rather than inherited from
# the shared parser, which has no default -- see archive_fetch.common_parser.
SEARCH_LIMIT = 20
# Handled by this wrapper, not the vendored CLI (folded in from upstream's
# standalone fastq-download-script skill). tests/test_archive_command_coverage.py
# uses this to tell "implemented here" apart from "implemented nowhere".
LOCAL_COMMANDS = ("download-script",)

DEMO_ACCESSION = "GSE30720"
DEMO_HTTP = _SKILL_DIR / "examples" / f"demo_{DEMO_ACCESSION}_http.json.gz"
DEMO_COMMANDS = ("metadata", "samples", "files", "metadata-table",
                 "samplesheet", "download-script")


def _build_parser():
    p = af.common_parser("geo_fetch.py", __doc__, COMMANDS)
    p.add_argument("--use-ncbi-credentials", action="store_true",
                   help="opt in to sending NCBI_EMAIL / NCBI_API_KEY to E-utilities")
    p.add_argument("--organism", help='search filter, e.g. "Arabidopsis thaliana"')
    p.add_argument("--type", dest="entry_type", help="entry type filter, e.g. gse, gds")
    p.add_argument("--matrix", action="store_true", help="download the series matrix")
    p.add_argument("--soft", action="store_true", help="download the SOFT family file")
    p.add_argument("--miniml", action="store_true", help="download the MINiML file")
    p.add_argument("--suppl", action="store_true", help="download supplementary files")
    p.add_argument("--assay", choices=["scrna", "bulk"],
                   help="samplesheet flavour: nf-core/scrnaseq or nf-core/rnaseq")
    p.add_argument("--strandedness",
                   choices=["auto", "forward", "reverse", "unstranded"], default="auto")
    p.add_argument("--group-by", default="sample_accession")
    p.add_argument("--local-dir")
    p.add_argument("--from-runtable", metavar="PATH")
    p.add_argument("--fastq-dir", metavar="DIR")
    p.add_argument("--fastq-naming", choices=["sra", "cellranger"])
    p.add_argument("--read-map", metavar="R1,R2")
    p.add_argument("--tool", choices=["curl", "wget"], default="curl",
                   help="download-script transfer tool; curl is the default "
                        "because it also retries HTTP 403 and aborts a stalled "
                        "transfer. Both resume a partial file.")
    p.add_argument("--no-slurm", action="store_true")
    p.add_argument("--partition")
    p.add_argument("--account")
    p.add_argument("--job-name", default="fastq_download")
    p.add_argument("--cpus", default="1")
    p.add_argument("--mem", default="4G")
    p.add_argument("--time", default="24:00:00")
    p.add_argument("--email")
    return p


def _cred(args) -> list[str]:
    return ["--use-ncbi-credentials"] if args.use_ncbi_credentials else []


def _to_upstream_argv(args, output_dir: Path) -> list[str]:
    cmd = args.command
    if cmd == "search":
        if not args.query:
            raise SystemExit("--command search needs --query")
        argv = ["search", args.query, "--limit", str(SEARCH_LIMIT if args.limit is None else args.limit)]
        if args.organism:
            argv += ["--organism", args.organism]
        if args.entry_type:
            argv += ["--type", args.entry_type]
        return argv + (["--json"] if args.json else []) + _cred(args)

    if not args.accession:
        raise SystemExit(f"--command {cmd} needs --accession")

    if cmd in ("metadata", "samples", "files"):
        return [cmd, args.accession] + (["--json"] if args.json else []) + _cred(args)
    if cmd == "download":
        argv = ["download", args.accession,
                "--out", str(af.resolve_out(args.out, output_dir, "downloads"))]
        for flag, on in (("--matrix", args.matrix), ("--soft", args.soft),
                         ("--miniml", args.miniml), ("--suppl", args.suppl)):
            if on:
                argv.append(flag)
        return argv + _cred(args)
    if cmd == "metadata-table":
        return ["metadata-table", args.accession,
                "--out", str(af.resolve_out(args.out, output_dir, "tables/metadata.tsv"))] + _cred(args)
    if cmd == "runtable":
        target = af.resolve_out(args.out, output_dir, "runtable")
        target.mkdir(parents=True, exist_ok=True)
        return ["runtable", args.accession, "--out", str(target)] + _cred(args)
    if cmd == "samplesheet":
        if not args.assay:
            raise SystemExit("--command samplesheet needs --assay scrna|bulk")
        argv = ["samplesheet", args.accession, "--assay", args.assay,
                "--strandedness", args.strandedness, "--group-by", args.group_by,
                "--out", str(af.resolve_out(args.out, output_dir, "samplesheet.csv"))]
        for flag, value in (("--local-dir", args.local_dir),
                            ("--from-runtable", args.from_runtable),
                            ("--fastq-dir", args.fastq_dir),
                            ("--fastq-naming", args.fastq_naming),
                            ("--read-map", args.read_map)):
            if value:
                argv += [flag, value]
        return argv + _cred(args)
    raise SystemExit(f"unsupported command: {cmd}")


def _run_download_script(args, output_dir: Path) -> str:
    """Emit a runnable download script from the samplesheet just written."""
    sheet = af.resolve_out(None, output_dir, "samplesheet.csv")
    if not sheet.exists():
        raise SystemExit(
            "download-script needs a samplesheet; run --command samplesheet first")
    slurm = None if args.no_slurm else SlurmOptions(
        job_name=args.job_name, partition=args.partition, account=args.account,
        cpus=args.cpus, mem=args.mem, time=args.time, email=args.email)
    path, n = write_download_script(
        urls_from_samplesheet(sheet), output_dir / "download_geo.sh",
        tool=args.tool, outdir="fastq", slurm=slurm)
    return (f"Wrote {path.name}: {n} download command(s) using {args.tool}.\n"
            "Nothing has been downloaded. Run it with --run, submit it with "
            "--submit, or execute it yourself.")


def _install_demo_transport() -> None:
    """Replay recorded responses by URL.

    GEO fans out across E-utilities, the GEO FTP listing, per-sample SOFT
    records and the ENA Portal, so the fixture is a URL->body map captured from
    one real run rather than a single payload. That keeps `--demo` on the real
    parsing path for all six commands.
    """
    recorded = json.loads(gzip.open(DEMO_HTTP, "rt", encoding="utf-8").read())

    def _offline_get(url: str, retries: int = 3) -> bytes:
        if url in recorded:
            return recorded[url].encode()
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
        args.assay = args.assay or "bulk"
        args.query = args.query or "Arabidopsis seedling transcriptome"

    accession = args.accession or DEMO_ACCESSION
    commands = list(DEMO_COMMANDS) if (args.demo and not args.command) else [args.command]

    sections: list[tuple[str, str]] = []
    for cmd in commands:
        args.command = cmd
        if cmd == "download-script":
            sections.append((cmd, _run_download_script(args, output_dir)))
            continue
        stdout, stderr = af.run_upstream(api, _to_upstream_argv(args, output_dir), output_dir)
        sections.append((cmd, stdout or stderr))

    report = af.write_report(
        output_dir,
        title=f"GEO report — {accession}",
        source_url=f"[Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession})",
        sections=sections,
    )
    written = [report]
    for rel in ("tables/metadata.tsv", "samplesheet.csv", "download_geo.sh"):
        if (output_dir / rel).exists():
            written.append(output_dir / rel)

    written.append(write_result_json(
        output_dir,
        skill=SKILL,
        version=VERSION,
        summary={"accession": accession, "commands": commands,
                 "demo": bool(args.demo),
                 "ncbi_credentials_sent": bool(args.use_ncbi_credentials)},
        data={"sections": dict(sections)},
        status="ok",
        ok=True,
    ))

    af.write_bundle(output_dir, script=Path(__file__).resolve(),
                    env_name="clawbio-geo-fetch", command=commands[0], written=written)
    print(f"Wrote {len(written)} artifact(s) to {output_dir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
