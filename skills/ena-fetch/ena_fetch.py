#!/usr/bin/env python3
"""ena-fetch — ClawBio entry point for the European Nucleotide Archive.

Ported from UKDRI/informatics_data_skills @ 7cc3e6e.
Copyright (c) 2026 UK Dementia Research Institute. Licensed MIT.

The archive logic lives in `ena_fetch_api.py`, vendored close to upstream. The
shared ClawBio machinery lives in `clawbio.common.archive_fetch`. This module
holds only what is specific to ENA: its subcommands, how they map onto the
upstream parser, and how `--demo` is served offline.

    python skills/ena-fetch/ena_fetch.py --demo --output /tmp/ena
    python skills/ena-fetch/ena_fetch.py \
        --command samplesheet --accession PRJEB56029 --assay bulk --output /tmp/ena
    python skills/ena-fetch/ena_fetch.py runs PRJEB56029 --output /tmp/ena
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
from clawbio.common.download_script import (  # noqa: E402
    SlurmOptions,
    urls_from_samplesheet,
    write_download_script,
)
from clawbio.common.report import write_result_json  # noqa: E402

api = af.load_sibling(_SKILL_DIR, "ena_fetch_api")

SKILL = "ena-fetch"
VERSION = "0.1.0"
COMMANDS = ("runs", "report", "fields", "search", "xml", "download",
            "metadata-table", "samplesheet", "download-script")

DEMO_ACCESSION = "PRJEB56029"
DEMO_FILEREPORT = _SKILL_DIR / "examples" / f"demo_{DEMO_ACCESSION}_filereport.tsv"
DEMO_SAMPLE_XML = _SKILL_DIR / "examples" / "demo_sample_xml.json"
DEMO_COMMANDS = ("runs", "metadata-table", "samplesheet", "download-script")


def _build_parser():
    p = af.common_parser("ena_fetch.py", __doc__, COMMANDS)
    p.add_argument("--result", default="read_run", help="ENA result type")
    p.add_argument("--fields", help="comma-separated field list")
    p.add_argument("--format", default="xml", help="xml, json, embl, fasta, text")
    p.add_argument("--submitted", action="store_true",
                   help="download submitted files rather than FASTQ")
    p.add_argument("--assay", choices=["scrna", "bulk"],
                   help="samplesheet flavour: nf-core/scrnaseq or nf-core/rnaseq")
    p.add_argument("--strandedness",
                   choices=["auto", "forward", "reverse", "unstranded"], default="auto")
    p.add_argument("--group-by", default="sample_accession")
    p.add_argument("--local-dir")
    p.add_argument("--fastq-dir", metavar="DIR")
    p.add_argument("--fastq-naming", choices=["sra", "cellranger"])
    p.add_argument("--read-map", metavar="R1,R2")
    # download-script (folded in from upstream's fastq-download-script skill)
    p.add_argument("--tool", choices=["wget", "curl"], default="wget")
    p.add_argument("--no-slurm", action="store_true")
    p.add_argument("--partition")
    p.add_argument("--account")
    p.add_argument("--job-name", default="fastq_download")
    p.add_argument("--cpus", default="1")
    p.add_argument("--mem", default="4G")
    p.add_argument("--time", default="24:00:00")
    p.add_argument("--email")
    return p


def _to_upstream_argv(args, output_dir: Path) -> list[str]:
    cmd = args.command
    if cmd == "fields":
        return ["fields", "--result", args.result]
    if cmd == "search":
        argv = ["search", "--result", args.result, "--limit", str(args.limit)]
        if args.query:
            argv += ["--query", args.query]
        if args.fields:
            argv += ["--fields", args.fields]
        if args.json:
            argv.append("--json")
        return argv

    if not args.accession:
        raise SystemExit(f"--command {cmd} needs --accession")

    if cmd == "runs":
        return ["runs", args.accession] + (["--json"] if args.json else [])
    if cmd == "xml":
        return ["xml", args.accession, "--format", args.format]
    if cmd == "report":
        argv = ["report", args.accession, "--result", args.result,
                "--limit", str(args.limit)]
        if args.fields:
            argv += ["--fields", args.fields]
        return argv + (["--json"] if args.json else [])
    if cmd == "download":
        argv = ["download", args.accession, "--result", args.result,
                "--out", str(af.resolve_out(args.out, output_dir, "downloads"))]
        return argv + (["--submitted"] if args.submitted else [])
    if cmd == "metadata-table":
        return ["metadata-table", args.accession,
                "--out", str(af.resolve_out(args.out, output_dir, "tables/metadata.tsv"))]
    if cmd == "samplesheet":
        if not args.assay:
            raise SystemExit("--command samplesheet needs --assay scrna|bulk")
        argv = ["samplesheet", args.accession, "--assay", args.assay,
                "--strandedness", args.strandedness, "--group-by", args.group_by,
                "--out", str(af.resolve_out(args.out, output_dir, "samplesheet.csv"))]
        for flag, value in (("--local-dir", args.local_dir),
                            ("--fastq-dir", args.fastq_dir),
                            ("--fastq-naming", args.fastq_naming),
                            ("--read-map", args.read_map)):
            if value:
                argv += [flag, value]
        return argv
    raise SystemExit(f"unsupported command: {cmd}")


def _run_download_script(args, output_dir: Path) -> str:
    """Emit a runnable download script from the samplesheet just written.

    Folded in from upstream's standalone `fastq-download-script` skill: the
    input is always another command's output, so it belongs here as a command
    rather than as a separate skill.
    """
    sheet = af.resolve_out(None, output_dir, "samplesheet.csv")
    if not sheet.exists():
        raise SystemExit(
            "download-script needs a samplesheet; run --command samplesheet first")
    slurm = None if args.no_slurm else SlurmOptions(
        job_name=args.job_name, partition=args.partition, account=args.account,
        cpus=args.cpus, mem=args.mem, time=args.time, email=args.email)
    path, n = write_download_script(
        urls_from_samplesheet(sheet), output_dir / "download_ena.sh",
        tool=args.tool, outdir="fastq", slurm=slurm)
    return (f"Wrote {path.name}: {n} download command(s) using {args.tool}.\n"
            "Nothing has been downloaded. Run it with --run, submit it with "
            "--submit, or execute it yourself.")


def _install_demo_transport() -> None:
    """Serve the vendored client's single HTTP entry point from fixtures.

    Every Portal call differs only by its `fields` parameter, so one stored
    superset TSV is projected down to whatever was asked for. That keeps
    `--demo` on the real parsing path instead of stubbing the commands out.
    """
    report_text = DEMO_FILEREPORT.read_text()
    header = report_text.splitlines()[0].split("\t")
    rows = api.parse_tsv(report_text)
    sample_xml = json.loads(DEMO_SAMPLE_XML.read_text())

    def _project(fields: str) -> bytes:
        wanted = [f for f in fields.split(",") if f in header] or header
        lines = ["\t".join(wanted)]
        lines += ["\t".join(r.get(f, "") for f in wanted) for r in rows]
        return ("\n".join(lines) + "\n").encode()

    def _offline_get(url: str, retries: int = 3) -> bytes:
        parsed = urllib.parse.urlparse(url)
        query = urllib.parse.parse_qs(parsed.query)
        if "/filereport" in parsed.path or "/search" in parsed.path:
            return _project(query.get("fields", [""])[0])
        if "/returnFields" in parsed.path:
            return b"columnId\tdescription\n" + b"\n".join(
                f"{h}\tdemo field".encode() for h in header)
        if "/xml/" in parsed.path:
            return sample_xml.get(parsed.path.rsplit("/", 1)[-1], "<SAMPLE_SET/>").encode()
        raise SystemExit(f"demo mode has no fixture for {url}")

    api.http_get = _offline_get


def main(argv: list[str] | None = None) -> int:
    argv = af.expand_positional(
        list(sys.argv[1:] if argv is None else argv), COMMANDS,
        query_commands=("search", "fields"))
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
        title=f"ENA report — {accession}",
        source_url=f"[European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/view/{accession})",
        sections=sections,
    )
    written = [report]
    for rel in ("tables/metadata.tsv", "samplesheet.csv", "download_ena.sh"):
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
                    env_name="clawbio-ena-fetch", command=commands[0], written=written)
    print(f"Wrote {len(written)} artifact(s) to {output_dir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
