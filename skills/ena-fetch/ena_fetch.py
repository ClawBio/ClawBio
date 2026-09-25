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
# Search hits shown in report.md. Explicit here rather than inherited from
# the shared parser, which has no default -- see archive_fetch.common_parser.
SEARCH_LIMIT = 20
# Handled by this wrapper, not the vendored CLI (folded in from upstream's
# standalone fastq-download-script skill). tests/test_archive_command_coverage.py
# uses this to tell "implemented here" apart from "implemented nowhere".
LOCAL_COMMANDS = ("download-script",)

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


def _to_upstream_argv(args, output_dir: Path) -> list[str]:
    cmd = args.command
    if cmd == "fields":
        return ["fields", "--result", args.result]
    if cmd == "search":
        argv = ["search", "--result", args.result,
                "--limit", str(SEARCH_LIMIT if args.limit is None else args.limit)]
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
        # No --limit unless the caller gave one: the vendored default is
        # 0 = no limit, and forwarding a wrapper default silently truncated a
        # 95-run report to 20 rows. `is not None` rather than truthiness --
        # --limit 0 is a deliberate request for everything, not an unset flag.
        argv = ["report", args.accession, "--result", args.result]
        if args.limit is not None:
            argv += ["--limit", str(args.limit)]
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
            "Nothing has been downloaded. Run it yourself with "
            f"`bash {path.name}`, or submit it with `sbatch {path.name}`.")


def _announce_download_size(args, output_dir: Path) -> None:
    """Print the transfer size before `download` blocks.

    Best effort: a failed size probe must never stop the download it exists to
    explain, so anything going wrong here is swallowed.
    """
    field = "submitted_bytes" if getattr(args, "submitted", False) else "fastq_bytes"
    try:
        text = api.portal_get(
            "filereport",
            {"accession": args.accession, "result": args.result,
             "fields": f"run_accession,{field}", "format": "tsv", "limit": "0"},
        )
        rows = api.parse_tsv(text)
    except Exception:
        return
    if rows:
        print(_download_preflight_line(args.accession, rows, field), file=sys.stderr)


def _truncation_warning(body: str, limit: int | None) -> str | None:
    """Flag a report whose row count exactly equals --limit.

    Correcting the default stops the accidental truncation; an explicit
    `--limit 50` on a 95-run study still returns 50 rows and still looks
    complete. Equality is the only signal available without a second query, so
    the wording says "probably" -- a study with exactly 50 runs is a false
    positive, which is the right trade against shipping a short table silently.
    """
    if not limit:  # None (unset) or 0 (explicitly unlimited)
        return None
    rows = [ln for ln in body.splitlines() if ln.strip()]
    if len(rows) - 1 != limit:  # minus the header
        return None
    return (f"[warning] report returned exactly {limit} row(s), the value of "
            f"--limit \u2014 the result is probably truncated. Re-run with "
            f"--limit 0 for the complete report.")


def _human_bytes(n: float) -> str:
    """Bytes as a size a human can judge a wait against."""
    for unit in ("B", "KB", "MB", "GB"):
        if abs(n) < 1024.0:
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024.0
    return f"{n:.1f} TB"


def _download_preflight_line(accession: str, rows: list, field: str) -> str:
    """Announce size before a blocking transfer.

    Printed from the wrapper, not the vendored command, because `run_upstream`
    redirects stdout and stderr into a StringIO: anything the download prints
    surfaces only after it returns, all at once, after the wait it was meant to
    explain. `warn_if_overwriting` reaches the terminal for the same reason --
    it runs outside the redirect.
    """
    total, count = 0, 0
    for row in rows:
        for value in filter(None, (row.get(field) or "").split(";")):
            count += 1
            try:
                total += int(value)
            except ValueError:  # ENA omits sizes on some records
                pass
    size = _human_bytes(total) if total else "size unknown"
    return (f"[download] {accession}: {count} file(s), ~{size} from "
            f"ftp.sra.ebi.ac.uk.\n"
            f"           One blocking transfer with no progress output; a slow "
            f"link looks the same as a hang. Ctrl-C is safe.")


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
        if cmd == "download":
            _announce_download_size(args, output_dir)
        stdout, stderr = af.run_upstream(api, _to_upstream_argv(args, output_dir), output_dir)
        body = stdout or stderr
        if cmd == "report":
            warning = _truncation_warning(body, args.limit)
            if warning:
                # Both channels on purpose: stderr so the caller sees it now,
                # and the section so report.md carries its own caveat rather
                # than relying on someone having watched the terminal.
                print(warning, file=sys.stderr)
                body = f"{warning}\n\n{body}"
        sections.append((cmd, body))

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
