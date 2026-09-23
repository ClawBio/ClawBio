#!/usr/bin/env python3
"""arrayexpress-fetch — ClawBio entry point for EMBL-EBI ArrayExpress.

Ported from UKDRI/informatics_data_skills @ 7cc3e6e.
Copyright (c) 2026 UK Dementia Research Institute. Licensed MIT.

The archive logic lives in `arrayexpress_fetch_api.py`, vendored close to
upstream. The shared ClawBio machinery lives in `clawbio.common.archive_fetch`.
This module holds only what is specific to ArrayExpress: its subcommands, how
they map onto the upstream parser, and how `--demo` is served offline.

    python skills/arrayexpress-fetch/arrayexpress_fetch.py --demo --output /tmp/ae
    python skills/arrayexpress-fetch/arrayexpress_fetch.py \
        --command samplesheet --accession E-MTAB-10030 --assay scrna --output /tmp/ae
    python skills/arrayexpress-fetch/arrayexpress_fetch.py metadata E-MTAB-10030 --output /tmp/ae
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

api = af.load_sibling(_SKILL_DIR, "arrayexpress_fetch_api")

SKILL = "arrayexpress-fetch"
VERSION = "0.1.0"
COMMANDS = ("metadata", "files", "sdrf", "download", "search",
            "metadata-table", "samplesheet", "download-script")
# Handled here rather than by the vendored CLI. `tests/test_archive_command_coverage.py`
# uses this to tell "implemented in the wrapper" apart from "implemented nowhere".
LOCAL_COMMANDS = ("download-script",)

DEMO_ACCESSION = "E-MTAB-10030"
DEMO_STUDY_FIXTURE = _SKILL_DIR / "examples" / f"demo_{DEMO_ACCESSION}.json"
DEMO_SDRF_FIXTURE = _SKILL_DIR / "examples" / f"demo_{DEMO_ACCESSION}.sdrf.txt"
DEMO_SEARCH_FIXTURE = _SKILL_DIR / "examples" / "demo_search.json"
# A demo run must write the whole non-optional Output Structure tree, so it
# exercises the SDRF and samplesheet paths rather than the cheapest command.
DEMO_COMMANDS = ("metadata", "files", "sdrf", "metadata-table", "samplesheet", "search")


def _build_parser():
    p = af.common_parser("arrayexpress_fetch.py", __doc__, COMMANDS)
    p.add_argument("--match", help="only files whose path contains this substring")
    p.add_argument("--magetab", action="store_true", help="download IDF + SDRF")
    p.add_argument("--processed", action="store_true", help="download processed matrices")
    p.add_argument("--raw", action="store_true", help="download raw data (FASTQ/CEL/BAM)")
    p.add_argument("--assay", choices=["scrna", "bulk"],
                   help="samplesheet flavour: scrna -> nf-core/scrnaseq, bulk -> nf-core/rnaseq")
    p.add_argument("--strandedness",
                   choices=["auto", "forward", "reverse", "unstranded"], default="auto",
                   help="strandedness column value (--assay bulk only)")
    p.add_argument("--local-dir", help="write local paths instead of FASTQ URLs")
    p.add_argument("--fastq-dir", metavar="DIR",
                   help="write paths into a flat fasterq-dump output directory")
    p.add_argument("--fastq-naming", choices=["sra", "cellranger"],
                   help="with --fastq-dir: how the FASTQ files are named")
    p.add_argument("--read-map", metavar="R1,R2",
                   help="declare which reads are the cDNA pair, 1-based (e.g. 3,4). "
                        "Required for 10x runs whose technical reads are separate files")
    # No --no-slurm/--run/--submit: this skill emits no download script. See
    # _refer_fastq_downloads_to_ena.
    return p


def _refer_fastq_downloads_to_ena(accession: str, output_dir: Path) -> None:
    """ArrayExpress brokers sequencing reads to ENA, so there is no FASTQ
    script to emit here.

    This exits rather than writing anything. It used to fall through to the
    vendored CLI, which rejected the unknown subcommand inside a redirected
    stderr -- exit 2, no message, and a stale report.md left on disk.

    `--command download` is unaffected: it still fetches whatever ArrayExpress
    itself hosts (IDF, SDRF, processed matrices, CEL, BAM).
    """
    raise SystemExit(
        "arrayexpress-fetch does not emit a FASTQ download script: ArrayExpress\n"
        "brokers sequencing reads to ENA and serves the bytes from there.\n"
        "\n"
        "Build the samplesheet here, then emit the script with ena-fetch against\n"
        f"this same output directory ({output_dir}):\n"
        "\n"
        f"    python skills/arrayexpress-fetch/arrayexpress_fetch.py \\\n"
        f"        --command samplesheet --accession {accession} --assay bulk \\\n"
        f"        --output {output_dir}\n"
        f"    python skills/ena-fetch/ena_fetch.py \\\n"
        f"        --command download-script --output {output_dir}\n"
        "\n"
        "ena-fetch reads samplesheet.csv and fetches whatever URLs it names, so\n"
        "it works whether the SDRF points at ftp.sra.ebi.ac.uk or at\n"
        "ftp.ebi.ac.uk. For runs not mirrored to ENA, or for 10x reads whose\n"
        "structure only fasterq-dump exposes reliably, use sra-tools directly:\n"
        "prefetch --option-file SRR_Acc_List.txt, then fasterq-dump.\n"
        "\n"
        "Use --command download for files ArrayExpress does host.")


def _to_upstream_argv(args, output_dir: Path) -> list[str]:
    """Translate the --command form into the vendored parser's argv."""
    cmd = args.command
    if cmd == "search":
        if not args.query:
            raise SystemExit("--command search needs --query")
        argv = ["search", args.query, "--limit", str(args.limit)]
        if args.json:
            argv.append("--json")
        return argv

    if not args.accession:
        raise SystemExit(f"--command {cmd} needs --accession")
    argv = [cmd, args.accession]

    if cmd in ("metadata", "files") and args.json:
        argv.append("--json")
    elif cmd == "download":
        for flag in ("magetab", "processed", "raw"):
            if getattr(args, flag):
                argv.append(f"--{flag}")
        if args.match:
            argv += ["--match", args.match]
        argv += ["--out", str(af.resolve_out(args.out, output_dir, "downloads"))]
    elif cmd == "metadata-table":
        argv += ["--out", str(af.resolve_out(args.out, output_dir, "tables/metadata.tsv"))]
    elif cmd == "samplesheet":
        if not args.assay:
            raise SystemExit("--command samplesheet needs --assay scrna|bulk")
        argv += ["--assay", args.assay,
                 "--out", str(af.resolve_out(args.out, output_dir, "samplesheet.csv"))]
        if args.assay == "bulk":
            argv += ["--strandedness", args.strandedness]
        for flag, value in (("--local-dir", args.local_dir),
                            ("--fastq-dir", args.fastq_dir),
                            ("--fastq-naming", args.fastq_naming),
                            ("--read-map", args.read_map)):
            if value:
                argv += [flag, value]
    return argv


def _install_demo_transport() -> None:
    """Serve the vendored module's HTTP layer from the committed fixtures.

    Patching `http_get` rather than the command functions keeps `--demo` on the
    same code path as a live run, so the demo exercises the real parsing.
    """
    study = DEMO_STUDY_FIXTURE.read_bytes()
    search = DEMO_SEARCH_FIXTURE.read_bytes()
    sdrf = DEMO_SDRF_FIXTURE.read_bytes()

    def _offline_get(url: str, retries: int = 3) -> bytes:
        if "/search" in url:
            return search
        if url.endswith(".sdrf.txt") or "sdrf" in url.rsplit("/", 1)[-1].lower():
            return sdrf
        if "/studies/" in url:
            return study
        if "/biosamples/" in url:
            return b"{}"
        if "ena/portal" in url:
            # The fixture's SDRF carries Comment[FASTQ_URI], so no ENA lookup
            # is needed; answer empty rather than reaching the network.
            return b""
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
        args.query = args.query or "microglia"
        args.assay = args.assay or "scrna"

    accession = args.accession or DEMO_ACCESSION
    commands = list(DEMO_COMMANDS) if args.demo else [args.command]

    sections: list[tuple[str, str]] = []
    for cmd in commands:
        args.command = cmd
        if cmd == "download-script":
            _refer_fastq_downloads_to_ena(accession, output_dir)
        stdout, stderr = af.run_upstream(api, _to_upstream_argv(args, output_dir), output_dir)
        sections.append((cmd, stdout or stderr))

    report = af.write_report(
        output_dir,
        title=f"ArrayExpress report — {accession}",
        source_url=f"[EMBL-EBI ArrayExpress](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{accession})",
        sections=sections,
    )
    written = [report]

    for rel in ("tables/metadata.tsv", "samplesheet.csv"):
        candidate = output_dir / rel
        if candidate.exists():
            written.append(candidate)

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
                    env_name="clawbio-arrayexpress-fetch",
                    command=commands[0], written=written)
    print(f"Wrote {len(written)} artifact(s) to {output_dir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
