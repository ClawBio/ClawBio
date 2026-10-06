#!/usr/bin/env python3
"""Run the upstream shotgun-metagenomics Snakemake workflow as a ClawBio skill.

This file validates inputs, renders the pipeline config, launches snakemake
against a pinned checkout of mobashirrahman/metagenomics-workflow, and turns the
pipeline's `final/` tables into `report.md` and `result.json`.

It reimplements no part of the analysis. Every number in the report is read from
a table the pipeline wrote or from the effective config; an unmeasured value is
reported as unmeasured rather than as zero, and a failed run reports the failure
rather than a partial result.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

# The repo root holds the `clawbio` package, which the reproducibility writers
# live in. The documented direct invocation puts only this directory on sys.path.
_PROJECT_ROOT = SKILL_DIR.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from magpipe_config import RunOptions, build_config, validate_config, write_config
from magpipe_errors import ErrorCode, SkillError
from magpipe_outputs import parse_outputs
from magpipe_preflight import (
    check_conda,
    check_output_dir,
    check_snakemake,
    default_cache_dir,
    network_requirements,
    resolve_pipeline_source,
)
from magpipe_report import write_report_md, write_reproducibility, write_result_json
from magpipe_runner import build_snakemake_command, run_snakemake
from magpipe_samplesheet import build_samplesheet
from magpipe_schemas import (
    CLI_ALIAS,
    DEFAULT_TIMEOUT_HOURS,
    PINNED_REF,
    PRESETS,
    SKILL_NAME,
)

# The configuration a demo run uses: upstream's, unchanged but for three paths.
DEMO_BASE_CONFIG = "workflow/config/config_test.yaml"
MAIN_BASE_CONFIG = "workflow/config/config.yaml"
DEFAULT_BINNERS = "metabat,maxbin,vamb"

# Upstream's own generator arguments, from docs/usage.md and test/run_test.sh.
DEMO_DATASET_ARGS = (
    "--samples", "2",
    "--genomes", "4",
    "--genome-length", "12000",
    "--pairs", "2000",
    "--seed", "1",
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="metagenome_mag_pipeline.py",
        description=(
            "Shotgun metagenomics QC, assembly, binning and MAG assessment. Wraps "
            "the pinned Snakemake workflow at "
            "github.com/mobashirrahman/metagenomics-workflow."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input",
        help="Samplesheet TSV, or a directory of FASTQ files to discover samples from",
    )
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Generate synthetic reads and run QC, assembly and reporting",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Validate the inputs and dry-run the workflow without executing it",
    )
    parser.add_argument(
        "--preset",
        choices=PRESETS,
        default="assembly",
        help="Which stages to enable (default: assembly)",
    )
    parser.add_argument("--cores", type=int, default=4, help="Cores for snakemake (default: 4)")
    parser.add_argument(
        "--mem-mb", type=int, default=8000, help="Total memory in MiB (default: 8000)"
    )
    parser.add_argument(
        "--coassemble", action="store_true", help="Pool reads by assembly group"
    )
    parser.add_argument(
        "--binners",
        default="metabat,maxbin,vamb",
        help="Comma-separated binning.binners (default: metabat,maxbin,vamb)",
    )
    parser.add_argument("--host-fasta", help="Decompressed host FASTA for read subtraction")
    parser.add_argument(
        "--skip-host-removal", action="store_true", help="Disable host read subtraction"
    )
    parser.add_argument("--checkm2-db", help="CheckM2 1.1.0 DIAMOND .dmnd database")
    parser.add_argument("--gtdbtk-db", help="Extracted GTDB-Tk reference directory")
    parser.add_argument(
        "--taxonomy-db", help="DIAMOND-indexed protein database whose subjects end in a taxid"
    )
    parser.add_argument("--taxdump", help="Unpacked NCBI taxdump directory")
    parser.add_argument("--kegg-db", help="DIAMOND KEGG-orthology database")
    parser.add_argument("--cog-db", help="DIAMOND COG database")
    parser.add_argument("--pfam-db", help="Pfam-A profile database")
    parser.add_argument("--pipeline-dir", help="Use an existing checkout of the workflow")
    parser.add_argument(
        "--conda-prefix",
        help="Where snakemake builds tool environments (default: <cache>/conda-envs)",
    )
    parser.add_argument(
        "--allow-remote-inputs",
        action="store_true",
        help="Permit samplesheet rows that name public SRA/ENA/DRA runs",
    )
    parser.add_argument(
        "--timeout-hours",
        type=float,
        default=DEFAULT_TIMEOUT_HOURS,
        help=f"Wall-clock cap in hours (default: {DEFAULT_TIMEOUT_HOURS:g})",
    )
    parser.add_argument(
        "--force", action="store_true", help="Write into a non-empty output directory"
    )
    return parser


def _options(args: argparse.Namespace) -> RunOptions:
    binners = tuple(part.strip() for part in args.binners.split(",") if part.strip())
    return RunOptions(
        preset=args.preset,
        cores=args.cores,
        mem_mb=args.mem_mb,
        coassemble=args.coassemble,
        host_fasta=Path(args.host_fasta) if args.host_fasta else None,
        skip_host_removal=args.skip_host_removal,
        allow_remote_inputs=args.allow_remote_inputs,
        binners=binners,
        checkm2_db=Path(args.checkm2_db) if args.checkm2_db else None,
        gtdbtk_db=Path(args.gtdbtk_db) if args.gtdbtk_db else None,
        taxonomy_db=Path(args.taxonomy_db) if args.taxonomy_db else None,
        taxdump=Path(args.taxdump) if args.taxdump else None,
        kegg_db=Path(args.kegg_db) if args.kegg_db else None,
        cog_db=Path(args.cog_db) if args.cog_db else None,
        pfam_db=Path(args.pfam_db) if args.pfam_db else None,
        demo=args.demo,
    )


def _demo_conflicts(args: argparse.Namespace) -> list[tuple[str, str]]:
    """Every flag this command gives the demo that the demo cannot honour.

    A demo run is upstream's own `workflow/config/config_test.yaml` with three
    paths moved: QC, assembly, QUAST, MultiQC and provenance, on synthetic reads
    that need no reference database. Anything that configures a stage it does not
    enable would be accepted and then silently not applied, leaving a user
    believing they ran a binned or annotated analysis when no such rule was
    enabled. The flags are refused instead, all of them named at once so the
    user edits their command once rather than one error per run.
    """
    conflicts: list[tuple[str, str]] = []
    if args.input:
        conflicts.append(("--input", str(args.input)))
    if args.preset != "assembly":
        conflicts.append(("--preset", args.preset))
    for flag, value in (
        ("--coassemble", args.coassemble),
        ("--host-fasta", args.host_fasta),
        ("--skip-host-removal", args.skip_host_removal),
        ("--checkm2-db", args.checkm2_db),
        ("--gtdbtk-db", args.gtdbtk_db),
        ("--taxonomy-db", args.taxonomy_db),
        ("--taxdump", args.taxdump),
        ("--kegg-db", args.kegg_db),
        ("--cog-db", args.cog_db),
        ("--pfam-db", args.pfam_db),
        ("--allow-remote-inputs", args.allow_remote_inputs),
        ("--binners", args.binners if args.binners != DEFAULT_BINNERS else None),
    ):
        if value:
            conflicts.append((flag, str(value) if isinstance(value, str) else ""))
    return conflicts


def _reject_demo_conflicts(args: argparse.Namespace) -> None:
    conflicts = _demo_conflicts(args)
    if not conflicts:
        return
    rendered = [f"{flag} {value}".strip() for flag, value in conflicts]
    raise SkillError(
        stage="config",
        error_code=ErrorCode.INVALID_CONFIG,
        message=(
            "--demo runs the pinned workflow's own test configuration "
            f"({DEMO_BASE_CONFIG}), which enables QC, assembly and reporting on "
            "synthetic reads, so it cannot honour: " + ", ".join(rendered)
        ),
        fix=(
            "Drop --demo to run those settings against your own reads with "
            "--input, or drop "
            + " ".join(rendered)
            + " to run the demo. --cores, --mem-mb, --timeout-hours, "
            "--pipeline-dir and --conda-prefix do apply to a demo run."
        ),
        details={
            # Bare flag names, so a caller can act on them without parsing prose.
            "conflicts": [flag for flag, _ in conflicts],
            "base_config": DEMO_BASE_CONFIG,
            "demo_config": "the demo enables assembly only; binning, MAG "
            "assessment, contig taxonomy, functional annotation and host "
            "decontamination are all off",
        },
    )


def _generate_demo_dataset(pipeline_dir: Path, data_dir: Path) -> Path:
    """Ask the pinned checkout's own generator for synthetic reads."""
    data_dir.mkdir(parents=True, exist_ok=True)
    command = [
        sys.executable,
        str(pipeline_dir / "test" / "make_test_dataset.py"),
        "--outdir",
        str(data_dir),
        "--dbdir",
        str(data_dir / "dbs"),
        *DEMO_DATASET_ARGS,
    ]
    result = subprocess.run(command, capture_output=True, text=True, timeout=900, check=False)
    if result.returncode != 0:
        raise SkillError(
            stage="samplesheet",
            error_code=ErrorCode.DEMO_REQUIRES_NETWORK,
            message="the upstream test data generator failed: "
            + (result.stderr.strip().splitlines() or ["(no stderr)"])[-1],
            fix="Check that Python can write to the output directory, then re-run "
            "--demo.",
            details={"command": command, "stderr": result.stderr[-4000:]},
        )
    sheet = data_dir / "samplesheet.tsv"
    if not sheet.is_file():
        raise SkillError(
            stage="samplesheet",
            error_code=ErrorCode.DEMO_REQUIRES_NETWORK,
            message=f"the upstream generator wrote no samplesheet at {sheet}",
            fix="Re-run --demo; if it persists, clone the workflow and pass "
            "--pipeline-dir.",
            details={"outdir": str(data_dir)},
        )
    return sheet


def _load_base_config(pipeline_dir: Path, *, demo: bool) -> dict:
    import yaml

    name = "config_test.yaml" if demo else "config.yaml"
    path = pipeline_dir / "workflow" / "config" / name
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as error:
        raise SkillError(
            stage="config",
            error_code=ErrorCode.INVALID_CONFIG,
            message=f"could not read the upstream configuration {path}: {error}",
            fix="Pass --pipeline-dir pointing at an intact checkout of the workflow.",
            details={"path": str(path)},
        ) from error
    if not isinstance(data, dict):
        raise SkillError(
            stage="config",
            error_code=ErrorCode.INVALID_CONFIG,
            message=f"{path} is not a configuration mapping",
            fix="Pass --pipeline-dir pointing at an intact checkout of the workflow.",
            details={"path": str(path)},
        )
    return data


def _reproducibility_files(output_dir: Path, results_dir: Path, config: dict) -> list[Path]:
    from magpipe_outputs import expected_outputs

    files = [
        output_dir / "report.md",
        output_dir / "result.json",
        output_dir / "pipeline" / "samplesheet.tsv",
        output_dir / "pipeline" / "config.yaml",
    ]
    for relative in expected_outputs(config):
        files.extend(sorted(results_dir.glob(relative)) or [results_dir / relative])
    return files


def run(args: argparse.Namespace) -> int:
    output_dir = Path(args.output).expanduser().resolve()
    pipeline_dir = output_dir / "pipeline"
    results_dir = output_dir / "results"
    warnings: list[str] = []

    check_output_dir(output_dir, force=args.force)
    pipeline_dir.mkdir(parents=True, exist_ok=True)

    # Before anything is fetched: a flag the demo cannot honour is the user's
    # own mistake to fix, and it should not cost them a clone to discover it.
    if args.demo:
        _reject_demo_conflicts(args)

    cache_dir = default_cache_dir()
    # The upstream workflow is public code: cloning it fetches nothing of the
    # user's, which is why this is unconditional rather than behind a flag.
    source = resolve_pipeline_source(
        Path(args.pipeline_dir) if args.pipeline_dir else None,
        cache_dir=cache_dir,
        ref=PINNED_REF,
        allow_fetch=True,
    )
    checkout = Path(source["path"])

    options = _options(args)
    if args.demo:
        generated = _generate_demo_dataset(checkout, pipeline_dir / "demo_data")
        samples = build_samplesheet(generated, pipeline_dir / "samplesheet.tsv")
    else:
        samples = build_samplesheet(Path(args.input), pipeline_dir / "samplesheet.tsv")

    base_config_name = DEMO_BASE_CONFIG if args.demo else MAIN_BASE_CONFIG
    config = build_config(
        _load_base_config(checkout, demo=args.demo),
        samples,
        options,
        samplesheet=pipeline_dir / "samplesheet.tsv",
        results_dir=results_dir,
        database_dir=output_dir / "db",
    )
    validate_config(config, samples)
    write_config(config, pipeline_dir / "config.yaml")

    check_snakemake()
    check_conda()
    conda_prefix = (
        Path(args.conda_prefix) if args.conda_prefix else cache_dir / "conda-envs"
    )
    requirements = network_requirements(config, samples)
    print(f"[{SKILL_NAME}] this run will use the network:", file=sys.stderr)
    for line in requirements:
        print(f"  - {line}", file=sys.stderr)

    argv, command = build_snakemake_command(
        pipeline_dir=checkout,
        config_path=pipeline_dir / "config.yaml",
        work_dir=pipeline_dir / "workdir",
        conda_prefix=conda_prefix,
        cores=args.cores,
        mem_mb=args.mem_mb,
        dry_run=args.check,
    )
    print(f"[{SKILL_NAME}] {' '.join(argv)}", file=sys.stderr)
    run_snakemake(
        argv,
        log_path=pipeline_dir / "snakemake.log",
        timeout_seconds=args.timeout_hours * 3600,
    )

    mode = "check" if args.check else ("demo" if args.demo else "run")
    status = "planned" if args.check else "ok"
    findings = None
    if not args.check:
        findings = parse_outputs(results_dir, config)

    report_source = dict(source)
    report_source["command"] = command
    # Which upstream config the effective config was derived from, so a report
    # reader can tell a demo run from a real one without guessing.
    report_source["base_config"] = base_config_name
    write_result_json(
        output_dir,
        status=status,
        mode=mode,
        findings=findings,
        config=config,
        pipeline_source=report_source,
        command=command,
        warnings=warnings,
        config_source={
            "base_config": base_config_name,
            "demo": bool(args.demo),
            "pipeline_commit": source.get("commit", "unknown"),
            "pipeline_ref": PINNED_REF,
        },
    )
    write_report_md(
        output_dir,
        status=status,
        mode=mode,
        findings=findings,
        config=config,
        pipeline_source=report_source,
        warnings=warnings,
    )
    write_reproducibility(
        output_dir,
        command=f"python clawbio.py run {CLI_ALIAS} " + _replay(args),
        files=_reproducibility_files(output_dir, results_dir, config),
    )
    for warning in warnings:
        print(f"[{SKILL_NAME}] warning: {warning}", file=sys.stderr)
    print(f"[{SKILL_NAME}] {status}: see {output_dir / 'report.md'}", file=sys.stderr)
    return 0


def _replay(args: argparse.Namespace) -> str:
    """The same run, as one copy-pasteable command line."""
    parts = []
    for flag, value in (
        ("--input", args.input),
        ("--output", args.output),
        ("--preset", args.preset),
        ("--cores", args.cores),
        ("--mem-mb", args.mem_mb),
        ("--binners", args.binners),
        ("--host-fasta", args.host_fasta),
        ("--checkm2-db", args.checkm2_db),
        ("--gtdbtk-db", args.gtdbtk_db),
        ("--taxonomy-db", args.taxonomy_db),
        ("--taxdump", args.taxdump),
        ("--kegg-db", args.kegg_db),
        ("--cog-db", args.cog_db),
        ("--pfam-db", args.pfam_db),
        ("--pipeline-dir", args.pipeline_dir),
        ("--conda-prefix", args.conda_prefix),
        ("--timeout-hours", args.timeout_hours),
    ):
        if value not in (None, ""):
            parts += [flag, str(value)]
    for flag, enabled in (
        ("--demo", args.demo),
        ("--check", args.check),
        ("--coassemble", args.coassemble),
        ("--skip-host-removal", args.skip_host_removal),
        ("--allow-remote-inputs", args.allow_remote_inputs),
        ("--force", args.force),
    ):
        if enabled:
            parts.append(flag)
    return " ".join(parts)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.input and not args.demo:
        error = SkillError(
            stage="input",
            error_code=ErrorCode.MISSING_INPUT,
            message="no input given",
            fix="Pass --input <samplesheet.tsv|reads_dir>, or run with --demo to "
            "generate synthetic reads.",
            details={"skill": SKILL_NAME, "alias": CLI_ALIAS},
        )
        return _fail(error, args)

    output_dir = Path(args.output).expanduser()
    try:
        return run(args)
    except SkillError as error:
        return _fail(error, args, output_dir=output_dir)
    except Exception as error:  # noqa: BLE001 - the CLI must not leak a traceback
        wrapped = SkillError(
            stage="unexpected",
            error_code=ErrorCode.UNEXPECTED_ERROR,
            message=f"{type(error).__name__}: {error}",
            fix=f"This is a bug in the {SKILL_NAME} skill. Re-run with --check to "
            "isolate which stage fails, and attach the traceback below.",
            details={"skill": SKILL_NAME, "alias": CLI_ALIAS},
        )
        import traceback

        wrapped.details["traceback"] = traceback.format_exc()[-4000:]
        return _fail(wrapped, args, output_dir=output_dir)


def _fail(error: SkillError, args: argparse.Namespace, output_dir: Path | None = None) -> int:
    print(json.dumps(error.to_dict(), indent=2), file=sys.stderr)
    directory = Path(output_dir if output_dir is not None else args.output).expanduser()
    if directory.is_dir():
        # A failed run still leaves a record: result.json must never claim success.
        try:
            source = {}
            if directory.is_dir():
                source = {"path": str(directory)}
            config = {}
            config_path = directory / "pipeline" / "config.yaml"
            if config_path.is_file():
                import yaml

                config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
            write_result_json(
                directory,
                status="error",
                mode="check" if args.check else ("demo" if args.demo else "run"),
                findings=None,
                config=config,
                pipeline_source=source,
                command="",
                warnings=[f"{error.error_code}: {error.message}"],
            )
        except Exception:  # noqa: BLE001 - never mask the original failure
            pass
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
