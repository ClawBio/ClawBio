#!/usr/bin/env python3
# nfcore_sarek_wrapper / scripts/preflight_check.py
"""Stage 2 of the pipeline.yaml DAG: run the wrapper's preflight rules.

Faithfully calls the ported ``preflight.run_preflight()`` (same directory,
byte-identical to ClawBio 1.x) with a params dict assembled from this skill's
DAG inputs plus the stage-1 samplesheet report. Every individual rule
(§5.1-§5.6 in the 1.x source: tool/mode pairing, reference-path existence,
flag compatibility, annotation cache, ...) is unmodified.

The one thing this script changes relative to 1.x's own ``--check`` path is
ORCHESTRATION, not RULES: 1.x's ``run_preflight()`` unconditionally probes for
Java >=17, Nextflow >=25.10.2, and the requested container backend (docker/
singularity/...) — appropriate for a wrapper about to launch Nextflow, but
this port never launches Nextflow (see SKILL.md), and this sandbox has neither
Java nor Nextflow installed. Set ``skip_runtime_checks=true`` (the sandbox/CI
path) to skip exactly those three probe calls and run every other rule for
real; set it to ``false`` (the real-run path) to call ``run_preflight()``
completely unmodified, which requires actual Java/Nextflow/backend binaries
on PATH and will correctly raise MISSING_JAVA / MISSING_NEXTFLOW /
MISSING_DOCKER (etc.) if they are absent — that is preflight doing its job,
not a bug.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

import preflight
import schemas
from command_builder import compose_profile
from errors import SkillError
from pipeline_source import resolve_pipeline_source


def _truthy(raw: str) -> bool:
    return str(raw).strip().lower() in ("1", "true", "yes", "on")


def _none_if_sentinel(raw: str) -> str | None:
    v = (raw or "").strip()
    return None if not v or v.lower() == "none" else v


def _tool_list(raw: str) -> list[str]:
    v = _none_if_sentinel(raw)
    return [t.strip().lower() for t in v.split(",") if t.strip()] if v else []


def _run_preflight_reduced(
    *,
    params: dict[str, Any],
    samplesheet: dict[str, Any],
    pipeline_source: dict[str, Any],
    output_dir: Path,
    repo_root: Path,
    allow_remote_inputs: bool,
) -> preflight.PreflightResult:
    """Same rule sequence as ``preflight.run_preflight()``, minus the three
    runtime-binary probes (_check_java / _check_nextflow /
    _check_backends_for_profile). Every call below is the identical private
    helper 1.x's own ``run_preflight()`` calls, in the same order — see
    preflight.py's ``run_preflight`` for the unabridged version this mirrors.
    """
    warnings_acc: list[str] = []
    notes_acc: list[str] = []

    preflight._check_demo_network(params)
    preflight._check_remote_inputs(params, samplesheet, allow_remote_inputs=allow_remote_inputs)

    profile = str(params.get("profile") or "")
    preflight._check_profile_string(profile)

    preflight._check_step(params)
    preflight._check_aligner(params)
    preflight._check_enums(params)
    warnings_acc.extend(preflight._check_genome_known(params))
    preflight._check_tools_required_for_step(params)
    preflight._check_tools_known(preflight._normalize_tools(params.get("tools")))

    # --- SKIPPED (sandbox mode): _check_backends_for_profile ---
    preflight._check_output_dir(output_dir, repo_root=repo_root, resume=False)
    warnings_acc.extend(preflight._check_macos_docker_tmp(profile, output_dir))
    preflight._check_email(params.get("email"))
    preflight._check_email(params.get("email_on_fail"), field_name="email_on_fail")
    preflight._check_pipeline_source(pipeline_source)

    analysis_mode = str(samplesheet.get("analysis_mode") or "germline")
    pairings = samplesheet.get("pairings") or []
    present_modes = {str(p.get("mode")) for p in pairings if p.get("mode")}
    if not present_modes:
        present_modes = {analysis_mode}
    rows_by_patient = preflight._extract_rows_by_patient(samplesheet)
    tools = preflight._normalize_tools(params.get("tools"))
    warnings_acc.extend(
        preflight._check_tools_against_pairing(
            params=params, tools=tools, present_modes=present_modes, rows_by_patient=rows_by_patient,
        )
    )

    warnings_acc.extend(preflight._check_flag_compatibility(params=params, tools=tools))
    preflight._check_reference_paths(params)
    warnings_acc.extend(preflight._check_reference_path_suffixes(params))
    warnings_acc.extend(preflight._check_annotation_cache(params=params, output_dir=output_dir))

    return preflight.PreflightResult(
        warnings=warnings_acc,
        notes=notes_acc,
        java_version="(skipped: skip_runtime_checks=true)",
        nextflow_version="(skipped: skip_runtime_checks=true)",
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run Sarek wrapper preflight checks")
    parser.add_argument("--samplesheet-report", required=True)
    parser.add_argument("--step", required=True)
    parser.add_argument("--tools", required=True)
    parser.add_argument("--skip-tools", required=True)
    parser.add_argument("--aligner", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--igenomes-ignore", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--wes", required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--demo", required=True)
    parser.add_argument("--arm", required=True)
    parser.add_argument("--gpu", required=True)
    parser.add_argument("--spark-profile", required=True)
    parser.add_argument("--mutect-profile", required=True)
    parser.add_argument("--allow-remote-inputs", required=True)
    parser.add_argument("--skip-runtime-checks", required=True)
    parser.add_argument("--extra-params-json", required=True)
    parser.add_argument("--pipeline-version", required=True)
    parser.add_argument("--pipeline-local", required=True)
    parser.add_argument("--sarek-output-dir", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    demo = _truthy(args.demo)
    arm = _truthy(args.arm)
    gpu = _truthy(args.gpu)
    spark = _truthy(args.spark_profile)
    mutect = _truthy(args.mutect_profile)
    skip_runtime = _truthy(args.skip_runtime_checks)
    allow_remote = _truthy(args.allow_remote_inputs)

    params: dict[str, Any] = {}
    extras_path = _none_if_sentinel(args.extra_params_json)
    if extras_path:
        params.update(json.loads(Path(extras_path).read_text(encoding="utf-8")))

    params["step"] = args.step
    params["tools"] = _tool_list(args.tools)
    params["skip_tools"] = _tool_list(args.skip_tools)
    if _none_if_sentinel(args.aligner):
        params["aligner"] = args.aligner
    if _none_if_sentinel(args.genome):
        params["genome"] = args.genome
    if _truthy(args.igenomes_ignore):
        params["igenomes_ignore"] = True
    if _none_if_sentinel(args.fasta):
        params["fasta"] = args.fasta
    params["wes"] = _truthy(args.wes)

    composed_profile = compose_profile(
        args.profile, demo=demo, arm=arm, gpu=gpu, spark=spark, mutect_profile=mutect,
    )
    params["profile"] = composed_profile

    samplesheet_report = json.loads(Path(args.samplesheet_report).read_text(encoding="utf-8"))

    pipeline_local = _none_if_sentinel(args.pipeline_local)
    source = resolve_pipeline_source(
        requested_version=args.pipeline_version,
        local_pipeline_dir=Path(pipeline_local).expanduser().resolve() if pipeline_local else None,
    )

    sarek_output_dir_raw = _none_if_sentinel(args.sarek_output_dir)
    # The real wrapper's OUTPUT_DIR_INSIDE_REPO / OUTPUT_DIR_NOT_EMPTY rules
    # govern the eventual Nextflow `--output` directory, which this port never
    # creates (it does not execute Nextflow). When the caller supplies a real
    # intended path, honor it and run the real check; otherwise use a
    # documented placeholder outside the repo so `_check_output_dir` still
    # exercises its repo-containment rule without requiring the caller to
    # pre-decide a real run's output location this far in advance.
    sarek_output_dir = (
        Path(sarek_output_dir_raw).expanduser().resolve()
        if sarek_output_dir_raw
        else Path("/tmp/clawbio-nfcore-sarek-wrapper-preflight-placeholder")
    )
    repo_root = schemas.PROJECT_ROOT

    try:
        if skip_runtime:
            result = _run_preflight_reduced(
                params=params,
                samplesheet=samplesheet_report,
                pipeline_source=source,
                output_dir=sarek_output_dir,
                repo_root=repo_root,
                allow_remote_inputs=allow_remote,
            )
        else:
            result = preflight.run_preflight(
                params=params,
                samplesheet=samplesheet_report,
                pipeline_source=source,
                output_dir=sarek_output_dir,
                repo_root=repo_root,
                resume=False,
                resume_manifest=None,
                allow_remote_inputs=allow_remote,
            )
    except SkillError as exc:
        print(json.dumps(exc.to_dict(), indent=2), file=sys.stderr)
        return 1

    out = {
        "ok": True,
        "skip_runtime_checks": skip_runtime,
        "warnings": result.warnings,
        "notes": result.notes,
        "java_version": result.java_version,
        "nextflow_version": result.nextflow_version,
        "composed_profile": composed_profile,
        "pipeline_source": source,
        "params_snapshot": params,
        "sarek_output_dir_checked": str(sarek_output_dir),
    }
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2, default=str) + "\n", encoding="utf-8")

    print(f"[preflight] passed (warnings: {len(result.warnings)}, skip_runtime_checks={skip_runtime})")
    for w in result.warnings:
        print(f"  · warning: {w}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
