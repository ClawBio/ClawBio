#!/usr/bin/env python3
"""DAG step: preflight.

Runs the pure-Python-validatable subset of the real, unmodified
``preflight.run_preflight`` (samplesheet cross-checks, reference/mutual-
exclusion/enum/threshold/contaminant/ribo/nextflow-config validation, resume
compatibility, GENCODE autodetection, and the accumulated warnings). Calls the
same private ``_check_*`` functions ``run_preflight`` calls, in the same
order -- this file does not reimplement any validation logic.

Deliberately NOT run here: ``_check_java``, ``_check_nextflow`` /
``_check_nextflow_presence``, and ``_check_profile`` (the docker/singularity/
conda/... daemon-readiness checks). Those three require real Java, Nextflow,
and a container runtime to be installed -- infrastructure this sandbox does
not have and that a real run genuinely needs (see SKILL.md Dependencies).
Folding a missing-binary check into a "gracefully degrades" story here would
misrepresent 1.x's audited contract, which hard-requires them even under
``--check``. This step's result JSON records that split explicitly via
``runtime_checked_here: false`` rather than silently skipping it.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from _isolated_imports import purge_foreign_bare_modules  # noqa: E402

purge_foreign_bare_modules("errors", "schemas", "preflight", "args_common")

import preflight as _pf  # noqa: E402
from args_common import build_args  # noqa: E402
from errors import SkillError  # noqa: E402

_RUNTIME_CHECK_NOTE = (
    "java/nextflow presence+version and the docker/singularity/conda/... daemon "
    "readiness gate (1.x preflight.py: _check_java, _check_nextflow[_presence], "
    "_check_profile) were NOT run by this DAG step. They require real "
    "infrastructure (Java 17+, Nextflow >=25.04.3, a container/env backend) this "
    "sandbox does not have, and are the genuine external prerequisite a real run "
    "needs beyond what this DAG validates -- see SKILL.md Dependencies. Every "
    "other 1.x preflight.py check (samplesheet cross-checks, reference/mutual-"
    "exclusion/enum/threshold/contaminant/ribo/nextflow-config validation, resume "
    "compatibility, GENCODE autodetection) ran for real above."
)


def _run_config_preflight(args, *, pipeline_source: dict, samplesheet_summary: dict) -> dict:
    """Same call sequence as preflight.run_preflight, minus the 3 runtime-binary checks."""
    warnings: list[str] = []
    _pf._check_demo_network(args)
    _pf._check_remote_inputs(args, samplesheet_summary)
    aligner_effective = _pf._check_supported_options(args)
    _pf._check_igenomes_genome(args, pipeline_source=pipeline_source, warnings=warnings)
    _pf._check_email(getattr(args, "email", None))
    _pf._check_email(getattr(args, "email_on_fail", None), field="email_on_fail")
    _pf._check_threshold_bounds(args)
    output_dir = Path(args.output).expanduser().resolve()
    _pf.check_output_dir_available(output_dir, resume=getattr(args, "resume", False))
    _pf._check_contaminant_paths(args)
    _pf._check_ribo_paths(args)
    _pf._check_nextflow_configs(args, warnings)
    _pf._check_acceleration_compatibility(args)
    _pf._check_contaminant_screening(args, warnings)
    _pf._check_mutual_exclusions(args)
    _pf._check_samplesheet_summary(args, samplesheet_summary=samplesheet_summary)
    bam_paths = samplesheet_summary.get("bam_paths", [])
    refs = (
        {}
        if (getattr(args, "demo", False) or getattr(args, "_noinput", False))
        else _pf._check_references(args, bam_reprocessing=bool(bam_paths), warnings=warnings)
    )
    _pf._check_skip_alignment_requirements(args, bam_paths=bam_paths)
    _pf._collect_bam_reprocessing_warning(bam_paths, aligner_effective, warnings)
    _pf._collect_index_aligner_warnings(args, warnings)
    gencode_autodetected = _pf._collect_gencode_autodetect_warning(args, warnings)
    _pf._check_resume_compatibility(args, output_dir=output_dir, pipeline_source=pipeline_source)
    _pf._collect_umi_warnings(args, warnings)
    _pf._collect_demo_warnings(args, warnings)
    _pf._collect_platform_warnings(args, output_dir, warnings)
    handoff_available = _pf._collect_alignment_warnings(args, warnings)

    return {
        "ok": True,
        "runtime_checked_here": False,
        "runtime_check_note": _RUNTIME_CHECK_NOTE,
        "pipeline_source": pipeline_source,
        "references": refs,
        "aligner_effective": aligner_effective,
        "handoff_available": handoff_available,
        "gencode_autodetected": gencode_autodetected,
        "warnings": warnings,
        "samplesheet": {
            "sample_count": samplesheet_summary.get("sample_count", 0),
            "sample_names": samplesheet_summary.get("sample_names", []),
            "unknown_columns": samplesheet_summary.get("unknown_columns", []),
            "strandedness_counts": samplesheet_summary.get("strandedness_counts", {}),
            "bam_count": len(samplesheet_summary.get("bam_paths", [])),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the pure-Python subset of nf-core/rnaseq wrapper preflight.")
    parser.add_argument("--output", required=True)
    parser.add_argument("--demo", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--aligner", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--prokaryotic", required=True)
    parser.add_argument("--allow-remote-inputs", required=True)
    parser.add_argument("--skip-alignment", required=True)
    parser.add_argument("--pipeline-source", required=True, help="pipeline_source.json from resolve_source")
    parser.add_argument("--samplesheet-summary", required=True, help="samplesheet_summary.json from prepare_samplesheet")
    parser.add_argument("--out", required=True, help="Path to write preflight_result.json")
    parsed = parser.parse_args()

    args = build_args(
        output=parsed.output,
        demo=parsed.demo,
        input=parsed.input,
        aligner=parsed.aligner,
        genome=parsed.genome,
        fasta=parsed.fasta,
        gtf=parsed.gtf,
        gff=parsed.gff,
        profile=parsed.profile,
        prokaryotic=parsed.prokaryotic,
        allow_remote_inputs=parsed.allow_remote_inputs,
        skip_alignment=parsed.skip_alignment,
    )
    pipeline_source = json.loads(Path(parsed.pipeline_source).read_text(encoding="utf-8"))
    samplesheet_summary = json.loads(Path(parsed.samplesheet_summary).read_text(encoding="utf-8"))

    out_path = Path(parsed.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        result = _run_config_preflight(args, pipeline_source=pipeline_source, samplesheet_summary=samplesheet_summary)
    except SkillError as exc:
        out_path.write_text(json.dumps(exc.to_dict(), indent=2) + "\n", encoding="utf-8")
        print(json.dumps(exc.to_dict(), indent=2), file=sys.stderr)
        return 1

    out_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"[preflight] passed (config-only; warnings: {len(result['warnings'])}); {_RUNTIME_CHECK_NOTE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
