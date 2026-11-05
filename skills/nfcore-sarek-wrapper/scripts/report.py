#!/usr/bin/env python3
# nfcore_sarek_wrapper / scripts/report.py
"""Stage 5 of the pipeline.yaml DAG: fan in all prior stages into a summary.

Not a port of any single 1.x module (1.x's ``reporting.py``/``provenance.py``
describe a completed Nextflow run's outputs, which this port never produces —
see SKILL.md). Writes ``report.md`` + ``result.json`` describing what WAS
validated/prepared and stating plainly that Nextflow was not executed.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Fan in DAG stages into report.md / result.json")
    parser.add_argument("--samplesheet-report", required=True)
    parser.add_argument("--preflight-report", required=True)
    parser.add_argument("--params-yaml", required=True)
    parser.add_argument("--command-json", required=True)
    parser.add_argument("--out-report", required=True)
    parser.add_argument("--out-result", required=True)
    args = parser.parse_args(argv)

    samplesheet = json.loads(Path(args.samplesheet_report).read_text(encoding="utf-8"))
    preflight_result = json.loads(Path(args.preflight_report).read_text(encoding="utf-8"))
    command_info = json.loads(Path(args.command_json).read_text(encoding="utf-8"))
    params_yaml_text = Path(args.params_yaml).read_text(encoding="utf-8")

    warnings = preflight_result.get("warnings") or []
    skip_runtime = bool(preflight_result.get("skip_runtime_checks"))

    lines: list[str] = []
    lines.append("# nfcore_sarek_wrapper — validation report")
    lines.append("")
    lines.append(
        "**Nextflow was NOT executed by this run.** This DAG validates a "
        "samplesheet, runs the wrapper's real preflight rules, and composes the "
        "exact nf-core/sarek 3.8.1 command a real run would execute — it stops "
        "there. See 'What a real run additionally needs' below."
    )
    lines.append("")
    lines.append("## Samplesheet")
    lines.append("")
    if samplesheet.get("mode") == "demo":
        lines.append("- Mode: `--demo` — samplesheet validation skipped (upstream `-profile test` supplies its own).")
    else:
        lines.append(f"- Samples: {samplesheet.get('sample_count')}")
        lines.append(f"- Patients: {samplesheet.get('patient_count')}")
        lines.append(f"- Analysis mode: `{samplesheet.get('analysis_mode')}`")
        pairings = samplesheet.get("pairings") or []
        if pairings:
            lines.append(f"- Pairings: {len(pairings)}")
    lines.append("")
    lines.append("## Preflight")
    lines.append("")
    lines.append(f"- Result: **passed**")
    lines.append(f"- Runtime probes (Java/Nextflow/backend): {'SKIPPED (sandbox mode)' if skip_runtime else 'checked'}")
    lines.append(f"- Java: {preflight_result.get('java_version') or 'n/a'}")
    lines.append(f"- Nextflow: {preflight_result.get('nextflow_version') or 'n/a'}")
    lines.append(f"- Composed profile: `{preflight_result.get('composed_profile')}`")
    if warnings:
        lines.append(f"- Warnings ({len(warnings)}):")
        for w in warnings:
            lines.append(f"  - {w}")
    else:
        lines.append("- Warnings: none")
    lines.append("")
    lines.append("## Composed Nextflow command")
    lines.append("")
    lines.append("```bash")
    lines.append(command_info.get("command_str", ""))
    lines.append("```")
    lines.append("")
    lines.append("## params.yaml")
    lines.append("")
    lines.append("```yaml")
    lines.append(params_yaml_text.rstrip("\n"))
    lines.append("```")
    lines.append("")
    lines.append("## What a real run additionally needs")
    lines.append("")
    lines.append("- Nextflow >=25.10.2 on PATH")
    lines.append("- Java >=17 on PATH")
    lines.append("- A container runtime matching the composed profile (Docker/Singularity/Apptainer/Podman/Conda/...)")
    lines.append("- Unless `--demo`: local reference genome data (iGenomes key or explicit `--fasta` + indices)")
    lines.append(
        "- Actually running `reproducibility/commands_preview.sh` (or the "
        "`nextflow_command.json` command list) from an environment with the above"
    )
    lines.append("")
    lines.append(
        "> ClawBio is a research and educational tool. It is not a medical device "
        "and does not provide clinical diagnoses. Consult a healthcare professional "
        "before making any medical decisions."
    )
    lines.append("")

    out_report = Path(args.out_report)
    out_report.parent.mkdir(parents=True, exist_ok=True)
    out_report.write_text("\n".join(lines), encoding="utf-8")

    result = {
        "status": "prepared_not_executed",
        "nextflow_executed": False,
        "samplesheet": samplesheet,
        "preflight": {
            "passed": True,
            "skip_runtime_checks": skip_runtime,
            "warnings": warnings,
            "java_version": preflight_result.get("java_version"),
            "nextflow_version": preflight_result.get("nextflow_version"),
            "composed_profile": preflight_result.get("composed_profile"),
        },
        "nextflow_command": command_info,
    }
    out_result = Path(args.out_result)
    out_result.parent.mkdir(parents=True, exist_ok=True)
    out_result.write_text(json.dumps(result, indent=2, default=str) + "\n", encoding="utf-8")

    print(f"[report] wrote {out_report} and {out_result}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
