#!/usr/bin/env python3
"""DAG step: report.

Fans in every upstream stage's output (pipeline source resolution, samplesheet
summary, config-only preflight result, params.yaml, and the composed-but-not-
run Nextflow command) into a single report.md + result.json, mirroring the
report.md/result.json pair the real 1.x wrapper writes -- but honestly framed
as a prepared-and-validated run, not a completed one, since no Nextflow
execution happened. This is new glue code (a fan-in write-up), not a port of
1.x's reporting.py, which composes its report from post-execution pipeline
outputs (count matrices, MultiQC, etc.) this DAG never produces.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description="Write report.md + result.json for the nf-core/rnaseq DAG-prepared run.")
    parser.add_argument("--output", required=True)
    parser.add_argument("--pipeline-source", required=True)
    parser.add_argument("--samplesheet-summary", required=True)
    parser.add_argument("--preflight-result", required=True)
    parser.add_argument("--params-json", required=True)
    parser.add_argument("--command-json", required=True)
    parser.add_argument("--out-report", required=True)
    parser.add_argument("--out-result", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output).expanduser().resolve()
    pipeline_source = json.loads(Path(args.pipeline_source).read_text(encoding="utf-8"))
    samplesheet_summary = json.loads(Path(args.samplesheet_summary).read_text(encoding="utf-8"))
    preflight_result = json.loads(Path(args.preflight_result).read_text(encoding="utf-8"))
    params = json.loads(Path(args.params_json).read_text(encoding="utf-8"))
    command = json.loads(Path(args.command_json).read_text(encoding="utf-8"))

    mode = samplesheet_summary.get("mode", "user")
    warnings = preflight_result.get("warnings", [])

    lines = [
        "# nf-core/rnaseq Wrapper -- ClawBio 2.0 DAG Report",
        "",
        "**This run was validated and prepared, but Nextflow was NOT executed.** "
        "The DAG covers preflight validation, samplesheet normalization, and "
        "params/command composition -- the pure-Python stages a real run also "
        "performs before handing off to Nextflow. Actually launching the pipeline "
        "requires Nextflow, a container/env runtime, and (for non-demo runs) "
        "reference genome data; see SKILL.md Dependencies.",
        "",
        "## Summary",
        f"- Mode: `{mode}`",
        f"- Aligner (effective): `{preflight_result.get('aligner_effective', '')}`",
        f"- Samples: `{samplesheet_summary.get('sample_count', 0)}`",
        f"- Pipeline source: `{pipeline_source.get('source_kind')}` -> `{pipeline_source.get('resolved_version')}`",
        f"- Handoff to rnaseq-de available: `{preflight_result.get('handoff_available')}`",
        f"- Preflight warnings: `{len(warnings)}`",
        "",
        "## Config-only preflight",
        "Ran every pure-Python 1.x preflight.py check (samplesheet cross-checks, "
        "reference/mutual-exclusion/enum/threshold/contaminant/ribo/nextflow-config "
        "validation, resume compatibility, GENCODE autodetection). Did NOT run "
        "java/nextflow presence+version or the container-runtime daemon-readiness "
        "gate -- see `runtime_check_note` in preflight_result.json.",
        "",
    ]
    if warnings:
        lines.append("### Warnings")
        for w in warnings:
            lines.append(f"- {w}")
        lines.append("")

    lines += [
        "## params.yaml (effective parameters)",
        "```yaml",
        *[f"{k}: {v}" for k, v in params.items()],
        "```",
        "",
        "## Prepared Nextflow command (NOT executed)",
        "```bash",
        f"cd {command.get('cwd')}",
        command.get("command_str", ""),
        "```",
        "",
        "## Next Steps",
        "1. Install Nextflow >=25.04.3, Java >=17, and a container/env runtime (Docker, "
        "Singularity, Apptainer, Podman, or Conda/Mamba).",
        "2. Run the command printed above (or `bash <output>/reproducibility/... `) to "
        "launch the real pipeline.",
        "3. Once counts land under `<output>/upstream/results/`, hand off to `rnaseq-de` "
        "for differential expression.",
        "",
    ]
    report_text = "\n".join(lines) + "\n"

    result_payload = {
        "ok": True,
        "status": "prepared_not_executed",
        "skill": "nfcore-rnaseq-wrapper",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": mode,
        "pipeline_source": pipeline_source,
        "samplesheet_summary": samplesheet_summary,
        "preflight_result": preflight_result,
        "params": params,
        "prepared_command": command,
    }

    out_report = Path(args.out_report)
    out_report.parent.mkdir(parents=True, exist_ok=True)
    out_report.write_text(report_text, encoding="utf-8")

    out_result = Path(args.out_result)
    out_result.parent.mkdir(parents=True, exist_ok=True)
    out_result.write_text(json.dumps(result_payload, indent=2) + "\n", encoding="utf-8")

    # Also mirror into the wrapper's own output tree (matches 1.x layout: report.md
    # and result.json live at the output-directory root).
    (output_dir / "report.md").write_text(report_text, encoding="utf-8")
    (output_dir / "result.json").write_text(json.dumps(result_payload, indent=2) + "\n", encoding="utf-8")

    print(f"[report] wrote {out_report} and {out_result} (mirrored under {output_dir})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
