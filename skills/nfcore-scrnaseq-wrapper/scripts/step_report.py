#!/usr/bin/env python3
"""clawbio step: report.

Fans in every prior step's artifact (pipeline source, samplesheet summary,
preflight report, effective params, composed command) into a single
`report.md` + `result.json`, mirroring `gwas_lookup`'s `merge -> report`
shape and `phylogenetics_builder`'s provenance fan-in. This is the DAG's
terminal artifact: it states plainly whether the composed command is safe to
run for real (`preflight.ok`), and if not, exactly which check failed and how
to fix it — but it never runs Nextflow itself.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description="Fan in all prior step artifacts into report.md + result.json.")
    parser.add_argument("--pipeline-source", required=True)
    parser.add_argument("--samplesheet-summary", required=True)
    parser.add_argument("--preflight", required=True)
    parser.add_argument("--params", required=True)
    parser.add_argument("--command", required=True)
    parser.add_argument("--out-report", required=True)
    parser.add_argument("--out-result", required=True)
    args = parser.parse_args()

    pipeline_source = json.loads(Path(args.pipeline_source).read_text(encoding="utf-8"))
    samplesheet_summary = json.loads(Path(args.samplesheet_summary).read_text(encoding="utf-8"))
    preflight = json.loads(Path(args.preflight).read_text(encoding="utf-8"))
    params = json.loads(Path(args.params).read_text(encoding="utf-8"))
    command = json.loads(Path(args.command).read_text(encoding="utf-8"))

    ok = bool(preflight.get("ok"))
    result = {
        "skill": "nfcore_scrnaseq_wrapper",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "nextflow_executed": False,
        "ready_to_run": ok,
        "pipeline_source": pipeline_source,
        "samplesheet_summary": samplesheet_summary,
        "preflight": preflight,
        "params": params,
        "command": command,
    }
    Path(args.out_result).write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")

    lines = []
    lines.append("# nfcore_scrnaseq_wrapper — run preparation report")
    lines.append("")
    lines.append(
        "**Nextflow was NOT executed.** This DAG validates and prepares a real "
        "nf-core/scrnaseq run; only the steps below actually happened."
    )
    lines.append("")
    lines.append("## Pipeline source")
    lines.append(f"- kind: `{pipeline_source.get('source_kind')}`")
    lines.append(f"- ref: `{pipeline_source.get('source_ref')}`")
    lines.append(f"- resolved version: `{pipeline_source.get('resolved_version')}`")
    lines.append("")
    lines.append("## Samplesheet")
    lines.append(f"- samples: {samplesheet_summary.get('sample_count')}")
    lines.append(f"- unknown columns omitted: {samplesheet_summary.get('unknown_columns')}")
    lines.append("")
    lines.append("## Preflight")
    if ok:
        lines.append("**PASSED** — every checked prerequisite (Java, Nextflow, backend, "
                      "references, protocol/preset compatibility, local-first input policy) "
                      "is satisfied. The composed command below is safe to run.")
    else:
        failure = preflight.get("failure", {})
        lines.append(f"**FAILED** at `{failure.get('error_code')}` ({failure.get('stage', 'preflight')}):")
        lines.append(f"> {failure.get('message')}")
        lines.append("")
        lines.append(f"**Fix:** {failure.get('fix')}")
    lines.append("")
    lines.append("## Effective params.yaml")
    lines.append("```yaml")
    for key in sorted(params):
        lines.append(f"{key}: {params[key]!r}")
    lines.append("```")
    lines.append("")
    lines.append("## Composed command (not executed)")
    lines.append("```bash")
    lines.append(command.get("command_str", ""))
    lines.append("```")
    lines.append("")
    lines.append("## Next steps")
    if ok:
        lines.append(
            "1. Install Nextflow >=25.04.0 and a container runtime (Docker/Singularity) "
            "if not already available.\n"
            "2. Run the composed command above from a shell.\n"
            "3. Point downstream ClawBio scRNA skills at the resulting `.h5ad` under "
            "`<output_dir>/upstream/results/**/mtx_conversions/`."
        )
    else:
        lines.append(
            "Resolve the preflight failure above, then re-run `clawbio start` with corrected "
            "inputs (a fresh run — this port has no in-place retry)."
        )

    Path(args.out_report).write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[report] ready_to_run={ok}")
    print(json.dumps({"ready_to_run": ok}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
