#!/usr/bin/env python3
"""Illumina Bridge: import DRAGEN export bundles into ClawBio."""

from __future__ import annotations

import argparse
import csv
import json
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.checksums import sha256_file
from clawbio.common.report import (
    generate_report_footer,
    generate_report_header,
    write_result_json,
)
from illumina_bundle import (
    BundleArtifacts,
    discover_bundle_artifacts,
    parse_qc_metrics,
    parse_sample_sheet,
    summarize_sample_sheet,
)
from illumina_providers import MetadataEnrichmentResult, build_metadata_provider


SKILL_VERSION = "0.1.0"
SKILL_NAME = "illumina-bridge"
SKILL_ALIAS = "illumina"
SKILL_DIR = Path(__file__).resolve().parent
DEMO_BUNDLE_DIR = SKILL_DIR / "demo_bundle"
DEFAULT_OUTPUT_DIR = "illumina_import"


def ensure_output_dir_ready(output_dir: Path) -> None:
    """Guard against silently overwriting a populated output directory."""

    output_dir.mkdir(parents=True, exist_ok=True)
    if any(output_dir.iterdir()):
        raise FileExistsError(
            f"Output directory already contains files: {output_dir}. "
            "Please choose a new directory."
        )


def merge_sample_metadata(
    sample_rows: list[dict[str, str]],
    metadata_result: MetadataEnrichmentResult,
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    """Merge optional ICA sample-level metadata into the sample manifest."""

    sample_lookup = {sample["sample_id"]: sample for sample in metadata_result.samples}
    merged_rows: list[dict[str, str]] = []
    matched = 0
    for row in sample_rows:
        enriched = dict(row)
        sample_meta = sample_lookup.get(row["sample_id"])
        if sample_meta:
            matched += 1
            enriched.update(
                {
                    "ica_sample_id": str(sample_meta.get("ica_sample_id", "")),
                    "ica_analysis_status": str(sample_meta.get("analysis_status", "")),
                    "ica_cohort": str(sample_meta.get("cohort", "")),
                    "ica_notes": str(sample_meta.get("notes", "")),
                }
            )
        else:
            enriched.update(
                {
                    "ica_sample_id": "",
                    "ica_analysis_status": "",
                    "ica_cohort": "",
                    "ica_notes": "",
                }
            )
        merged_rows.append(enriched)

    return merged_rows, {
        "samples_in_bundle": len(sample_rows),
        "samples_enriched": matched,
        "samples_unmatched": len(sample_rows) - matched,
    }


def build_downstream_routing_hints(
    *,
    vcf_path: Path,
    sample_count: int,
) -> list[dict[str, str]]:
    """Provide conservative suggestions for next ClawBio steps."""

    cohort_reason = (
        "Use the exported cohort VCF to quantify representation and heterozygosity."
        if sample_count > 1
        else "Use the exported single-sample VCF for follow-up interpretation workflows."
    )
    return [
        {
            "skill": "equity",
            "reason": cohort_reason,
            "example_command": f"python clawbio.py run equity --input {shlex.quote(str(vcf_path))} --output <dir>",
        },
        {
            "skill": "gwas",
            "reason": "Look up high-interest rsIDs from the imported VCF across external variant databases.",
            "example_command": "python clawbio.py run gwas --rsid <rsid> --output <dir>",
        },
        {
            "skill": "clinpgx",
            "reason": "Follow up candidate pharmacogenes or drug labels after DRAGEN variant review.",
            "example_command": "python clawbio.py run clinpgx --gene <gene_symbol> --output <dir>",
        },
    ]


def build_summary_and_data(
    *,
    bundle: BundleArtifacts,
    sample_rows: list[dict[str, str]],
    sample_summary: dict[str, Any],
    qc_summary: dict[str, Any],
    metadata_result: MetadataEnrichmentResult,
    metadata_merge: dict[str, Any],
    downstream_hints: list[dict[str, str]],
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Assemble deterministic result payloads for result.json."""

    summary = {
        "platform": "illumina",
        "source_type": "dragen_export",
        "sample_count": sample_summary["sample_count"],
        "metadata_status": metadata_result.status,
        "has_qc_metrics": True,
        "bundle_dir": str(bundle.bundle_dir),
    }
    data = {
        "platform": "illumina",
        "source_type": "dragen_export",
        "artifacts": {
            "bundle_dir": str(bundle.bundle_dir),
            "vcf": str(bundle.vcf_path),
            "qc_metrics": str(bundle.qc_path),
            "sample_sheet": str(bundle.sample_sheet_path),
        },
        "samples": sample_rows,
        "sample_summary": sample_summary,
        "qc_summary": qc_summary,
        "metadata_enrichment": {
            **metadata_result.to_dict(),
            "merge": metadata_merge,
        },
        "downstream_routing_hints": downstream_hints,
    }
    return summary, data


def write_sample_manifest(output_dir: Path, sample_rows: list[dict[str, str]]) -> Path:
    """Write the normalized sample manifest table."""

    tables_dir = output_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = tables_dir / "sample_manifest.csv"
    fieldnames = [
        "sample_id",
        "sample_name",
        "sample_project",
        "sample_type",
        "pair_id",
        "sample_feature",
        "lane",
        "index",
        "index2",
        "index_id",
        "project_name",
        "library_name",
        "library_prep_kit_name",
        "index_adapter_kit_name",
        "description",
        "ica_sample_id",
        "ica_analysis_status",
        "ica_cohort",
        "ica_notes",
    ]
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in sample_rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})
    return manifest_path


def write_reproducibility_bundle(
    output_dir: Path,
    *,
    bundle: BundleArtifacts,
    metadata_provider: str,
    metadata_status: str,
    ica_project_id: str | None,
    ica_run_id: str | None,
    demo: bool,
) -> Path:
    """Write commands, environment metadata, and checksums."""

    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)

    commands = [
        "#!/usr/bin/env bash",
        "# Reproduce this Illumina bundle import",
        f"# Date: {datetime.now(timezone.utc).isoformat()}",
    ]
    parts = [
        f"python clawbio.py run {SKILL_ALIAS}",
        "--demo" if demo else f"--input {shlex.quote(str(bundle.bundle_dir))}",
        f"--output {shlex.quote(str(output_dir))}",
    ]
    if metadata_provider:
        parts.append(f"--metadata-provider {metadata_provider}")
    if ica_project_id:
        parts.append(f"--ica-project-id {shlex.quote(ica_project_id)}")
    if ica_run_id:
        parts.append(f"--ica-run-id {shlex.quote(ica_run_id)}")
    if metadata_status == "mocked-demo":
        parts.insert(0, "ILLUMINA_ICA_API_KEY=")
    commands.append(" ".join(parts))
    (repro_dir / "commands.sh").write_text("\n".join(commands) + "\n", encoding="utf-8")

    env_yaml = [
        f"skill: {SKILL_NAME}",
        f"version: {SKILL_VERSION}",
        f"metadata_provider: {metadata_provider}",
        f"metadata_status: {metadata_status}",
        f"bundle_dir: {bundle.bundle_dir}",
        f"vcf: {bundle.vcf_path}",
        f"sample_sheet: {bundle.sample_sheet_path}",
        f"qc_metrics: {bundle.qc_path}",
    ]
    if ica_project_id:
        env_yaml.append(f"ica_project_id: {ica_project_id}")
    if ica_run_id:
        env_yaml.append(f"ica_run_id: {ica_run_id}")
    (repro_dir / "environment.yml").write_text("\n".join(env_yaml) + "\n", encoding="utf-8")

    checksum_lines = []
    for path in (bundle.vcf_path, bundle.sample_sheet_path, bundle.qc_path):
        checksum_lines.append(f"{sha256_file(path)}  {path.name}")
    (repro_dir / "checksums.sha256").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")
    return repro_dir


def write_markdown_report(
    output_dir: Path,
    *,
    bundle: BundleArtifacts,
    sample_rows: list[dict[str, str]],
    sample_summary: dict[str, Any],
    qc_summary: dict[str, Any],
    metadata_result: MetadataEnrichmentResult,
    downstream_hints: list[dict[str, str]],
) -> Path:
    """Write the main report.md artifact."""

    header = generate_report_header(
        title="Illumina / DRAGEN Import Report",
        skill_name=SKILL_NAME,
        input_files=[bundle.vcf_path, bundle.sample_sheet_path, bundle.qc_path],
        extra_metadata={
            "Platform": "Illumina",
            "Source type": "DRAGEN export",
            "Metadata provider": metadata_result.provider,
            "Metadata status": metadata_result.status,
        },
    )

    sample_lines = [
        "| Sample ID | Sample Name | Type | Project | Lane |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in sample_rows:
        sample_lines.append(
            f"| {row['sample_id']} | {row['sample_name'] or '-'} | "
            f"{row.get('sample_type') or '-'} | {row['sample_project'] or row.get('project_name') or '-'} | "
            f"{row['lane'] or '-'} |"
        )

    qc_lines = []
    for key in (
        "run_id",
        "instrument",
        "analysis_software",
        "workflow_version",
        "yield_gb",
        "percent_q30",
        "percent_q30_read1",
        "percent_q30_read2",
        "pf_reads",
        "percent_pf_reads",
        "cluster_density_k_mm2",
        "reported_sample_count",
        "completed_samples",
    ):
        if key in qc_summary:
            qc_lines.append(f"- **{key}**: {qc_summary[key]}")
    if not qc_lines:
        qc_lines.append("- No normalized QC metrics were available.")

    metadata_lines = [
        f"- **Provider**: {metadata_result.provider}",
        f"- **Status**: {metadata_result.status}",
    ]
    if metadata_result.project:
        metadata_lines.append(
            f"- **Project**: {metadata_result.project.get('name') or metadata_result.project.get('id')}"
        )
    if metadata_result.run:
        metadata_lines.append(
            f"- **Analysis**: {metadata_result.run.get('name') or metadata_result.run.get('id')}"
        )
        metadata_lines.append(
            f"- **Analysis status**: {metadata_result.run.get('status') or 'unknown'}"
        )
    for warning in metadata_result.warnings:
        metadata_lines.append(f"- **Warning**: {warning}")

    next_steps = [
        f"- **{hint['skill']}**: {hint['reason']} Example: `{hint['example_command']}`"
        for hint in downstream_hints
    ]

    body = "\n".join(
        [
            "## Import Summary",
            "",
            f"- **Bundle directory**: `{bundle.bundle_dir}`",
            f"- **Sample count**: {sample_summary['sample_count']}",
            f"- **Sample projects**: {', '.join(sample_summary['sample_projects']) or 'None'}",
            f"- **Sample types**: {', '.join(sample_summary.get('sample_types', [])) or 'Unknown'}",
            "",
            "## Discovered Artifacts",
            "",
            f"- **VCF**: `{bundle.vcf_path}`",
            f"- **SampleSheet**: `{bundle.sample_sheet_path}`",
            f"- **QC metrics**: `{bundle.qc_path}`",
            "",
            "## Sample Manifest",
            "",
            *sample_lines,
            "",
            "## QC Summary",
            "",
            *qc_lines,
            "",
            "## Metadata Enrichment",
            "",
            *metadata_lines,
            "",
            "## Recommended Next Steps",
            "",
            *next_steps,
            "",
            "## Privacy Boundary",
            "",
            "- This adapter only imports exported DRAGEN artifacts plus optional metadata-only ICA context.",
            "- No genomic payload is uploaded by this skill.",
            "",
            generate_report_footer(),
        ]
    )

    report_path = output_dir / "report.md"
    report_path.write_text(header + body, encoding="utf-8")
    return report_path


def import_bundle(
    *,
    bundle_dir: Path,
    output_dir: Path,
    vcf_override: str | Path | None = None,
    qc_override: str | Path | None = None,
    sample_sheet_override: str | Path | None = None,
    metadata_provider_name: str = "none",
    ica_project_id: str | None = None,
    ica_run_id: str | None = None,
    allow_mock_metadata: bool = False,
    check_output_dir: bool = True,
) -> dict[str, Any]:
    """Import an Illumina/DRAGEN bundle and materialize ClawBio artifacts."""

    if check_output_dir:
        ensure_output_dir_ready(output_dir)
    else:
        output_dir.mkdir(parents=True, exist_ok=True)

    bundle = discover_bundle_artifacts(
        bundle_dir,
        vcf_override=vcf_override,
        qc_override=qc_override,
        sample_sheet_override=sample_sheet_override,
    )
    sample_rows = parse_sample_sheet(bundle.sample_sheet_path)
    sample_summary = summarize_sample_sheet(sample_rows)
    qc_summary = parse_qc_metrics(bundle.qc_path)

    provider = build_metadata_provider(metadata_provider_name)
    metadata_result = provider.enrich(
        bundle_dir=bundle.bundle_dir,
        project_id=ica_project_id,
        run_id=ica_run_id,
        allow_mock=allow_mock_metadata,
    )

    merged_rows, metadata_merge = merge_sample_metadata(sample_rows, metadata_result)
    downstream_hints = build_downstream_routing_hints(
        vcf_path=bundle.vcf_path,
        sample_count=sample_summary["sample_count"],
    )
    summary, data = build_summary_and_data(
        bundle=bundle,
        sample_rows=merged_rows,
        sample_summary=sample_summary,
        qc_summary=qc_summary,
        metadata_result=metadata_result,
        metadata_merge=metadata_merge,
        downstream_hints=downstream_hints,
    )

    write_sample_manifest(output_dir, merged_rows)
    write_reproducibility_bundle(
        output_dir,
        bundle=bundle,
        metadata_provider=metadata_provider_name,
        metadata_status=metadata_result.status,
        ica_project_id=ica_project_id,
        ica_run_id=ica_run_id,
        demo=allow_mock_metadata and bundle.bundle_dir == DEMO_BUNDLE_DIR,
    )
    write_markdown_report(
        output_dir,
        bundle=bundle,
        sample_rows=merged_rows,
        sample_summary=sample_summary,
        qc_summary=qc_summary,
        metadata_result=metadata_result,
        downstream_hints=downstream_hints,
    )
    write_result_json(
        output_dir=output_dir,
        skill=SKILL_ALIAS,
        version=SKILL_VERSION,
        summary=summary,
        data=data,
        input_checksum=sha256_file(bundle.vcf_path),
    )

    return {
        "summary": summary,
        "data": data,
        "output_dir": str(output_dir),
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """CLI argument parser."""

    parser = argparse.ArgumentParser(
        description="Import Illumina / DRAGEN export bundles into ClawBio."
    )
    parser.add_argument("--input", help="Bundle directory containing VCF, SampleSheet, and QC metrics")
    parser.add_argument("--output", "-o", default=DEFAULT_OUTPUT_DIR, help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run against the bundled synthetic DRAGEN demo")
    parser.add_argument("--vcf", help="Explicit VCF path override")
    parser.add_argument("--qc", help="Explicit QC metrics path override")
    parser.add_argument("--sample-sheet", dest="sample_sheet", help="Explicit SampleSheet path override")
    parser.add_argument("--metadata-provider", choices=["none", "ica"], default="none")
    parser.add_argument("--ica-project-id", help="ICA project ID for metadata enrichment")
    parser.add_argument(
        "--ica-run-id",
        help="ICA analysis ID for metadata enrichment (not a sequencing run ID)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """CLI entrypoint."""

    args = parse_args(argv)
    if not args.demo and not args.input:
        print("ERROR: Use --input <bundle_dir> or --demo.", file=sys.stderr)
        return 1

    bundle_dir = DEMO_BUNDLE_DIR if args.demo else Path(args.input).expanduser().resolve()
    output_dir = Path(args.output).expanduser().resolve()

    try:
        result = import_bundle(
            bundle_dir=bundle_dir,
            output_dir=output_dir,
            vcf_override=args.vcf,
            qc_override=args.qc,
            sample_sheet_override=args.sample_sheet,
            metadata_provider_name=args.metadata_provider,
            ica_project_id=args.ica_project_id,
            ica_run_id=args.ica_run_id,
            allow_mock_metadata=args.demo,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print("Illumina Bridge import complete")
    print(f"  Output: {result['output_dir']}")
    print(f"  Samples: {result['summary']['sample_count']}")
    print(f"  Metadata status: {result['summary']['metadata_status']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
