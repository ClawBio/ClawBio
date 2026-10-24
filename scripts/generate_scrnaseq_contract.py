#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import urllib.request
from pathlib import Path
from typing import Any

_SKILL_DIR = (
    Path(__file__).resolve().parent.parent / "skills" / "nfcore-scrnaseq-wrapper"
)


def fetch_schema(version: str) -> dict[str, Any]:
    url = f"https://raw.githubusercontent.com/nf-core/scrnaseq/{version}/nextflow_schema.json"
    print(f"Fetching {url}...")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as response:
        return json.loads(response.read().decode("utf-8"))


def parse_params(schema: dict[str, Any]) -> dict[str, dict[str, Any]]:
    params = {}
    defs = schema.get("$defs", schema.get("definitions", {}))

    for section_name, section in defs.items():
        properties = section.get("properties", {})
        for param_name, meta in properties.items():
            param_def = {}
            if "type" in meta:
                param_def["type"] = meta["type"]
            if "default" in meta:
                param_def["default"] = meta["default"]
            if "enum" in meta:
                param_def["enum"] = meta["enum"]
            if "pattern" in meta:
                param_def["pattern"] = meta["pattern"]
            if meta.get("deprecated"):
                param_def["deprecated"] = True

            # The schema might also define required at the section level
            if param_name in section.get("required", []):
                param_def["required"] = True

            params[param_name] = param_def

    # Manual default value overrides for parameters whose defaults are defined
    # in nextflow.config rather than nextflow_schema.json (so they match the live run).
    if "save_align_intermeds" in params:
        params["save_align_intermeds"]["default"] = True

    return params


def format_params_dict(params: dict[str, dict[str, Any]]) -> str:
    lines = ["OFFICIAL_PARAMS: dict[str, dict[str, object]] = {"]
    for name, meta in params.items():
        line = f'    "{name}": {repr(meta)},'
        lines.append(line)
    lines.append("}")
    return "\n".join(lines)


def generate_contract_code(version: str, params_str: str) -> str:
    # Use the static template for the rest of the file
    version_sanitized = version.replace(".", "_")
    template = f'''from __future__ import annotations

NFCORE_SCRNASEQ_VERSION = "{version}"

# Verbatim STAR_ALIGN ext.args from nf-core/scrnaseq {version} conf/modules.config.
# The macOS Docker workaround appends ``--outTmpDir`` to this exact base; pinning
# it here (one source of the {version} fact) means the override can never silently
# drop an upstream flag. ``test_pinned_star_args_match_sibling_checkout_if_present``
# cross-checks this against a real checkout when one is available.
STAR_ALIGN_BASE_EXT_ARGS = (
    "--readFilesCommand zcat --runDirPerm All_RWX --outWigType bedGraph "
    "--twopassMode Basic --outSAMtype BAM SortedByCoordinate "
    "--limitBAMsortRAM ${{task.memory.toBytes()}}"
)

# Resource ceilings copied from conf/test.config (with the 1 h cap raised to 4 h
# for emulation overhead). These are appropriate ONLY for the small nf-core test
# profile / ``--demo`` runs and MUST NOT be applied to real datasets — a human
# STAR index needs far more than 15 GB.
MACOS_DEMO_RESOURCE_LIMITS = {{"cpus": 4, "memory": "15.GB", "time": "4.h"}}

ALIGNER_OUTPUT_DIRS = {{
    "simpleaf": "simpleaf",
    "star": "star",
    "kallisto": "kallisto",
    "cellranger": "cellranger",
    "cellrangerarc": "cellrangerarc",
    "cellrangermulti": "cellrangermulti",
}}
COMMON_REQUIRED_OUTPUTS = ["pipeline_info"]
MULTIQC_REQUIRED_OUTPUT = "multiqc/multiqc_report.html"

# The three Cell Ranger presets (Cell Ranger is not on bioconda; docker/singularity
# only). Centralised so preflight and output validation share one definition.
CELLRANGER_FAMILY_PRESETS = frozenset(
    {{"cellranger", "cellrangerarc", "cellrangermulti"}}
)
# Aligners that publish a top-level fastqc/ tree (hard gate). 4.1.0 runs FASTQC on
# the shared ch_fastq before aligner branching, so EVERY aligner (incl. the Cell
# Ranger family) is FastQC'd unless --skip-fastqc (audit H-02).
FASTQC_GATED_ALIGNERS = frozenset(ALIGNER_OUTPUT_DIRS)

{params_str}

WRAPPER_SUPPORTED_UPSTREAM_PARAMS = {{
    "input",
    "outdir",
    "email",
    "multiqc_title",
    "barcode_whitelist",
    "aligner",
    "protocol",
    "skip_multiqc",
    "skip_fastqc",
    "skip_cellbender",
    "genome",
    "fasta",
    "igenomes_ignore",
    "transcript_fasta",
    "gtf",
    "save_reference",
    "save_align_intermeds",
    "igenomes_base",
    "txp2gene",
    "simpleaf_index",
    "simpleaf_umi_resolution",
    "star_index",
    "star_ignore_sjdbgtf",
    "seq_center",
    "star_feature",
    "kallisto_index",
    "kb_t1c",
    "kb_t2c",
    "kb_workflow",
    "cellranger_index",
    "skip_cellranger_renaming",
    "motifs",
    "cellrangerarc_config",
    "cellrangerarc_reference",
    "cellranger_vdj_index",
    "skip_cellrangermulti_vdjref",
    "gex_frna_probe_set",
    "gex_target_panel",
    "gex_cmo_set",
    "fb_reference",
    "vdj_inner_enrichment_primers",
    "gex_barcode_sample_assignment",
    "cellranger_multi_barcodes",
    "email_on_fail",
    "multiqc_config",
    "multiqc_logo",
    "multiqc_methods_description",
    "publish_dir_mode",
    "trace_report_suffix",
    "monochrome_logs",
    "skip_emptydrops",
}}

INTENTIONALLY_UNSUPPORTED_PARAMS = {{
    "custom_config_version",
    "custom_config_base",
    "config_profile_name",
    "config_profile_description",
    "config_profile_contact",
    "config_profile_url",
    "version",
    "plaintext_email",
    "max_multiqc_email_size",
    "hook_url",
    "validate_params",
    "pipelines_testdata_base_path",
    "help",
    "help_full",
    "show_hidden",
}}
DEPRECATED_PARAMS = {{name for name, meta in OFFICIAL_PARAMS.items() if meta.get("deprecated")}}
WRAPPER_DEPRECATED_ALIAS_PARAMS = {{"skip_emptydrops"}}

# Per-parameter rationale for every intentionally unsupported official parameter.
INTENTIONALLY_UNSUPPORTED_REASONS = {{
    "custom_config_version": "institutional nf-core/configs metadata; the wrapper pins its own config policy",
    "custom_config_base": "institutional nf-core/configs metadata; the wrapper pins its own config policy",
    "config_profile_name": "institutional profile metadata, set by site configs not the wrapper",
    "config_profile_description": "institutional profile metadata, set by site configs not the wrapper",
    "config_profile_contact": "institutional profile metadata, set by site configs not the wrapper",
    "config_profile_url": "institutional profile metadata, set by site configs not the wrapper",
    "version": "interactive version flag; the wrapper pins the pipeline version explicitly",
    "plaintext_email": "email-delivery tuning, outside the wrapper's local-first scope",
    "max_multiqc_email_size": "email-delivery tuning, outside the wrapper's local-first scope",
    "hook_url": "external notification webhook (Slack/Teams); not local-first",
    "validate_params": "disabling nf-schema validation would weaken the wrapper's fixed validation policy",
    "pipelines_testdata_base_path": "test-data source override; the wrapper owns demo/test wiring",
    "help": "interactive help flag handled by the wrapper's own --help",
    "help_full": "interactive help flag handled by the wrapper's own --help",
    "show_hidden": "interactive help flag handled by the wrapper's own --help",
}}

# ── Protocol routing matrix — faithful mirror of assets/protocols.json @ {version} ─
# Source of truth (verbatim semantics):
#   https://github.com/nf-core/scrnaseq/{version}/assets/protocols.json
# Keys are wrapper *presets* (standard == simpleaf); values are the protocol
# tokens that aligner's block defines, NORMALISED to match
# preflight._normalize_protocol_token (lowercase, separators stripped).
#
# This is the single source of truth for protocol rules the wrapper enforces
# ahead of upstream: which presets accept `auto`, which accept `smartseq`, and
# which presets may pass unknown/custom protocol strings through. cellrangermulti
# is intentionally ABSENT: that path is configured by the multi samplesheet, not
# the --protocol value.
PROTOCOLS_JSON_{version_sanitized} = {{
    "standard": ("10xv1", "10xv2", "10xv3", "10xv4", "dropseq"),  # simpleaf/alevin
    "star": ("10xv1", "10xv2", "10xv3", "10xv4", "dropseq", "smartseq"),
    "kallisto": ("10xv1", "10xv2", "10xv3", "10xv4", "dropseq", "smartseq"),
    "cellranger": ("auto", "10xv1", "10xv2", "10xv3", "10xv4"),
    "cellrangerarc": ("auto",),
}}
KNOWN_PROTOCOL_TOKENS = frozenset(
    token for tokens in PROTOCOLS_JSON_{version_sanitized}.values() for token in tokens
)
PRESETS_SUPPORTING_SMARTSEQ_PROTOCOL = frozenset(
    preset for preset, tokens in PROTOCOLS_JSON_{version_sanitized}.items() if "smartseq" in tokens
)
# Presets that may forward an unknown/custom protocol string (non-Cell-Ranger).
PRESETS_SUPPORTING_CUSTOM_PROTOCOL = (
    frozenset(PROTOCOLS_JSON_{version_sanitized}) - CELLRANGER_FAMILY_PRESETS
)

POLICY_SOURCE_NFCORE_DOCS = "nfcore_scrnaseq_{version_sanitized}_docs"
POLICY_SOURCE_CLAWBIO = "clawbio_wrapper"
'''
    return template


def main():
    parser = argparse.ArgumentParser(
        description="Generate nfcore-scrnaseq-wrapper contract from nf-core schema."
    )
    parser.add_argument(
        "--version",
        default="4.1.0",
        help="nf-core/scrnaseq version to fetch (e.g. 4.1.0)",
    )
    parser.add_argument(
        "--output",
        help="Output file path (default: overwrites the contract file in the skill directory)",
    )
    args = parser.parse_args()

    schema = fetch_schema(args.version)
    params = parse_params(schema)
    params_str = format_params_dict(params)

    code = generate_contract_code(args.version, params_str)

    version_sanitized = args.version.replace(".", "_")
    default_out = _SKILL_DIR / f"nfcore_{version_sanitized}_contract.py"
    output_path = Path(args.output) if args.output else default_out

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as f:
        f.write(code)

    print(f"Contract file written to {output_path}")
    print("NOTE: You may need to run `black` or `ruff format` on the output file.")
    print(
        "NOTE: If upgrading versions, ensure WRAPPER_SUPPORTED_UPSTREAM_PARAMS covers any new parameters."
    )


if __name__ == "__main__":
    main()
