from __future__ import annotations

NFCORE_SCRNASEQ_VERSION = "4.1.0"

# Verbatim STAR_ALIGN ext.args from nf-core/scrnaseq 4.1.0 conf/modules.config.
# The macOS Docker workaround appends ``--outTmpDir`` to this exact base; pinning
# it here (one source of the 4.1.0 fact) means the override can never silently
# drop an upstream flag. ``test_pinned_star_args_match_sibling_checkout_if_present``
# cross-checks this against a real checkout when one is available.
STAR_ALIGN_BASE_EXT_ARGS = (
    "--readFilesCommand zcat --runDirPerm All_RWX --outWigType bedGraph "
    "--twopassMode Basic --outSAMtype BAM SortedByCoordinate "
    "--limitBAMsortRAM ${task.memory.toBytes()}"
)

# Resource ceilings copied from conf/test.config (with the 1 h cap raised to 4 h
# for emulation overhead). These are appropriate ONLY for the small nf-core test
# profile / ``--demo`` runs and MUST NOT be applied to real datasets — a human
# STAR index needs far more than 15 GB.
MACOS_DEMO_RESOURCE_LIMITS = {"cpus": 4, "memory": "15.GB", "time": "4.h"}

ALIGNER_OUTPUT_DIRS = {
    "simpleaf": "simpleaf",
    "star": "star",
    "kallisto": "kallisto",
    "cellranger": "cellranger",
    "cellrangerarc": "cellrangerarc",
    "cellrangermulti": "cellrangermulti",
}
COMMON_REQUIRED_OUTPUTS = ["pipeline_info"]
MULTIQC_REQUIRED_OUTPUT = "multiqc/multiqc_report.html"

# The three Cell Ranger presets. Cell Ranger is NOT distributed via
# bioconda/biocontainers (10x Genomics licensing), so it is shipped only in the
# docker/singularity containers nf-core builds for it; a conda/mamba profile
# cannot resolve it. Centralised here so preflight and output validation share
# one definition (audit F-01/F-02).
CELLRANGER_FAMILY_PRESETS = frozenset(
    {"cellranger", "cellrangerarc", "cellrangermulti"}
)
# Aligners for which nf-core/scrnaseq 4.1.0 publishes a top-level ``fastqc/`` tree
# on the raw input reads, so a missing FastQC report is a genuine failure. This is
# EVERY aligner: 4.1.0 runs FASTQC on the shared ``ch_fastq`` channel before any
# aligner-specific branching (``if (!params.skip_fastqc) FASTQC_CHECK(ch_fastq)`` in
# workflows/scrnaseq.nf), and the Cell Ranger family — including cellrangermulti —
# consumes that same channel, so its reads are FastQC'd too. output.md confirms it:
# "FastQC is applied to all aligners' input reads" (results/fastqc). Treating the
# Cell Ranger family as exempt would let a partially-failed run pass silently
# (audit H-02 corrects the earlier F-02 leniency). ``--skip-fastqc`` disables it.
FASTQC_GATED_ALIGNERS = frozenset(ALIGNER_OUTPUT_DIRS)

OFFICIAL_PARAMS: dict[str, dict[str, object]] = {
    "input": {"type": "string", "pattern": "^\\S+\\.csv$", "required": True},
    "outdir": {"type": "string", "required": True},
    "email": {
        "type": "string",
        "pattern": "^([a-zA-Z0-9_\\-\\.]+)@([a-zA-Z0-9_\\-\\.]+)\\.([a-zA-Z]{2,5})$",
    },
    "multiqc_title": {"type": "string"},
    "barcode_whitelist": {"type": "string"},
    "aligner": {
        "type": "string",
        "default": "simpleaf",
        "enum": [
            "kallisto",
            "star",
            "simpleaf",
            "cellranger",
            "cellrangerarc",
            "cellrangermulti",
        ],
    },
    "protocol": {"type": "string", "default": "auto"},
    "skip_multiqc": {"type": "boolean"},
    "skip_fastqc": {"type": "boolean"},
    "skip_cellbender": {"type": "boolean"},
    "skip_emptydrops": {"type": "boolean", "deprecated": True},
    "genome": {"type": "string"},
    "fasta": {"type": "string", "pattern": "^\\S+\\.fn?a(sta)?(\\.gz)?$"},
    "igenomes_ignore": {"type": "boolean"},
    "transcript_fasta": {"type": "string"},
    "gtf": {"type": "string"},
    "save_reference": {"type": "boolean"},
    "save_align_intermeds": {"type": "boolean"},
    "igenomes_base": {"type": "string", "default": "s3://ngi-igenomes/igenomes/"},
    "txp2gene": {"type": "string"},
    "simpleaf_index": {"type": "string"},
    "simpleaf_umi_resolution": {
        "type": "string",
        "default": "cr-like",
        "enum": [
            "cr-like",
            "cr-like-em",
            "parsimony",
            "parsimony-em",
            "parsimony-gene",
            "parsimony-gene-em",
        ],
    },
    "star_index": {"type": "string"},
    "star_ignore_sjdbgtf": {"type": "string"},
    "seq_center": {"type": "string"},
    "star_feature": {
        "type": "string",
        "default": "Gene",
        "enum": ["Gene", "GeneFull", "Gene Velocyto"],
    },
    "kallisto_index": {"type": "string"},
    "kb_t1c": {"type": "string"},
    "kb_t2c": {"type": "string"},
    "kb_workflow": {
        "type": "string",
        "default": "standard",
        "enum": ["standard", "lamanno", "nac"],
    },
    "cellranger_index": {"type": "string"},
    "skip_cellranger_renaming": {"type": "boolean"},
    "motifs": {"type": "string"},
    "cellrangerarc_config": {"type": "string"},
    "cellrangerarc_reference": {"type": "string"},
    "cellranger_vdj_index": {"type": "string"},
    "skip_cellrangermulti_vdjref": {"type": "boolean"},
    "gex_frna_probe_set": {"type": "string"},
    "gex_target_panel": {"type": "string"},
    "gex_cmo_set": {"type": "string"},
    "fb_reference": {"type": "string"},
    "vdj_inner_enrichment_primers": {"type": "string"},
    "gex_barcode_sample_assignment": {"type": "string"},
    "cellranger_multi_barcodes": {"type": "string"},
    "custom_config_version": {"type": "string", "default": "master"},
    "custom_config_base": {
        "type": "string",
        "default": "https://raw.githubusercontent.com/nf-core/configs/master",
    },
    "config_profile_name": {"type": "string"},
    "config_profile_description": {"type": "string"},
    "config_profile_contact": {"type": "string"},
    "config_profile_url": {"type": "string"},
    "version": {"type": "boolean"},
    "publish_dir_mode": {
        "type": "string",
        "default": "copy",
        "enum": ["symlink", "rellink", "link", "copy", "copyNoFollow", "move"],
    },
    "email_on_fail": {
        "type": "string",
        "pattern": "^([a-zA-Z0-9_\\-\\.]+)@([a-zA-Z0-9_\\-\\.]+)\\.([a-zA-Z]{2,5})$",
    },
    "plaintext_email": {"type": "boolean"},
    "max_multiqc_email_size": {
        "type": "string",
        "default": "25.MB",
        "pattern": "^\\d+(\\.\\d+)?\\.?\\s*(K|M|G|T)?B$",
    },
    "monochrome_logs": {"type": "boolean"},
    "hook_url": {"type": "string"},
    "multiqc_config": {"type": "string"},
    "multiqc_logo": {"type": "string"},
    "multiqc_methods_description": {"type": "string"},
    "validate_params": {"type": "boolean", "default": True},
    "pipelines_testdata_base_path": {
        "type": "string",
        "default": "https://raw.githubusercontent.com/nf-core/test-datasets/",
    },
    "trace_report_suffix": {"type": "string"},
    "help": {"type": ["boolean", "string"]},
    "help_full": {"type": "boolean"},
    "show_hidden": {"type": "boolean"},
}

WRAPPER_SUPPORTED_UPSTREAM_PARAMS = {
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
}

INTENTIONALLY_UNSUPPORTED_PARAMS = {
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
}
DEPRECATED_PARAMS = {
    name for name, meta in OFFICIAL_PARAMS.items() if meta.get("deprecated")
}
WRAPPER_DEPRECATED_ALIAS_PARAMS = {"skip_emptydrops"}

# Per-parameter rationale for every intentionally unsupported official parameter,
# so the CLI/--check surface can explain *why* each is excluded rather than just
# listing it (audit finding #7). Categories: institutional config metadata,
# interactive help/version flags, and notification/validation knobs that would
# weaken the wrapper's fixed reproducibility/validation policy.
INTENTIONALLY_UNSUPPORTED_REASONS = {
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
}

# ── Protocol routing matrix — faithful mirror of assets/protocols.json @ 4.1.0 ─
# Source of truth (verbatim semantics):
#   https://github.com/nf-core/scrnaseq/4.1.0/assets/protocols.json
# Keys are wrapper *presets* (standard == simpleaf); values are the protocol
# tokens that aligner's block defines, NORMALISED to match
# preflight._normalize_protocol_token (lowercase, separators stripped).
#
# This is the single source of truth for protocol rules the wrapper enforces
# ahead of upstream: which presets accept `auto`, which accept `smartseq`, and
# which presets may pass unknown/custom protocol strings through. cellrangermulti
# is intentionally ABSENT: that path is configured by the multi samplesheet, not
# the --protocol value.
PROTOCOLS_JSON_4_1_0 = {
    "standard": ("10xv1", "10xv2", "10xv3", "10xv4", "dropseq"),  # simpleaf/alevin
    "star": ("10xv1", "10xv2", "10xv3", "10xv4", "dropseq", "smartseq"),
    "kallisto": ("10xv1", "10xv2", "10xv3", "10xv4", "dropseq", "smartseq"),
    "cellranger": ("auto", "10xv1", "10xv2", "10xv3", "10xv4"),
    "cellrangerarc": ("auto",),
}
KNOWN_PROTOCOL_TOKENS = frozenset(
    token for tokens in PROTOCOLS_JSON_4_1_0.values() for token in tokens
)
PRESETS_SUPPORTING_AUTO_PROTOCOL = frozenset(
    preset for preset, tokens in PROTOCOLS_JSON_4_1_0.items() if "auto" in tokens
)
PRESETS_SUPPORTING_SMARTSEQ_PROTOCOL = frozenset(
    preset for preset, tokens in PROTOCOLS_JSON_4_1_0.items() if "smartseq" in tokens
)
# Map-routed presets with no `auto` fallback must be given an explicit protocol.
# (cellranger/cellrangerarc default to auto; cellrangermulti is samplesheet-driven
# and therefore absent from the matrix, so it is never forced.)
PRESETS_REQUIRING_EXPLICIT_PROTOCOL = (
    frozenset(PROTOCOLS_JSON_4_1_0) - PRESETS_SUPPORTING_AUTO_PROTOCOL
)
# Presets that may forward an unknown/custom protocol string straight to the
# aligner. The 4.1.0 usage docs document custom values only for the non-Cell
# Ranger aligners (simpleaf/star/kallisto). Derived from the routing matrix minus
# the Cell Ranger family so this can never drift from CELLRANGER_FAMILY_PRESETS
# (audit F-8). Conceptually distinct from PRESETS_REQUIRING_EXPLICIT_PROTOCOL even
# though both resolve to {standard, star, kallisto} in 4.1.0.
PRESETS_SUPPORTING_CUSTOM_PROTOCOL = (
    frozenset(PROTOCOLS_JSON_4_1_0) - CELLRANGER_FAMILY_PRESETS
)

POLICY_SOURCE_NFCORE_DOCS = "nfcore_scrnaseq_4_1_0_docs"
POLICY_SOURCE_CLAWBIO = "clawbio_wrapper"
