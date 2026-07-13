from __future__ import annotations

import json
import re
import shlex
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

from clawbio.common.portable_commands import build_portable_commands_sh
from clawbio.common.report import generate_report_footer, generate_report_header, write_result_json
from clawbio.common.textio import write_text_lf, write_text_lf_atomic

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))


sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules

purge_foreign_bare_modules("schemas")

from schemas import DEFAULT_TIMEOUT_HOURS, DEFAULT_TRIMMER, SKILL_ALIAS, SKILL_DIR, SKILL_NAME, SKILL_VERSION

_REPRO_PATH_FLAGS = (
    "fasta",
    "gtf",
    "gff",
    "transcript_fasta",
    "additional_fasta",
    "gene_bed",
    "splicesites",
    "star_index",
    "rsem_index",
    "hisat2_index",
    "bowtie2_index",
    "salmon_index",
    "kallisto_index",
    "sortmerna_index",
    "multiqc_config",
    "multiqc_logo",
    "multiqc_methods_description",
    "ribo_database_manifest",
    "kraken_db",
    "sylph_db",
    "sylph_taxonomy",
    "bbsplit_fasta_list",
    "bbsplit_index",
    "igenomes_base",
)
_REPRO_MULTI_PATH_FLAGS = frozenset({"sylph_db", "sylph_taxonomy"})
_OPTIONAL_VALUE_FLAGS = (
    "pseudo_aligner",
    "trimmer",
    "genome",
    "ribo_removal_tool",
    "umi_dedup_tool",
    "umitools_bc_pattern",
    "umitools_bc_pattern2",
    "umitools_umi_separator",
    "umitools_extract_method",
    "umi_discard_read",
    "umitools_grouping_method",
    "stranded_threshold",
    "unstranded_threshold",
    "min_mapped_reads",
    "pseudo_aligner_kmer_size",
    "extra_trimgalore_args",
    "extra_fastp_args",
    "extra_fqlint_args",
    "extra_salmon_quant_args",
    "extra_kallisto_quant_args",
    "min_trimmed_reads",
    "seq_center",
    "seq_platform",
    "salmon_quant_libtype",
    "kallisto_quant_fraglen",
    "kallisto_quant_fraglen_sd",
    "featurecounts_feature_type",
    "featurecounts_group_type",
    "gtf_extra_attributes",
    "gtf_group_features",
    "rseqc_modules",
    "extra_star_align_args",
    "extra_bowtie2_align_args",
    "hisat2_build_memory",
    "gpu_container_options",
    "email",
    "multiqc_title",
    "contaminant_screening",
    "contaminant_screening_input",
    "bracken_precision",
    "rsem_extra_args",
    "metadata",
    "formula",
    "contrast",
    "downstream_output",
    "publish_dir_mode",
    "email_on_fail",
)
_BOOLEAN_FLAGS = (
    "skip_quantification_merge",
    "skip_trimming",
    "skip_alignment",
    "skip_pseudo_alignment",
    "skip_markduplicates",
    "skip_bigwig",
    "skip_stringtie",
    "skip_fastqc",
    "skip_dupradar",
    "skip_qualimap",
    "skip_rseqc",
    "skip_biotype_qc",
    "skip_deseq2_qc",
    "skip_multiqc",
    "skip_qc",
    "skip_linting",
    "skip_gtf_filter",
    "skip_gtf_transcript_filter",
    "save_reference",
    "save_trimmed",
    "save_align_intermeds",
    "save_unaligned",
    "save_merged_fastq",
    "save_non_ribo_reads",
    "save_umi_intermeds",
    "save_bbsplit_reads",
    "remove_ribo_rna",
    "with_umi",
    "skip_umi_extract",
    "umitools_dedup_stats",
    "umitools_dedup_primary_only",
    "gencode",
    "star_ignore_sjdbgtf",
    "bam_csi_index",
    "stringtie_ignore_gtf",
    "gffread_transcript_fasta",
    "arm",
    "prokaryotic",
    "rapid_quant",
    "enable_preseq",
    "save_kraken_assignments",
    "save_kraken_unassigned",
    "use_rustqc",
    "use_parabricks_star",
    "use_sentieon_star",
    "use_gpu_ribodetector",
    "run_downstream",
    "skip_downstream",
    "skip_bbsplit",
    "allow_pipeline_version_override",
    # Recorded so a run that opted into remote inputs replays faithfully: rnaseq's
    # commands.sh re-invokes the wrapper (which re-runs preflight), so without this
    # the replay would fail REMOTE_INPUT_NOT_ALLOWED even though the bundle's
    # samplesheet/params.yaml already carry the remote URIs. Only emitted when the
    # user actually passed --allow-remote-inputs, like every other flag here.
    "allow_remote_inputs",
)

_PORTABILITY_NOTICE = """\

# ── Portability notice ────────────────────────────────────────────────────────
# Input paths in samplesheet.valid.csv are machine-specific. Local inputs are
# absolute (required by Nextflow): FASTQ paths for fresh runs, BAM paths for
# reprocessing runs. A run started with --allow-remote-inputs instead carries
# remote URIs (http(s)://, s3://, gs://, ...) that resolve from anywhere and need
# no remapping. Before replaying on a different machine:
#
#   1. Remap any local FASTQ/BAM paths in samplesheet.valid.csv (skip for remote URIs):
#        python3 reproducibility/remap_paths.py --old /original/prefix --new /new/prefix
#
#   2. Remap reference/index paths in commands.sh (if references moved):
#        python3 reproducibility/remap_paths.py --refs-old /original/refs --refs-new /new/refs
#
#   3. Update the --output path above if the output directory changed:
#        python3 reproducibility/remap_paths.py --output-dir /new/output/dir
#
#   4. Verify everything:
#        python3 reproducibility/remap_paths.py --verify
#
# If ClawBio is installed at a non-standard path on this machine:
#   CLAWBIO_REPO=/path/to/ClawBio bash reproducibility/commands.sh
"""

_CLAWBIO_REPO_FALLBACK = (
    'if [[ ! -d "$REPO_ROOT/skills" ]]; then\n'
    '  echo "ERROR: Could not locate repo root (no skills/ directory found)" >&2\n'
    '  exit 1\n'
    'fi'
)
# The reproducibility bundle always lives OUTSIDE the ClawBio checkout (the wrapper
# forbids --output inside the repo), so the template's "walk up to find skills/" never
# succeeds and the bundle must be told where the checkout is. Bake the generating
# checkout as the CLAWBIO_REPO default so same-machine replay works with no manual
# setup; a different machine overrides CLAWBIO_REPO. (Sarek/scRNA-seq do not need this —
# they replay Nextflow directly and never re-invoke the wrapper.)
_REPO_ROOT = _SKILL_DIR.parent.parent


def _clawbio_repo_fallback_patched() -> str:
    repo = _REPO_ROOT.as_posix()
    return (
        'if [[ ! -d "$REPO_ROOT/skills" ]]; then\n'
        f'  REPO_ROOT="${{CLAWBIO_REPO:-{repo}}}"\n'
        '  if [[ ! -d "$REPO_ROOT/skills" ]]; then\n'
        '    echo "ERROR: Could not locate the ClawBio checkout (no skills/ directory)." >&2\n'
        '    echo "Set CLAWBIO_REPO to your ClawBio checkout and re-run:" >&2\n'
        '    echo "  CLAWBIO_REPO=/path/to/ClawBio bash commands.sh" >&2\n'
        '    exit 1\n'
        '  fi\n'
        'fi'
    )


# The shared template header claims replay works "from anywhere inside the repository
# clone", but this bundle lives outside the clone. Correct the wording to match the
# baked CLAWBIO_REPO default above.
_ANCHOR_HEADER_REPLAY_HINT = "# from anywhere inside the repository clone."
_ANCHOR_HEADER_REPLAY_HINT_PATCHED = (
    "# from any directory. This bundle lives outside the ClawBio checkout; the script\n"
    "# locates the checkout it was generated from automatically (set CLAWBIO_REPO to\n"
    "# replay from a different checkout)."
)


def write_report(
    output_dir: Path,
    *,
    args,
    pipeline_source: dict[str, object],
    preflight_result: dict[str, object],
    parsed_outputs: dict[str, object],
    command_str: str,
    post_run_warnings: list[str] | None = None,
) -> Path:
    lines = build_report_lines(
        output_dir,
        args=args,
        pipeline_source=pipeline_source,
        preflight_result=preflight_result,
        parsed_outputs=parsed_outputs,
        command_str=command_str,
        post_run_warnings=post_run_warnings,
    )
    report_path = output_dir / "report.md"
    write_text_lf(report_path, "\n".join(lines))
    return report_path


def build_report_lines(
    output_dir: Path,
    *,
    args,
    pipeline_source: dict[str, object],
    preflight_result: dict[str, object],
    parsed_outputs: dict[str, object],
    command_str: str,
    post_run_warnings: list[str] | None = None,
) -> list[str]:
    header = generate_report_header(
        "nf-core/rnaseq Wrapper Report",
        SKILL_NAME,
        SKILL_VERSION,
        extra_metadata={
            "Aligner": str(getattr(args, "aligner", "")),
            "Profile": str(getattr(args, "profile", "")),
            "Pipeline source": str(pipeline_source.get("source_kind", "")),
            "Pipeline ref": str(pipeline_source.get("resolved_version", "")),
        },
    )
    warnings = preflight_result.get("warnings", [])
    warning_text = "none" if not warnings else ", ".join(str(w) for w in warnings)
    commands_sh = output_dir / "reproducibility" / "commands.sh"
    lines: list[str] = [
        header,
        "## Summary",
        "",
        f"- Aligner: `{getattr(args, 'aligner', '')}`",
        f"- Pseudo-aligner: `{getattr(args, 'pseudo_aligner', None) or 'none'}`",
        f"- Trimmer: `{getattr(args, 'trimmer', '')}`",
        f"- Samples: `{parsed_outputs.get('samples_detected', 0)}`",
        f"- Strandedness mode: `{_strandedness_summary(preflight_result)}`",
        "",
        "## Preflight",
        "",
        f"- Java: `{preflight_result.get('java', {}).get('version', '')}`",
        f"- Nextflow: `{preflight_result.get('nextflow', {}).get('version', '')}`",
        f"- Backend profile: `{getattr(args, 'profile', '')}`",
        f"- Warnings: {warning_text}",
        "",
        "## Outputs",
        "",
        f"- Preferred counts TSV: `{_value_or(parsed_outputs.get('preferred_counts_tsv', ''), 'not available')}`",
        f"- SummarizedExperiment RDS: `{_value_or(parsed_outputs.get('rds_file', ''), 'not available')}`",
        f"- TPM TSV: `{_value_or(parsed_outputs.get('tpm_tsv', ''), 'not available')}`",
        f"- MultiQC report: `{_value_or(parsed_outputs.get('multiqc_report', ''), 'not found')}`",
        f"- Pipeline info: `{_value_or(parsed_outputs.get('pipeline_info_dir', ''), 'not found')}`",
        f"- Per-sample mode: `{str(bool(parsed_outputs.get('skip_quantification_merge', False))).lower()}`",
        f"- HISAT2 no-quant mode: `{str(bool(parsed_outputs.get('hisat2_no_quant', False))).lower()}`",
        "",
        "## Reproducibility",
        "",
        f"- Command: `{command_str}`",
        f"- Repro bundle: `{output_dir / 'reproducibility'}`",
        f"- Replay: `CLAWBIO_REPO=<repo_root> bash {shlex.quote(str(commands_sh))}`",
        "",
        *_build_handoff_lines(parsed_outputs, args=args),
    ]
    if post_run_warnings:
        lines += [
            "## Post-run Warnings",
            "",
            *[f"- {w}" for w in post_run_warnings],
            "",
        ]
    lines.append(generate_report_footer())
    return lines


def _build_handoff_lines(parsed_outputs: dict[str, object] | str, *, args=None) -> list[str]:
    if isinstance(parsed_outputs, str):
        parsed_outputs = {"preferred_counts_tsv": parsed_outputs, "handoff_available": bool(parsed_outputs)}
    counts = str(parsed_outputs.get("preferred_counts_tsv", ""))
    if parsed_outputs.get("handoff_available") and counts:
        quoted_counts = shlex.quote(counts)
        metadata = getattr(args, "metadata", None) if args is not None else None
        formula = getattr(args, "formula", None) if args is not None else None
        contrast = getattr(args, "contrast", None) if args is not None else None
        downstream_output = getattr(args, "downstream_output", None) if args is not None else None
        metadata_arg = shlex.quote(str(metadata)) if metadata else "<your_metadata.csv>"
        formula_arg = shlex.quote(str(formula)) if formula else '"~ batch + condition"'
        contrast_arg = shlex.quote(str(contrast)) if contrast else '"condition,treated,control"'
        output_arg = shlex.quote(str(downstream_output)) if downstream_output else "<dir>"
        return [
            "## Next Steps",
            "",
            "```bash",
            '"${CLAWBIO_REPO}/clawbio.py" run rnaseq \\',
            f"  --counts {quoted_counts} \\",
            f"  --metadata {metadata_arg} \\",
            f"  --formula {formula_arg} \\",
            f"  --contrast {contrast_arg} \\",
            f"  --output {output_arg}",
            "```",
            "",
        ]
    if parsed_outputs.get("hisat2_no_quant"):
        detail = "HISAT2 mode does not produce merged quantification counts for rnaseq downstream analysis."
    elif parsed_outputs.get("skip_quantification_merge"):
        detail = "Per-sample quantification mode does not produce the merged count matrix required by rnaseq downstream analysis."
    else:
        detail = "No preferred merged count matrix was found for rnaseq downstream analysis."
    return ["## Next Steps", "", f"- {detail}", ""]


def write_repro_commands(output_dir: Path, *, args) -> None:
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)
    command_args = build_repro_command_args(output_dir, args=args)
    # Build the whole script in memory and write it ONCE, atomically. An in-place
    # --resume replay re-invokes the wrapper, which regenerates the very commands.sh
    # that bash is still executing; a truncate-and-rewrite (or several) would corrupt
    # bash's mid-run read. Composing in memory also removes the previous read/patch/
    # rewrite cycles.
    content = build_portable_commands_sh(
        skill_name=SKILL_NAME,
        # Phase 4 owns creating this orchestrator entrypoint; Phase 3 records
        # the planned replay target so generated bundles remain stable.
        script_name="nfcore_rnaseq_wrapper.py",
        args=command_args,
        generated_at=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
    )
    content = _apply_repo_fallback(content)
    content = _apply_idempotent_resume(content, output_dir)
    if not getattr(args, "demo", False):
        content += _PORTABILITY_NOTICE
    write_text_lf_atomic(repro_dir / "commands.sh", content)
    _write_remap_script(repro_dir)


def build_repro_command_args(output_dir: Path, *, args) -> dict[str, str | None]:
    command_args: dict[str, str | None] = {"--output": output_dir.as_posix()}
    # Mirror F7 (params_builder): a self-contained test profile owns the aligner, so
    # only record --aligner for a real run or when the user explicitly chose it.
    # Otherwise replaying commands.sh would pass --aligner and diverge from the
    # aligner-less params.yaml the original run wrote.
    if not getattr(args, "_noinput", False) or getattr(args, "_aligner_explicit", True):
        command_args["--aligner"] = getattr(args, "aligner", "")
    command_args["--profile"] = getattr(args, "profile", "")
    command_args["--pipeline-version"] = getattr(args, "pipeline_version", "")
    if getattr(args, "pipeline_local", None):
        command_args["--pipeline-local"] = Path(getattr(args, "pipeline_local")).expanduser().resolve().as_posix()
    if getattr(args, "demo", False):
        command_args["--demo"] = None
    elif not getattr(args, "_noinput", False):
        # _noinput: self-contained test profile supplies params.input; write neither --input nor --demo.
        # For normal runs, reference the bundled samplesheet copy so the replay command is
        # portable: commands.sh lives in reproducibility/ and ${SCRIPT_DIR} resolves to it.
        command_args["--input"] = "${SCRIPT_DIR}/samplesheet.valid.csv"
    for field in _OPTIONAL_VALUE_FLAGS:
        value = getattr(args, field, None)
        if value is None or value == "":
            continue
        # Omit the default --trimmer so commands.sh stays consistent with params.yaml
        # ("omit = trust upstream default", mirroring params_builder._add_aligner_params)
        # and never records a flag the user did not choose.
        if field == "trimmer" and value == DEFAULT_TRIMMER:
            continue
        command_args[f"--{field.replace('_', '-')}"] = str(value)
    # --timeout-hours has a non-None default, so it cannot ride the optional-value
    # loop above: only record it when the run overrode the pinned default, keeping
    # the replay command minimal yet faithful.
    timeout_hours = getattr(args, "timeout_hours", None)
    if timeout_hours is not None and float(timeout_hours) != float(DEFAULT_TIMEOUT_HOURS):
        command_args["--timeout-hours"] = str(timeout_hours)
    for field in _BOOLEAN_FLAGS:
        if getattr(args, field, False):
            command_args[f"--{field.replace('_', '-')}"] = None
    if getattr(args, "resume", False):
        command_args["--resume"] = None
    _deseq2_vst = getattr(args, "deseq2_vst", None)
    if _deseq2_vst is False:
        command_args["--no-deseq2-vst"] = None
    elif _deseq2_vst is True:
        command_args["--deseq2-vst"] = None
    for field in _REPRO_PATH_FLAGS:
        value = getattr(args, field, None)
        if value:
            command_args[f"--{field.replace('_', '-')}"] = _serialize_repro_path(field, str(value))
    user_configs = getattr(args, "nextflow_config", None) or []
    if user_configs:
        command_args["--nextflow-config"] = _stage_user_nextflow_configs(output_dir, user_configs)
    return command_args


def _safe_config_basename(name: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._")
    return safe or "nextflow.config"


def _stage_user_nextflow_configs(output_dir: Path, config_paths: list[str]) -> list[str]:
    """Copy each user ``-c``/``--nextflow-config`` file into
    ``reproducibility/nextflow_configs/`` and return ``${SCRIPT_DIR}``-relative
    references, so ``commands.sh`` replays the exact configs portably instead of baking
    a host-specific absolute path — which broke replay when the config lived outside the
    output dir (the argv rewriter had to fall back to an absolute path or ``<EDIT_ME>``).
    ``${SCRIPT_DIR}`` resolves to the bundle's ``reproducibility/`` dir at replay, so the
    copied configs travel with the bundle. Remote URIs and paths whose source file is
    missing fall back to the resolved value unchanged. Mirrors
    nfcore-scrnaseq-wrapper._copy_user_nextflow_configs (adapted to rnaseq's
    wrapper-re-invoking commands.sh, which anchors on ``${SCRIPT_DIR}``)."""
    config_dir = output_dir / "reproducibility" / "nextflow_configs"
    staged: list[str] = []
    for index, raw in enumerate(config_paths, start=1):
        raw = str(raw)
        if "://" in raw:
            staged.append(raw)  # remote config URI — not a local file to copy
            continue
        source = Path(raw).expanduser().resolve()
        if not source.is_file():
            staged.append(source.as_posix())  # defensive: nothing to copy in
            continue
        config_dir.mkdir(parents=True, exist_ok=True)
        # Idempotent re-bundling: an in-place --resume replay re-invokes the wrapper with
        # --nextflow-config already pointing at THIS bundle's staged copy
        # (${SCRIPT_DIR}/nextflow_configs/config_NN_*.config resolves inside config_dir).
        # Reference it in place instead of copying it again under a new config_NN_ prefix,
        # which would accumulate config_01_config_01_... on successive replays.
        if source.parent == config_dir.resolve():
            staged.append(f"${{SCRIPT_DIR}}/nextflow_configs/{source.name}")
            continue
        destination = config_dir / f"config_{index:02d}_{_safe_config_basename(source.name)}"
        shutil.copyfile(source, destination)
        staged.append(f"${{SCRIPT_DIR}}/nextflow_configs/{destination.name}")
    return staged


def _serialize_repro_path(field: str, value: str) -> str:
    if field in _REPRO_MULTI_PATH_FLAGS:
        return ",".join(_serialize_one_repro_path(item.strip()) for item in value.split(",") if item.strip())
    return _serialize_one_repro_path(value)


def _serialize_one_repro_path(value: str) -> str:
    if "://" in value:
        return value
    return Path(value).expanduser().resolve().as_posix()


def write_result(
    output_dir: Path,
    *,
    args,
    pipeline_source: dict[str, object],
    parsed_outputs: dict[str, object],
    command_str: str,
    post_run_warnings: list[str] | None = None,
) -> Path:
    output_artifacts = build_output_artifacts(parsed_outputs)
    return write_result_json(
        output_dir,
        skill=SKILL_ALIAS,
        version=SKILL_VERSION,
        # Minimal shared cross-wrapper contract: a successful run is ok/"ok".
        status="ok",
        ok=True,
        summary=build_result_summary(
            args=args,
            pipeline_source=pipeline_source,
            parsed_outputs=parsed_outputs,
            output_artifacts=output_artifacts,
            post_run_warnings=post_run_warnings,
        ),
        data=build_result_data(parsed_outputs=parsed_outputs, output_artifacts=output_artifacts, command_str=command_str),
    )


def build_output_artifacts(parsed_outputs: dict[str, object]) -> dict[str, object]:
    return dict(parsed_outputs)


def build_result_summary(
    *,
    args,
    pipeline_source: dict[str, object],
    parsed_outputs: dict[str, object],
    output_artifacts: dict[str, object],
    post_run_warnings: list[str] | None = None,
) -> dict[str, object]:
    return {
        "aligner": getattr(args, "aligner", ""),
        "aligner_effective": parsed_outputs.get("aligner_effective", ""),
        "pseudo_aligner": getattr(args, "pseudo_aligner", None),
        "pipeline_source_kind": pipeline_source.get("source_kind", ""),
        "pipeline_version_or_commit": pipeline_source.get("resolved_version", ""),
        "profile": getattr(args, "profile", ""),
        "resume_used": bool(getattr(args, "resume", False)),
        "preferred_counts_tsv": parsed_outputs.get("preferred_counts_tsv", ""),
        "rds_file": parsed_outputs.get("rds_file", ""),
        "multiqc_report": parsed_outputs.get("multiqc_report", ""),
        "handoff_available": parsed_outputs.get("handoff_available", False),
        "samples_detected": parsed_outputs.get("samples_detected", 0),
        "skip_quantification_merge": parsed_outputs.get("skip_quantification_merge", False),
        "hisat2_no_quant": parsed_outputs.get("hisat2_no_quant", False),
        "output_artifacts": output_artifacts,
        "post_run_warnings": list(post_run_warnings) if post_run_warnings else [],
    }


def build_result_data(
    *,
    parsed_outputs: dict[str, object],
    output_artifacts: dict[str, object],
    command_str: str,
) -> dict[str, object]:
    return {
        "canonical_skill_name": SKILL_NAME,
        "cli_alias": SKILL_ALIAS,
        "command": command_str,
        "output_artifacts": output_artifacts,
        "outputs": parsed_outputs,
    }


def write_check_result(output_dir: Path, payload: dict[str, object]) -> Path:
    path = output_dir / "check_result.json"
    write_text_lf(path, json.dumps(payload, indent=2))
    return path


def _apply_repo_fallback(content: str) -> str:
    """Bake the CLAWBIO_REPO default and correct the misleading replay-hint wording."""
    if _CLAWBIO_REPO_FALLBACK in content:
        content = content.replace(_CLAWBIO_REPO_FALLBACK, _clawbio_repo_fallback_patched())
    content = content.replace(_ANCHOR_HEADER_REPLAY_HINT, _ANCHOR_HEADER_REPLAY_HINT_PATCHED)
    return content


# The portable interpreter (`${PYTHON:-python3}`) is emitted by the shared
# clawbio/common/portable_commands template, so no per-skill rewrite is needed here.
_REPLAY_INVOCATION_LINE = '"${PYTHON:-python3}" "$SKILL_SCRIPT" \\'


def _apply_idempotent_resume(content: str, output_dir: Path) -> str:
    """Make an in-place replay idempotent.

    Unlike the Sarek/scRNA-seq bundles (which replay Nextflow directly, and Nextflow
    tolerates a populated output directory), the RNA-seq bundle re-invokes the wrapper,
    whose preflight rejects a non-empty --output with OUTPUT_DIR_NOT_EMPTY. So a plain
    re-run of commands.sh in the same directory fails. This injects a guard that adds
    `--resume` when the target output directory already holds a completed run of this
    bundle (reproducibility/manifest.json present); a fresh or `remap_paths.py
    --output-dir`-relocated output directory has no manifest and runs normally.

    This applies to `--demo` bundles too: a successful demo run populates its own
    output dir, so its bundle needs the guard just as much as a real run's. `-resume`
    is orthogonal to `-profile test` (nf-core documents no incompatibility) and the
    demo run's work tree and Nextflow cache both persist under the output dir.

    Skipped only when the command already carries `--resume` (the run always resumes).
    """
    if _REPLAY_INVOCATION_LINE not in content:
        return content
    if "\n    --resume" in content:
        return content
    manifest = f"{output_dir.as_posix()}/reproducibility/manifest.json"
    guard = (
        "# Idempotent replay: resume in place if this output directory already holds a\n"
        "# completed run of this bundle (reproducibility/manifest.json). A fresh or\n"
        "# remapped output directory (see the remap_paths.py --output-dir note below)\n"
        "# has no manifest and runs normally.\n"
        'CLAWBIO_RESUME=""\n'
        f'if [[ -f "{manifest}" ]]; then\n'
        '  CLAWBIO_RESUME="--resume"\n'
        "fi\n"
        '"${PYTHON:-python3}" "$SKILL_SCRIPT" $CLAWBIO_RESUME \\'
    )
    return content.replace(_REPLAY_INVOCATION_LINE, guard, 1)


def _write_remap_script(repro_dir: Path) -> None:
    repro_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(SKILL_DIR / "remap_paths.py", repro_dir / "remap_paths.py")


def _value_or(value: object, fallback: str) -> str:
    return str(value) if value else fallback


def _strandedness_summary(preflight_result: dict[str, object]) -> str:
    samplesheet = preflight_result.get("samplesheet", {})
    counts = samplesheet.get("strandedness_counts", {}) if isinstance(samplesheet, dict) else {}
    if not counts:
        return "not recorded"
    return ", ".join(f"{k}: {v}" for k, v in sorted(counts.items()))
