from __future__ import annotations

import math
import re
import json
import sys
from pathlib import Path

from clawbio.common.nfcore.parameters import serialize_params_yaml as _serialize_params
from clawbio.common.nfcore.parameters import write_params_yaml as _write_params
from clawbio.common.nfcore.paths import resolve_path_or_uri


_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))
sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules

purge_foreign_bare_modules("errors", "schemas")

from errors import ErrorCode, SkillError
from schemas import (
    ALL_REFERENCE_PATH_FIELDS,
    GENOME_REFERENCE_FIELDS,
    PRESET_ALIGNERS,
)

_WHITESPACE_RE = re.compile(r"\s")

_SKIP_FLAGS = (
    "skip_fastqc",
    "skip_multiqc",
    "skip_cellranger_renaming",
    "skip_cellrangermulti_vdjref",
)

# Every local path written to params.yaml (genome references + auxiliary files).
# Only GENOME_REFERENCE_FIELDS suppress iGenomes; see _add_reference_path_params.
_REFERENCE_PATH_FIELDS = ALL_REFERENCE_PATH_FIELDS


def _posix(value: str) -> str:
    """Resolve a user-supplied path to absolute and convert to forward-slash notation.

    Nextflow (Java) accepts forward slashes on all platforms, including Windows.
    Using as_posix() avoids YAML escape-sequence ambiguity with backslashes and
    ensures relative paths are anchored before Nextflow changes its working directory.
    """
    return Path(value).expanduser().resolve().as_posix()


def _posix_or_uri(value: str) -> str:
    resolved = resolve_path_or_uri(str(value), base_dir=Path.cwd())
    return resolved if isinstance(resolved, str) else resolved.as_posix()


def _validate_native_scalar_params(params: dict[str, object]) -> None:
    from clawbio.common.nfcore.parameters import validate_parameter_scalars

    contract = json.loads((_SKILL_DIR / "reproducibility/upstream_scalar_contract.json").read_text())
    try:
        validate_parameter_scalars(params, contract["parameters"])
    except ValueError as exc:
        raise SkillError(
            stage="validation", error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message=f"Invalid native scrnaseq parameter: {exc}",
            fix="Use the scalar types and numeric bounds in the pinned nf-core/scrnaseq schema.",
            details={"error": str(exc)},
        ) from exc


def build_params_file(
    args, *, normalized_samplesheet: Path, output_dir: Path
) -> tuple[Path, dict[str, object]]:
    params = build_effective_params(
        args, normalized_samplesheet=normalized_samplesheet, output_dir=output_dir
    )
    params_path = write_params_yaml(params, output_dir=output_dir)
    return params_path, params


def build_effective_params(
    args, *, normalized_samplesheet: Path, output_dir: Path
) -> dict[str, object]:
    params = _build_base_params(
        args, normalized_samplesheet=normalized_samplesheet, output_dir=output_dir
    )
    # Fully hermetic demo: the upstream `test` profile owns EVERY pipeline parameter
    # (input, references, protocol, and all QC/skip/tuning/save/reporting knobs).
    # Because a -params-file value overrides profile config in Nextflow, writing any
    # of them would silently alter the test run. Only the wrapper-forced essentials
    # remain — outdir, aligner=star, igenomes_ignore, skip_cellbender — all already
    # set by _build_base_params, so demo returns here without emitting anything else
    # (audit H-01).
    if getattr(args, "demo", False):
        return params
    _add_input_metadata_params(params, args)
    _add_generic_nfcore_params(params, args)
    _add_skip_params(params, args)
    _add_aligner_tuning_params(params, args)
    _add_save_flags(params, args)
    _add_symbolic_reference_params(params, args)
    _add_reference_path_params(params, args)
    _validate_native_scalar_params(params)
    return params


def _canonical_protocol(value: str) -> str:
    # Upstream protocols.json lookups are case-sensitive. Preflight accepts
    # these aliases, so emit the corresponding key rather than a custom token.
    if not isinstance(value, str):
        return value  # scalar validation reports the type error
    normalized = re.sub(r"[\s_-]+", "", value.lower())
    if re.fullmatch(r"10xv[1-4]", normalized):
        return normalized.upper()
    if normalized in {"auto", "dropseq", "smartseq"}:
        return normalized
    return value


def _build_base_params(
    args, *, normalized_samplesheet: Path, output_dir: Path
) -> dict[str, object]:
    params: dict[str, object] = {
        # Relative to cwd: the live run executes from output_dir and replay's
        # commands.sh cds into output_dir, so this resolves to the same place on
        # every machine — keeping params.yaml free of absolute prefixes (also
        # making resume checksums portable). Forward slashes are safe on all platforms.
        "outdir": "upstream/results",
        "aligner": PRESET_ALIGNERS[args.preset],
    }
    # skip_cellbender defaults to false upstream, so only emit it when requested
    # (keeps params.yaml to genuine overrides, matching the other skip flags).
    if _skip_cellbender_enabled(args):
        params["skip_cellbender"] = True
    if not args.demo:
        params["input"] = _schema_safe_input_path(
            normalized_samplesheet, output_dir=output_dir
        )
    else:
        # nf-schema 4.x validates igenomes_base (an S3 URL) even when the test
        # profile provides explicit fasta/gtf. DNS failure aborts before any task
        # runs. igenomes_ignore suppresses that validation when iGenomes is unused.
        params["igenomes_ignore"] = True
    # protocol is supplied by the test profile in demo mode (see hermetic-demo
    # rationale in build_effective_params); never override it from the CLI.
    if args.protocol and not args.demo:
        params["protocol"] = _canonical_protocol(args.protocol)
    return params


def _skip_cellbender_enabled(args) -> bool:
    return bool(
        getattr(args, "skip_cellbender", False)
        or getattr(args, "skip_emptydrops", False)
    )


def _schema_safe_input_path(normalized_samplesheet: Path, *, output_dir: Path) -> str:
    """Return an nf-core schema-compatible path to the normalized samplesheet.

    nf-core/scrnaseq 4.1.0 validates --input with ``^\\S+\\.csv$``. The wrapper
    stores the normalized samplesheet under ``output/reproducibility`` and runs
    Nextflow from ``output`` so this can be a stable, whitespace-free relative
    path even when the repository or output directory contains spaces.
    """
    normalized_samplesheet = normalized_samplesheet.resolve()
    output_dir = output_dir.resolve()
    try:
        input_path = normalized_samplesheet.relative_to(output_dir).as_posix()
    except ValueError:
        input_path = normalized_samplesheet.as_posix()
    if _WHITESPACE_RE.search(input_path):
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="The normalized samplesheet path is not compatible with the nf-core/scrnaseq schema.",
            fix=(
                "Use an output directory layout where the relative path to "
                "reproducibility/samplesheet.valid.csv contains no whitespace."
            ),
            details={"input_path": input_path},
        )
    return input_path


def _add_input_metadata_params(params: dict[str, object], args) -> None:
    # Input/output metadata — only written when provided.
    if getattr(args, "email", None):
        params["email"] = args.email
    if getattr(args, "multiqc_title", None):
        params["multiqc_title"] = args.multiqc_title


def _add_generic_nfcore_params(params: dict[str, object], args) -> None:
    if getattr(args, "email_on_fail", None):
        params["email_on_fail"] = args.email_on_fail
    for param_name in ("multiqc_config", "multiqc_logo", "multiqc_methods_description"):
        value = getattr(args, param_name, None)
        if value:
            params[param_name] = _posix_or_uri(value)
    if getattr(args, "publish_dir_mode", None):
        params["publish_dir_mode"] = args.publish_dir_mode
    if getattr(args, "trace_report_suffix", None):
        params["trace_report_suffix"] = args.trace_report_suffix
    if getattr(args, "monochrome_logs", False):
        params["monochrome_logs"] = True


def _add_skip_params(params: dict[str, object], args) -> None:
    # Skip flags — only written when True to keep params.yaml clean.
    for flag_name in _SKIP_FLAGS:
        if getattr(args, flag_name, False):
            params[flag_name] = True


def _add_aligner_tuning_params(params: dict[str, object], args) -> None:
    # STARsolo extras.
    if getattr(args, "star_ignore_sjdbgtf", False):
        # Schema type is 'string' (not boolean) — write "true" to pass nf-core JSON schema validation.
        # Groovy evaluates any non-empty string as truthy, so this correctly disables sjdbGTF.
        params["star_ignore_sjdbgtf"] = "true"
    if getattr(args, "seq_center", None):
        params["seq_center"] = args.seq_center

    # Aligner-specific tuning — only written when explicitly set.
    if getattr(args, "star_feature", None):
        params["star_feature"] = args.star_feature
    if getattr(args, "simpleaf_umi_resolution", None):
        params["simpleaf_umi_resolution"] = args.simpleaf_umi_resolution
    if getattr(args, "kb_workflow", None):
        params["kb_workflow"] = args.kb_workflow


def _add_symbolic_reference_params(params: dict[str, object], args) -> None:
    # iGenomes shortcut and CellRanger ARC reference — symbolic names, not paths.
    if getattr(args, "genome", None):
        params["genome"] = args.genome
    if getattr(args, "cellrangerarc_reference", None):
        params["cellrangerarc_reference"] = args.cellrangerarc_reference
    # igenomes_base is a base location for iGenomes (often an s3:// URL or a local
    # mirror). Pass it through verbatim — it must NOT be resolved to an absolute
    # local path (that would corrupt s3:// URLs), and it does not trigger
    # igenomes_ignore (setting it means iGenomes IS used, from a custom base).
    if getattr(args, "igenomes_base", None):
        params["igenomes_base"] = args.igenomes_base
    if getattr(args, "igenomes_ignore", False):
        params["igenomes_ignore"] = True


def _add_save_flags(params: dict[str, object], args) -> None:
    if getattr(args, "save_reference", False):
        params["save_reference"] = True
    save_align_intermeds = getattr(args, "save_align_intermeds", None)
    if save_align_intermeds is not None:
        params["save_align_intermeds"] = bool(save_align_intermeds)


def _add_reference_path_params(params: dict[str, object], args) -> None:
    # Local file paths — resolved to absolute POSIX before writing so Nextflow
    # can locate them regardless of its own working directory at runtime. Remote
    # URIs are preserved verbatim (see the _posix_or_uri call below).
    explicit_genome_refs = []
    for param_name in _REFERENCE_PATH_FIELDS:
        value = getattr(args, param_name, None)
        if value:
            # Preserve remote URIs verbatim; only local paths are resolved.
            params[param_name] = _posix_or_uri(value)
            if param_name in GENOME_REFERENCE_FIELDS:
                explicit_genome_refs.append(param_name)
    # Only a real GENOME reference (e.g. a prebuilt star_index, or fasta+gtf)
    # means iGenomes is unused, so suppress nf-schema DNS validation of the
    # default igenomes_base S3 URL. Auxiliary files (barcode whitelists, CMO/
    # probe/feature sets, primers, multi barcode samplesheets) are NOT genome
    # references and must never trigger igenomes_ignore — doing so would silently
    # break a `--genome <shortcut>` run that also supplies one of those files.
    if explicit_genome_refs:
        params.setdefault("igenomes_ignore", True)


def write_params_yaml(params: dict[str, object], *, output_dir: Path) -> Path:
    return _write_params(params, output_dir=output_dir, pipeline_name="nf-core/scrnaseq")


def serialize_params_yaml(params: dict[str, object], *, line_width: float = math.inf) -> str:
    return _serialize_params(params, pipeline_name="nf-core/scrnaseq", line_width=line_width)
