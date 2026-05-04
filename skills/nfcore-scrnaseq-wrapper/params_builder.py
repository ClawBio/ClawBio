from __future__ import annotations

import re
import sys
from pathlib import Path

import yaml

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))

from errors import ErrorCode, SkillError
from schemas import PRESET_ALIGNERS

_WHITESPACE_RE = re.compile(r"\s")

_SKIP_FLAGS = (
    "skip_fastqc",
    "skip_multiqc",
    "skip_cellranger_renaming",
    "skip_cellrangermulti_vdjref",
)

_REFERENCE_PATH_FIELDS = (
    "fasta",
    "gtf",
    "transcript_fasta",
    "txp2gene",
    "simpleaf_index",
    "kallisto_index",
    "star_index",
    "cellranger_index",
    "barcode_whitelist",
    "kb_t1c",
    "kb_t2c",
    "motifs",
    "cellrangerarc_config",
    "cellranger_vdj_index",
    "gex_frna_probe_set",
    "gex_target_panel",
    "gex_cmo_set",
    "fb_reference",
    "vdj_inner_enrichment_primers",
    "gex_barcode_sample_assignment",
    "cellranger_multi_barcodes",
)


def _posix(value: str) -> str:
    """Resolve a user-supplied path to absolute and convert to forward-slash notation.

    Nextflow (Java) accepts forward slashes on all platforms, including Windows.
    Using as_posix() avoids YAML escape-sequence ambiguity with backslashes and
    ensures relative paths are anchored before Nextflow changes its working directory.
    """
    return Path(value).expanduser().resolve().as_posix()


def build_params_file(args, *, normalized_samplesheet: Path, output_dir: Path) -> tuple[Path, dict[str, object]]:
    params = build_effective_params(args, normalized_samplesheet=normalized_samplesheet, output_dir=output_dir)
    params_path = write_params_yaml(params, output_dir=output_dir)
    return params_path, params


def build_effective_params(args, *, normalized_samplesheet: Path, output_dir: Path) -> dict[str, object]:
    params = _build_base_params(args, normalized_samplesheet=normalized_samplesheet, output_dir=output_dir)
    _add_input_metadata_params(params, args)
    _add_skip_params(params, args)
    _add_aligner_tuning_params(params, args)
    _add_symbolic_reference_params(params, args)
    _add_save_flags(params, args)
    _add_reference_path_params(params, args)
    return params


def _build_base_params(args, *, normalized_samplesheet: Path, output_dir: Path) -> dict[str, object]:
    upstream_outdir = output_dir / "upstream" / "results"
    params: dict[str, object] = {
        # Forward slashes: safe on all platforms and unambiguous in YAML.
        "outdir": upstream_outdir.as_posix(),
        "aligner": PRESET_ALIGNERS[args.preset],
        "skip_cellbender": _skip_cellbender_enabled(args),
    }
    if not args.demo:
        params["input"] = _schema_safe_input_path(normalized_samplesheet, output_dir=output_dir)
    else:
        # nf-schema 4.x validates igenomes_base (an S3 URL) even when the test
        # profile provides explicit fasta/gtf. DNS failure aborts before any task
        # runs. igenomes_ignore suppresses that validation when iGenomes is unused.
        params["igenomes_ignore"] = True
    if args.protocol:
        params["protocol"] = args.protocol
    return params


def _skip_cellbender_enabled(args) -> bool:
    return bool(getattr(args, "skip_cellbender", False) or getattr(args, "skip_emptydrops", False))


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


def _add_save_flags(params: dict[str, object], args) -> None:
    if getattr(args, "save_reference", False):
        params["save_reference"] = True
    if getattr(args, "save_align_intermeds", False):
        params["save_align_intermeds"] = True


def _add_reference_path_params(params: dict[str, object], args) -> None:
    # All file paths — resolved to absolute POSIX before writing so Nextflow
    # can locate them regardless of its own working directory at runtime.
    explicit_refs = []
    for param_name in _REFERENCE_PATH_FIELDS:
        value = getattr(args, param_name, None)
        if value:
            params[param_name] = _posix(value)
            explicit_refs.append(param_name)
    # When fasta/gtf are provided explicitly, iGenomes is unused. Suppress
    # nf-schema DNS validation of the default igenomes_base S3 URL.
    if explicit_refs and "fasta" in explicit_refs:
        params.setdefault("igenomes_ignore", True)


def write_params_yaml(params: dict[str, object], *, output_dir: Path) -> Path:
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)
    params_path = repro_dir / "params.yaml"
    params_path.write_text(serialize_params_yaml(params), encoding="utf-8")
    return params_path


def serialize_params_yaml(params: dict[str, object]) -> str:
    return yaml.dump(params, allow_unicode=True, sort_keys=False)
