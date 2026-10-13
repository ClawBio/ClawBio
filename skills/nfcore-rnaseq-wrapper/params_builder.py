from __future__ import annotations

from clawbio.common.nfcore.parameters import serialize_params_yaml as _serialize_params, write_params_yaml as _write_params

import sys
from pathlib import Path

import yaml

from clawbio.common.textio import write_text_lf

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))
sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules

purge_foreign_bare_modules("errors", "schemas")

from errors import ErrorCode, SkillError
from schemas import DEFAULT_TRIMMER as _DEFAULT_TRIMMER
from schemas import WHITESPACE_RE as _WHITESPACE_RE

_REFERENCE_PATH_FIELDS = (
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
)

_CONTAMINANT_PATH_FIELDS = (
    "kraken_db",
    "bbsplit_fasta_list",
    "bbsplit_index",
)
# sylph_db and sylph_taxonomy accept comma-separated path lists; each item
# must be expanded independently without mangling URI schemes.
_CONTAMINANT_MULTI_PATH_FIELDS = (
    "sylph_db",
    "sylph_taxonomy",
)

_MULTIQC_PATH_FIELDS = (
    "multiqc_config",
    "multiqc_logo",
    "multiqc_methods_description",
)

_SKIP_FLAGS = (
    "skip_trimming",
    "skip_alignment",
    "skip_pseudo_alignment",
    "skip_quantification_merge",
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
)

_SAVE_FLAGS = (
    "save_reference",
    "save_trimmed",
    "save_align_intermeds",
    "save_unaligned",
    "save_merged_fastq",
    "save_non_ribo_reads",
    "save_umi_intermeds",
    "save_kraken_assignments",
    "save_kraken_unassigned",
    "save_bbsplit_reads",
)

_NUMERIC_FIELDS = (
    "stranded_threshold",
    "unstranded_threshold",
    "min_mapped_reads",
    "pseudo_aligner_kmer_size",
    "min_trimmed_reads",
    "kallisto_quant_fraglen",
    "kallisto_quant_fraglen_sd",
)
_STRING_FIELDS = (
    "extra_trimgalore_args",
    "extra_fastp_args",
    "extra_fqlint_args",
    "extra_salmon_quant_args",
    "extra_kallisto_quant_args",
    "seq_center",
    "seq_platform",
    "salmon_quant_libtype",
    "featurecounts_feature_type",
    "gtf_extra_attributes",
    "gtf_group_features",
    "rseqc_modules",
    "email",
    "email_on_fail",
    "multiqc_title",
    "extra_star_align_args",
    "extra_bowtie2_align_args",
    "hisat2_build_memory",
    "gpu_container_options",
    "publish_dir_mode",
)
# featurecounts_group_type is kept out of _STRING_FIELDS so we can guard it
# against --gencode (upstream auto-flips to gene_type when --gencode is set).

_UMI_FIELDS = (
    "umi_dedup_tool",
    "umitools_bc_pattern",
    "umitools_bc_pattern2",
    "umitools_umi_separator",
    "umitools_extract_method",
    "umi_discard_read",
    "umitools_grouping_method",
)

_BOOLEAN_FIELDS = (
    "arm",
    # prokaryotic is a real hidden boolean parameter in the nf-core/rnaseq 3.26.0 schema.
    # Writing it to params.yaml ensures full reproducibility: replaying the run with only
    # params.yaml (without -profile prokaryotic) will still configure STAR correctly.
    "prokaryotic",
    "bam_csi_index",
    "stringtie_ignore_gtf",
    "gffread_transcript_fasta",
    "umitools_dedup_stats",
    "umitools_dedup_primary_only",
    "use_rustqc",
    "use_parabricks_star",
    "use_sentieon_star",
    "use_gpu_ribodetector",
)


def build_params_file(
    args, *, normalized_samplesheet: Path, output_dir: Path, gencode_autodetected: bool = False
) -> tuple[Path, dict[str, object]]:
    params = build_effective_params(
        args,
        normalized_samplesheet=normalized_samplesheet,
        output_dir=output_dir,
        gencode_autodetected=gencode_autodetected,
    )
    params_path = write_params_yaml(params, output_dir=output_dir)
    return params_path, params


def build_effective_params(
    args, *, normalized_samplesheet: Path, output_dir: Path, gencode_autodetected: bool = False
) -> dict[str, object]:
    params = _build_base_params(args, normalized_samplesheet=normalized_samplesheet, output_dir=output_dir)
    # Fully-hermetic profiles: --demo AND any self-contained nf-core test profile
    # (--profile test/test_full/test_prokaryotic/…, flagged as _noinput) own EVERY
    # pipeline parameter (input, references, and all QC/skip/tuning/contaminant/ribo/
    # umi/reporting knobs). Because a -params-file value overrides profile config in
    # Nextflow, writing any of them would silently alter the profile's hermetic run.
    # Only the wrapper-forced essentials remain — outdir, and aligner when explicit —
    # already set by _build_base_params, so we return here without emitting anything
    # else (audit OBS-3/F1; mirrors the scrnaseq wrapper). _apply_noinput_overrides
    # also clears the user's reference args upstream; this is defence-in-depth so a
    # tuning flag can never leak into a test-profile params.yaml.
    if getattr(args, "demo", False) or getattr(args, "_noinput", False):
        return params
    _add_aligner_params(params, args)
    _add_reference_params(params, args)
    _add_contaminant_params(params, args)
    _add_ribo_params(params, args)
    _add_umi_params(params, args)
    _add_flags(params, args)
    _add_tuning_params(params, args, gencode_autodetected=gencode_autodetected)
    return params


def _build_base_params(args, *, normalized_samplesheet: Path, output_dir: Path) -> dict[str, object]:
    params: dict[str, object] = {
        # Relative so nf-core's ^\S+$ schema validator accepts paths with spaces.
        # Nextflow runs with cwd=output_dir, so "upstream/results" resolves correctly.
        "outdir": "upstream/results",
    }
    # Self-contained nf-core test profiles (_noinput) ship their own aligner in the
    # profile config. A -params-file value overrides profile config, so only write
    # aligner for a real run, or when the user explicitly chose --aligner. Defaults
    # to writing it (back-compat) when _aligner_explicit was never recorded.
    if not getattr(args, "_noinput", False) or getattr(args, "_aligner_explicit", True):
        params["aligner"] = args.aligner
    if getattr(args, "demo", False):
        params["igenomes_ignore"] = True
    elif not getattr(args, "_noinput", False):
        params["input"] = _schema_safe_input_path(normalized_samplesheet, output_dir=output_dir)
    # _noinput: profile provides its own params.input; the wrapper writes neither
    # input nor igenomes_ignore — the profile config controls both.
    return params


def _schema_safe_input_path(normalized_samplesheet: Path, *, output_dir: Path) -> str:
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
            message="The normalized samplesheet path is not compatible with the nf-core/rnaseq schema.",
            fix=(
                "Use an output directory layout where the relative path to "
                "reproducibility/samplesheet.valid.csv contains no whitespace."
            ),
            details={"input_path": input_path},
        )
    return input_path


def _add_aligner_params(params: dict[str, object], args) -> None:
    if getattr(args, "pseudo_aligner", None):
        params["pseudo_aligner"] = args.pseudo_aligner
    if getattr(args, "trimmer", None) and args.trimmer != _DEFAULT_TRIMMER:
        params["trimmer"] = args.trimmer
    if getattr(args, "star_ignore_sjdbgtf", False):
        params["star_ignore_sjdbgtf"] = True


def _add_reference_params(params: dict[str, object], args) -> None:
    if getattr(args, "genome", None):
        params["genome"] = args.genome
    if getattr(args, "igenomes_base", None):
        params["igenomes_base"] = _posix_or_uri(args.igenomes_base)
    explicit_ref_names: list[str] = []
    for field in _REFERENCE_PATH_FIELDS:
        value = getattr(args, field, None)
        if value:
            params[field] = _posix_or_uri(value)
            explicit_ref_names.append(field)
    # Ignore the iGenomes bundle whenever the user supplies their own reference
    # (any explicit ref) without a named --genome — including the fully-prebuilt
    # index route (no --fasta) — so a dated bundled annotation/index is never
    # silently mixed in. --genome keeps iGenomes active; demo forces ignore.
    if getattr(args, "demo", False) or (not getattr(args, "genome", None) and explicit_ref_names):
        params.setdefault("igenomes_ignore", True)


def _add_contaminant_params(params: dict[str, object], args) -> None:
    for field in ("contaminant_screening", "contaminant_screening_input", "bracken_precision"):
        value = getattr(args, field, None)
        if value:
            params[field] = value
    bbsplit_enabled = False
    for field in _CONTAMINANT_PATH_FIELDS:
        value = getattr(args, field, None)
        if value:
            params[field] = _posix_or_uri(value)
            if field in {"bbsplit_fasta_list", "bbsplit_index"}:
                bbsplit_enabled = True
    for field in _CONTAMINANT_MULTI_PATH_FIELDS:
        value = getattr(args, field, None)
        if value:
            params[field] = _posix_multi(value)
    if getattr(args, "skip_bbsplit", False):
        params["skip_bbsplit"] = True
    elif bbsplit_enabled:
        params["skip_bbsplit"] = False


def _posix_or_uri(value: str) -> str:
    """Expand to absolute POSIX path, but pass URI schemes (s3://, http://) through unchanged."""
    if "://" in value:
        return value
    return Path(value).expanduser().resolve().as_posix()


def _posix_multi(value: str) -> str:
    """Expand each comma-separated path individually; pass URI schemes through unchanged."""
    items = [item.strip() for item in value.split(",") if item.strip()]
    return ",".join(_posix_or_uri(item) for item in items)


def _add_ribo_params(params: dict[str, object], args) -> None:
    if getattr(args, "remove_ribo_rna", False):
        params["remove_ribo_rna"] = True
    if getattr(args, "ribo_removal_tool", None):
        params["ribo_removal_tool"] = args.ribo_removal_tool
    if getattr(args, "ribo_database_manifest", None):
        params["ribo_database_manifest"] = _posix_or_uri(args.ribo_database_manifest)


def _add_umi_params(params: dict[str, object], args) -> None:
    if getattr(args, "with_umi", False):
        params["with_umi"] = True
    if getattr(args, "skip_umi_extract", False):
        params["skip_umi_extract"] = True
    for field in _UMI_FIELDS:
        value = getattr(args, field, None)
        if value:
            params[field] = value


def _add_flags(params: dict[str, object], args) -> None:
    if getattr(args, "enable_preseq", False):
        params["skip_preseq"] = False
    for flag in (*_SKIP_FLAGS, *_SAVE_FLAGS, *_BOOLEAN_FIELDS):
        if getattr(args, flag, False):
            params[flag] = True


def _add_tuning_params(params: dict[str, object], args, *, gencode_autodetected: bool = False) -> None:
    gencode = bool(getattr(args, "gencode", False) or gencode_autodetected)
    if gencode:
        params["gencode"] = True
    # deseq2_vst uses a tri-state sentinel (None / True / False).
    # None → omit from params.yaml → pipeline applies its upstream default (true = VST).
    # If nf-core ever changes that upstream default, this sentinel must be revisited.
    if getattr(args, "deseq2_vst", None) is False:
        params["deseq2_vst"] = False
    elif getattr(args, "deseq2_vst", None) is True:
        params["deseq2_vst"] = True
    for field in _NUMERIC_FIELDS:
        value = getattr(args, field, None)
        if value is not None:
            params[field] = value
    for field in _STRING_FIELDS:
        value = getattr(args, field, None)
        if value:
            params[field] = value
    for field in _MULTIQC_PATH_FIELDS:
        value = getattr(args, field, None)
        if value:
            params[field] = _posix_or_uri(value)
    # featurecounts_group_type is only written when --gencode is NOT set.
    # When --gencode is set, upstream auto-flips it to "gene_type" — writing
    # it explicitly would override the auto-flip and break biotype QC.
    if not gencode:
        value = getattr(args, "featurecounts_group_type", None)
        if value:
            params["featurecounts_group_type"] = value


def write_params_yaml(params: dict[str, object], *, output_dir: Path) -> Path:
    return _write_params(params, output_dir=output_dir, pipeline_name="nf-core/rnaseq")


def serialize_params_yaml(params: dict[str, object]) -> str:
    return _serialize_params(params, pipeline_name="nf-core/rnaseq")
