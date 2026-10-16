# ruff: noqa: E402
from __future__ import annotations

import re
import json
import sys
from pathlib import Path

from clawbio.common.nfcore.parameters import serialize_params_yaml as _serialize_params
from clawbio.common.nfcore.parameters import write_params_yaml as _write_params
from clawbio.common.nfcore.paths import resolve_path_or_uri

import yaml

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))


sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules

purge_foreign_bare_modules("errors", "schemas")

from errors import ErrorCode, SkillError


_WHITESPACE_RE = re.compile(r"\s")

# Reference path fields that should be POSIX-resolved when present.
_REFERENCE_PATH_FIELDS = (
    "fasta", "fasta_fai", "dict",
    "bwa", "bwamem2", "dragmap",
    "dbsnp", "dbsnp_tbi",
    "known_indels", "known_indels_tbi",
    "known_snps", "known_snps_tbi",
    "germline_resource", "germline_resource_tbi",
    "pon", "pon_tbi",
    "intervals",
    "ascat_alleles", "ascat_loci", "ascat_loci_gc", "ascat_loci_rt",
    "chr_dir", "mappability",
    "msisensor2_models", "msisensorpro_scan",
    "ngscheckmate_bed",
    "sentieon_dnascope_model",
    "snpeff_cache", "vep_cache",
    "cnvkit_reference",
    "cf_chrom_len",
    # Read filtering (BBSplit) — path to genome TSV / pre-built index.
    "bbsplit_fasta_list", "bbsplit_index",
)

# Annotation plugin paths (kept separate so we can validate them with their flags)
_ANNOTATION_PATH_FIELDS = (
    "dbnsfp", "dbnsfp_tbi",
    "mastermind_file",
    "phenotypes_file", "phenotypes_file_tbi",
    "spliceai_snv", "spliceai_snv_tbi",
    "spliceai_indel", "spliceai_indel_tbi",
    "bcftools_annotations", "bcftools_annotations_tbi",
    "bcftools_columns", "bcftools_header_lines",
    "condel_config",
    "snpsift_databases",
    "varlociraptor_scenario_tumor_only",
    "varlociraptor_scenario_somatic",
    "varlociraptor_scenario_germline",
)

_MULTIQC_PATH_FIELDS = ("multiqc_config", "multiqc_logo", "multiqc_methods_description")
_PROTECTED_EXTRA_PARAMS = {"input", "input_restart", "outdir"}
_PARAMS_FILE_PATH_FIELDS = (
    "input",
    "input_restart",
    "outdir_cache",
    "igenomes_base",
    *_REFERENCE_PATH_FIELDS,
    *_ANNOTATION_PATH_FIELDS,
    *_MULTIQC_PATH_FIELDS,
)

# Flags that should be omitted unless the user explicitly opted in.
# Sarek uses --skip_tools string instead — no individual skip_* flags.
_SKIP_FLAGS = ()

_SAVE_FLAGS = (
    "save_mapped", "save_output_as_bam", "save_reference",
    "save_trimmed", "save_split_fastqs",
    "save_bbsplit_reads",
)

# NOTE: `arm` is intentionally NOT here. Unlike nf-core/rnaseq (where `arm` is a
# real hidden boolean parameter), nf-core/sarek 3.8.1 has no `arm` parameter — ARM
# support is the `arm64` *profile* only. `--arm` therefore composes the profile via
# command_builder and is tracked as wrapper run-state (provenance manifest + resume
# drift), but must never be written to params.yaml or the schema validator rejects it.
_BOOLEAN_FLAGS = (
    "trim_fastq", "trim_nextseq",
    "umi_in_read_header",
    "sentieon_consensus",
    "wes", "no_intervals",
    "joint_germline", "joint_mutect2",
    "only_paired_variant_calling",
    "ignore_soft_clipped_bases",
    "cf_contamination_adjustment",
    "filter_vcfs", "normalize_vcfs",
    "snv_consensus_calling", "concatenate_vcfs",
    "vep_include_fasta", "vep_condel", "vep_dbnsfp",
    "vep_loftee", "vep_mastermind", "vep_phenotypes",
    "vep_spliceai", "vep_spliceregion",
    "build_only_index", "download_cache",
    "igenomes_ignore",
    "mastermind_mutations", "mastermind_var_iden", "mastermind_url",
)
# Note: save_reference lives only in _SAVE_FLAGS to avoid a redundant double-write.

_NUMERIC_FIELDS = (
    "split_fastq", "nucleotides_per_second",
    "clip_r1", "clip_r2",
    "three_prime_clip_r1", "three_prime_clip_r2",
    "length_required",
    "umi_length", "umi_base_skip",
    "markduplicates_pixel_distance",
    "ascat_min_base_qual", "ascat_min_counts", "ascat_min_map_qual",
    "ascat_ploidy", "ascat_purity",
    "cf_coeff", "cf_contamination", "cf_minqual", "cf_mincov",
    "cf_window",
    "consensus_min_count", "varlociraptor_chunk_size",
)

_STRING_FIELDS = (
    "tools", "skip_tools",
    "dbsnp_vqsr", "known_indels_vqsr", "known_snps_vqsr",
    "umi_read_structure", "group_by_umi_strategy",
    "umi_location", "umi_tag",
    "use_gatk_spark",
    "ascat_genome",
    "sentieon_haplotyper_emit_mode", "sentieon_dnascope_emit_mode",
    "sentieon_dnascope_pcr_indel_model", "gatk_pcr_indel_model",
    # freebayes_filter is `string` upstream — a vcflib/vcffilter expression
    # (default "30"), not an integer. See https://nf-co.re/sarek/3.8.1/parameters/ .
    "freebayes_filter",
    # cf_ploidy is `string` upstream (default "2") — accepts a comma-separated
    # list (e.g. "2,3,4") for Control-FREEC to try multiple ploidies.
    "cf_ploidy",
    "bcftools_filter_criteria",
    "vep_custom_args", "vep_version", "vep_out_format",
    "vep_cache_version", "vep_genome", "vep_species",
    "snpeff_db",
    "dbnsfp_consequence", "dbnsfp_fields",
    "phenotypes_include_types",
    "seq_center", "seq_platform",
    "email", "email_on_fail",
    "multiqc_title",
    "publish_dir_mode",
    "hook_url",
    "max_multiqc_email_size",
    "trace_report_suffix",
)

def _validate_native_scalar_params(params: dict[str, object], *, patterns: bool = True) -> None:
    from clawbio.common.nfcore.parameters import validate_against_contract

    contract = json.loads((_SKILL_DIR / "reproducibility/upstream_scalar_contract.json").read_text())
    try:
        validate_against_contract(params, contract, patterns=patterns)
    except ValueError as exc:
        raise SkillError(
            stage="validation", error_code=ErrorCode.INVALID_FLAG_COMBINATION,
            message=f"Invalid native sarek parameter: {exc}",
            fix="Use the scalar types and numeric bounds in the pinned nf-core/sarek schema.",
            details={"error": str(exc)},
        ) from exc


def load_run_params(args) -> dict[str, object]:
    params = _load_user_params_file(getattr(args, "params_file", None))
    if getattr(args, "demo", False):
        for name in (
            *_REFERENCE_PATH_FIELDS, *_ANNOTATION_PATH_FIELDS,
            "input", "input_restart", "genome", "igenomes_base", "igenomes_ignore",
            "dbsnp_vqsr", "known_indels_vqsr", "known_snps_vqsr",
        ):
            params.pop(name, None)
    return params


def build_params_file(
    args, *, normalized_samplesheet: Path, output_dir: Path
) -> tuple[Path, dict[str, object]]:
    params = build_effective_params(
        args,
        normalized_samplesheet=normalized_samplesheet,
        output_dir=output_dir,
    )
    params_path = write_params_yaml(params, output_dir=output_dir)
    return params_path, params


def build_effective_params(
    args, *, normalized_samplesheet: Path, output_dir: Path
) -> dict[str, object]:
    params = load_run_params(args)
    # Sarek 3.8.1 assigns params.input_restart = retrieveInput(...) at runtime,
    # so a user-provided input_restart value is overwritten before consumption.
    # A restart sheet validated by the wrapper is intentionally sent as input.
    params.pop("input_restart", None)
    params.update(_build_base_params(args, normalized_samplesheet=normalized_samplesheet, output_dir=output_dir))
    _add_aligner_params(params, args)
    _add_reference_params(params, args)
    _add_annotation_params(params, args)
    _add_flags(params, args)
    _add_tuning_params(params, args)
    _add_extra_params(params, args)
    _validate_native_scalar_params(params)
    return params


def _load_user_params_file(params_file: str | None) -> dict[str, object]:
    if not params_file:
        return {}
    path = Path(params_file).expanduser()
    try:
        from clawbio.common.nfcore.parameters import load_unique_yaml

        loaded = load_unique_yaml(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.MISSING_INPUT,
            message="--params-file could not be read.",
            fix="Provide a readable YAML params file or remove --params-file.",
            details={"params_file": str(path), "error": str(exc)},
        )
    except (ValueError, yaml.YAMLError) as exc:
        raise SkillError(
            stage="validation", error_code=ErrorCode.INVALID_FLAG_COMBINATION,
            message="--params-file is not valid, unambiguous UTF-8 YAML.",
            fix="Correct YAML syntax and remove duplicate keys.",
            details={"params_file": str(path), "error": str(exc)},
        ) from exc
    if loaded is None:
        return {}
    if not isinstance(loaded, dict):
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_FLAG_COMBINATION,
            message="--params-file must contain a YAML mapping.",
            fix="Use a standard nf-core/Sarek params YAML file with key/value pairs.",
            details={"params_file": str(path)},
        )
    params = dict(loaded)
    _validate_native_scalar_params(params, patterns=False)
    # Nextflow is launched from output_dir, not from the user's invocation
    # directory. Preserve native relative-path meaning by resolving only fields
    # that Sarek consumes as file/directory resources before changing cwd.
    for field in _PARAMS_FILE_PATH_FIELDS:
        value = params.get(field)
        if isinstance(value, str) and value.strip() and value.strip().lower() != "false":
            params[field] = _posix_or_uri(value)
    return params


def _build_base_params(args, *, normalized_samplesheet: Path, output_dir: Path) -> dict[str, object]:
    params: dict[str, object] = {
        # Relative so nf-core's ^\S+$ schema validator accepts paths with spaces.
        # Nextflow runs with cwd=output_dir, so "upstream/results" resolves correctly.
        "outdir": "upstream/results",
    }
    step = getattr(args, "step", None)
    if step:
        params["step"] = step
    # --demo runs the upstream `test` profile, whose test.config already defines
    # input, genome, igenomes_base, tools and split_fastq. The wrapper must NOT
    # write input or igenomes_ignore in that case — doing so would override the
    # profile's matched test reference and break the run. Same for self-contained
    # profiles (_noinput).
    if getattr(args, "demo", False) or getattr(args, "_noinput", False):
        pass
    elif getattr(args, "build_only_index", False) and not getattr(args, "input", None):
        # Build-only flow: cache download may be layered on, but is not input-free itself.
        params["input"] = False
    else:
        params["input"] = _schema_safe_input_path(normalized_samplesheet, output_dir=output_dir)
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
            message="The normalized samplesheet path is not compatible with the nf-core/sarek schema.",
            fix=(
                "Use an output directory layout where the relative path to "
                "reproducibility/samplesheet.valid.csv contains no whitespace."
            ),
            details={"input_path": input_path},
        )
    return input_path


def _add_aligner_params(params: dict[str, object], args) -> None:
    aligner = getattr(args, "aligner", None)
    if aligner:
        # Written even when it equals the pipeline default so the reproducibility
        # replay records the exact aligner used (matches the rnaseq wrapper).
        params["aligner"] = aligner
    use_gatk_spark = getattr(args, "use_gatk_spark", None)
    if use_gatk_spark:
        params["use_gatk_spark"] = use_gatk_spark


def _add_reference_params(params: dict[str, object], args) -> None:
    if getattr(args, "genome", None):
        params["genome"] = args.genome
    if getattr(args, "igenomes_base", None):
        params["igenomes_base"] = _posix_or_uri(args.igenomes_base)
    for field in _REFERENCE_PATH_FIELDS:
        value = getattr(args, field, None)
        if value not in (None, ""):
            params[field] = _reference_value(value)


def _add_annotation_params(params: dict[str, object], args) -> None:
    for field in _ANNOTATION_PATH_FIELDS:
        value = getattr(args, field, None)
        if value not in (None, ""):
            params[field] = _reference_value(value)


def _add_flags(params: dict[str, object], args) -> None:
    for flag in (*_SKIP_FLAGS, *_SAVE_FLAGS, *_BOOLEAN_FLAGS):
        if getattr(args, flag, False):
            params[flag] = True


# Minimum bounds from nextflow_schema.json. split_fastq is special: 0 (off) or >=250.
_NUMERIC_MIN = {
    "umi_length": 1, "umi_base_skip": 0, "length_required": 1,
    "clip_r1": 0, "clip_r2": 0, "three_prime_clip_r1": 0, "three_prime_clip_r2": 0,
    "nucleotides_per_second": 1,
    "ascat_min_base_qual": 0, "ascat_min_counts": 0, "ascat_min_map_qual": 0,
    "ascat_ploidy": 0, "cf_contamination": 0, "cf_minqual": 0, "cf_mincov": 0,
    "consensus_min_count": 1,
    "varlociraptor_chunk_size": 1,
}


def _validate_numeric_constraints(args) -> None:
    """Enforce the schema's numeric minimums before launching Nextflow."""
    def _num(field):
        value = getattr(args, field, None)
        if value in (None, ""):
            return None
        try:
            return float(value)
        except (TypeError, ValueError):
            return None  # non-numeric strings (e.g. cf_ploidy lists) are not checked here

    sf = _num("split_fastq")
    if sf is not None and sf != 0 and sf < 250:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_FLAG_COMBINATION,
            message=f"--split_fastq must be 0 (off) or >= 250 (got {getattr(args, 'split_fastq', None)}).",
            fix="Set --split-fastq 0 to disable splitting, or a value >= 250.",
            details={"split_fastq": getattr(args, "split_fastq", None)},
        )
    for field, minimum in _NUMERIC_MIN.items():
        v = _num(field)
        if v is not None and v < minimum:
            raise SkillError(
                stage="validation",
                error_code=ErrorCode.INVALID_FLAG_COMBINATION,
                message=f"--{field.replace('_', '-')} must be >= {minimum} (got {getattr(args, field, None)}).",
                fix=f"Provide a value >= {minimum}.",
                details={field: getattr(args, field, None), "minimum": minimum},
            )


def _add_tuning_params(params: dict[str, object], args) -> None:
    _validate_numeric_constraints(args)
    for field in _NUMERIC_FIELDS:
        value = getattr(args, field, None)
        if value is not None:
            params[field] = value
    for field in _STRING_FIELDS:
        value = getattr(args, field, None)
        if value:
            params[field] = _schema_safe_string(field, str(value))
    for field in _MULTIQC_PATH_FIELDS:
        value = getattr(args, field, None)
        if value:
            # MultiQC config/logo/methods may be a local file or a remote URL;
            # preserve URIs verbatim (parity with the genome reference fields).
            params[field] = _posix_or_uri(value)
    outdir_cache = getattr(args, "outdir_cache", None)
    if outdir_cache:
        # Upstream declares this as directory-path; resolve local values from
        # the user's invocation cwd before Nextflow is launched in output_dir.
        params["outdir_cache"] = _posix_or_uri(outdir_cache)


def _add_extra_params(params: dict[str, object], args) -> None:
    for key, value in (getattr(args, "_extras", None) or {}).items():
        if key in _PROTECTED_EXTRA_PARAMS:
            raise SkillError(
                stage="validation",
                error_code=ErrorCode.INVALID_FLAG_COMBINATION,
                message=f"--extra-param cannot override wrapper-managed parameter '{key}'.",
                fix="Use --input/--input-restart for samplesheets and --output for the wrapper output root.",
                details={"key": key},
            )
        # This is the escape hatch for native Sarek params not exposed as a
        # first-class wrapper flag; treat it as an explicit CLI override of a
        # reusable base --params-file.
        # `_merge_extra_params` has already typed each official exposed or
        # generic parameter against the upstream schema. Leave unknown native
        # escape-hatch values as text rather than guessing a type.
        coerced = value
        if key in {*_REFERENCE_PATH_FIELDS, *_ANNOTATION_PATH_FIELDS}:
            params[key] = _reference_value(coerced)
        elif key in {*_MULTIQC_PATH_FIELDS, "outdir_cache"} and coerced is not False:
            params[key] = _posix_or_uri(str(coerced))
        else:
            params[key] = coerced


def _schema_safe_string(field: str, value: str) -> str:
    """Canonicalise string params whose official schema is comma-list based."""
    if field in {"tools", "skip_tools", "use_gatk_spark"}:
        return ",".join(part.strip() for part in re.split(r"[,\s]+", value) if part.strip())
    if field in {"sentieon_haplotyper_emit_mode", "sentieon_dnascope_emit_mode"}:
        return ",".join(part.strip() for part in value.split(",") if part.strip())
    return value


def _posix_or_uri(value: str) -> str:
    resolved = resolve_path_or_uri(str(value), base_dir=Path.cwd())
    return resolved if isinstance(resolved, str) else resolved.as_posix()


def _reference_value(value: object) -> str | bool:
    """Preserve Sarek's documented `false` sentinel for reference overrides."""
    if value is False or (isinstance(value, str) and value.strip().lower() == "false"):
        return False
    return _posix_or_uri(str(value))


def write_params_yaml(params: dict[str, object], *, output_dir: Path) -> Path:
    return _write_params(params, output_dir=output_dir, pipeline_name="nf-core/sarek")


def serialize_params_yaml(params: dict[str, object]) -> str:
    return _serialize_params(params, pipeline_name="nf-core/sarek")
