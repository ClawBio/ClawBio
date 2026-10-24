from __future__ import annotations

from argparse import Namespace
from pathlib import Path
import sys

import pytest
import yaml

_SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_SKILL_DIR))


def _purge_foreign_modules(*names: str) -> None:
    for name in names:
        module = sys.modules.get(name)
        module_file = Path(getattr(module, "__file__", "") or "")
        if module is not None and _SKILL_DIR not in module_file.parents and module_file != _SKILL_DIR / f"{name}.py":
            sys.modules.pop(name, None)


def _purge_local_modules(*names: str) -> None:
    for name in names:
        module = sys.modules.get(name)
        module_file = Path(getattr(module, "__file__", "") or "")
        if module is not None and (module_file == _SKILL_DIR / f"{name}.py" or _SKILL_DIR in module_file.parents):
            sys.modules.pop(name, None)


def _remove_skill_dir_from_sys_path() -> None:
    while str(_SKILL_DIR) in sys.path:
        sys.path.remove(str(_SKILL_DIR))


_purge_foreign_modules("errors", "schemas", "params_builder")

from params_builder import build_effective_params, build_params_file

_purge_local_modules("errors", "schemas", "params_builder")
_remove_skill_dir_from_sys_path()


def _args(**overrides) -> Namespace:
    defaults = dict(
        demo=False,
        aligner="star_salmon",
        pseudo_aligner=None,
        trimmer="trimgalore",
        remove_ribo_rna=False,
        ribo_removal_tool=None,
        pseudo_aligner_kmer_size=None,
        extra_trimgalore_args=None,
        extra_fastp_args=None,
        extra_salmon_quant_args=None,
        extra_kallisto_quant_args=None,
        min_trimmed_reads=None,
        seq_center=None,
        seq_platform=None,
        salmon_quant_libtype=None,
        kallisto_quant_fraglen=None,
        kallisto_quant_fraglen_sd=None,
        bam_csi_index=False,
        stringtie_ignore_gtf=False,
        featurecounts_group_type=None,
        umitools_extract_method=None,
        umi_discard_read=None,
        umitools_grouping_method=None,
        with_umi=False,
        umi_dedup_tool=None,
        umitools_bc_pattern=None,
        umitools_bc_pattern2=None,
        umitools_umi_separator=None,
        skip_umi_extract=False,
        genome=None,
        fasta=None,
        gtf=None,
        gff=None,
        transcript_fasta=None,
        additional_fasta=None,
        gene_bed=None,
        splicesites=None,
        star_index=None,
        rsem_index=None,
        hisat2_index=None,
        bowtie2_index=None,
        salmon_index=None,
        kallisto_index=None,
        star_ignore_sjdbgtf=False,
        gencode=False,
        extra_star_align_args=None,
        extra_bowtie2_align_args=None,
        skip_trimming=False,
        skip_alignment=False,
        skip_pseudo_alignment=False,
        skip_quantification_merge=False,
        skip_markduplicates=False,
        skip_bigwig=False,
        skip_stringtie=False,
        skip_fastqc=False,
        skip_dupradar=False,
        skip_qualimap=False,
        skip_rseqc=False,
        skip_biotype_qc=False,
        skip_deseq2_qc=False,
        skip_multiqc=False,
        enable_preseq=False,
        skip_qc=False,
        save_reference=False,
        save_trimmed=False,
        save_align_intermeds=False,
        save_unaligned=False,
        save_merged_fastq=False,
        save_non_ribo_reads=False,
        save_umi_intermeds=False,
        stranded_threshold=None,
        unstranded_threshold=None,
        min_mapped_reads=None,
        featurecounts_feature_type=None,
        gtf_extra_attributes=None,
        gtf_group_features=None,
        deseq2_vst=None,
        rseqc_modules=None,
        email=None,
        multiqc_title=None,
        multiqc_config=None,
        multiqc_logo=None,
        prokaryotic=False,
        contaminant_screening=None,
        contaminant_screening_input=None,
        kraken_db=None,
        bracken_precision=None,
        sylph_db=None,
        sylph_taxonomy=None,
        bbsplit_fasta_list=None,
        bbsplit_index=None,
        save_kraken_assignments=False,
        save_kraken_unassigned=False,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


def _samplesheet(output_dir: Path) -> Path:
    path = output_dir / "reproducibility" / "samplesheet.valid.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("sample,fastq_1,strandedness\n", encoding="utf-8")
    return path


def _load_params(tmp_path: Path, args: Namespace | None = None, samplesheet: Path | None = None) -> dict[str, object]:
    args = args or _args()
    samplesheet = samplesheet or _samplesheet(tmp_path)
    params_path, payload = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    assert payload == yaml.safe_load(params_path.read_text(encoding="utf-8"))
    return payload


# ── aligner ──────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("aligner", ["star_salmon", "hisat2"])
def test_writes_requested_aligner(tmp_path, aligner):
    params = _load_params(tmp_path, _args(aligner=aligner))
    assert params["aligner"] == aligner


# ── outdir contract ───────────────────────────────────────────────────────────


def test_outdir_is_relative_upstream_results(tmp_path):
    # outdir is kept relative so nf-core's ^\S+$ validator accepts paths with
    # spaces (common on macOS). Nextflow resolves it against its CWD=output_dir.
    params = _load_params(tmp_path)
    assert params["outdir"] == "upstream/results"
    assert not Path(params["outdir"]).is_absolute()


# ── demo mode ────────────────────────────────────────────────────────────────


def test_demo_omits_input_and_sets_igenomes_ignore(tmp_path):
    params = _load_params(tmp_path, _args(demo=True))
    assert "input" not in params
    assert params["igenomes_ignore"] is True


# ── passthrough params present when set ──────────────────────────────────────


def test_reviewed_public_params_are_not_silently_dropped(tmp_path):
    params = _load_params(
        tmp_path,
        _args(
            pseudo_aligner_kmer_size=31,
            extra_trimgalore_args="--clip_R1 5",
            extra_fastp_args="--cut_front",
            extra_salmon_quant_args="--validateMappings",
            extra_kallisto_quant_args="--bias",
            min_trimmed_reads=10,
            seq_center="ClawBio",
            seq_platform="illumina",
            salmon_quant_libtype="A",
            kallisto_quant_fraglen=200,
            kallisto_quant_fraglen_sd=20,
            bam_csi_index=True,
            stringtie_ignore_gtf=True,
            featurecounts_group_type="gene_id",
            umitools_extract_method="regex",
            umi_discard_read=2,
            umitools_grouping_method="directional",
        ),
    )

    assert params["pseudo_aligner_kmer_size"] == 31
    assert params["extra_trimgalore_args"] == "--clip_R1 5"
    assert params["extra_fastp_args"] == "--cut_front"
    assert params["extra_salmon_quant_args"] == "--validateMappings"
    assert params["extra_kallisto_quant_args"] == "--bias"
    assert params["min_trimmed_reads"] == 10
    assert params["seq_center"] == "ClawBio"
    assert params["seq_platform"] == "illumina"
    assert params["salmon_quant_libtype"] == "A"
    assert params["kallisto_quant_fraglen"] == 200
    assert params["kallisto_quant_fraglen_sd"] == 20
    assert params["bam_csi_index"] is True
    assert params["stringtie_ignore_gtf"] is True
    assert params["featurecounts_group_type"] == "gene_id"
    assert params["umitools_extract_method"] == "regex"
    assert params["umi_discard_read"] == 2
    assert params["umitools_grouping_method"] == "directional"


# ── passthrough params absent when unset ─────────────────────────────────────


def test_reviewed_public_params_are_omitted_when_unset_or_false(tmp_path):
    params = _load_params(tmp_path)
    for key in (
        "pseudo_aligner_kmer_size",
        "extra_trimgalore_args",
        "extra_fastp_args",
        "min_trimmed_reads",
        "seq_center",
        "seq_platform",
        "salmon_quant_libtype",
        "kallisto_quant_fraglen",
        "kallisto_quant_fraglen_sd",
        "bam_csi_index",
        "stringtie_ignore_gtf",
        "featurecounts_group_type",
        "umitools_extract_method",
        "umi_discard_read",
        "umitools_grouping_method",
    ):
        assert key not in params


# ── reference handling ────────────────────────────────────────────────────────


def test_genome_shortcut_is_symbolic(tmp_path):
    params = _load_params(tmp_path, _args(genome="GRCh38"))
    assert params["genome"] == "GRCh38"


def test_reference_paths_are_absolute_posix(tmp_path):
    """Spot-check three representative reference flags (sequence, annotation, index)."""
    for flag in ("fasta", "gtf", "salmon_index"):
        filename = {"fasta": "ref.fa", "gtf": "ref.gtf", "salmon_index": "index"}[flag]
        ref = tmp_path / flag / filename
        ref.parent.mkdir(exist_ok=True)
        ref.write_text("x", encoding="utf-8")
        params = _load_params(tmp_path / flag, _args(**{flag: str(ref)}))
        assert params[flag] == ref.resolve().as_posix()
        assert "\\" not in params[flag]


def test_explicit_refs_set_igenomes_ignore(tmp_path):
    """Providing explicit fasta (or fasta+gtf) must inject igenomes_ignore=True."""
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    gtf = tmp_path / "genes.gtf"
    gtf.write_text("", encoding="utf-8")
    params = _load_params(tmp_path, _args(fasta=str(fasta), gtf=str(gtf)))
    assert params["igenomes_ignore"] is True


def test_igenomes_ignore_not_set_when_using_genome_shortcut(tmp_path):
    params = _load_params(tmp_path, _args(genome="GRCh38", fasta=None, gtf=None))
    assert "igenomes_ignore" not in params
    assert params.get("genome") == "GRCh38"


# ── star / gencode options ────────────────────────────────────────────────────


def test_gencode_is_written_without_featurecounts_group_type_override(tmp_path):
    params = _load_params(tmp_path, _args(gencode=True))
    assert params["gencode"] is True
    assert "featurecounts_group_type" not in params


def test_extra_aligner_args_are_written(tmp_path):
    params = _load_params(tmp_path, _args(extra_star_align_args="--limitBAMsortRAM 1", extra_bowtie2_align_args="-k 10"))
    assert params["extra_star_align_args"] == "--limitBAMsortRAM 1"
    assert params["extra_bowtie2_align_args"] == "-k 10"


# ── skip flags ────────────────────────────────────────────────────────────────


def test_skip_flags_written_only_when_true(tmp_path):
    """All skip_* flags must be absent when False and present when True."""
    skip_flags = [
        "skip_trimming",
        "skip_alignment",
        "skip_markduplicates",
        "skip_bigwig",
        "skip_fastqc",
        "skip_multiqc",
        "skip_qc",
    ]
    # None set — all absent
    params_false = _load_params(tmp_path / "false", _args(**{f: False for f in skip_flags}))
    for flag in skip_flags:
        assert flag not in params_false, f"{flag} should be absent when False"
    # All set — all present
    params_true = _load_params(tmp_path / "true", _args(**{f: True for f in skip_flags}))
    for flag in skip_flags:
        assert params_true[flag] is True, f"{flag} should be True when set"


# ── save flags ────────────────────────────────────────────────────────────────


def test_save_flags_written_only_when_true(tmp_path):
    """All save_* flags must be absent when False and present when True."""
    save_flags = [
        "save_reference",
        "save_trimmed",
        "save_align_intermeds",
        "save_unaligned",
        "save_merged_fastq",
        "save_non_ribo_reads",
        "save_umi_intermeds",
    ]
    params_false = _load_params(tmp_path / "false", _args(**{f: False for f in save_flags}))
    for flag in save_flags:
        assert flag not in params_false, f"{flag} should be absent when False"
    params_true = _load_params(tmp_path / "true", _args(**{f: True for f in save_flags}))
    for flag in save_flags:
        assert params_true[flag] is True, f"{flag} should be True when set"


# ── contaminant screening ─────────────────────────────────────────────────────


def test_contaminant_screening_params_are_written_with_upstream_names(tmp_path):
    params = _load_params(
        tmp_path,
        _args(
            contaminant_screening="kraken2_bracken",
            contaminant_screening_input="trimmed",
            kraken_db=str(tmp_path / "kraken"),
            bracken_precision="G",
            sylph_db=str(tmp_path / "sylph.db"),
            sylph_taxonomy=str(tmp_path / "taxonomy.tsv"),
            save_kraken_assignments=True,
            save_kraken_unassigned=True,
        ),
    )
    assert params["contaminant_screening"] == "kraken2_bracken"
    assert params["contaminant_screening_input"] == "trimmed"
    assert params["kraken_db"].endswith("kraken")
    assert params["bracken_precision"] == "G"
    assert params["sylph_db"].endswith("sylph.db")
    assert params["sylph_taxonomy"].endswith("taxonomy.tsv")
    assert params["save_kraken_assignments"] is True
    assert params["save_kraken_unassigned"] is True


# ── gencode_autodetected kwarg ────────────────────────────────────────────────


def test_build_effective_params_gencode_autodetected_flag(tmp_path):
    """build_effective_params honours gencode_autodetected kwarg without touching args."""
    samplesheet = _samplesheet(tmp_path / "a")
    args = _args(gencode=False)
    assert not hasattr(args, "_gencode_autodetected")

    # gencode_autodetected=True → params["gencode"] must be True
    params_with = build_effective_params(
        args,
        normalized_samplesheet=samplesheet,
        output_dir=tmp_path / "a",
        gencode_autodetected=True,
    )
    assert params_with.get("gencode") is True

    # gencode_autodetected=False → "gencode" must NOT be in params
    samplesheet2 = _samplesheet(tmp_path / "b")
    params_without = build_effective_params(
        args,
        normalized_samplesheet=samplesheet2,
        output_dir=tmp_path / "b",
        gencode_autodetected=False,
    )
    assert "gencode" not in params_without

    # args must never be mutated by either call
    assert not hasattr(args, "_gencode_autodetected"), (
        "build_effective_params must not write _gencode_autodetected onto args"
    )


# ── bbsplit ───────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("field", ["bbsplit_fasta_list"])
def test_bbsplit_index_or_fasta_enables_bbsplit(tmp_path, field):
    value = tmp_path / field
    params = _load_params(tmp_path, _args(**{field: str(value)}))
    assert params[field].endswith(field)
    assert params["skip_bbsplit"] is False


# ── deseq2_vst tri-state ──────────────────────────────────────────────────────


def test_deseq2_vst_default_omitted(tmp_path):
    params = _load_params(tmp_path, _args())
    assert "deseq2_vst" not in params


@pytest.mark.parametrize("value", [False, True])
def test_deseq2_vst_explicit_values_written(tmp_path, value):
    params = _load_params(tmp_path, _args(deseq2_vst=value))
    assert params["deseq2_vst"] is value


# ── Tier 2 passthrough flags ──────────────────────────────────────────────────


def test_umi_params_are_written_with_mixed_types(tmp_path):
    params = _load_params(
        tmp_path,
        _args(
            with_umi=True,
            umi_dedup_tool="umitools",
            umitools_bc_pattern="NNNN",
            umitools_bc_pattern2="CCCC",
            umitools_umi_separator=":",
            skip_umi_extract=True,
            save_umi_intermeds=True,
        ),
    )
    assert params["with_umi"] is True
    assert params["umi_dedup_tool"] == "umitools"
    assert params["umitools_bc_pattern"] == "NNNN"
    assert params["umitools_bc_pattern2"] == "CCCC"
    assert params["umitools_umi_separator"] == ":"
    assert params["skip_umi_extract"] is True
    assert params["save_umi_intermeds"] is True


# ── salmon / kallisto tuning ──────────────────────────────────────────────────


# ── input metadata / labels ───────────────────────────────────────────────────


# ── remote reference URIs passthrough ────────────────────────────────────────


@pytest.mark.parametrize("scheme,field", [
    ("s3", "fasta"),
    ("gs", "gtf"),
])
def test_remote_ref_uri_passed_through_unchanged(tmp_path, scheme, field):
    """params_builder must not mangle remote URIs through Path().resolve().

    Path(...).resolve() would turn 'https://example.org/genome.fa' into an
    absolute POSIX path like '/https:/example.org/genome.fa', breaking the Nextflow run.
    Reference fields must use _posix_or_uri() which passes URIs unchanged.
    """
    suffix = "gtf" if field == "gtf" else "fa"
    uri = f"{scheme}://bucket.example.org/refs/file.{suffix}"
    params = _load_params(tmp_path, _args(**{field: uri}))
    assert params.get(field) == uri, (
        f"params['{field}'] = {params.get(field)!r}; expected URI {uri!r} to be passed through unchanged"
    )


def test_remote_ribo_manifest_uri_passed_through_unchanged(tmp_path):
    uri = "s3://bucket.example.org/ribo/ribodetector_manifest.txt"
    params = _load_params(tmp_path, _args(ribo_database_manifest=uri))
    assert params["ribo_database_manifest"] == uri


@pytest.mark.parametrize("field", ["kraken_db", "bbsplit_fasta_list", "bbsplit_index"])
def test_remote_contaminant_path_uri_passed_through_unchanged(tmp_path, field):
    uri = f"s3://bucket.example.org/contaminants/{field}"
    params = _load_params(tmp_path, _args(**{field: uri}))
    assert params[field] == uri


@pytest.mark.parametrize("field", ["multiqc_config", "multiqc_logo", "multiqc_methods_description"])
def test_remote_multiqc_path_uri_passed_through_unchanged(tmp_path, field):
    uri = f"https://example.org/multiqc/{field}.yaml"
    params = _load_params(tmp_path, _args(**{field: uri}))
    assert params[field] == uri


# ── Audit L3: the default trimmer is centralised in schemas (single source) ───


def test_default_trimmer_centralised_in_schemas():
    _SKILL_DIR_LOCAL = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(_SKILL_DIR_LOCAL))
    try:
        import schemas as _schemas
        import params_builder as _pb
    finally:
        while str(_SKILL_DIR_LOCAL) in sys.path:
            sys.path.remove(str(_SKILL_DIR_LOCAL))
    assert _schemas.DEFAULT_TRIMMER == "trimgalore"
    assert _schemas.DEFAULT_TRIMMER in _schemas.SUPPORTED_TRIMMERS
    # params_builder must consume the shared constant, not a private literal.
    assert _pb._DEFAULT_TRIMMER == _schemas.DEFAULT_TRIMMER


# ── F2: igenomes_ignore on a fully-prebuilt reference route ───────────────────


def test_igenomes_ignore_set_for_prebuilt_index_route_without_fasta(tmp_path):
    """Prebuilt indices + annotation with neither --genome nor --fasta must still
    set igenomes_ignore=True so no iGenomes bundle is silently consulted."""
    params = _load_params(
        tmp_path,
        _args(
            aligner="star_salmon",
            star_index="/refs/star",
            salmon_index="/refs/salmon",
            gtf="/refs/genes.gtf",
            genome=None,
            fasta=None,
        ),
    )
    assert params.get("igenomes_ignore") is True


# ── F7: self-contained test profiles own the aligner ──────────────────────────


def test_noinput_profile_omits_aligner_when_not_user_chosen(tmp_path):
    """A self-contained nf-core test profile owns the aligner; the wrapper must not
    override it via -params-file unless the user explicitly chose one."""
    params = _load_params(tmp_path, _args(aligner="star_salmon", _noinput=True, _aligner_explicit=False))
    assert "aligner" not in params


def test_noinput_profile_writes_aligner_when_user_chose_it(tmp_path):
    """When the user explicitly passes --aligner, honour it even for test profiles."""
    params = _load_params(tmp_path, _args(aligner="star_rsem", _noinput=True, _aligner_explicit=True))
    assert params["aligner"] == "star_rsem"
