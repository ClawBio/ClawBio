from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pytest

_SKILL_DIR = Path(__file__).resolve().parent.parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))


def _purge_foreign_bare_modules(*names: str) -> None:
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


_purge_foreign_bare_modules("reporting", "schemas")

from reporting import _REPO_ROOT, _REPRO_PATH_FLAGS, build_repro_command_args, write_report, write_repro_commands, write_result
from schemas import DEFAULT_TIMEOUT_HOURS

_purge_local_modules("reporting", "schemas")
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))


def _args(tmp_path: Path, **kwargs) -> argparse.Namespace:
    defaults = dict(
        input=str(tmp_path / "samplesheet.csv"),
        output=str(tmp_path / "out"),
        aligner="star_salmon",
        pseudo_aligner=None,
        trimmer="trimgalore",
        profile="docker",
        pipeline_version="3.26.0",
        pipeline_local=None,
        resume=False,
        demo=False,
        check=False,
        run_downstream=False,
        metadata=None,
        formula=None,
        contrast=None,
        downstream_output=None,
        skip_quantification_merge=False,
        skip_trimming=False,
        skip_alignment=False,
        skip_pseudo_alignment=False,
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
        skip_qc=False,
        save_reference=False,
        save_trimmed=False,
        save_align_intermeds=False,
        save_unaligned=False,
        save_merged_fastq=False,
        save_non_ribo_reads=False,
        save_umi_intermeds=False,
        gencode=False,
        deseq2_vst=None,
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
        genome=None,
        remove_ribo_rna=False,
        ribo_removal_tool=None,
        with_umi=False,
        skip_umi_extract=False,
        umi_dedup_tool=None,
        umitools_bc_pattern=None,
        umitools_bc_pattern2=None,
        umitools_umi_separator=None,
        umitools_extract_method=None,
        umi_discard_read=None,
        umitools_grouping_method=None,
        stranded_threshold=None,
        unstranded_threshold=None,
        min_mapped_reads=None,
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
        featurecounts_feature_type=None,
        featurecounts_group_type=None,
        gtf_extra_attributes=None,
        gtf_group_features=None,
        rseqc_modules=None,
        extra_star_align_args=None,
        extra_bowtie2_align_args=None,
        star_ignore_sjdbgtf=False,
        bam_csi_index=False,
        stringtie_ignore_gtf=False,
        email=None,
        multiqc_title=None,
        multiqc_config=None,
        multiqc_logo=None,
        enable_preseq=False,
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
        rsem_extra_args=None,
        skip_downstream=False,
        nextflow_config=None,
        skip_bbsplit=False,
        igenomes_base=None,
        publish_dir_mode=None,
        email_on_fail=None,
    )
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


def _source() -> dict:
    return {"source_kind": "remote_tag", "source_ref": "3.26.0", "resolved_version": "3.26.0"}


def _preflight() -> dict:
    return {
        "java": {"version": "21.0.1"},
        "nextflow": {"version": "25.04.3"},
        "warnings": [],
    }


def _parsed(**kwargs) -> dict:
    defaults = dict(
        aligner_effective="star_salmon",
        preferred_counts_tsv="/out/star_salmon/salmon.merged.gene_counts_length_scaled.tsv",
        raw_counts_tsv="/out/star_salmon/salmon.merged.gene_counts.tsv",
        tpm_tsv="/out/star_salmon/salmon.merged.gene_tpm.tsv",
        rds_file="/out/star_salmon/salmon.merged.gene.SummarizedExperiment.rds",
        tx2gene_augmented="/out/star_salmon/salmon.merged.tx2gene_augmented.tsv",
        per_sample_quant_dirs=[],
        rsem_genes_results=[],
        multiqc_report="/out/multiqc/star_salmon/multiqc_report.html",
        pipeline_info_dir="/out/pipeline_info",
        samples_detected=2,
        handoff_available=True,
        hisat2_no_quant=False,
        skip_quantification_merge=False,
    )
    defaults.update(kwargs)
    return defaults


def test_write_report_creates_rnaseq_report(tmp_path):
    path = write_report(tmp_path, args=_args(tmp_path), pipeline_source=_source(), preflight_result=_preflight(), parsed_outputs=_parsed(), command_str="nextflow run nf-core/rnaseq")
    assert path == tmp_path / "report.md"
    assert "# nf-core/rnaseq Wrapper Report" in path.read_text(encoding="utf-8")


def test_report_header_includes_skill_version(tmp_path):
    """report.md must show the skill version, consistent with result.json/manifest.json/
    catalog.json (all 0.1.0). An empty Version field is an inconsistency."""
    from schemas import SKILL_VERSION
    write_report(tmp_path, args=_args(tmp_path), pipeline_source=_source(), preflight_result=_preflight(), parsed_outputs=_parsed(), command_str="cmd")
    content = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert f"**Version**: {SKILL_VERSION}" in content
    assert "**Version**: \n" not in content, "report Version field must not be empty"


def test_report_summary_contains_bulk_rnaseq_fields(tmp_path):
    write_report(tmp_path, args=_args(tmp_path, pseudo_aligner="salmon", trimmer="fastp"), pipeline_source=_source(), preflight_result=_preflight(), parsed_outputs=_parsed(samples_detected=3), command_str="cmd")
    content = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "Pseudo-aligner: `salmon`" in content
    assert "Trimmer: `fastp`" in content
    assert "Samples: `3`" in content


def test_report_next_steps_emit_rnaseq_de_template(tmp_path):
    write_report(tmp_path, args=_args(tmp_path), pipeline_source=_source(), preflight_result=_preflight(), parsed_outputs=_parsed(), command_str="cmd")
    content = (tmp_path / "report.md").read_text(encoding="utf-8")
    # Must use CLAWBIO_REPO anchor so the command works from any directory.
    assert '"${CLAWBIO_REPO}/clawbio.py" run rnaseq' in content
    assert "python clawbio.py run rnaseq" not in content
    assert "--counts /out/star_salmon/salmon.merged.gene_counts_length_scaled.tsv" in content
    assert "<your_metadata.csv>" in content
    assert "--formula" in content
    assert "--contrast" in content


def test_report_replay_line_includes_clawbio_repo(tmp_path):
    """The Replay line in report.md must include the CLAWBIO_REPO prefix so users
    can copy-paste it directly — commands.sh requires the env var to be set."""
    write_report(tmp_path, args=_args(tmp_path), pipeline_source=_source(), preflight_result=_preflight(), parsed_outputs=_parsed(), command_str="cmd")
    content = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "CLAWBIO_REPO=" in content
    assert "bash" in content
    assert "commands.sh" in content


def test_report_handoff_quotes_counts_path_with_spaces(tmp_path):
    write_report(
        tmp_path,
        args=_args(tmp_path),
        pipeline_source=_source(),
        preflight_result=_preflight(),
        parsed_outputs=_parsed(preferred_counts_tsv="/tmp/run with spaces/counts.tsv"),
        command_str="nextflow run nf-core/rnaseq",
    )
    content = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "--counts '/tmp/run with spaces/counts.tsv' \\" in content


def test_report_next_steps_explain_hisat2_no_handoff(tmp_path):
    parsed = _parsed(preferred_counts_tsv="", handoff_available=False, hisat2_no_quant=True)
    write_report(tmp_path, args=_args(tmp_path, aligner="hisat2"), pipeline_source=_source(), preflight_result=_preflight(), parsed_outputs=parsed, command_str="cmd")
    content = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "HISAT2" in content
    assert "rnaseq" in content


def test_write_result_summary_contains_handoff_keys(tmp_path):
    result_path = write_result(tmp_path, args=_args(tmp_path), pipeline_source=_source(), parsed_outputs=_parsed(), command_str="cmd")
    payload = json.loads(result_path.read_text(encoding="utf-8"))
    summary = payload["summary"]
    assert summary["preferred_counts_tsv"] == _parsed()["preferred_counts_tsv"]
    assert summary["handoff_available"] is True
    assert summary["hisat2_no_quant"] is False
    assert summary["skip_quantification_merge"] is False


def test_write_result_has_shared_status_ok_contract(tmp_path):
    """Issue 5: a successful result.json must carry the shared status/ok
    discriminators (consistent with nfcore-sarek/scrnaseq)."""
    result_path = write_result(tmp_path, args=_args(tmp_path), pipeline_source=_source(), parsed_outputs=_parsed(), command_str="cmd")
    payload = json.loads(result_path.read_text(encoding="utf-8"))
    assert payload["status"] == "ok"
    assert payload["ok"] is True


def test_write_result_data_contains_full_outputs(tmp_path):
    parsed = _parsed(samples_detected=7)
    result_path = write_result(tmp_path, args=_args(tmp_path), pipeline_source=_source(), parsed_outputs=parsed, command_str="cmd")
    payload = json.loads(result_path.read_text(encoding="utf-8"))
    assert payload["data"]["canonical_skill_name"] == "nfcore-rnaseq-wrapper"
    assert payload["data"]["cli_alias"] == "rnaseq-pipeline"
    assert payload["data"]["outputs"] == parsed
    assert payload["data"]["output_artifacts"] == parsed


def test_write_repro_commands_uses_rnaseq_script_name(tmp_path):
    write_repro_commands(tmp_path, args=_args(tmp_path))
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "nfcore_rnaseq_wrapper.py" in content
    assert "nfcore_scrnaseq_wrapper.py" not in content


def test_write_repro_commands_appends_portability_notice_for_real_run(tmp_path):
    write_repro_commands(tmp_path, args=_args(tmp_path, demo=False))
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "Portability notice" in content
    assert "remap_paths.py" in content


def test_write_repro_commands_skips_portability_notice_for_demo(tmp_path):
    write_repro_commands(tmp_path, args=_args(tmp_path, demo=True))
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "Portability notice" not in content


def test_portability_notice_does_not_claim_remote_inputs_are_absolute(tmp_path):
    """A run started with --allow-remote-inputs carries remote URIs, not absolute
    local paths. The portability notice must not make a blanket claim that every
    input path in the samplesheet is absolute, and it must acknowledge the remote
    case so the remap guidance is not misapplied to URLs (which resolve anywhere).
    """
    write_repro_commands(tmp_path, args=_args(tmp_path, demo=False, allow_remote_inputs=True))
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "Portability notice" in content
    # The unconditional "are absolute" assertion is false for remote inputs.
    assert "in samplesheet.valid.csv are absolute" not in content
    # The notice must cover the remote-input case explicitly.
    assert "--allow-remote-inputs" in content
    assert "remote" in content.lower()


def test_build_repro_command_args_includes_rnaseq_core_flags(tmp_path):
    command_args = build_repro_command_args(tmp_path, args=_args(tmp_path, aligner="star_rsem", profile="singularity"))
    assert command_args["--output"] == tmp_path.as_posix()
    assert command_args["--aligner"] == "star_rsem"
    assert command_args["--profile"] == "singularity"
    assert command_args["--pipeline-version"] == "3.26.0"


def test_repro_commands_sh_uses_portable_python_interpreter(tmp_path):
    """commands.sh must not hardcode a bare `python` — modern macOS and many Linux
    distros ship only `python3`, so a bare `python` replay fails with
    'python: command not found'. The bundle must use a portable interpreter
    (`${PYTHON:-python3}`) so replay works on any compatible OS / another folder."""
    write_repro_commands(tmp_path, args=_args(tmp_path))
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert '"${PYTHON:-python3}" "$SKILL_SCRIPT"' in content
    assert "\npython \"$SKILL_SCRIPT\"" not in content, "bare `python` is not portable"


def test_repro_commands_sh_has_idempotent_resume_guard(tmp_path):
    """Re-running commands.sh in place must be idempotent. The RNA-seq bundle
    re-invokes the wrapper (which re-runs preflight and rejects a non-empty output
    dir, unlike the Sarek/scRNA-seq bundles that replay Nextflow directly). The
    replay therefore guards with --resume when the target output dir already holds a
    completed run of this bundle (reproducibility/manifest.json present); a fresh or
    remapped output dir has no manifest and runs normally."""
    write_repro_commands(tmp_path, args=_args(tmp_path))
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    manifest = f"{tmp_path.as_posix()}/reproducibility/manifest.json"
    assert f'if [[ -f "{manifest}" ]]; then' in content
    assert 'CLAWBIO_RESUME="--resume"' in content
    assert '"${PYTHON:-python3}" "$SKILL_SCRIPT" $CLAWBIO_RESUME' in content


def test_repro_commands_sh_has_idempotent_resume_guard_for_demo(tmp_path):
    """A demo bundle must be replayable in place, exactly like a real-run bundle.

    A successful `--demo` run leaves its own output dir populated, so replaying
    `bash reproducibility/commands.sh` re-invokes the wrapper against a non-empty
    --output and preflight rejects it with OUTPUT_DIR_NOT_EMPTY. Nextflow's -resume
    is orthogonal to -profile test (nf-core documents no incompatibility), and the
    demo run's work dir (<output>/upstream/work) and Nextflow cache (<output>/.nextflow)
    both persist in place — so the demo bundle gets the same guard as any other."""
    write_repro_commands(tmp_path, args=_args(tmp_path, demo=True))
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    manifest = f"{tmp_path.as_posix()}/reproducibility/manifest.json"
    assert f'if [[ -f "{manifest}" ]]; then' in content
    assert 'CLAWBIO_RESUME="--resume"' in content
    assert '"${PYTHON:-python3}" "$SKILL_SCRIPT" $CLAWBIO_RESUME' in content
    assert "\n    --demo" in content


def test_repro_commands_sh_resume_guard_skipped_when_run_used_resume(tmp_path):
    """If the original run already used --resume the command always resumes, so the
    conditional guard is redundant and omitted."""
    write_repro_commands(tmp_path, args=_args(tmp_path, resume=True))
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "CLAWBIO_RESUME" not in content
    assert "\n    --resume" in content


def test_write_repro_commands_replaces_commands_sh_atomically(tmp_path):
    """An in-place --resume replay regenerates commands.sh while bash is still executing
    it. The write must be atomic (os.replace to a fresh inode), so a reader that already
    has the file open keeps seeing the ORIGINAL bytes rather than a truncated/rewritten
    file. Regenerating with different content must not disturb the open reader."""
    write_repro_commands(tmp_path, args=_args(tmp_path, aligner="star_salmon"))
    commands_sh = tmp_path / "reproducibility" / "commands.sh"
    original = commands_sh.read_text(encoding="utf-8")
    with commands_sh.open("r", encoding="utf-8") as held:
        head = held.read(20)
        # Regenerate in place with different content (simulates the wrapper re-invoked
        # by the very commands.sh under execution).
        write_repro_commands(tmp_path, args=_args(tmp_path, aligner="star_rsem"))
        rest = held.read()
    assert head + rest == original, "open reader saw a truncated/rewritten file — not atomic"
    assert "star_rsem" not in (head + rest)
    # The on-disk file is the regenerated version.
    assert "star_rsem" in commands_sh.read_text(encoding="utf-8")


def test_repro_commands_sh_bakes_clawbio_repo_default(tmp_path):
    """The bundle always lives OUTSIDE the ClawBio checkout, so walking up from the
    script never finds skills/. commands.sh must bake the generating checkout as the
    CLAWBIO_REPO default so a same-machine replay needs no manual CLAWBIO_REPO, while
    remaining overridable. The misleading "inside the repository clone" hint is gone."""
    write_repro_commands(tmp_path, args=_args(tmp_path))
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert f'REPO_ROOT="${{CLAWBIO_REPO:-{_REPO_ROOT.as_posix()}}}"' in content
    assert "inside the repository clone" not in content


def test_build_repro_command_args_omits_default_trimmer(tmp_path):
    """commands.sh must not record the default --trimmer (trimgalore) — it would
    diverge from params.yaml, which omits defaults ('omit = trust upstream default'),
    and record a flag the user never chose. Mirrors params_builder._add_aligner_params."""
    command_args = build_repro_command_args(tmp_path, args=_args(tmp_path))  # default trimmer
    assert "--trimmer" not in command_args


def test_build_repro_command_args_records_non_default_trimmer(tmp_path):
    """A user-chosen non-default trimmer must be recorded for a faithful replay."""
    command_args = build_repro_command_args(tmp_path, args=_args(tmp_path, trimmer="fastp"))
    assert command_args["--trimmer"] == "fastp"


def test_build_repro_command_args_includes_reference_paths(tmp_path):
    fasta = tmp_path / "ref.fa"
    gtf = tmp_path / "genes.gtf"
    command_args = build_repro_command_args(tmp_path, args=_args(tmp_path, fasta=str(fasta), gtf=str(gtf)))
    assert command_args["--fasta"] == fasta.resolve().as_posix()
    assert command_args["--gtf"] == gtf.resolve().as_posix()


def test_build_repro_command_args_records_allow_remote_inputs(tmp_path):
    """A run that opted into remote inputs must replay with --allow-remote-inputs.
    rnaseq's commands.sh re-invokes the wrapper (which re-runs preflight), so without
    the flag the replay fails REMOTE_INPUT_NOT_ALLOWED even though the samplesheet/
    params.yaml in the bundle carry the remote URIs. (sarek/scrnaseq replay Nextflow
    directly and never hit the gate, so this only affects rnaseq.)"""
    command_args = build_repro_command_args(tmp_path, args=_args(tmp_path, allow_remote_inputs=True))
    assert "--allow-remote-inputs" in command_args
    assert command_args["--allow-remote-inputs"] is None


def test_build_repro_command_args_omits_allow_remote_inputs_when_unset(tmp_path):
    """The flag is recorded only when the user actually opted in, mirroring every other
    boolean flag — a local-only run's bundle must not silently enable remote fetches."""
    command_args = build_repro_command_args(tmp_path, args=_args(tmp_path, allow_remote_inputs=False))
    assert "--allow-remote-inputs" not in command_args


def test_build_repro_command_args_demo_omits_input(tmp_path):
    command_args = build_repro_command_args(tmp_path, args=_args(tmp_path, demo=True))
    assert "--demo" in command_args
    assert "--input" not in command_args


def test_build_repro_command_args_real_run_uses_bundle_samplesheet(tmp_path):
    command_args = build_repro_command_args(tmp_path, args=_args(tmp_path, input=str(tmp_path / "ss.csv")))
    assert command_args["--input"] == "${SCRIPT_DIR}/samplesheet.valid.csv"


# ── P1: reporting fixes for _noinput, remote refs, --nextflow-config ──────────


def test_build_repro_command_args_noinput_includes_profile(tmp_path):
    """commands.sh for self-contained profiles must preserve --profile so replay is correct."""
    args = _args(tmp_path, input=None, demo=False, profile="docker,test_full")
    args._noinput = True
    result = build_repro_command_args(tmp_path, args=args)
    assert result.get("--profile") == "docker,test_full"


def test_build_repro_command_args_noinput_omits_aligner_when_not_explicit(tmp_path):
    """Mirror F7: a self-contained profile owns the aligner, so commands.sh must NOT
    bake in --aligner — replaying it would flip the run to an explicit aligner and
    diverge from the original (aligner-less) params.yaml."""
    args = _args(tmp_path, input=None, demo=False, profile="docker,test_full", aligner="star_salmon")
    args._noinput = True
    args._aligner_explicit = False
    result = build_repro_command_args(tmp_path, args=args)
    assert "--aligner" not in result


def test_build_repro_command_args_noinput_keeps_aligner_when_explicit(tmp_path):
    """If the user explicitly chose --aligner with a test profile, the replay keeps it."""
    args = _args(tmp_path, input=None, demo=False, profile="docker,test_full", aligner="star_rsem")
    args._noinput = True
    args._aligner_explicit = True
    result = build_repro_command_args(tmp_path, args=args)
    assert result["--aligner"] == "star_rsem"


@pytest.mark.parametrize("scheme,field", [
    ("https", "fasta"),
])
def test_build_repro_command_args_remote_refs_not_mangled(tmp_path, scheme, field):
    """Remote URI refs in commands.sh must be written unchanged, not resolved through Path.

    Path(value).resolve() turns 'https://example.org/genome.fa' into
    '/.../https:/example.org/genome.fa', breaking replay.  _REPRO_PATH_FLAGS
    fields must use URI-aware serialization.
    """
    uri = f"{scheme}://bucket.example.org/refs/file.fa"
    args = _args(tmp_path, **{field: uri})
    result = build_repro_command_args(tmp_path, args=args)
    flag = f"--{field.replace('_', '-')}"
    assert result.get(flag) == uri, (
        f"commands.sh must write remote URI {uri!r} unchanged; got {result.get(flag)!r}"
    )


def test_build_repro_command_args_sylph_multi_paths_resolved_individually(tmp_path, monkeypatch):
    db1 = tmp_path / "sylph" / "db1"
    db2 = tmp_path / "sylph" / "db2"
    db1.mkdir(parents=True)
    db2.mkdir(parents=True)
    monkeypatch.chdir(tmp_path)
    args = _args(tmp_path, sylph_db="sylph/db1,sylph/db2")
    result = build_repro_command_args(tmp_path, args=args)
    assert result["--sylph-db"] == f"{db1.resolve().as_posix()},{db2.resolve().as_posix()}"


def test_build_repro_command_args_nextflow_config_included(tmp_path):
    """--nextflow-config configs are copied into the bundle and replayed via
    ${SCRIPT_DIR}, so a run with --nextflow-config hpc.config --nextflow-config
    rsem.config replays portably (not via host-specific absolute paths)."""
    cfg1 = tmp_path / "hpc.config"
    cfg2 = tmp_path / "rsem.config"
    cfg1.write_text("process { cpus = 2 }\n", encoding="utf-8")
    cfg2.write_text("params { x = 1 }\n", encoding="utf-8")
    args = _args(tmp_path, nextflow_config=[str(cfg1), str(cfg2)])
    result = build_repro_command_args(tmp_path, args=args)
    assert "--nextflow-config" in result, "--nextflow-config must appear in repro command args"
    configs = result["--nextflow-config"]
    assert configs == [
        "${SCRIPT_DIR}/nextflow_configs/config_01_hpc.config",
        "${SCRIPT_DIR}/nextflow_configs/config_02_rsem.config",
    ]
    # The files were actually copied into the bundle, verbatim.
    staged_dir = tmp_path / "reproducibility" / "nextflow_configs"
    assert (staged_dir / "config_01_hpc.config").read_text() == "process { cpus = 2 }\n"
    assert (staged_dir / "config_02_rsem.config").read_text() == "params { x = 1 }\n"


def test_build_repro_command_args_stages_config_outside_output_dir(tmp_path):
    """The reported bug: a -c config living OUTSIDE the output directory must still
    replay out-of-the-box. It is copied into reproducibility/nextflow_configs/ and
    referenced via ${SCRIPT_DIR}, so commands.sh bakes no absolute host path and no
    <EDIT_ME> placeholder."""
    external = tmp_path / "elsewhere" / "limits.config"
    external.parent.mkdir()
    external.write_text("process { resourceLimits = [ memory: '12.GB' ] }\n", encoding="utf-8")
    output_dir = tmp_path / "run"
    output_dir.mkdir()
    result = build_repro_command_args(output_dir, args=_args(output_dir, nextflow_config=[str(external)]))
    configs = result["--nextflow-config"]
    assert configs == ["${SCRIPT_DIR}/nextflow_configs/config_01_limits.config"]
    for value in configs:
        assert "<EDIT_ME>" not in value
        assert not value.startswith("/"), "must not bake a host-specific absolute path"
    copied = output_dir / "reproducibility" / "nextflow_configs" / "config_01_limits.config"
    assert copied.read_text() == "process { resourceLimits = [ memory: '12.GB' ] }\n"


def test_build_repro_command_args_config_restage_is_idempotent(tmp_path):
    """On an in-place --resume replay, commands.sh re-invokes the wrapper with
    --nextflow-config pointing at the ALREADY-staged
    ${SCRIPT_DIR}/nextflow_configs/config_01_*.config (a file inside THIS bundle).
    Re-staging must reference it in place, not copy it again under a new config_NN_
    prefix — otherwise successive replays accumulate config_01_config_01_..."""
    external = tmp_path / "elsewhere" / "limits.config"
    external.parent.mkdir()
    external.write_text("process { resourceLimits = [ memory: '12.GB' ] }\n", encoding="utf-8")
    output_dir = tmp_path / "run"
    output_dir.mkdir()
    first = build_repro_command_args(output_dir, args=_args(output_dir, nextflow_config=[str(external)]))
    assert first["--nextflow-config"] == ["${SCRIPT_DIR}/nextflow_configs/config_01_limits.config"]
    staged = output_dir / "reproducibility" / "nextflow_configs" / "config_01_limits.config"
    assert staged.is_file()
    # Replay: the wrapper is re-invoked with the already-staged (absolute) path.
    second = build_repro_command_args(output_dir, args=_args(output_dir, nextflow_config=[str(staged)]))
    assert second["--nextflow-config"] == ["${SCRIPT_DIR}/nextflow_configs/config_01_limits.config"]
    config_dir = output_dir / "reproducibility" / "nextflow_configs"
    names = sorted(p.name for p in config_dir.iterdir())
    assert names == ["config_01_limits.config"], f"double-prefixed copy created: {names}"


def test_build_repro_command_args_remote_config_uri_passed_through(tmp_path):
    """A remote config URI is not a local file, so it is replayed verbatim (not copied)."""
    result = build_repro_command_args(
        tmp_path, args=_args(tmp_path, nextflow_config=["https://example.org/hpc.config"])
    )
    assert result["--nextflow-config"] == ["https://example.org/hpc.config"]


def test_commands_sh_replays_every_run_affecting_wrapper_flag(tmp_path):
    """Replay fidelity: every wrapper flag that affects the run must be emitted by
    build_repro_command_args, so a newly added flag cannot be silently dropped from
    commands.sh and make the reproducibility bundle an unfaithful replay."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "rnaseq_wrapper_repro_parity", _SKILL_DIR / "nfcore_rnaseq_wrapper.py"
    )
    wrapper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(wrapper)
    parser = wrapper.build_parser()

    ns = argparse.Namespace()
    wrapper_flags: set[str] = set()
    for action in parser._actions:
        if not action.option_strings or "-h" in action.option_strings:
            continue
        wrapper_flags.update(action.option_strings)
        dest = action.dest
        if action.nargs == 0:
            setattr(ns, dest, True)
        elif action.choices:
            setattr(ns, dest, sorted(str(c) for c in action.choices)[0])
        elif action.type is int:
            setattr(ns, dest, 7)
        elif action.type is float:
            setattr(ns, dest, 0.7)
        else:
            setattr(ns, dest, f"VAL_{dest}")
    ns.demo = False
    ns._noinput = False
    ns._aligner_explicit = True
    ns.timeout_hours = 99.0  # != pinned default so it is recorded
    ns.deseq2_vst = True
    ns.nextflow_config = ["c.config"]
    ns.resume = True

    emitted = set(build_repro_command_args(tmp_path, args=ns).keys())
    # Handled specially / intentionally not replayed verbatim:
    #  --input/--output → emitted in bundle-relative form by dedicated logic;
    #  --demo → off in this real-run scenario; --check → means "do not run";
    #  --no-deseq2-vst → alternate store_false form of --deseq2-vst (mutually exclusive);
    #  -v/--verbose/--no-banner → launcher-only logging/display, not run-affecting;
    #  -c/--config → aliases of --nextflow-config, replayed under that canonical name.
    special = {
        "--input",
        "--output",
        "--demo",
        "--check",
        "--no-deseq2-vst",
        "-v",
        "--verbose",
        "--no-banner",
        "-c",
        "--config",
        "--work-dir",
        "--allow-remote-inputs",
    }
    missing = wrapper_flags - emitted - special
    assert missing == set(), f"wrapper flags not replayed in commands.sh: {sorted(missing)}"


def test_repro_path_flags_are_all_remappable():
    """Every reference/path flag reporting writes into commands.sh must be one that
    remap_paths.py can rewrite (_REFERENCE_FLAGS). Otherwise moving the bundle to a
    new machine silently leaves that path pointing at the old location."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "rnaseq_remap_paths_parity", _SKILL_DIR / "remap_paths.py"
    )
    remap = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(remap)
    # reporting stores fields with underscores; commands.sh / remap use dashes.
    reporting_dash = {field.replace("_", "-") for field in _REPRO_PATH_FLAGS}
    unremappable = reporting_dash - set(remap._REFERENCE_FLAGS)
    assert unremappable == set(), (
        f"reporting writes path flags remap_paths.py cannot rewrite: {sorted(unremappable)}"
    )


# ── P1: four flags missing from build_repro_command_args ──────────────────────


def test_build_repro_command_args_includes_igenomes_base(tmp_path):
    base = tmp_path / "igenomes"
    result = build_repro_command_args(tmp_path, args=_args(tmp_path, igenomes_base=str(base)))
    assert "--igenomes-base" in result
    assert str(base.resolve()) in result["--igenomes-base"]


# ── Audit F-12: --timeout-hours / --allow-pipeline-version-override in replay ──


def test_repro_command_includes_non_default_timeout_hours(tmp_path):
    result = build_repro_command_args(tmp_path, args=_args(tmp_path, timeout_hours=48.0))
    assert result.get("--timeout-hours") == "48.0"


def test_repro_command_omits_default_timeout_hours(tmp_path):
    result = build_repro_command_args(tmp_path, args=_args(tmp_path, timeout_hours=float(DEFAULT_TIMEOUT_HOURS)))
    assert "--timeout-hours" not in result


def test_repro_command_records_pipeline_version_override(tmp_path):
    result = build_repro_command_args(tmp_path, args=_args(tmp_path, allow_pipeline_version_override=True))
    assert "--allow-pipeline-version-override" in result


def test_repro_command_omits_version_override_when_unset(tmp_path):
    result = build_repro_command_args(tmp_path, args=_args(tmp_path))
    assert "--allow-pipeline-version-override" not in result
    assert "--timeout-hours" not in result


def test_strandedness_summary_reads_from_strandedness_counts(tmp_path):
    preflight = _preflight()
    preflight["samplesheet"] = {"strandedness_counts": {"forward": 2, "reverse": 1}}
    write_report(tmp_path, args=_args(tmp_path), pipeline_source=_source(), preflight_result=preflight, parsed_outputs=_parsed(), command_str="cmd")
    content = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "forward: 2" in content
    assert "reverse: 1" in content


# ── Task 3: Verify _noinput guard ensures test-profile runs omit --input ──────


def test_repro_command_omits_input_for_noinput_runs(tmp_path):
    """build_repro_command_args must omit --input when args._noinput is True.

    Test-profile runs (test_full, test_prokaryotic, etc.) set _noinput=True.
    The replay command must not include --input so re-running the profile
    uses the bundled test data, not an external samplesheet.
    """
    args = _args(tmp_path, demo=False)
    args._noinput = True  # simulates test-profile run

    command_args = build_repro_command_args(tmp_path, args=args)

    assert "--input" not in command_args, (
        "--input must be absent from replay command for _noinput runs; "
        f"got keys: {list(command_args.keys())}"
    )
    assert "--demo" not in command_args, "--demo must not be set for _noinput runs"


# ── Step 3: post_run_warnings surfaced in report.md and result.json ──────────

def test_write_report_renders_post_run_warnings_section(tmp_path):
    """When post_run_warnings is non-empty, report.md must contain a
    '## Post-run Warnings' section listing each warning."""
    warnings = ["Provenance bundle could not be written (OSError: no space left)."]
    write_report(
        tmp_path,
        args=_args(tmp_path),
        pipeline_source=_source(),
        preflight_result=_preflight(),
        parsed_outputs=_parsed(),
        command_str="cmd",
        post_run_warnings=warnings,
    )
    content = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "## Post-run Warnings" in content, "report.md must include a Post-run Warnings section"
    assert "Provenance bundle" in content or "no space left" in content, (
        "report.md must include the warning text"
    )


def test_write_result_includes_post_run_warnings_in_json(tmp_path):
    """When post_run_warnings is non-empty, result.json summary must include
    a 'post_run_warnings' key listing the warnings."""
    warnings = ["Provenance bundle could not be written (RuntimeError: disk full)."]
    write_result(
        tmp_path,
        args=_args(tmp_path),
        pipeline_source=_source(),
        parsed_outputs=_parsed(),
        command_str="cmd",
        post_run_warnings=warnings,
    )
    result = json.loads((tmp_path / "result.json").read_text(encoding="utf-8"))
    stored = result.get("summary", {}).get("post_run_warnings", None)
    assert stored is not None, "result.json summary must have 'post_run_warnings' key"
    assert len(stored) == 1
    assert "disk full" in stored[0] or "Provenance" in stored[0]


def test_replay_line_is_shell_safe_for_apostrophe_path(tmp_path):
    """report.md replay line must survive output dirs containing single quotes."""
    tricky_out = tmp_path / "it's" / "out"
    tricky_out.mkdir(parents=True)
    write_report(
        tricky_out,
        args=_args(tmp_path, output=str(tricky_out)),
        pipeline_source=_source(),
        preflight_result=_preflight(),
        parsed_outputs=_parsed(),
        command_str="nextflow run nf-core/rnaseq",
    )
    content = (tricky_out / "report.md").read_text(encoding="utf-8")
    commands_sh_path = str(tricky_out / "reproducibility" / "commands.sh")
    # The old broken form must not appear verbatim.
    assert f"bash '{commands_sh_path}'" not in content
    # The replay line must still reference commands.sh somehow.
    assert "commands.sh" in content
