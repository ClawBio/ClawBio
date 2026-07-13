from __future__ import annotations

from argparse import Namespace
import json
from pathlib import Path
import sys

import pytest

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


_purge_foreign_modules("errors", "schemas", "preflight")

import preflight
from errors import ErrorCode, SkillError
from schemas import SUPPORTED_IGENOMES_NAMES

_purge_local_modules("errors", "schemas", "preflight")
_remove_skill_dir_from_sys_path()


_PIPELINE_SOURCE = {
    "source_kind": "local_checkout",
    "source_ref": "/repo/rnaseq",
    "resolved_version": "abc123",
    "branch": "main",
    "dirty": False,
}
_DEFAULT_REMOTE_PIPELINE_SOURCE = {
    "source_kind": "remote_repo",
    "source_ref": "nf-core/rnaseq",
    "resolved_version": "3.26.0",
}


def _touch(path: Path, text: str = "x") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _args(tmp_path: Path, **overrides) -> Namespace:
    fasta = _touch(tmp_path / "refs" / "genome.fa", ">chr1\nACGT\n")
    gtf = _touch(tmp_path / "refs" / "genes.gtf", "chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"g1\";\n")
    defaults = dict(
        output=str(tmp_path / "out"),
        resume=False,
        demo=False,
        allow_remote_inputs=False,
        aligner="star_salmon",
        pseudo_aligner=None,
        trimmer="trimgalore",
        remove_ribo_rna=False,
        ribo_removal_tool=None,
        umi_dedup_tool=None,
        with_umi=False,
        umitools_bc_pattern=None,
        skip_umi_extract=False,
        profile="docker",
        email=None,
        genome=None,
        fasta=str(fasta),
        gtf=str(gtf),
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
        sortmerna_index=None,
        extra_star_align_args=None,
        extra_bowtie2_align_args=None,
        skip_alignment=False,
        run_downstream=False,
        prokaryotic=False,
        contaminant_screening=None,
        contaminant_screening_input=None,
        bbsplit_fasta_list=None,
        bbsplit_index=None,
        skip_bbsplit=False,
        save_bbsplit_reads=False,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


def _samplesheet(tmp_path: Path, **overrides) -> dict[str, object]:
    fastq = _touch(tmp_path / "reads" / "sample_R1.fastq.gz")
    summary = {
        "sample_count": 1,
        "sample_names": ["sampleA"],
        "fastq_paths": [fastq],
        "bam_paths": [],
        "strandedness_counts": {"auto": 1},
        "unknown_columns": [],
    }
    summary.update(overrides)
    return summary


def _mock_env(monkeypatch):
    # Pin the platform so host-OS-specific warnings (e.g. the macOS/VirtioFS /tmp
    # warning) do not make these preflight assertions non-deterministic across
    # runners. Tests that exercise darwin behaviour override this afterwards.
    monkeypatch.setattr(preflight.sys, "platform", "linux")
    monkeypatch.setattr(preflight, "_check_java", lambda: {"path": "/usr/bin/java", "version": "17.0.8"})
    monkeypatch.setattr(preflight, "_check_nextflow", lambda: {"path": "/usr/bin/nextflow", "version": "25.04.3"})
    monkeypatch.setattr(preflight, "_check_profile", lambda profile: {"profile": profile, "backend_path": "/usr/bin/docker", "backend_ready": True})


def _run(args: Namespace, samplesheet: dict[str, object], pipeline_source: dict[str, object] | None = None):
    return preflight.run_preflight(args, pipeline_source=pipeline_source or _PIPELINE_SOURCE, samplesheet_summary=samplesheet)


def test_preflight_happy_path(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    result = _run(_args(tmp_path), _samplesheet(tmp_path))
    assert result["ok"] is True
    assert result["references"]["fasta"].endswith("genome.fa")
    assert result["samplesheet"]["sample_count"] == 1
    assert result["warnings"] == []


def test_bam_reprocessing_emits_aligner_match_reminder(tmp_path, monkeypatch):
    """nf-core cannot mix quantifier types between BAM generation and reprocessing;
    the wrapper can't infer a BAM's origin, so preflight must warn whenever BAM
    columns are present (the default --aligner star_salmon silently mismatches RSEM BAMs)."""
    _mock_env(monkeypatch)
    bam = _touch(tmp_path / "bams" / "sampleA.markdup.sorted.bam")
    summary = _samplesheet(tmp_path, bam_paths=[bam])
    result = _run(_args(tmp_path, skip_alignment=True), summary)
    assert any("BAM reprocessing" in w and "aligner" in w.lower() for w in result["warnings"])


def test_no_bam_reprocessing_reminder_for_fastq_runs(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    result = _run(_args(tmp_path), _samplesheet(tmp_path))
    assert not any("BAM reprocessing" in w for w in result["warnings"])


def test_rejects_invalid_aligner_before_environment_checks(tmp_path, monkeypatch):
    monkeypatch.setattr(preflight, "_check_java", lambda: (_ for _ in ()).throw(AssertionError("java should not run")))
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, aligner="bad"), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_ALIGNER


def test_parse_version_tuple_java(monkeypatch=None):
    assert preflight._parse_version_tuple('openjdk version "17.0.8"') == (17, 0, 8)


def test_java_version_too_old(monkeypatch):
    monkeypatch.setattr(preflight.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(preflight, "_command_output", lambda args: 'openjdk version "11.0.1"')
    with pytest.raises(SkillError) as exc:
        preflight._check_java()
    assert exc.value.error_code == ErrorCode.JAVA_VERSION_TOO_OLD


def test_nextflow_version_too_old(monkeypatch):
    monkeypatch.setattr(preflight.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(preflight, "_command_output", lambda args: "Nextflow version 24.10.0")
    with pytest.raises(SkillError) as exc:
        preflight._check_nextflow()
    assert exc.value.error_code == ErrorCode.NEXTFLOW_VERSION_TOO_OLD


def test_nextflow_version_preserves_zero_padded_month(monkeypatch):
    """A detected version like '26.04.3' must be reported verbatim, never
    reconstructed from the comparison tuple as '26.4.3' (not a real release)."""
    monkeypatch.setattr(preflight.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(preflight, "_command_output", lambda args: "Nextflow version 26.04.3 build 5928")
    assert preflight._check_nextflow()["version"] == "26.04.3"


def test_java_version_preserves_exact_string(monkeypatch):
    monkeypatch.setattr(preflight.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(preflight, "_command_output", lambda args: 'openjdk version "17.0.10" 2024-01-16')
    assert preflight._check_java()["version"] == "17.0.10"


def test_check_mode_does_not_execute_nextflow_version(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(preflight, "_check_java", lambda: {"path": "/usr/bin/java", "version": "17"})
    monkeypatch.setattr(preflight, "_check_profile", lambda profile: {"profile": profile, "backend_path": "", "backend_ready": True})
    monkeypatch.setattr(preflight.shutil, "which", lambda name: f"/usr/bin/{name}")

    def fake_command_output(cmd):
        calls.append(cmd)
        if cmd[0] == "nextflow":
            raise AssertionError("check mode must not execute nextflow")
        return ""

    monkeypatch.setattr(preflight, "_command_output", fake_command_output)
    result = preflight.run_preflight(
        _args(tmp_path, check=True, demo=True, profile="test"),
        pipeline_source={"source_kind": "remote", "source_ref": "3.26.0"},
        samplesheet_summary={"sample_count": 0, "fastq_paths": []},
    )
    assert result["nextflow"]["version_checked"] is False
    assert not any(cmd and cmd[0] == "nextflow" for cmd in calls)


def test_output_dir_not_empty_rejected_without_resume(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    out = tmp_path / "out"
    out.mkdir()
    (out / "result.json").write_text("{}", encoding="utf-8")
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, output=str(out)), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.OUTPUT_DIR_NOT_EMPTY


def test_output_dir_incomplete_prior_run_gives_actionable_error(tmp_path, monkeypatch):
    """A dir from a prior FAILED run (result.json present, no manifest) cannot be
    resumed (no manifest) yet blocks a fresh run (not empty). Preflight must detect
    this and tell the user precisely what to do, instead of the generic message (F8)."""
    _mock_env(monkeypatch)
    out = tmp_path / "out"
    out.mkdir()
    (out / "result.json").write_text('{"ok": false}', encoding="utf-8")
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, output=str(out)), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.OUTPUT_DIR_NOT_EMPTY
    assert exc.value.details.get("prior_run_incomplete") is True
    assert "delete" in exc.value.fix.lower()


def test_output_dir_with_only_reproducibility_bundle_is_accepted(tmp_path):
    # The reproducibility/ bundle is entirely wrapper-generated (commands.sh,
    # checksums.sha256, environment.yml, remap_paths.py, provenance JSON, …), never
    # user data, so an output dir containing only it (e.g. a copied bundle) must NOT
    # trip OUTPUT_DIR_NOT_EMPTY. Otherwise a bundle-only dir cannot be re-run without
    # --resume (Catch-22), and the wrappers diverge — nfcore-sarek already accepts it.
    out = tmp_path / "out"
    repro = out / "reproducibility"
    repro.mkdir(parents=True)
    for f in (
        "params.yaml", "samplesheet.valid.csv", "manifest.json", "commands.sh",
        "checksums.sha256", "environment.yml", "remap_paths.py", "inputs.json", "run.json",
    ):
        (repro / f).write_text("x", encoding="utf-8")
    # Must not raise for a fresh (non-resume) invocation.
    preflight.check_output_dir_available(out, resume=False)


def test_output_dir_inside_repo_is_rejected(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    inside = _SKILL_DIR / "tmp-test-output"
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, output=str(inside)), _samplesheet(tmp_path))
    # A dedicated, self-describing code (matching nfcore-sarek/scrnaseq) rather than
    # the generic OUTPUT_DIR_NOT_WRITABLE, which wrongly implies a permissions issue.
    assert exc.value.error_code == ErrorCode.OUTPUT_DIR_INSIDE_REPO


def test_reference_accepts_genome_shortcut(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    result = _run(_args(tmp_path, genome="GRCh38", fasta=None, gtf=None), _samplesheet(tmp_path))
    assert result["references"]["genome"] == "GRCh38"


def test_supported_igenomes_names_match_audited_3260_keys():
    assert SUPPORTED_IGENOMES_NAMES == frozenset(
        {
            "GRCh37",
            "GRCh38",
            "CHM13",
            "GRCm38",
            "TAIR10",
            "EB2",
            "UMD3.1",
            "WBcel235",
            "CanFam3.1",
            "GRCz10",
            "BDGP6",
            "EquCab2",
            "EB1",
            "Galgal4",
            "Gm01",
            "Mmul_1",
            "IRGSP-1.0",
            "CHIMP2.1.4",
            "Rnor_5.0",
            "Rnor_6.0",
            "R64-1-1",
            "EF2",
            "Sbi1",
            "Sscrofa10.2",
            "AGPv3",
            "hg38",
            "hg19",
            "mm10",
            "bosTau8",
            "ce10",
            "canFam3",
            "danRer10",
            "dm6",
            "equCab2",
            "galGal4",
            "panTro4",
            "rn6",
            "sacCer3",
            "susScr3",
        }
    )


def test_unknown_genome_warns_but_does_not_reject(tmp_path, monkeypatch):
    """Non-iGenomes names should warn but allow execution so custom genome catalogues work."""
    _mock_env(monkeypatch)
    result = _run(
        _args(tmp_path, genome="grch37", fasta=None, gtf=None),
        _samplesheet(tmp_path),
        pipeline_source=_DEFAULT_REMOTE_PIPELINE_SOURCE,
    )
    assert result["ok"] is True
    assert any("grch37" in w or "iGenomes" in w for w in result["warnings"])


def test_gencode_gtf_autodetect_sets_marker_and_warns(tmp_path, monkeypatch):
    args = _args(tmp_path, gencode=False)
    Path(args.gtf).write_text(
        '#comment\nchr1\tHAVANA\tgene\t1\t10\t.\t+\t.\tgene_id "ENSG1"; gene_type "protein_coding"; havana_gene "OTTHUMG1";\n',
        encoding="utf-8",
    )
    _mock_env(monkeypatch)

    result = _run(args, _samplesheet(tmp_path))

    assert result["gencode_autodetected"] is True
    assert not hasattr(args, "_gencode_autodetected")
    assert any("gencode" in warning.lower() for warning in result["warnings"])


def test_reference_accepts_fasta_gff_group(tmp_path, monkeypatch):
    gff = _touch(tmp_path / "refs" / "genes.gff")
    _mock_env(monkeypatch)
    result = _run(_args(tmp_path, gtf=None, gff=str(gff)), _samplesheet(tmp_path))
    assert result["references"]["gff"].endswith("genes.gff")


def test_reference_rejects_missing_group(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, fasta=None, gtf=None, gff=None), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.MISSING_REFERENCE


def test_reference_rejects_genome_with_explicit_reference(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, genome="GRCh38"), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.CONFLICTING_REFERENCES


def test_reference_rejects_genome_combined_with_index_shortcuts(tmp_path, monkeypatch):
    field = "star_index"
    index_path = _touch(tmp_path / "refs" / field / "index")
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(
            _args(tmp_path, genome="GRCh38", fasta=None, gtf=None, **{field: str(index_path)}),
            _samplesheet(tmp_path),
        )
    assert exc.value.error_code == ErrorCode.CONFLICTING_REFERENCES
    assert field in exc.value.details.get("explicit_set", [])


def test_reference_path_not_found_uses_specific_error_code(tmp_path, monkeypatch):
    field = "fasta"
    _mock_env(monkeypatch)
    args = _args(tmp_path, **{field: str(tmp_path / "missing.fa")})
    with pytest.raises(SkillError) as exc:
        _run(args, _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.REFERENCE_PATH_NOT_FOUND
    assert exc.value.details["field"] == field


def test_transcript_and_additional_fasta_are_mutually_exclusive(tmp_path, monkeypatch):
    tf = _touch(tmp_path / "refs" / "tx.fa")
    af = _touch(tmp_path / "refs" / "additional.fa")
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, transcript_fasta=str(tf), additional_fasta=str(af)), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.CONFLICTING_FASTA_ARGS


def test_skip_alignment_requires_pseudo_aligner(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, skip_alignment=True), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION


def test_skip_alignment_and_skip_pseudo_alignment_is_preprocessing_only_mode(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    result = _run(
        _args(tmp_path, skip_alignment=True, skip_pseudo_alignment=True),
        _samplesheet(tmp_path),
    )
    assert result["ok"] is True


def test_with_umi_requires_pattern_unless_extract_is_skipped(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, with_umi=True), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION

    assert _run(_args(tmp_path / "ok", with_umi=True, skip_umi_extract=True), _samplesheet(tmp_path / "ok"))["ok"] is True


def test_ribo_database_manifest_must_exist_when_set(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, ribo_database_manifest=str(tmp_path / "missing_manifest.txt")), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.REFERENCE_PATH_NOT_FOUND
    assert exc.value.details["field"] == "ribo_database_manifest"


def test_rsem_extra_args_emits_no_effect_warning(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    result = _run(_args(tmp_path, rsem_extra_args="--estimate-rspd"), _samplesheet(tmp_path))
    assert any("extra_rsem_quant_args" in warning for warning in result["warnings"])


def test_hisat2_downstream_emits_warning_and_disables_handoff(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    result = _run(_args(tmp_path, aligner="hisat2", run_downstream=True), _samplesheet(tmp_path))
    assert result["handoff_available"] is False
    assert any("hisat2" in warning.lower() for warning in result["warnings"])


def test_samplesheet_fastq_paths_are_revalidated(tmp_path, monkeypatch):
    bad = _touch(tmp_path / "reads" / "sample R1.fastq.gz")
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path), _samplesheet(tmp_path, fastq_paths=[bad]))
    assert exc.value.error_code == ErrorCode.INVALID_FASTQ


def test_samplesheet_remote_fastq_uri_not_existence_checked(tmp_path, monkeypatch):
    """With --allow-remote-inputs, a remote FASTQ URI passes preflight (existence is
    deferred to Nextflow). Without the flag it is rejected (see the gate tests)."""
    _mock_env(monkeypatch)
    result = _run(
        _args(tmp_path, allow_remote_inputs=True),
        _samplesheet(tmp_path, fastq_paths=["s3://bucket.example.org/reads/sample_R1.fastq.gz"]),
    )
    assert result["ok"] is True


def test_bam_reprocessing_requires_skip_alignment(tmp_path, monkeypatch):
    bam = _touch(tmp_path / "sample.bam")
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path), _samplesheet(tmp_path, bam_paths=[bam]))
    assert exc.value.error_code == ErrorCode.BAM_REPROCESSING_INCOMPLETE


def test_bam_reprocessing_skip_alignment_no_reference_is_allowed(tmp_path, monkeypatch):
    """Official nf-core/rnaseq 3.26.0 docs show --skip_alignment with BAMs and NO reference.

    The pipeline handles reference-less BAM reprocessing natively (e.g. using
    transcript info from the BAM header). The wrapper must not block this.
    See: https://nf-co.re/rnaseq/3.26.0/docs/usage/ (BAM input for reprocessing).
    """
    bam = _touch(tmp_path / "sample.bam")
    _mock_env(monkeypatch)
    # Explicitly clear all reference fields that _args provides by default
    result = _run(
        _args(tmp_path, skip_alignment=True, fasta=None, gtf=None, genome=None),
        _samplesheet(tmp_path, bam_paths=[bam]),
    )
    assert result["samplesheet"]["bam_count"] == 1


def test_remote_bam_uri_skip_alignment_not_existence_checked(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    result = _run(
        _args(tmp_path, skip_alignment=True, fasta=None, gtf=None, genome=None, allow_remote_inputs=True),
        _samplesheet(tmp_path, bam_paths=["s3://bucket.example.org/bams/sample.genome.bam"]),
    )
    assert result["samplesheet"]["bam_count"] == 1


def test_bam_paths_must_exist(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(
            _args(tmp_path, skip_alignment=True, pseudo_aligner="salmon"),
            _samplesheet(tmp_path, bam_paths=[tmp_path / "missing.bam"]),
        )
    assert exc.value.error_code == ErrorCode.REFERENCE_PATH_NOT_FOUND


def test_resume_requires_manifest(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, resume=True), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_RESUME_STATE


def test_resume_rejects_incompatible_manifest_fields(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    out = tmp_path / "out"
    repro = out / "reproducibility"
    repro.mkdir(parents=True)
    manifest = {
        "aligner": "hisat2",  # differs from star_salmon in args
        "pseudo_aligner": None,
        "profile": "docker",
        "prokaryotic": False,
        "pipeline_source": {"source_kind": "local_checkout", "resolved_version": "abc123"},
    }
    (repro / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, resume=True, output=str(out)), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_RESUME_STATE


def test_resume_accepts_matching_manifest(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    out = tmp_path / "out"
    repro = out / "reproducibility"
    repro.mkdir(parents=True)
    (repro / "manifest.json").write_text(
        json.dumps(
            {
                "aligner": "star_salmon",
                "pseudo_aligner": None,
                "profile": "docker",
                "prokaryotic": False,
                "pipeline_source": {"source_kind": "local_checkout", "resolved_version": "abc123"},
            }
        ),
        encoding="utf-8",
    )
    assert _run(_args(tmp_path, resume=True, output=str(out)), _samplesheet(tmp_path))["ok"] is True


def _write_manifest(out: Path, manifest: dict) -> None:
    repro = out / "reproducibility"
    repro.mkdir(parents=True, exist_ok=True)
    (repro / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")


_MATCHING_MANIFEST = {
    "aligner": "star_salmon",
    "pseudo_aligner": None,
    "profile": "docker",
    "prokaryotic": False,
    "pipeline_source": {"source_kind": "local_checkout", "resolved_version": "abc123"},
}


def test_resume_rejects_arm_flag_drift_user_added(tmp_path, monkeypatch):
    """Adding --arm on a resume from a manifest that ran without --arm must fail.

    arm changes the generated .nextflow_macos_docker.config (presence of
    --platform linux/amd64). Letting the run continue would silently swap
    container architecture between the original and resumed steps.
    """
    _mock_env(monkeypatch)
    out = tmp_path / "out"
    _write_manifest(out, {**_MATCHING_MANIFEST, "arm": False})
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, resume=True, output=str(out), arm=True), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_RESUME_STATE
    assert exc.value.details["field"] == "arm"


def test_resume_rejects_arm_flag_drift_user_removed(tmp_path, monkeypatch):
    """Removing --arm on a resume from a manifest that ran with --arm must fail."""
    _mock_env(monkeypatch)
    out = tmp_path / "out"
    _write_manifest(out, {**_MATCHING_MANIFEST, "arm": True})
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, resume=True, output=str(out), arm=False), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_RESUME_STATE
    assert exc.value.details["field"] == "arm"


def test_resume_accepts_matching_arm_flag(tmp_path, monkeypatch):
    """Resuming with --arm against a manifest where arm=true must succeed."""
    _mock_env(monkeypatch)
    out = tmp_path / "out"
    _write_manifest(out, {**_MATCHING_MANIFEST, "arm": True})
    assert _run(_args(tmp_path, resume=True, output=str(out), arm=True), _samplesheet(tmp_path))["ok"] is True


def test_resume_legacy_manifest_without_arm_treats_as_false(tmp_path, monkeypatch):
    """Legacy manifests (pre-arm key) must resume successfully when the new run
    also has arm disabled. Backward compatibility for runs created before
    the arm field was added to the manifest schema.
    """
    _mock_env(monkeypatch)
    out = tmp_path / "out"
    # _MATCHING_MANIFEST intentionally omits the 'arm' key — same as old manifests.
    _write_manifest(out, _MATCHING_MANIFEST)
    assert _run(_args(tmp_path, resume=True, output=str(out), arm=False), _samplesheet(tmp_path))["ok"] is True


def test_resume_legacy_manifest_without_arm_rejects_arm_true(tmp_path, monkeypatch):
    """Legacy manifests without arm must reject resumes that explicitly add --arm."""
    _mock_env(monkeypatch)
    out = tmp_path / "out"
    _write_manifest(out, _MATCHING_MANIFEST)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, resume=True, output=str(out), arm=True), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_RESUME_STATE
    assert exc.value.details["field"] == "arm"


def test_resume_rejects_demo_drift_against_real_manifest(tmp_path, monkeypatch):
    """Resuming a real run's output dir with --demo must fail with a clear error.

    --demo composes the upstream `test` profile, which supplies its own samplesheet
    AND its own bundled references. Attaching that to a real run's manifest would
    swap both inputs and references underneath the resume, so the demo flag is part
    of the resume compatibility contract — exactly like aligner/profile/arm.
    """
    _mock_env(monkeypatch)
    out = tmp_path / "out"
    _write_manifest(out, {**_MATCHING_MANIFEST, "demo": False})
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, resume=True, output=str(out), demo=True), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_RESUME_STATE
    assert exc.value.details["field"] == "demo"


def test_resume_rejects_real_run_against_demo_manifest(tmp_path, monkeypatch):
    """The mirror case: a real run must not resume a demo output dir."""
    _mock_env(monkeypatch)
    out = tmp_path / "out"
    _write_manifest(out, {**_MATCHING_MANIFEST, "demo": True})
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, resume=True, output=str(out), demo=False), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_RESUME_STATE
    assert exc.value.details["field"] == "demo"


def test_resume_accepts_matching_demo_flag(tmp_path, monkeypatch):
    """Replaying a demo bundle in place (bash reproducibility/commands.sh) resumes the
    demo run, so a demo manifest resumed with --demo must be accepted."""
    _mock_env(monkeypatch)
    out = tmp_path / "out"
    _write_manifest(out, {**_MATCHING_MANIFEST, "demo": True})
    assert _run(_args(tmp_path, resume=True, output=str(out), demo=True), _samplesheet(tmp_path))["ok"] is True


def test_resume_legacy_manifest_without_demo_treats_as_false(tmp_path, monkeypatch):
    """Manifests written before the demo key existed must still resume a real run."""
    _mock_env(monkeypatch)
    out = tmp_path / "out"
    # _MATCHING_MANIFEST intentionally omits the 'demo' key — same as old manifests.
    _write_manifest(out, _MATCHING_MANIFEST)
    assert _run(_args(tmp_path, resume=True, output=str(out), demo=False), _samplesheet(tmp_path))["ok"] is True


def test_resume_params_checksum_rejects_manifest_checksum_mismatch(tmp_path):
    payload = {"outdir": "/tmp/out/upstream/results", "aligner": "star_salmon"}
    repro = tmp_path / "reproducibility"
    repro.mkdir()
    (repro / "manifest.json").write_text(json.dumps({"params_checksum": "sha256:not-the-current-payload"}), encoding="utf-8")

    with pytest.raises(SkillError) as exc:
        preflight.check_resume_params_checksum(payload, tmp_path)

    assert exc.value.error_code == ErrorCode.INVALID_RESUME_STATE
    assert exc.value.details["field"] == "params_checksum"


def test_demo_skips_input_and_reference_checks_and_warns_on_overrides(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    result = _run(
        _args(tmp_path, demo=True, aligner="hisat2", fasta=str(tmp_path / "missing.fa"), gtf=None),
        {"sample_count": 0, "sample_names": [], "fastq_paths": [], "bam_paths": [], "strandedness_counts": {}, "unknown_columns": []},
    )
    assert result["references"] == {}
    assert result["aligner_effective"] == "star_salmon"
    assert any("demo" in warning.lower() for warning in result["warnings"])


def test_macos_docker_tmp_warning(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    monkeypatch.setattr(preflight.sys, "platform", "darwin")
    # Bypass check_output_dir_available so /tmp path validation doesn't fail on
    # macOS (sandbox/permissions prevent creating /tmp subdirs in CI/test env).
    monkeypatch.setattr(preflight, "check_output_dir_available", lambda *a, **kw: None)
    result = _run(_args(tmp_path, output="/tmp/clawbio-rnaseq-test"), _samplesheet(tmp_path))
    tmp_warnings = [w for w in result["warnings"] if "/tmp" in w]
    assert tmp_warnings, result["warnings"]
    msg = " ".join(tmp_warnings)
    # Must describe the real failure mode (a hard '.command.run: No such file or
    # directory' because Colima does not share /tmp into its VM) and give the
    # actionable fix — parity with nfcore-sarek / nfcore-scrnaseq, not the old
    # misleading "may be slow / VirtioFS" wording.
    assert "No such file or directory" in msg
    assert "home" in msg.lower()
    assert "slow" not in msg.lower()


def test_macos_docker_tmp_warning_resolves_private_tmp(tmp_path, monkeypatch):
    """Detection is resolve-based (is_under_tmp), so the macOS-canonical
    /private/tmp is caught just like /tmp."""
    _mock_env(monkeypatch)
    monkeypatch.setattr(preflight.sys, "platform", "darwin")
    monkeypatch.setattr(preflight, "check_output_dir_available", lambda *a, **kw: None)
    result = _run(_args(tmp_path, output="/private/tmp/clawbio-rnaseq-test"), _samplesheet(tmp_path))
    assert any("No such file or directory" in w for w in result["warnings"])


def test_macos_docker_non_tmp_output_has_no_tmp_warning(tmp_path, monkeypatch):
    """A normal (non-/tmp) output directory must not trigger the /tmp guard."""
    _mock_env(monkeypatch)
    monkeypatch.setattr(preflight.sys, "platform", "darwin")
    # pytest's tmp_path lives under /tmp on Linux CI but under /private/var/... on
    # macOS, so the real is_under_tmp would fire the guard on Linux and make this
    # test (which asserts the guard does NOT fire for a non-/tmp output) host-
    # dependent. Force the "not under /tmp" precondition deterministically; the
    # under-/tmp branch is covered by its own test.
    monkeypatch.setattr(preflight, "is_under_tmp", lambda _p: False)
    out = tmp_path / "out"
    result = _run(_args(tmp_path, output=str(out)), _samplesheet(tmp_path))
    assert not any("No such file or directory" in w for w in result["warnings"])


def test_demo_under_nxf_offline_is_rejected(tmp_path, monkeypatch):
    """--demo runs nf-core's remote `-profile test`; under NXF_OFFLINE it cannot
    fetch its test data, so preflight must fail fast with a clear, actionable
    DEMO_REQUIRES_NETWORK error instead of a cryptic Nextflow 'does not exist'."""
    _mock_env(monkeypatch)
    monkeypatch.setenv("NXF_OFFLINE", "true")
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, demo=True), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.DEMO_REQUIRES_NETWORK


def test_demo_without_nxf_offline_not_blocked(tmp_path, monkeypatch):
    """A demo with network available (no NXF_OFFLINE) must not trip the guard."""
    _mock_env(monkeypatch)
    monkeypatch.delenv("NXF_OFFLINE", raising=False)
    _run(_args(tmp_path, demo=True), _samplesheet(tmp_path))


def test_nxf_offline_real_run_not_blocked(tmp_path, monkeypatch):
    """A real (non-demo) run is fully local, so NXF_OFFLINE must not block it —
    the guard is specific to --demo's remote test profile."""
    _mock_env(monkeypatch)
    monkeypatch.setenv("NXF_OFFLINE", "true")
    _run(_args(tmp_path, demo=False), _samplesheet(tmp_path))


def test_email_on_fail_invalid_rejected(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    with pytest.raises(SkillError):
        _run(_args(tmp_path, email_on_fail="bad address@example.org"), _samplesheet(tmp_path))


def test_min_trimmed_reads_negative_rejected(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, min_trimmed_reads=-1), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION


def test_kallisto_quant_fraglen_below_minimum_rejected(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, kallisto_quant_fraglen=0), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION


def test_min_mapped_reads_above_maximum_rejected(tmp_path, monkeypatch):
    """min_mapped_reads > 100 is not a valid percentage and must be rejected."""
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, min_mapped_reads=150), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION


# ── Audit F9: pseudo_aligner_kmer_size must be an odd integer in [1, 31] ───────
# Both Salmon and Kallisto represent the index k-mer in a 64-bit machine word, so
# the k-mer must be odd and ≤ 31 (kallisto/salmon hard cap). An out-of-range or
# even value is a guaranteed indexing crash; the wrapper must catch it early like
# every other numeric tuning field rather than let the pipeline abort late.


def test_pseudo_aligner_kmer_size_below_minimum_rejected(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, pseudo_aligner_kmer_size=0), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION
    assert exc.value.details.get("field") == "pseudo_aligner_kmer_size"


def test_pseudo_aligner_kmer_size_above_maximum_rejected(tmp_path, monkeypatch):
    """k > 31 cannot be represented in Salmon/Kallisto's 64-bit k-mer word."""
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, pseudo_aligner_kmer_size=32), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION
    assert exc.value.details.get("maximum") == 31


def test_pseudo_aligner_kmer_size_even_rejected(tmp_path, monkeypatch):
    """Both pseudo-aligners require an odd k-mer; an even value always fails indexing."""
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, pseudo_aligner_kmer_size=30), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION


def test_pseudo_aligner_kmer_size_valid_odd_accepted(tmp_path, monkeypatch):
    """A valid odd k-mer in [1, 31] must pass preflight unchanged."""
    _mock_env(monkeypatch)
    result = _run(_args(tmp_path, pseudo_aligner_kmer_size=25), _samplesheet(tmp_path))
    assert result["ok"] is True


def test_composite_profile_passes_preflight(tmp_path, monkeypatch):
    """Composite profiles (comma-separated) must not raise INVALID_PROFILE."""
    _mock_env(monkeypatch)
    result = _run(_args(tmp_path, profile="docker,prokaryotic"), _samplesheet(tmp_path))
    assert result["ok"] is True


def test_remote_fasta_url_passes_preflight(tmp_path, monkeypatch):
    """--fasta <URI> must not raise REFERENCE_PATH_NOT_FOUND."""
    _mock_env(monkeypatch)
    monkeypatch.setattr(preflight, "check_output_dir_available", lambda *a, **kw: None)
    result = _run(
        _args(tmp_path, fasta="s3://example.org/genome.fa", gtf="s3://example.org/genes.gtf", allow_remote_inputs=True),
        _samplesheet(tmp_path),
    )
    assert result["ok"] is True, f"Remote s3:// refs should pass preflight with --allow-remote-inputs; got: {result}"


def test_local_ref_path_not_found_still_raises(tmp_path, monkeypatch):
    """A local path that does not exist must still raise REFERENCE_PATH_NOT_FOUND."""
    _mock_env(monkeypatch)
    monkeypatch.setattr(preflight, "check_output_dir_available", lambda *a, **kw: None)
    with pytest.raises(preflight.SkillError) as exc_info:
        _run(
            _args(tmp_path, fasta="/nonexistent/genome.fa", gtf="/nonexistent/genes.gtf"),
            _samplesheet(tmp_path),
        )
    assert exc_info.value.error_code == preflight.ErrorCode.REFERENCE_PATH_NOT_FOUND


_BACKEND_MISSING_ERROR_CODES = {
    "conda": {ErrorCode.MISSING_CONDA},
    "mamba": {ErrorCode.MISSING_CONDA},
    "singularity": {ErrorCode.MISSING_SINGULARITY},
    "apptainer": {ErrorCode.MISSING_SINGULARITY},
    "podman": {ErrorCode.MISSING_PODMAN},
    "shifter": {ErrorCode.MISSING_HPC_RUNTIME},
    "charliecloud": {ErrorCode.MISSING_HPC_RUNTIME},
}


def test_profile_backend_missing_binary_is_rejected(monkeypatch):
    monkeypatch.setattr(preflight.shutil, "which", lambda name: None)
    monkeypatch.setattr(preflight.subprocess, "run", lambda *a, **k: Namespace(returncode=0, stdout="", stderr=""))
    with pytest.raises(SkillError) as exc:
        preflight._check_profile("singularity")
    assert exc.value.error_code in _BACKEND_MISSING_ERROR_CODES["singularity"]



# ── _collect_reference_values ────────────────────────────────────────────────


def test_collect_reference_values_omits_empty_fields():
    # When only --fasta and --gtf are provided, the returned dict must contain
    # ONLY those two fields (plus "genome" if non-empty). Empty-string fields
    # must NOT be included — they cause _reference_checksums to hash the CWD.
    args = Namespace(genome=None, fasta="/ref/genome.fa", gtf="/ref/genes.gtf",
                     gff=None, transcript_fasta=None, additional_fasta=None,
                     gene_bed=None, splicesites=None, star_index=None,
                     rsem_index=None, hisat2_index=None, bowtie2_index=None,
                     salmon_index=None, kallisto_index=None)
    result = preflight._collect_reference_values(args)
    assert "fasta" in result and result["fasta"] == "/ref/genome.fa"
    assert "gtf" in result and result["gtf"] == "/ref/genes.gtf"
    for key in ("genome", "gff", "transcript_fasta", "additional_fasta", "gene_bed",
                "splicesites", "star_index", "rsem_index", "hisat2_index",
                "bowtie2_index", "salmon_index", "kallisto_index"):
        assert key not in result, (
            f"Empty-string field '{key}' must not appear in returned dict"
        )


# ── Schema-compliance validations (nf-core/rnaseq 3.26.0) ─────────────────────
# Red tests: document expected behaviour; turn green once implementation lands.


def test_email_long_tld_rejected(tmp_path, monkeypatch):
    """user@example.technology (11-char TLD) must be rejected per nf-core email schema."""
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, email="user@example.technology"), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION
    assert "email" in exc.value.details.get("field", "")


def test_genome_identifier_with_spaces_rejected(tmp_path, monkeypatch):
    """--genome 'bad genome' must be rejected — schema pattern ^[a-zA-Z0-9_\\-\\.]+$ forbids spaces."""
    _mock_env(monkeypatch)
    args = _args(tmp_path, genome="bad genome", fasta=None, gtf=None)
    with pytest.raises(SkillError) as exc:
        _run(args, _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION
    assert "genome" in str(exc.value.details).lower()


def test_fasta_wrong_extension_rejected(tmp_path, monkeypatch):
    """genome.txt must be rejected — schema forbids non-FASTA extensions."""
    _mock_env(monkeypatch)
    txt = tmp_path / "refs" / "genome.txt"
    txt.parent.mkdir(parents=True, exist_ok=True)
    txt.write_text(">chr1\nACGT\n", encoding="utf-8")
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, fasta=str(txt)), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION
    assert "fasta" in exc.value.details.get("field", "")


def test_hisat2_build_memory_bad_format_rejected(tmp_path, monkeypatch):
    """'nonsense' must be rejected per nf-core hisat2_build_memory schema pattern."""
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, hisat2_build_memory="nonsense"), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION
    assert "hisat2" in str(exc.value.details).lower()


def test_umitools_umi_separator_whitespace_rejected(tmp_path, monkeypatch):
    """A whitespace separator ' ' must be rejected — schema pattern ^\\S+$ forbids whitespace."""
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, umitools_umi_separator=" "), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION
    assert "umitools" in str(exc.value.details).lower()


def test_umitools_umi_separator_multichar_rejected(tmp_path, monkeypatch):
    """'::' (2 chars) must be rejected — schema has maxLength: 1, single character only."""
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, umitools_umi_separator="::"), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION
    assert "umitools" in str(exc.value.details).lower()


# ── GTF/GFF/BED extension validation ─────────────────────────────────────────

def test_gtf_wrong_extension_rejected(tmp_path, monkeypatch):
    """genes.txt for --gtf must be rejected — schema requires .gtf or .gtf.gz extension."""
    _mock_env(monkeypatch)
    f = _touch(tmp_path / "refs" / "genes.txt")
    fasta = _touch(tmp_path / "refs" / "genome.fa")
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, fasta=str(fasta), gtf=str(f)), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION
    assert "gtf" in exc.value.details.get("field", "")


def test_gff_wrong_extension_rejected(tmp_path, monkeypatch):
    """genes.txt for --gff must be rejected — schema: ^\\S+\\.gff3?(\\.gz)?$"""
    _mock_env(monkeypatch)
    f = _touch(tmp_path / "refs" / "genes.txt")
    fasta = _touch(tmp_path / "refs" / "genome.fa")
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, fasta=str(fasta), gtf=None, gff=str(f)), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION
    assert "gff" in exc.value.details.get("field", "")


def test_gene_bed_wrong_extension_rejected(tmp_path, monkeypatch):
    """genes.txt for --gene-bed must be rejected — schema requires .bed or .bed.gz extension."""
    _mock_env(monkeypatch)
    f = _touch(tmp_path / "refs" / "genes.txt")
    fasta = _touch(tmp_path / "refs" / "genome.fa")
    gtf = _touch(tmp_path / "refs" / "genes.gtf")
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, fasta=str(fasta), gtf=str(gtf), gene_bed=str(f)), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION
    assert "gene_bed" in exc.value.details.get("field", "")


# ---------------------------------------------------------------------------
# Enum validation — umitools_extract_method, umitools_grouping_method,
# salmon_quant_libtype, publish_dir_mode
# (nf-core/rnaseq 3.26.0 nextflow_schema.json enum constraints)
# ---------------------------------------------------------------------------

def test_umitools_extract_method_invalid_rejected(tmp_path, monkeypatch):
    """umitools_extract_method must be 'string' or 'regex' (schema enum)."""
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, umitools_extract_method="pattern"), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION


def test_umitools_grouping_method_invalid_rejected(tmp_path, monkeypatch):
    """umitools_grouping_method must be one of the schema enum values."""
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, umitools_grouping_method="random"), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION


def test_salmon_quant_libtype_invalid_rejected(tmp_path, monkeypatch):
    """salmon_quant_libtype must be one of the 16 schema enum values."""
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, salmon_quant_libtype="INVALID"), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION


def test_publish_dir_mode_invalid_rejected(tmp_path, monkeypatch):
    """publish_dir_mode must be one of the schema enum values."""
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, publish_dir_mode="invalid_mode"), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION


def test_sortmerna_index_nonexistent_path_rejected(tmp_path, monkeypatch):
    """sortmerna_index must be existence-checked in preflight like all other index paths.

    nf-core/rnaseq 3.26.0 schema: sortmerna_index is a string path parameter.
    A non-existent path must raise REFERENCE_PATH_NOT_FOUND so the user gets
    a clear error before Nextflow launches.
    """
    _mock_env(monkeypatch)
    bad_path = str(tmp_path / "nonexistent_sortmerna_idx")
    args = _args(tmp_path, remove_ribo_rna=True, ribo_removal_tool="sortmerna",
                 sortmerna_index=bad_path)
    with pytest.raises(SkillError) as exc:
        _run(args, _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.REFERENCE_PATH_NOT_FOUND
    assert exc.value.details.get("field") == "sortmerna_index"


def test_sortmerna_index_compatible_with_genome(tmp_path, monkeypatch):
    """--sortmerna-index must NOT conflict with --genome.

    sortmerna_index is an rRNA filter database, not a genome reference.
    nf-core/rnaseq allows --genome GRCh38 together with --sortmerna-index /db.
    If the wrapper raises CONFLICTING_REFERENCES here, it is broken.
    """
    _mock_env(monkeypatch)
    db_dir = tmp_path / "sortmerna_db"
    db_dir.mkdir()
    # fasta/gtf=None: using --genome, so explicit refs must be absent.
    args = _args(tmp_path, genome="GRCh38", fasta=None, gtf=None,
                 remove_ribo_rna=True, ribo_removal_tool="sortmerna",
                 sortmerna_index=str(db_dir))
    # Must not raise — genome + sortmerna_index is a valid combination.
    _run(args, _samplesheet(tmp_path))


def test_extra_star_align_args_rejected_for_star_rsem(tmp_path, monkeypatch):
    """--extra-star-align-args is only valid for star_salmon per the nf-core schema.

    The fix message must reference ext.args (not ext.args2, which does not exist
    in the RSEM_CALCULATEEXPRESSION module) and must mention --star-options.
    """
    _mock_env(monkeypatch)
    args = _args(tmp_path, aligner="star_rsem", extra_star_align_args="--outSAMtype BAM SortedByCoordinate")
    with pytest.raises(SkillError) as exc_info:
        _run(args, _samplesheet(tmp_path))
    err = exc_info.value
    assert err.error_code == ErrorCode.INCOMPATIBLE_ALIGNER_ARGS
    fix = err.fix or ""
    assert "ext.args2" not in fix, "fix must not reference ext.args2 (that key does not exist in RSEM_CALCULATEEXPRESSION)"
    assert "star-options" in fix or "star_options" in fix, "fix should guide users to --star-options"


def test_extra_star_align_args_rejected_for_hisat2(tmp_path, monkeypatch):
    """--extra-star-align-args must also be rejected for --aligner hisat2."""
    _mock_env(monkeypatch)
    args = _args(tmp_path, aligner="hisat2", extra_star_align_args="--some-star-flag")
    with pytest.raises(SkillError) as exc_info:
        _run(args, _samplesheet(tmp_path))
    err = exc_info.value
    assert err.error_code == ErrorCode.INCOMPATIBLE_ALIGNER_ARGS


# ── Audit follow-up: F-02 timeout bounds ──────────────────────────────────────


def test_timeout_hours_zero_accepted_disables_cap(tmp_path, monkeypatch):
    """--timeout-hours 0 is accepted: it disables the wall-clock cap (parity with
    nfcore-scrnaseq/sarek). Only negative values are rejected at preflight."""
    _mock_env(monkeypatch)
    result = _run(_args(tmp_path, timeout_hours=0), _samplesheet(tmp_path))
    assert result["ok"] is True


def test_timeout_hours_negative_rejected(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, timeout_hours=-3), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION


def test_timeout_hours_positive_accepted(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    result = _run(_args(tmp_path, timeout_hours=48), _samplesheet(tmp_path))
    assert result["ok"] is True


# ── Audit follow-up: F-03 transcriptome-only pseudo route ─────────────────────


def test_transcriptome_only_pseudo_route_is_accepted(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    tx = _touch(tmp_path / "refs" / "transcriptome.fa", ">t1\nACGT\n")
    args = _args(
        tmp_path,
        fasta=None,
        transcript_fasta=str(tx),
        pseudo_aligner="salmon",
        skip_alignment=True,
    )
    result = _run(args, _samplesheet(tmp_path))
    assert result["ok"] is True


def test_transcriptome_only_without_annotation_still_rejected(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    tx = _touch(tmp_path / "refs" / "transcriptome.fa", ">t1\nACGT\n")
    args = _args(
        tmp_path,
        fasta=None,
        gtf=None,
        transcript_fasta=str(tx),
        pseudo_aligner="salmon",
        skip_alignment=True,
    )
    with pytest.raises(SkillError) as exc:
        _run(args, _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.MISSING_REFERENCE


def test_transcript_source_without_skip_alignment_is_rejected(tmp_path, monkeypatch):
    # A genome aligner (default star_salmon) still runs, so a genome FASTA is required.
    _mock_env(monkeypatch)
    tx = _touch(tmp_path / "refs" / "transcriptome.fa", ">t1\nACGT\n")
    args = _args(tmp_path, fasta=None, transcript_fasta=str(tx), pseudo_aligner="salmon")
    with pytest.raises(SkillError) as exc:
        _run(args, _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.MISSING_REFERENCE


# ── Audit follow-up: F-04 rseqc-modules validation ────────────────────────────


def test_invalid_rseqc_module_rejected(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    args = _args(tmp_path, rseqc_modules="bam_stat,infer_experimnt")
    with pytest.raises(SkillError) as exc:
        _run(args, _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION
    assert "infer_experimnt" in exc.value.details["unknown"]


def test_valid_rseqc_modules_including_tin_accepted(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    args = _args(tmp_path, rseqc_modules="bam_stat,tin")
    result = _run(args, _samplesheet(tmp_path))
    assert result["ok"] is True


# ── Audit follow-up: F-10 contaminant-screening tool↔database cross-check ──────


def test_kraken2_screening_without_db_rejected(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    args = _args(tmp_path, contaminant_screening="kraken2", kraken_db=None)
    with pytest.raises(SkillError) as exc:
        _run(args, _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION
    assert exc.value.details["field"] == "kraken_db"


def test_kraken2_bracken_screening_without_db_rejected(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    args = _args(tmp_path, contaminant_screening="kraken2_bracken", kraken_db=None)
    with pytest.raises(SkillError) as exc:
        _run(args, _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION


def test_sylph_screening_without_db_rejected(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    args = _args(tmp_path, contaminant_screening="sylph", sylph_db=None)
    with pytest.raises(SkillError) as exc:
        _run(args, _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION
    assert exc.value.details["field"] == "sylph_db"


def test_sylph_screening_remote_db_uri_accepted(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    args = _args(tmp_path, contaminant_screening="sylph", sylph_db="s3://bucket.example.org/sylph/db.syldb")
    result = _run(args, _samplesheet(tmp_path))
    assert result["ok"] is True


def test_kraken2_screening_with_db_accepted(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    db = tmp_path / "kraken_db"
    db.mkdir()
    args = _args(tmp_path, contaminant_screening="kraken2", kraken_db=str(db))
    result = _run(args, _samplesheet(tmp_path))
    assert result["ok"] is True


def test_bracken_precision_without_bracken_tool_warns(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    db = tmp_path / "kraken_db"
    db.mkdir()
    # bracken_precision only applies to kraken2_bracken; with plain kraken2 it is inert.
    args = _args(tmp_path, contaminant_screening="kraken2", kraken_db=str(db), bracken_precision="G")
    result = _run(args, _samplesheet(tmp_path))
    assert any("bracken" in str(w).lower() for w in result["warnings"])


def test_nextflow_config_with_params_override_rejected(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    cfg = tmp_path / "override.config"
    cfg.write_text("params.aligner = 'hisat2'\n", encoding="utf-8")
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, nextflow_config=[str(cfg)]), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION
    assert exc.value.details["field"] == "nextflow_config"


def test_nextflow_config_with_params_block_override_rejected(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    cfg = tmp_path / "override-block.config"
    cfg.write_text("params {\n  aligner = 'hisat2'\n}\n", encoding="utf-8")
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, nextflow_config=[str(cfg)]), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION
    assert exc.value.details["field"] == "nextflow_config"


def test_nextflow_config_custom_genome_catalogue_accepted(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    cfg = tmp_path / "my_genomes.config"
    cfg.write_text(
        "params.genomes {\n"
        "  'MY_GENOME' {\n"
        "    fasta = '/refs/my.fa'\n"
        "    gtf = '/refs/my.gtf'\n"
        "  }\n"
        "}\n",
        encoding="utf-8",
    )
    result = _run(
        _args(tmp_path, genome="MY_GENOME", fasta=None, gtf=None, nextflow_config=[str(cfg)]),
        _samplesheet(tmp_path),
        pipeline_source=_DEFAULT_REMOTE_PIPELINE_SOURCE,
    )
    assert result["ok"] is True


def test_nextflow_config_without_params_override_accepted(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    cfg = tmp_path / "hpc.config"
    cfg.write_text("process.executor = 'slurm'\n", encoding="utf-8")
    result = _run(_args(tmp_path, nextflow_config=[str(cfg)]), _samplesheet(tmp_path))
    assert result["ok"] is True


def test_parabricks_star_requires_star_salmon_aligner(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, aligner="hisat2", use_parabricks_star=True), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION
    assert exc.value.details["field"] == "use_parabricks_star"


def test_sentieon_star_requires_star_aligner(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, aligner="bowtie2_salmon", use_sentieon_star=True), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION
    assert exc.value.details["field"] == "use_sentieon_star"


def test_gpu_ribodetector_requires_matching_ribo_tool(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, use_gpu_ribodetector=True, ribo_removal_tool="sortmerna"), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION
    assert exc.value.details["field"] == "use_gpu_ribodetector"


# ── Audit follow-up: F-11 UMI options set without --with-umi warn (silent no-op)


def test_umi_options_without_with_umi_warn(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    args = _args(tmp_path, with_umi=False, umitools_bc_pattern="NNNNNN")
    result = _run(args, _samplesheet(tmp_path))
    assert result["ok"] is True
    assert any("with-umi" in str(w).lower() for w in result["warnings"])


def test_umi_options_with_with_umi_do_not_warn(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    args = _args(tmp_path, with_umi=True, umitools_bc_pattern="NNNNNN")
    result = _run(args, _samplesheet(tmp_path))
    assert not any("with-umi" in str(w).lower() for w in result["warnings"])


def test_no_umi_options_no_umi_warning(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    result = _run(_args(tmp_path), _samplesheet(tmp_path))
    assert not any("with-umi" in str(w).lower() for w in result["warnings"])


# ── Audit follow-up: F-12 prebuilt-index reference completeness (no --fasta) ───


def _index_dir(tmp_path: Path, name: str) -> str:
    d = tmp_path / "refs" / name
    d.mkdir(parents=True, exist_ok=True)
    return str(d)


def test_star_salmon_prebuilt_index_plus_transcript_fasta_accepted(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    tx = _touch(tmp_path / "refs" / "transcriptome.fa", ">t1\nACGT\n")
    args = _args(
        tmp_path,
        fasta=None,
        star_index=_index_dir(tmp_path, "star"),
        transcript_fasta=str(tx),
    )
    result = _run(args, _samplesheet(tmp_path))
    assert result["ok"] is True


def test_star_salmon_prebuilt_index_plus_salmon_index_accepted(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    args = _args(
        tmp_path,
        fasta=None,
        star_index=_index_dir(tmp_path, "star"),
        salmon_index=_index_dir(tmp_path, "salmon"),
    )
    result = _run(args, _samplesheet(tmp_path))
    assert result["ok"] is True


def test_star_salmon_index_without_transcript_source_rejected(tmp_path, monkeypatch):
    # star_salmon still needs a transcript source (fasta/transcript_fasta/salmon_index)
    # for Salmon quantification; a STAR index + GTF alone is incomplete.
    _mock_env(monkeypatch)
    args = _args(tmp_path, fasta=None, star_index=_index_dir(tmp_path, "star"))
    with pytest.raises(SkillError) as exc:
        _run(args, _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.MISSING_REFERENCE


def test_hisat2_prebuilt_index_plus_gtf_accepted(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    args = _args(tmp_path, aligner="hisat2", fasta=None, hisat2_index=_index_dir(tmp_path, "hisat2"))
    result = _run(args, _samplesheet(tmp_path))
    assert result["ok"] is True


def test_bowtie2_salmon_prebuilt_index_plus_salmon_index_accepted(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    args = _args(
        tmp_path,
        aligner="bowtie2_salmon",
        fasta=None,
        bowtie2_index=_index_dir(tmp_path, "bt2"),
        salmon_index=_index_dir(tmp_path, "salmon"),
    )
    result = _run(args, _samplesheet(tmp_path))
    assert result["ok"] is True


def test_star_rsem_rsem_index_accepted(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    args = _args(tmp_path, aligner="star_rsem", fasta=None, rsem_index=_index_dir(tmp_path, "rsem"))
    result = _run(args, _samplesheet(tmp_path))
    assert result["ok"] is True


def test_star_rsem_star_index_without_rsem_source_rejected(tmp_path, monkeypatch):
    # RSEM quantification needs --rsem-index or --fasta; a bare STAR index is not enough.
    _mock_env(monkeypatch)
    args = _args(tmp_path, aligner="star_rsem", fasta=None, star_index=_index_dir(tmp_path, "star"))
    with pytest.raises(SkillError) as exc:
        _run(args, _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.MISSING_REFERENCE


def test_prebuilt_index_without_annotation_rejected(tmp_path, monkeypatch):
    # Annotation (GTF/GFF) is always required for quantification.
    _mock_env(monkeypatch)
    tx = _touch(tmp_path / "refs" / "transcriptome.fa", ">t1\nACGT\n")
    args = _args(
        tmp_path,
        fasta=None,
        gtf=None,
        star_index=_index_dir(tmp_path, "star"),
        transcript_fasta=str(tx),
    )
    with pytest.raises(SkillError) as exc:
        _run(args, _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.MISSING_REFERENCE


# ── Audit follow-up F-06: FASTA extension regex matches real extensions only ──


def test_fasta_regex_accepts_canonical_extensions():
    for name in ("genome.fa", "genome.fasta", "genome.fna", "g.fa.gz", "g.fasta.gz", "g.fna.gz"):
        assert preflight._FASTA_EXT_RE.match(name), name


def test_fasta_regex_rejects_bogus_extensions():
    for name in ("genome.fnasta", "genome.fnasta.gz", "genome.txt", "genome.faa"):
        assert not preflight._FASTA_EXT_RE.match(name), name


# ── Audit follow-up F-1: reference paths that resolve into a directory with ──
# whitespace fail the nf-core ^\S+ schema pattern at runtime. Catch them at
# preflight with a precise, dedicated error (mirrors the samplesheet input guard)
# instead of letting Nextflow abort late with a misleading message.


def test_reference_path_with_whitespace_in_directory_rejected(tmp_path):
    spaced_dir = tmp_path / "my refs"
    fasta = _touch(spaced_dir / "genome.fa", ">chr1\nACGT\n")
    gtf = _touch(spaced_dir / "genes.gtf", "chr1\tsrc\tgene\t1\t4\t.\t+\t.\tg\n")
    with pytest.raises(SkillError) as exc:
        preflight._check_reference_paths_exist({"fasta": str(fasta), "gtf": str(gtf)})
    assert exc.value.error_code == ErrorCode.REFERENCE_PATH_HAS_WHITESPACE
    assert exc.value.details["field"] == "fasta"


def test_reference_path_whitespace_check_precedes_extension_check(tmp_path):
    # An absolute path containing whitespace must report the whitespace problem,
    # not a misleading "FASTA extension" error (the ^\S+ extension regex also fails
    # on whitespace, so ordering matters).
    spaced = _touch(tmp_path / "ref dir" / "genome.fa", ">c\nA\n")
    with pytest.raises(SkillError) as exc:
        preflight._check_reference_paths_exist({"fasta": str(spaced)})
    assert exc.value.error_code == ErrorCode.REFERENCE_PATH_HAS_WHITESPACE


def test_whitespace_free_reference_path_still_accepted(tmp_path):
    fasta = _touch(tmp_path / "refs" / "genome.fa", ">c\nA\n")
    gtf = _touch(tmp_path / "refs" / "genes.gtf", "chr1\tsrc\tgene\t1\t4\t.\t+\t.\tg\n")
    # Must not raise.
    preflight._check_reference_paths_exist({"fasta": str(fasta), "gtf": str(gtf)})


# ── Audit follow-up F-5: --check mode intentionally skips the Nextflow version ──
# gate (presence only). That leniency must be surfaced as a warning so the run is
# not silently assumed to have a compatible Nextflow.


def test_check_mode_warns_that_nextflow_version_unverified(tmp_path, monkeypatch):
    monkeypatch.setattr(preflight, "_check_java", lambda: {"path": "/usr/bin/java", "version": "17"})
    monkeypatch.setattr(preflight, "_check_profile", lambda profile: {"profile": profile, "backend_path": "", "backend_ready": True})
    monkeypatch.setattr(preflight.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(preflight, "_command_output", lambda cmd: "")
    result = preflight.run_preflight(
        _args(tmp_path, check=True, demo=True, profile="test"),
        pipeline_source={"source_kind": "remote", "source_ref": "3.26.0"},
        samplesheet_summary={"sample_count": 0, "fastq_paths": []},
    )
    assert result["nextflow"]["version_checked"] is False
    assert any("version" in w.lower() and "check" in w.lower() for w in result["warnings"])


def test_real_run_does_not_emit_version_unverified_warning(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    result = _run(_args(tmp_path), _samplesheet(tmp_path))
    assert not any("not verified" in w.lower() for w in result["warnings"])


# ── Audit M1: --genome allows additive annotation/transcriptome overrides ─────
# nf-core/rnaseq supports iGenomes + an annotation/transcriptome override
# (e.g. --genome GRCh38 --gtf custom.gtf, or --additional-fasta spike_ins.fa).
# Only a second *genome sequence* source (--fasta) or a genome *index*
# (--star-index/--rsem-index/--hisat2-index/--bowtie2-index) genuinely conflicts.


def test_genome_with_additional_fasta_accepted(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    af = _touch(tmp_path / "refs" / "spikein.fa", ">ercc\nACGT\n")
    result = _run(
        _args(tmp_path, genome="GRCh38", fasta=None, gtf=None, additional_fasta=str(af)),
        _samplesheet(tmp_path),
    )
    assert result["references"]["genome"] == "GRCh38"
    assert result["references"]["additional_fasta"].endswith("spikein.fa")


def test_genome_with_gtf_override_accepted(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    gtf = _touch(tmp_path / "refs" / "custom.gtf", "chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"g1\";\n")
    result = _run(
        _args(tmp_path, genome="GRCh38", fasta=None, gtf=str(gtf)),
        _samplesheet(tmp_path),
    )
    assert result["references"]["genome"] == "GRCh38"
    assert result["references"]["gtf"].endswith("custom.gtf")


def test_genome_with_transcript_fasta_accepted(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    tf = _touch(tmp_path / "refs" / "tx.fa", ">tx1\nACGT\n")
    result = _run(
        _args(tmp_path, genome="GRCh38", fasta=None, gtf=None, transcript_fasta=str(tf)),
        _samplesheet(tmp_path),
    )
    assert result["references"]["genome"] == "GRCh38"


def test_genome_with_explicit_fasta_still_rejected(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, genome="GRCh38", gtf=None), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.CONFLICTING_REFERENCES
    assert "fasta" in exc.value.details.get("explicit_set", [])


def test_genome_with_star_index_still_rejected(tmp_path, monkeypatch):
    field = "star_index"
    index_path = _touch(tmp_path / "refs" / field / "SAindex")
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        _run(
            _args(tmp_path, genome="GRCh38", fasta=None, gtf=None, **{field: str(index_path.parent)}),
            _samplesheet(tmp_path),
        )
    assert exc.value.error_code == ErrorCode.CONFLICTING_REFERENCES
    assert field in exc.value.details.get("explicit_set", [])


# ── Audit M2: --nextflow-config params-override filter is no longer evadable ───


def test_nextflow_config_params_subscript_override_rejected(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    cfg = tmp_path / "subscript.config"
    cfg.write_text("params['aligner'] = 'hisat2'\n", encoding="utf-8")
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, nextflow_config=[str(cfg)]), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION
    assert exc.value.details["field"] == "nextflow_config"


def test_nextflow_config_params_put_override_rejected(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    cfg = tmp_path / "put.config"
    cfg.write_text("params.put('aligner', 'hisat2')\n", encoding="utf-8")
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, nextflow_config=[str(cfg)]), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION


def test_nextflow_config_include_local_params_override_rejected(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    inner = tmp_path / "inner.config"
    inner.write_text("params.aligner = 'hisat2'\n", encoding="utf-8")
    outer = tmp_path / "outer.config"
    outer.write_text("includeConfig 'inner.config'\n", encoding="utf-8")
    with pytest.raises(SkillError) as exc:
        _run(_args(tmp_path, nextflow_config=[str(outer)]), _samplesheet(tmp_path))
    assert exc.value.error_code == ErrorCode.INVALID_PRESET_CONFIGURATION


def test_nextflow_config_include_clean_local_accepted(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    inner = tmp_path / "inner_clean.config"
    inner.write_text("process.cpus = 8\n", encoding="utf-8")
    outer = tmp_path / "outer_clean.config"
    outer.write_text("includeConfig 'inner_clean.config'\n", encoding="utf-8")
    result = _run(_args(tmp_path, nextflow_config=[str(outer)]), _samplesheet(tmp_path))
    assert result["ok"] is True


def test_nextflow_config_include_remote_warns_not_blocks(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    cfg = tmp_path / "remote_include.config"
    cfg.write_text("includeConfig 'https://example.org/foo.config'\n", encoding="utf-8")
    result = _run(_args(tmp_path, nextflow_config=[str(cfg)]), _samplesheet(tmp_path))
    assert result["ok"] is True
    assert any("includeConfig" in w and "audit" in w.lower() for w in result["warnings"])


def test_nextflow_config_include_cycle_is_safe(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    a = tmp_path / "a.config"
    b = tmp_path / "b.config"
    a.write_text("includeConfig 'b.config'\n", encoding="utf-8")
    b.write_text("includeConfig 'a.config'\n", encoding="utf-8")
    result = _run(_args(tmp_path, nextflow_config=[str(a)]), _samplesheet(tmp_path))
    assert result["ok"] is True


# ── Audit L1: preflight shares the single-source FASTQ regex (no local shadow) ─


def test_preflight_does_not_shadow_shared_fastq_regex():
    src = (_SKILL_DIR / "preflight.py").read_text(encoding="utf-8")
    assert "_FASTQ_BASENAME_RE = re.compile" not in src, (
        "preflight must use the shared schemas.FASTQ_BASENAME_RE, not a local redefinition"
    )
    assert "FASTQ_BASENAME_RE as _FASTQ_BASENAME_RE" in src


# ── Audit L2: macOS Docker + STAR index build warns about the memory ceiling ──


def test_macos_docker_star_build_without_index_warns(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    monkeypatch.setattr(preflight.sys, "platform", "darwin")
    result = _run(
        _args(tmp_path, profile="docker", aligner="star_salmon"),  # fasta set, no star_index
        _samplesheet(tmp_path),
    )
    assert any("memory" in w.lower() and "index" in w.lower() for w in result["warnings"])


def test_macos_docker_star_with_prebuilt_index_no_memory_warning(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    monkeypatch.setattr(preflight.sys, "platform", "darwin")
    star_dir = tmp_path / "refs" / "star"
    star_dir.mkdir(parents=True)
    result = _run(
        _args(tmp_path, profile="docker", aligner="star_salmon", fasta=None, gtf=None,
              genome="GRCh38", star_index=None),  # genome route, no local fasta build
        _samplesheet(tmp_path),
    )
    assert not any("memory" in w.lower() and "STAR" in w for w in result["warnings"])


def test_macos_docker_non_star_aligner_no_memory_warning(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    monkeypatch.setattr(preflight.sys, "platform", "darwin")
    result = _run(
        _args(tmp_path, profile="docker", aligner="hisat2"),
        _samplesheet(tmp_path),
    )
    assert not any("STAR index" in w for w in result["warnings"])


# ── Audit follow-up F2: --gtf + --gff together must match upstream behaviour ───
# nf-core/rnaseq docs: "If --gff is provided … the latter [GTF] will be used if
# both are provided." The wrapper therefore keeps --gtf, drops --gff, and warns —
# it must NOT reject a configuration the pipeline accepts.


def test_gtf_and_gff_together_prefers_gtf_with_warning(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    gff = _touch(tmp_path / "refs" / "genes.gff3", "##gff-version 3\n")
    args = _args(tmp_path, gff=str(gff))  # default _args sets fasta+gtf
    result = _run(args, _samplesheet(tmp_path))
    assert result["references"]["gtf"].endswith("genes.gtf")
    assert "gff" not in result["references"]
    assert args.gff is None, "--gff must be dropped so it is not written to params.yaml"
    assert any("gff" in w.lower() and "gtf" in w.lower() for w in result["warnings"])


def test_gtf_and_gff_together_with_genome_prefers_gtf_with_warning(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    gtf = _touch(tmp_path / "refs" / "a.gtf", "x\n")
    gff = _touch(tmp_path / "refs" / "a.gff3", "##gff-version 3\n")
    args = _args(tmp_path, genome="GRCh38", fasta=None, gtf=str(gtf), gff=str(gff))
    result = _run(args, _samplesheet(tmp_path))
    assert result["references"]["genome"] == "GRCh38"
    assert result["references"]["gtf"].endswith("a.gtf")
    assert "gff" not in result["references"]
    assert args.gff is None
    assert any("gff" in w.lower() for w in result["warnings"])


def test_genome_with_only_gff_accepted(tmp_path, monkeypatch):
    _mock_env(monkeypatch)
    gff = _touch(tmp_path / "refs" / "only.gff3", "##gff-version 3\n")
    result = _run(_args(tmp_path, genome="GRCh38", fasta=None, gtf=None, gff=str(gff)), _samplesheet(tmp_path))
    assert result["references"]["genome"] == "GRCh38"


# ── Audit F3: warn when an index flag the chosen aligner cannot consume is given ─


def test_warns_when_salmon_index_supplied_for_star_rsem(tmp_path, monkeypatch):
    """star_rsem quantifies with RSEM, not Salmon — a --salmon-index is silently
    unused, so the wrapper must warn rather than let the user believe it took effect."""
    _mock_env(monkeypatch)
    salmon_index = _touch(tmp_path / "refs" / "salmon_idx" / "info.json")
    args = _args(tmp_path, aligner="star_rsem", salmon_index=str(salmon_index))
    result = _run(args, _samplesheet(tmp_path))
    assert any("salmon-index" in w and "star_rsem" in w for w in result["warnings"])


def test_no_index_mismatch_warning_when_index_matches_aligner(tmp_path, monkeypatch):
    """star_salmon DOES consume a --salmon-index — no mismatch warning should fire."""
    _mock_env(monkeypatch)
    salmon_index = _touch(tmp_path / "refs" / "salmon_idx" / "info.json")
    args = _args(tmp_path, aligner="star_salmon", salmon_index=str(salmon_index))
    result = _run(args, _samplesheet(tmp_path))
    assert not any("does not use" in w.lower() for w in result["warnings"])


def test_kallisto_index_relevant_for_kallisto_pseudo_aligner(tmp_path, monkeypatch):
    """A --kallisto-index is consumed when --pseudo-aligner kallisto runs, so no warning."""
    _mock_env(monkeypatch)
    kallisto_index = _touch(tmp_path / "refs" / "kallisto_idx" / "index")
    args = _args(tmp_path, aligner="star_salmon", pseudo_aligner="kallisto", kallisto_index=str(kallisto_index))
    result = _run(args, _samplesheet(tmp_path))
    assert not any("does not use" in w.lower() for w in result["warnings"])


# ── F4: GENCODE autodetect handles an uppercase .GZ gzip suffix ───────────────


def test_gtf_has_gencode_markers_detects_uppercase_gz_suffix(tmp_path):
    """A gzipped GTF whose name ends in uppercase .GZ must still be opened as gzip
    (case-insensitive suffix), so GENCODE markers are detected, not silently missed."""
    import gzip

    gtf = tmp_path / "genes.gtf.GZ"
    with gzip.open(gtf, "wt", encoding="utf-8") as handle:
        handle.write('chr1\tHAVANA\tgene\t1\t10\t.\t+\t.\tgene_id "ENSG1"; gene_type "protein_coding";\n')
    assert preflight._gtf_has_gencode_markers(gtf) is True


# --- --allow-remote-inputs gate (local-first by default) ---------------------

def test_remote_fastq_rejected_by_default():
    from argparse import Namespace as _NS
    with pytest.raises(SkillError) as exc:
        preflight._check_remote_inputs(
            _NS(allow_remote_inputs=False), {"fastq_paths": ["s3://bucket/n_R1.fastq.gz"], "bam_paths": []}
        )
    assert exc.value.error_code == ErrorCode.REMOTE_INPUT_NOT_ALLOWED


def test_remote_reference_rejected_by_default():
    from argparse import Namespace as _NS
    with pytest.raises(SkillError) as exc:
        preflight._check_remote_inputs(_NS(allow_remote_inputs=False, fasta="gs://bucket/genome.fa"), {})
    assert exc.value.error_code == ErrorCode.REMOTE_INPUT_NOT_ALLOWED


def test_remote_inputs_allowed_with_flag_warns(capsys):
    from argparse import Namespace as _NS
    preflight._check_remote_inputs(_NS(allow_remote_inputs=True), {"fastq_paths": ["https://x/n_R1.fastq.gz"]})
    assert "--allow-remote-inputs" in capsys.readouterr().err


def test_local_inputs_pass_the_remote_gate():
    from argparse import Namespace as _NS
    preflight._check_remote_inputs(_NS(allow_remote_inputs=False), {"fastq_paths": ["/d/n_R1.fastq.gz"]})
