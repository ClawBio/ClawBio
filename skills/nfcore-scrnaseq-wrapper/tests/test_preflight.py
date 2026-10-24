from __future__ import annotations

from argparse import Namespace
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import preflight
from errors import SkillError
from nfcore_4_1_0_contract import POLICY_SOURCE_NFCORE_DOCS

_PIPELINE_SOURCE = {
    "source_kind": "local_checkout",
    "source_ref": "scrnaseq-checkout",
    "resolved_version": "abc123",
    "branch": "main",
    "dirty": False,
}


def _mock_env(monkeypatch):
    monkeypatch.setattr(
        preflight, "_check_java", lambda: {"path": "/usr/bin/java", "version": "17.0.8"}
    )
    monkeypatch.setattr(
        preflight,
        "_check_nextflow",
        lambda: {"path": "/usr/bin/nextflow", "version": "25.04.0"},
    )
    monkeypatch.setattr(
        preflight,
        "_check_profile",
        lambda profile: {
            "profile": profile,
            "backend_path": "/usr/bin/docker",
            "backend_ready": True,
        },
    )


def _args(tmp_path: Path) -> Namespace:
    return Namespace(
        output=str(tmp_path / "out"),
        resume=False,
        preset="star",
        protocol="10XV3",
        profile="docker",
        demo=False,
        email=None,
        email_on_fail=None,
        genome=None,
        fasta=str(tmp_path / "genome.fa"),
        gtf=str(tmp_path / "genes.gtf"),
        transcript_fasta=None,
        txp2gene=None,
        simpleaf_index=None,
        kallisto_index=None,
        star_index=None,
        cellranger_index=None,
        barcode_whitelist=None,
        cellrangerarc_reference=None,
        cellrangerarc_config=None,
        motifs=None,
        kb_t1c=None,
        kb_t2c=None,
        cellranger_vdj_index=None,
        skip_cellrangermulti_vdjref=False,
        gex_frna_probe_set=None,
        gex_target_panel=None,
        gex_cmo_set=None,
        fb_reference=None,
        vdj_inner_enrichment_primers=None,
        gex_barcode_sample_assignment=None,
        cellranger_multi_barcodes=None,
    )


def test_reference_paths_accept_remote_uris():
    """Remote reference URIs (s3://, gs://, https://) are passed through: Nextflow
    resolves them at runtime, so existence is not pre-checked (parity with
    nfcore-rnaseq/sarek and the nf-core schema)."""
    preflight._check_reference_paths_exist(
        {
            "fasta": "s3://bucket/genome.fa",
            "gtf": "https://example.org/genes.gtf",
            "star_index": "gs://bucket/star/",
        }
    )


# --- --allow-remote-inputs gate (local-first by default) ---------------------

def test_remote_fastq_rejected_by_default():
    """Remote samplesheet inputs are rejected unless --allow-remote-inputs is set."""
    with pytest.raises(SkillError) as exc:
        preflight._check_remote_inputs(
            Namespace(allow_remote_inputs=False),
            {"fastq_paths": ["s3://bucket/sampleA_R1.fastq.gz"]},
        )
    assert exc.value.error_code == "REMOTE_INPUT_NOT_ALLOWED"
    assert "s3://bucket/sampleA_R1.fastq.gz" in exc.value.details["remote_paths"]


def test_remote_reference_rejected_by_default():
    """Remote reference/index paths are rejected unless --allow-remote-inputs is set."""
    with pytest.raises(SkillError) as exc:
        preflight._check_remote_inputs(
            Namespace(allow_remote_inputs=False, fasta="https://example.org/genome.fa"),
            {"fastq_paths": []},
        )
    assert exc.value.error_code == "REMOTE_INPUT_NOT_ALLOWED"


def test_remote_inputs_allowed_with_flag_warns(capsys):
    """With --allow-remote-inputs, remote paths pass and a network-fetch WARNING is emitted."""
    preflight._check_remote_inputs(
        Namespace(allow_remote_inputs=True),
        {"fastq_paths": ["gs://bucket/sampleA_R1.fastq.gz"]},
    )
    err = capsys.readouterr().err
    assert "--allow-remote-inputs" in err and "gs://bucket/sampleA_R1.fastq.gz" in err


def test_local_inputs_pass_the_gate():
    """Local paths never trigger the remote-input gate."""
    preflight._check_remote_inputs(
        Namespace(allow_remote_inputs=False),
        {"fastq_paths": ["/data/sampleA_R1.fastq.gz"]},
    )


def test_preflight_happy_path(tmp_path, monkeypatch):
    args = _args(tmp_path)
    Path(args.fasta).write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(args.gtf).write_text(
        'chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id "g1";\n', encoding="utf-8"
    )
    _mock_env(monkeypatch)
    result = preflight.run_preflight(
        args,
        pipeline_source=_PIPELINE_SOURCE,
        samplesheet_summary={"sample_count": 1, "unknown_columns": []},
    )
    assert result["ok"] is True


def test_demo_under_nxf_offline_is_rejected(tmp_path, monkeypatch):
    """--demo runs nf-core's remote `-profile test`; under NXF_OFFLINE it cannot
    fetch its test data, so preflight must fail fast with a clear DEMO_REQUIRES_NETWORK
    error instead of a cryptic Nextflow 'does not exist'."""
    args = _args(tmp_path)
    args.demo = True
    _mock_env(monkeypatch)
    monkeypatch.setenv("NXF_OFFLINE", "true")
    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={"sample_count": 1, "unknown_columns": []},
        )
    assert exc.value.error_code == "DEMO_REQUIRES_NETWORK"


def test_demo_without_nxf_offline_not_blocked(tmp_path, monkeypatch):
    """A demo with network available (no NXF_OFFLINE) must not trip the guard."""
    args = _args(tmp_path)
    args.demo = True
    _mock_env(monkeypatch)
    monkeypatch.delenv("NXF_OFFLINE", raising=False)
    result = preflight.run_preflight(
        args,
        pipeline_source=_PIPELINE_SOURCE,
        samplesheet_summary={"sample_count": 1, "unknown_columns": []},
    )
    assert result["ok"] is True


def test_nxf_offline_real_run_not_blocked(tmp_path, monkeypatch):
    """A real (non-demo) run is fully local, so NXF_OFFLINE must not block it."""
    args = _args(tmp_path)
    Path(args.fasta).write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(args.gtf).write_text(
        'chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id "g1";\n', encoding="utf-8"
    )
    _mock_env(monkeypatch)
    monkeypatch.setenv("NXF_OFFLINE", "true")
    result = preflight.run_preflight(
        args,
        pipeline_source=_PIPELINE_SOURCE,
        samplesheet_summary={"sample_count": 1, "unknown_columns": []},
    )
    assert result["ok"] is True


def test_profile_accepts_safe_institutional_executor_component(monkeypatch):
    monkeypatch.setattr(
        preflight,
        "_check_singularity_compatible_profile",
        lambda profile: {
            "profile": profile,
            "backend_path": "/usr/bin/singularity",
            "backend_ready": True,
        },
    )
    result = preflight._check_profile("singularity,slurm")

    assert result["profile"] == "singularity,slurm"
    assert result["backend_ready"] is True
    assert result["components"] == ["singularity", "slurm"]


def test_profile_rejects_unsafe_unknown_component():
    with pytest.raises(SkillError) as exc:
        preflight._check_profile("docker,bad profile")

    assert exc.value.error_code == "INVALID_PROFILE"
    assert exc.value.details["invalid_components"] == ["bad profile"]


@pytest.mark.parametrize(
    "checker,binary,code",
    [
        ("_check_docker_profile", "docker", "DOCKER_NOT_RUNNING"),
        ("_check_podman_profile", "podman", "PODMAN_NOT_RUNNING"),
    ],
)
def test_backend_daemon_down_message_matches_sibling_wording(monkeypatch, checker, binary, code):
    """When the backend is installed but its daemon/socket is down, the message must
    use the same "<binary> is installed but its service is not available." wording as
    the nfcore-rnaseq / nfcore-sarek wrappers (cross-wrapper consistency), and be
    internally consistent between docker and podman."""
    monkeypatch.setattr(preflight.shutil, "which", lambda name: f"/usr/bin/{name}")

    class _Info:
        returncode = 1

    monkeypatch.setattr(preflight.subprocess, "run", lambda *a, **k: _Info())
    with pytest.raises(SkillError) as exc:
        getattr(preflight, checker)(binary)
    assert exc.value.error_code == code
    assert exc.value.message == f"{binary.capitalize()} is installed but its service is not available."


def test_preflight_rejects_missing_nextflow_config(tmp_path, monkeypatch):
    args = _args(tmp_path)
    args.extra_config = [str(tmp_path / "missing.config")]
    Path(args.fasta).write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(args.gtf).write_text(
        'chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id "g1";\n', encoding="utf-8"
    )
    _mock_env(monkeypatch)

    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={"sample_count": 1, "unknown_columns": []},
        )

    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["field"] == "extra_config"


def test_fasta_path_with_whitespace_is_left_to_nf_schema(tmp_path, monkeypatch):
    args = _args(tmp_path)
    fasta_dir = tmp_path / "refs with spaces"
    fasta_dir.mkdir()
    fasta = fasta_dir / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    gtf.write_text('chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id "g1";\n', encoding="utf-8")
    args.fasta = str(fasta)
    args.gtf = str(gtf)
    _mock_env(monkeypatch)

    # The fasta pattern ^\S+\.fn?a(sta)?(\.gz)?$ covers the whole path, so nf-schema
    # refuses it when Nextflow starts.
    preflight.run_preflight(
        args,
        pipeline_source=_PIPELINE_SOURCE,
        samplesheet_summary={"sample_count": 1, "unknown_columns": []},
    )


def test_preflight_rejects_email_not_matching_nfcore_schema(tmp_path, monkeypatch):
    args = _args(tmp_path)
    args.email = "bad address@example.org"
    _mock_env(monkeypatch)

    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={"sample_count": 1, "unknown_columns": []},
        )

    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["field"] == "email"


@pytest.mark.parametrize("preset", ["cellranger", "cellrangerarc", "cellrangermulti"])
def test_missing_protocol_preserves_cellranger_auto_detection(
    tmp_path, monkeypatch, preset
):
    args = _args(tmp_path)
    args.preset = preset
    args.protocol = None
    args.genome = "GRCh38"
    args.fasta = None
    args.gtf = None
    _mock_env(monkeypatch)

    result = preflight.run_preflight(
        args,
        pipeline_source=_PIPELINE_SOURCE,
        samplesheet_summary={
            "sample_count": 1,
            "unknown_columns": [],
            "feature_types": ["gex"],
        },
    )

    assert result["ok"] is True


@pytest.mark.parametrize("preset", ["star", "kallisto"])
def test_scrnaseq_accepts_smartseq_protocol_for_documented_presets(
    tmp_path, monkeypatch, preset
):
    args = _args(tmp_path)
    args.preset = preset
    args.protocol = "smartseq"
    Path(args.fasta).write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(args.gtf).write_text(
        'chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id "g1";\n', encoding="utf-8"
    )
    _mock_env(monkeypatch)

    result = preflight.run_preflight(
        args,
        pipeline_source=_PIPELINE_SOURCE,
        samplesheet_summary={
            "sample_count": 1,
            "unknown_columns": [],
            "feature_types": ["gex"],
        },
    )

    assert result["ok"] is True


@pytest.mark.parametrize(
    "preset", ["standard", "cellranger", "cellrangerarc", "cellrangermulti"]
)
def test_scrnaseq_rejects_smartseq_protocol_for_unsupported_presets(
    tmp_path, monkeypatch, preset
):
    args = _args(tmp_path)
    args.preset = preset
    args.protocol = "smartseq"
    if preset in {"cellranger", "cellrangerarc", "cellrangermulti"}:
        args.genome = "GRCh38"
        args.fasta = None
        args.gtf = None
    _mock_env(monkeypatch)

    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={"sample_count": 1, "unknown_columns": []},
        )

    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert "Smart-seq" in exc.value.message
    assert exc.value.details["policy_source"] == POLICY_SOURCE_NFCORE_DOCS


@pytest.mark.parametrize(
    ("preset", "protocol"),
    [
        ("cellranger", "dropseq"),
        ("cellranger", "customchemistry"),
        ("cellrangerarc", "10XV3"),
        ("cellrangerarc", "customchemistry"),
    ],
)
def test_cellranger_family_rejects_protocols_outside_nfcore_matrix(
    tmp_path, monkeypatch, preset, protocol
):
    args = _args(tmp_path)
    args.preset = preset
    args.protocol = protocol
    args.genome = "GRCh38"
    args.fasta = None
    args.gtf = None
    _mock_env(monkeypatch)

    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={
                "sample_count": 1,
                "unknown_columns": [],
                "feature_types": ["gex"],
            },
        )

    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["preset"] == preset
    assert exc.value.details["protocol"] == protocol
    assert exc.value.details["policy_source"] == POLICY_SOURCE_NFCORE_DOCS


@pytest.mark.parametrize("preset", ["standard", "star", "kallisto"])
def test_custom_protocol_passthrough_is_kept_for_documented_aligners(
    tmp_path, monkeypatch, preset
):
    args = _args(tmp_path)
    args.preset = preset
    args.protocol = "customchemistry"
    Path(args.fasta).write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(args.gtf).write_text(
        'chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id "g1";\n', encoding="utf-8"
    )
    _mock_env(monkeypatch)

    result = preflight.run_preflight(
        args,
        pipeline_source=_PIPELINE_SOURCE,
        samplesheet_summary={
            "sample_count": 1,
            "unknown_columns": [],
            "feature_types": ["gex"],
        },
    )

    assert result["ok"] is True


def test_output_dir_with_only_reproducibility_bundle_is_accepted(tmp_path):
    # The reproducibility/ bundle is entirely wrapper-generated (commands.sh,
    # checksums.sha256, environment.yml, remap_paths.py, resource_limits.config,
    # provenance JSON, …), never user data, so an output dir containing only it must
    # NOT trip OUTPUT_DIR_NOT_EMPTY — a bundle-only dir must be re-runnable, and the
    # three wrappers must agree (nfcore-sarek already accepts it).
    out = tmp_path / "out"
    repro = out / "reproducibility"
    repro.mkdir(parents=True)
    for f in (
        "params.yaml", "samplesheet.valid.csv", "manifest.json", "commands.sh",
        "checksums.sha256", "environment.yml", "remap_paths.py", "resource_limits.config",
        "inputs.json",
    ):
        (repro / f).write_text("x", encoding="utf-8")
    preflight.check_output_dir_available(out, resume=False)


def test_output_dir_not_empty_raises(tmp_path, monkeypatch):
    out = tmp_path / "out"
    out.mkdir()
    (out / "result.json").write_text("{}", encoding="utf-8")
    args = _args(tmp_path)
    args.output = str(out)
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    gtf.write_text('chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id "g1";\n', encoding="utf-8")
    args.fasta = str(fasta)
    args.gtf = str(gtf)
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={"sample_count": 1, "unknown_columns": []},
        )
    assert exc.value.error_code == "OUTPUT_DIR_NOT_EMPTY"


def test_genome_shortcut_satisfies_reference_requirement(tmp_path, monkeypatch):
    """--genome (iGenomes shortcut) must satisfy the reference requirement for any preset."""
    args = _args(tmp_path)
    args.genome = "GRCh38"
    args.fasta = None
    args.gtf = None
    _mock_env(monkeypatch)
    result = preflight.run_preflight(
        args,
        pipeline_source=_PIPELINE_SOURCE,
        samplesheet_summary={"sample_count": 1, "unknown_columns": []},
    )
    assert result["ok"] is True


def test_cellrangermulti_antibody_capture_requires_feature_barcode_reference(
    tmp_path, monkeypatch
):
    args = _args(tmp_path)
    args.preset = "cellrangermulti"
    args.genome = "GRCh38"
    args.fasta = None
    args.gtf = None
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={
                "sample_count": 1,
                "unknown_columns": [],
                "feature_types": ["ab"],
            },
        )
    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["missing_field"] == "fb_reference"


def test_cellrangermulti_ffpe_probe_set_requires_barcodes_samplesheet(
    tmp_path, monkeypatch
):
    args = _args(tmp_path)
    args.preset = "cellrangermulti"
    args.genome = "GRCh38"
    args.fasta = None
    args.gtf = None
    probe_set = tmp_path / "probe_set.csv"
    probe_set.write_text("gene_id\n", encoding="utf-8")
    args.gex_frna_probe_set = str(probe_set)
    _mock_env(monkeypatch)

    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={
                "sample_count": 1,
                "unknown_columns": [],
                "feature_types": ["gex"],
            },
        )

    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["missing_field"] == "cellranger_multi_barcodes"
    assert exc.value.details["field"] == "gex_frna_probe_set"


def test_cellrangermulti_rejects_ffpe_and_cmo_combination(tmp_path, monkeypatch):
    args = _args(tmp_path)
    args.preset = "cellrangermulti"
    args.genome = "GRCh38"
    args.fasta = None
    args.gtf = None
    probe_set = tmp_path / "probe_set.csv"
    barcodes = tmp_path / "barcodes.csv"
    probe_set.write_text("gene_id\n", encoding="utf-8")
    barcodes.write_text("sample,cmo_ids\nsampleA,CMO301\n", encoding="utf-8")
    args.gex_frna_probe_set = str(probe_set)
    args.cellranger_multi_barcodes = str(barcodes)
    _mock_env(monkeypatch)

    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={
                "sample_count": 1,
                "unknown_columns": [],
                "feature_types": ["gex", "cmo"],
            },
        )

    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["conflicting_fields"] == ["gex_frna_probe_set", "cmo"]
    assert exc.value.details["policy_source"] == POLICY_SOURCE_NFCORE_DOCS


def _write_multi_barcodes(path: Path, body_rows: str) -> None:
    """Write a cellranger_multi_barcodes CSV with the documented 4.1.0 header."""
    path.write_text(
        "sample,multiplexed_sample_id,probe_barcode_ids,cmo_ids,ocm_ids,description\n"
        + body_rows,
        encoding="utf-8",
    )


def test_cellrangermulti_rejects_mixed_modes_in_one_barcodes_row(tmp_path, monkeypatch):
    """nf-core/scrnaseq usage docs: the --cellranger-multi-barcodes samplesheet
    encodes the multiplexing mode per sample via probe_barcode_ids (FFPE), cmo_ids
    (CMO) and ocm_ids (OCM), which are mutually exclusive. A single sample mixing
    FFPE and OCM must fail fast in preflight (audit F-2)."""
    args = _args(tmp_path)
    args.preset = "cellrangermulti"
    args.genome = "GRCh38"
    args.fasta = None
    args.gtf = None
    barcodes = tmp_path / "barcodes.csv"
    _write_multi_barcodes(barcodes, "S1,S1_a,BC001,,OB1,desc\n")  # FFPE + OCM
    args.cellranger_multi_barcodes = str(barcodes)
    _mock_env(monkeypatch)

    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={
                "sample_count": 1,
                "unknown_columns": [],
                "feature_types": ["gex"],
            },
        )

    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["conflicting_modes"] == ["FFPE", "OCM"]
    assert exc.value.details["sample"] == "S1"
    assert exc.value.details["policy_source"] == POLICY_SOURCE_NFCORE_DOCS


def test_cellrangermulti_rejects_mixed_modes_across_sample_rows(tmp_path, monkeypatch):
    """Two rows of the same physical sample using different modes (CMO then OCM) is
    rejected — 'using more than one for a single sample' (audit F-2)."""
    args = _args(tmp_path)
    args.preset = "cellrangermulti"
    args.genome = "GRCh38"
    args.fasta = None
    args.gtf = None
    barcodes = tmp_path / "barcodes.csv"
    _write_multi_barcodes(barcodes, "S1,S1_a,,CMO301,,desc\nS1,S1_b,,,OB2,desc\n")
    args.cellranger_multi_barcodes = str(barcodes)
    _mock_env(monkeypatch)

    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={
                "sample_count": 1,
                "unknown_columns": [],
                "feature_types": ["gex"],
            },
        )

    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["conflicting_modes"] == ["CMO", "OCM"]


def test_cellrangermulti_allows_single_mode_ocm_barcodes(tmp_path, monkeypatch):
    """A barcodes samplesheet using a single mode (OCM only, across rows) must NOT
    be rejected (audit F-2)."""
    args = _args(tmp_path)
    args.preset = "cellrangermulti"
    args.genome = "GRCh38"
    args.fasta = None
    args.gtf = None
    barcodes = tmp_path / "barcodes.csv"
    _write_multi_barcodes(barcodes, "S1,S1_a,,,OB1,desc\nS1,S1_b,,,OB2,desc\n")
    args.cellranger_multi_barcodes = str(barcodes)
    _mock_env(monkeypatch)

    result = preflight.run_preflight(
        args,
        pipeline_source=_PIPELINE_SOURCE,
        samplesheet_summary={
            "sample_count": 1,
            "unknown_columns": [],
            "feature_types": ["gex"],
        },
    )

    assert result["ok"] is True


def test_cellrangermulti_barcodes_check_tolerates_unparseable_csv(tmp_path, monkeypatch):
    """A barcodes CSV the wrapper cannot parse locally is deferred to upstream
    nf-schema rather than raising a spurious wrapper error (audit F-2)."""
    args = _args(tmp_path)
    args.preset = "cellrangermulti"
    args.genome = "GRCh38"
    args.fasta = None
    args.gtf = None
    barcodes = tmp_path / "barcodes.csv"
    # No recognised mode columns at all → nothing to validate, must not raise.
    barcodes.write_text("sample,multiplexed_sample_id\nS1,S1_a\n", encoding="utf-8")
    args.cellranger_multi_barcodes = str(barcodes)
    _mock_env(monkeypatch)

    result = preflight.run_preflight(
        args,
        pipeline_source=_PIPELINE_SOURCE,
        samplesheet_summary={
            "sample_count": 1,
            "unknown_columns": [],
            "feature_types": ["gex"],
        },
    )

    assert result["ok"] is True


def test_genome_and_fasta_mutually_exclusive(tmp_path, monkeypatch):
    """--genome and --fasta together must raise CONFLICTING_REFERENCES before Nextflow runs."""
    args = _args(tmp_path)
    args.genome = "GRCh38"
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    args.fasta = str(fasta)
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={"sample_count": 1, "unknown_columns": []},
        )
    assert exc.value.error_code == "CONFLICTING_REFERENCES"


# ── mamba profile ─────────────────────────────────────────────────────────────


# ── podman profile ────────────────────────────────────────────────────────────


# ── HPC profiles (shifter, charliecloud) ─────────────────────────────────────


# ── macOS + Docker + /tmp warning ─────────────────────────────────────────────


def test_seq_center_with_non_star_preset_warns(tmp_path, monkeypatch, capsys):
    """seq_center only affects STARsolo BAM read groups in nf-core/scrnaseq 4.1.0;
    forwarding it under another aligner is a harmless no-op, so the wrapper warns
    rather than silently ignoring it (audit F-6)."""
    args = _args(tmp_path)
    args.preset = "kallisto"
    args.seq_center = "SequencingCenterX"
    Path(args.fasta).write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(args.gtf).write_text(
        'chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id "g1";\n', encoding="utf-8"
    )
    _mock_env(monkeypatch)

    preflight.run_preflight(
        args,
        pipeline_source=_PIPELINE_SOURCE,
        samplesheet_summary={"sample_count": 1, "unknown_columns": []},
    )

    err = capsys.readouterr().err
    assert "seq_center" in err and "STARsolo" in err


def test_seq_center_with_star_preset_does_not_warn(tmp_path, monkeypatch, capsys):
    args = _args(tmp_path)  # preset defaults to star
    args.seq_center = "SequencingCenterX"
    Path(args.fasta).write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(args.gtf).write_text(
        'chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id "g1";\n', encoding="utf-8"
    )
    _mock_env(monkeypatch)

    preflight.run_preflight(
        args,
        pipeline_source=_PIPELINE_SOURCE,
        samplesheet_summary={"sample_count": 1, "unknown_columns": []},
    )

    assert "seq_center" not in capsys.readouterr().err


# nf-core/scrnaseq 4.1.0's nextflow_schema.json declares every reference and auxiliary
# path the wrapper exposes `exists: true` except motifs and cellrangerarc_config, and
# gives fasta the ^\S+\.fn?a(sta)?(\.gz)?$ pattern: nf-schema refuses a missing path or
# a wrong FASTA name when Nextflow starts, before any task runs.
_SCHEMA_CHECKED_PATHS = (
    "fasta", "gtf", "transcript_fasta", "txp2gene", "simpleaf_index", "kallisto_index",
    "star_index", "cellranger_index", "barcode_whitelist", "kb_t1c", "kb_t2c",
    "cellranger_vdj_index", "gex_frna_probe_set", "gex_target_panel", "gex_cmo_set",
    "fb_reference", "vdj_inner_enrichment_primers", "gex_barcode_sample_assignment",
    "cellranger_multi_barcodes",
)


@pytest.mark.parametrize("field", _SCHEMA_CHECKED_PATHS)
def test_paths_the_schema_checks_are_left_to_nf_schema(tmp_path, field):
    preflight._check_reference_paths_exist({field: str(tmp_path / "missing" / field)})


def test_a_fasta_name_the_schema_refuses_is_left_to_nf_schema(tmp_path, monkeypatch):
    args = _args(tmp_path)
    for name in ("genome.txt", "genes.gtf"):
        (tmp_path / name).write_text("x", encoding="utf-8")
    args.fasta = str(tmp_path / "genome.txt")
    _mock_env(monkeypatch)
    preflight.run_preflight(args, pipeline_source=_PIPELINE_SOURCE,
                            samplesheet_summary={"sample_count": 1, "unknown_columns": []})


@pytest.mark.parametrize("field", ["motifs", "cellrangerarc_config"])
def test_paths_the_schema_does_not_check_must_exist(tmp_path, field):
    """No `exists` for these in 4.1.0: workflows/scrnaseq.nf reads motifs with a plain
    file(), so a missing local file would fail only inside the ARC task."""
    with pytest.raises(SkillError) as exc:
        preflight._check_reference_paths_exist({field: str(tmp_path / f"missing_{field}")})
    assert exc.value.error_code == "MISSING_REFERENCE"
    assert exc.value.details["field"] == field


# The minimum versions are Nextflow's to enforce: each pinned release declares
# manifest.nextflowVersion with `!`, so Nextflow refuses an older engine when the run
# starts, and the Nextflow launcher refuses a Java outside 17-26 before that. The wrapper
# checks only that `java` and `nextflow` run and records their versions for the bundle.
@pytest.mark.parametrize("java_text,nextflow_text", [
    ('openjdk version "11.0.1"', "Nextflow version 24.10.0 build 5928"),
    ('openjdk version "16.0.2"', "Nextflow version 25.04.2 build 5947"),
])
def test_old_java_and_nextflow_are_left_to_nextflow(monkeypatch, java_text, nextflow_text):
    import preflight as p
    monkeypatch.setattr(p.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(p, "_command_output", lambda args: java_text if args[0] == "java" else nextflow_text)
    assert p._check_java()["version"] == java_text.split('"')[1]
    assert p._check_nextflow()["version"] == nextflow_text.split()[2]


# nf-core/scrnaseq 4.1.0 refuses protocol `auto` for every aligner but the Cell Ranger
# family itself, at the start of workflows/scrnaseq.nf ("Only cellranger supports
# `protocol = 'auto'`"), and `auto` is its default, so a run without --protocol lands
# there too. The wrapper keeps no copy of that rule; params_builder writes an `auto`
# alias in the exact form the pipeline compares.
@pytest.mark.parametrize("preset", ["standard", "star", "kallisto"])
@pytest.mark.parametrize("protocol", [None, "auto", "AUTO", "Auto"])
def test_auto_and_missing_protocol_are_left_to_the_pipeline(tmp_path, monkeypatch, preset, protocol):
    import params_builder

    args = _args(tmp_path)
    args.preset = preset
    args.protocol = protocol
    _mock_env(monkeypatch)
    result = preflight.run_preflight(
        args,
        pipeline_source=_PIPELINE_SOURCE,
        samplesheet_summary={"sample_count": 1, "unknown_columns": []},
    )
    assert result["ok"] is True
    if protocol is not None:
        assert params_builder._canonical_protocol(protocol) == "auto"
