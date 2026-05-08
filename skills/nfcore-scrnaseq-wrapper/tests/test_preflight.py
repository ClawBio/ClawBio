from __future__ import annotations

from argparse import Namespace
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import preflight
from errors import SkillError

_PIPELINE_SOURCE = {
    "source_kind": "local_checkout",
    "source_ref": "scrnaseq-checkout",
    "resolved_version": "abc123",
    "branch": "main",
    "dirty": False,
}


def _mock_env(monkeypatch):
    monkeypatch.setattr(preflight, "_check_java", lambda: {"path": "/usr/bin/java", "version": "17.0.8"})
    monkeypatch.setattr(preflight, "_check_nextflow", lambda: {"path": "/usr/bin/nextflow", "version": "25.04.0"})
    monkeypatch.setattr(preflight, "_check_profile", lambda profile: {"profile": profile, "backend_path": "/usr/bin/docker", "backend_ready": True})


def _args(tmp_path: Path) -> Namespace:
    return Namespace(
        output=str(tmp_path / "out"),
        resume=False,
        preset="star",
        protocol="10XV3",
        profile="docker",
        demo=False,
        email=None,
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


def test_preflight_happy_path(tmp_path, monkeypatch):
    args = _args(tmp_path)
    Path(args.fasta).write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(args.gtf).write_text("chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"g1\";\n", encoding="utf-8")
    _mock_env(monkeypatch)
    result = preflight.run_preflight(
        args,
        pipeline_source=_PIPELINE_SOURCE,
        samplesheet_summary={"sample_count": 1, "unknown_columns": []},
    )
    assert result["ok"] is True


def test_preflight_rejects_missing_reference(tmp_path, monkeypatch):
    args = _args(tmp_path)
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={"sample_count": 1, "unknown_columns": []},
        )
    assert exc.value.error_code == "MISSING_REFERENCE"


def test_fasta_path_with_whitespace_fails_before_nextflow_schema(tmp_path, monkeypatch):
    args = _args(tmp_path)
    fasta_dir = tmp_path / "refs with spaces"
    fasta_dir.mkdir()
    fasta = fasta_dir / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    gtf.write_text("chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"g1\";\n", encoding="utf-8")
    args.fasta = str(fasta)
    args.gtf = str(gtf)
    _mock_env(monkeypatch)

    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={"sample_count": 1, "unknown_columns": []},
        )

    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["field"] == "fasta"


def test_preflight_rejects_fasta_extension_not_matching_nfcore_schema(tmp_path, monkeypatch):
    args = _args(tmp_path)
    fasta = tmp_path / "genome.txt"
    gtf = tmp_path / "genes.gtf"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    gtf.write_text("chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"g1\";\n", encoding="utf-8")
    args.fasta = str(fasta)
    args.gtf = str(gtf)
    _mock_env(monkeypatch)

    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={"sample_count": 1, "unknown_columns": []},
        )

    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["field"] == "fasta"
    assert exc.value.details["schema_pattern"] == r"^\S+\.fn?a(sta)?(\.gz)?$"


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


def test_standard_preset_requires_explicit_non_auto_protocol(tmp_path, monkeypatch):
    args = _args(tmp_path)
    args.preset = "standard"
    args.protocol = None
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={"sample_count": 1, "unknown_columns": []},
        )
    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["protocol"] == "auto"


def test_standard_simpleaf_rejects_smartseq_protocol(tmp_path, monkeypatch):
    args = _args(tmp_path)
    args.preset = "standard"
    args.protocol = "smartseq"
    _mock_env(monkeypatch)

    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={"sample_count": 1, "unknown_columns": []},
        )

    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details == {"preset": "standard", "protocol": "smartseq"}


def test_star_accepts_smartseq_protocol_supported_by_upstream(tmp_path, monkeypatch):
    args = _args(tmp_path)
    args.preset = "star"
    args.protocol = "smartseq"
    Path(args.fasta).write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(args.gtf).write_text("chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"g1\";\n", encoding="utf-8")
    _mock_env(monkeypatch)

    result = preflight.run_preflight(
        args,
        pipeline_source=_PIPELINE_SOURCE,
        samplesheet_summary={"sample_count": 1, "unknown_columns": []},
    )

    assert result["ok"] is True


def test_cellrangerarc_accepts_explicit_protocol(tmp_path, monkeypatch):
    """cellrangerarc must not hard-block non-auto protocols — the pipeline only warns and
    passes them verbatim to the aligner (Utils.groovy:17-18)."""
    args = _args(tmp_path)
    args.preset = "cellrangerarc"
    args.protocol = "10XV3"
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


def test_cellrangerarc_accepts_auto_protocol(tmp_path, monkeypatch):
    """cellrangerarc must accept protocol=auto — this was the only allowed value before
    F4 fix and must remain valid after the hard-block was removed."""
    args = _args(tmp_path)
    args.preset = "cellrangerarc"
    args.protocol = "auto"
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


def test_cellranger_accepts_known_10x_protocol(tmp_path, monkeypatch):
    """cellranger passes preflight with a standard 10x protocol."""
    args = _args(tmp_path)
    args.preset = "cellranger"
    args.protocol = "10XV3"
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


def test_cellranger_accepts_auto_protocol(tmp_path, monkeypatch):
    """cellranger must accept protocol=auto — the pipeline auto-detects chemistry.
    This was previously hard-blocked before F2 fix (Utils.groovy:17-18)."""
    args = _args(tmp_path)
    args.preset = "cellranger"
    args.protocol = "auto"
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


def test_cellranger_accepts_custom_protocol_string(tmp_path, monkeypatch):
    """cellranger must not hard-block unrecognized protocols — the pipeline only warns and
    passes them verbatim to the aligner (Utils.groovy:17-18)."""
    args = _args(tmp_path)
    args.preset = "cellranger"
    args.protocol = "dropseq"
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



def test_missing_conda_uses_correct_error_code(monkeypatch):
    import shutil
    monkeypatch.setattr(shutil, "which", lambda name: None)
    with pytest.raises(SkillError) as exc:
        preflight._check_profile("conda")
    assert exc.value.error_code == "MISSING_CONDA"


def test_missing_singularity_uses_correct_error_code(monkeypatch):
    import shutil
    monkeypatch.setattr(shutil, "which", lambda name: None)
    with pytest.raises(SkillError) as exc:
        preflight._check_profile("singularity")
    assert exc.value.error_code == "MISSING_SINGULARITY"


def test_missing_apptainer_uses_correct_error_code(monkeypatch):
    import shutil
    monkeypatch.setattr(shutil, "which", lambda name: None)
    with pytest.raises(SkillError) as exc:
        preflight._check_profile("apptainer")
    assert exc.value.error_code == "MISSING_SINGULARITY"


def test_singularity_profile_falls_back_to_apptainer_binary(monkeypatch):
    """--profile singularity must accept an apptainer binary (they are compatible)."""
    import shutil
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/apptainer" if name == "apptainer" else None)
    result = preflight._check_profile("singularity")
    assert result["backend_ready"] is True
    assert result["backend_path"] == "/usr/bin/apptainer"


def test_apptainer_profile_falls_back_to_singularity_binary(monkeypatch):
    """--profile apptainer must accept a singularity binary (they are compatible)."""
    import shutil
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/singularity" if name == "singularity" else None)
    result = preflight._check_profile("apptainer")
    assert result["backend_ready"] is True
    assert result["backend_path"] == "/usr/bin/singularity"


def test_singularity_error_mentions_both_binaries(monkeypatch):
    """Error details must list both singularity and apptainer so the user knows what to install."""
    import shutil
    monkeypatch.setattr(shutil, "which", lambda name: None)
    with pytest.raises(SkillError) as exc:
        preflight._check_profile("singularity")
    tried = exc.value.details.get("tried", [])
    assert "singularity" in tried
    assert "apptainer" in tried


def test_output_dir_not_empty_raises(tmp_path, monkeypatch):
    out = tmp_path / "out"
    out.mkdir()
    (out / "result.json").write_text("{}", encoding="utf-8")
    args = _args(tmp_path)
    args.output = str(out)
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    gtf.write_text("chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"g1\";\n", encoding="utf-8")
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


def test_parse_version_tuple_handles_java_format():
    assert preflight._parse_version_tuple('openjdk version "17.0.8" 2023-07-18') == (17, 0, 8)


def test_parse_version_tuple_handles_nextflow_format():
    assert preflight._parse_version_tuple("      N E X T F L O W\n      version 25.04.0 build 5940") == (25, 4, 0)


def test_parse_version_tuple_returns_empty_for_no_version():
    assert preflight._parse_version_tuple("no version here at all") == ()


def test_unsupported_profile_raises_invalid_profile(monkeypatch):
    import shutil
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/docker")
    with pytest.raises(SkillError) as exc:
        preflight._check_profile("unsupported_backend")
    assert exc.value.error_code == "INVALID_PROFILE"


def test_command_output_returns_empty_on_timeout(monkeypatch):
    import subprocess as sp
    def _raise(*a, **kw):
        raise sp.TimeoutExpired(cmd=["java"], timeout=10)
    monkeypatch.setattr(sp, "run", _raise)
    assert preflight._command_output(["java", "-version"]) == ""


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


def test_cellrangerarc_preset_accepted_by_preflight(tmp_path, monkeypatch):
    """cellrangerarc passes preflight when its reference and ARC-specific files are provided."""
    args = _args(tmp_path)
    args.preset = "cellrangerarc"
    args.protocol = None
    args.genome = "GRCh38"
    args.fasta = None
    args.gtf = None
    arc_config = tmp_path / "arc_config.json"
    motifs = tmp_path / "motifs.pfm"
    arc_config.write_text("{}", encoding="utf-8")
    motifs.write_text("MOTIF\n", encoding="utf-8")
    args.cellrangerarc_config = str(arc_config)
    args.motifs = str(motifs)
    args.cellrangerarc_reference = "GRCh38-2024-A"
    _mock_env(monkeypatch)
    result = preflight.run_preflight(
        args,
        pipeline_source=_PIPELINE_SOURCE,
        samplesheet_summary={"sample_count": 1, "unknown_columns": []},
    )
    assert result["ok"] is True


def test_cellrangermulti_preset_accepted_by_preflight(tmp_path, monkeypatch):
    """cellrangermulti must be a recognised preset that passes preflight when reference is provided."""
    args = _args(tmp_path)
    args.preset = "cellrangermulti"
    args.genome = "GRCh38"
    args.fasta = None
    args.gtf = None
    barcodes = tmp_path / "multi_barcodes.csv"
    barcodes.write_text("sample,barcode\n", encoding="utf-8")
    args.cellranger_multi_barcodes = str(barcodes)
    _mock_env(monkeypatch)
    result = preflight.run_preflight(
        args,
        pipeline_source=_PIPELINE_SOURCE,
        samplesheet_summary={"sample_count": 1, "unknown_columns": []},
    )
    assert result["ok"] is True


def test_cellrangerarc_reference_alone_does_not_satisfy_requirement(tmp_path, monkeypatch):
    """cellrangerarc_reference is symbolic metadata, not a reference/index by itself."""
    args = _args(tmp_path)
    args.preset = "cellrangerarc"
    args.protocol = None
    args.genome = None
    args.fasta = None
    args.gtf = None
    args.cellrangerarc_reference = "GRCh38-2024-A"
    arc_config = tmp_path / "arc_config.json"
    motifs = tmp_path / "motifs.pfm"
    arc_config.write_text("{}", encoding="utf-8")
    motifs.write_text("MOTIF\n", encoding="utf-8")
    args.cellrangerarc_config = str(arc_config)
    args.motifs = str(motifs)
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={"sample_count": 1, "unknown_columns": []},
        )
    assert exc.value.error_code == "MISSING_REFERENCE"


def test_cellrangerarc_missing_config_fails_preflight(tmp_path, monkeypatch):
    """ARC-specific files must be validated before launching Nextflow."""
    args = _args(tmp_path)
    args.preset = "cellrangerarc"
    args.protocol = None
    args.genome = "GRCh38"
    args.fasta = None
    args.gtf = None
    args.cellrangerarc_reference = "GRCh38-2024-A"
    motifs = tmp_path / "motifs.pfm"
    motifs.write_text("MOTIF\n", encoding="utf-8")
    args.motifs = str(motifs)
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={"sample_count": 1, "unknown_columns": []},
        )
    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["missing_field"] == "cellrangerarc_config"


def test_cellrangerarc_missing_reference_name_fails_preflight(tmp_path, monkeypatch):
    """ARC build mode needs the symbolic reference name used by the ARC config."""
    args = _args(tmp_path)
    args.preset = "cellrangerarc"
    args.protocol = None
    args.genome = "GRCh38"
    args.fasta = None
    args.gtf = None
    arc_config = tmp_path / "arc_config.json"
    motifs = tmp_path / "motifs.pfm"
    arc_config.write_text("{}", encoding="utf-8")
    motifs.write_text("MOTIF\n", encoding="utf-8")
    args.cellrangerarc_config = str(arc_config)
    args.motifs = str(motifs)
    _mock_env(monkeypatch)
    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={"sample_count": 1, "unknown_columns": []},
        )
    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["missing_field"] == "cellrangerarc_reference"


def test_cellrangerarc_motifs_are_optional_when_building_reference(tmp_path, monkeypatch):
    """nf-core/scrnaseq treats motifs as optional for CellRanger ARC mkref."""
    args = _args(tmp_path)
    args.preset = "cellrangerarc"
    args.protocol = None
    args.genome = "GRCh38"
    args.fasta = None
    args.gtf = None
    args.cellrangerarc_reference = "GRCh38-2024-A"
    arc_config = tmp_path / "arc_config.json"
    arc_config.write_text("{}", encoding="utf-8")
    args.cellrangerarc_config = str(arc_config)
    # motifs intentionally omitted
    _mock_env(monkeypatch)
    result = preflight.run_preflight(
        args,
        pipeline_source=_PIPELINE_SOURCE,
        samplesheet_summary={"sample_count": 1, "unknown_columns": []},
    )
    assert result["ok"] is True


def test_cellrangerarc_prebuilt_index_skips_build_files(tmp_path, monkeypatch):
    """A prebuilt ARC/CellRanger index should not require motifs/config build inputs."""
    args = _args(tmp_path)
    args.preset = "cellrangerarc"
    args.protocol = None
    args.fasta = None
    args.gtf = None
    index = tmp_path / "cellranger_arc_index"
    index.mkdir()
    args.cellranger_index = str(index)
    _mock_env(monkeypatch)
    result = preflight.run_preflight(
        args,
        pipeline_source=_PIPELINE_SOURCE,
        samplesheet_summary={"sample_count": 1, "unknown_columns": []},
    )
    assert result["ok"] is True


def test_cellrangermulti_without_multiplexed_cmo_barcodes_passes_preflight(tmp_path, monkeypatch):
    args = _args(tmp_path)
    args.preset = "cellrangermulti"
    args.genome = "GRCh38"
    args.fasta = None
    args.gtf = None
    _mock_env(monkeypatch)
    result = preflight.run_preflight(
        args,
        pipeline_source=_PIPELINE_SOURCE,
        samplesheet_summary={"sample_count": 1, "unknown_columns": [], "feature_types": ["gex", "vdj"]},
    )
    assert result["ok"] is True


def test_cellrangermulti_cmo_requires_barcodes(tmp_path, monkeypatch):
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
            samplesheet_summary={"sample_count": 1, "unknown_columns": [], "feature_types": ["cmo"]},
        )
    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["missing_field"] == "cellranger_multi_barcodes"


def test_cellrangermulti_antibody_capture_requires_feature_barcode_reference(tmp_path, monkeypatch):
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
            samplesheet_summary={"sample_count": 1, "unknown_columns": [], "feature_types": ["ab"]},
        )
    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["missing_field"] == "fb_reference"


def test_cellrangermulti_vdj_rows_reject_skip_vdjref_without_prebuilt_index(tmp_path, monkeypatch):
    args = _args(tmp_path)
    args.preset = "cellrangermulti"
    args.genome = "GRCh38"
    args.fasta = None
    args.gtf = None
    args.skip_cellrangermulti_vdjref = True
    _mock_env(monkeypatch)

    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={"sample_count": 1, "unknown_columns": [], "feature_types": ["gex", "vdj"]},
        )

    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["feature_type"] == "vdj"
    assert exc.value.details["missing_field"] == "cellranger_vdj_index"


def test_cellrangermulti_rejects_cmo_and_ffpe_probe_set_together(tmp_path, monkeypatch):
    args = _args(tmp_path)
    args.preset = "cellrangermulti"
    args.genome = "GRCh38"
    args.fasta = None
    args.gtf = None
    probe_set = tmp_path / "probe_set.csv"
    barcodes = tmp_path / "multi_barcodes.csv"
    probe_set.write_text("gene_id\n", encoding="utf-8")
    barcodes.write_text("sample,barcode\n", encoding="utf-8")
    args.gex_frna_probe_set = str(probe_set)
    args.cellranger_multi_barcodes = str(barcodes)
    _mock_env(monkeypatch)

    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            args,
            pipeline_source=_PIPELINE_SOURCE,
            samplesheet_summary={"sample_count": 1, "unknown_columns": [], "feature_types": ["gex", "cmo"]},
        )

    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["active_modes"] == ["cmo", "ffpe"]


def test_cellrangermulti_ffpe_probe_set_requires_barcodes_samplesheet(tmp_path, monkeypatch):
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
            samplesheet_summary={"sample_count": 1, "unknown_columns": [], "feature_types": ["gex"]},
        )

    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["missing_field"] == "cellranger_multi_barcodes"
    assert exc.value.details["field"] == "gex_frna_probe_set"


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


def test_docker_info_timeout_raises_docker_not_running(monkeypatch):
    import shutil
    import subprocess as sp
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/docker" if name == "docker" else None)
    def _raise(*a, **kw):
        raise sp.TimeoutExpired(cmd=["docker"], timeout=30)
    monkeypatch.setattr(sp, "run", _raise)
    with pytest.raises(SkillError) as exc:
        preflight._check_profile("docker")
    assert exc.value.error_code == "DOCKER_NOT_RUNNING"


def test_check_output_dir_available_allows_params_yaml(tmp_path):
    """reproducibility/params.yaml must not trigger OUTPUT_DIR_NOT_EMPTY."""
    from preflight import check_output_dir_available
    repro = tmp_path / "reproducibility"
    repro.mkdir()
    (repro / "params.yaml").write_text("aligner: star\n", encoding="utf-8")
    check_output_dir_available(tmp_path, resume=False)


def test_check_output_dir_available_allows_demo_samplesheet(tmp_path):
    """reproducibility/samplesheet.demo.csv must not trigger OUTPUT_DIR_NOT_EMPTY."""
    from preflight import check_output_dir_available
    repro = tmp_path / "reproducibility"
    repro.mkdir()
    (repro / "samplesheet.demo.csv").write_text("sample,fastq_1\n", encoding="utf-8")
    check_output_dir_available(tmp_path, resume=False)


# ---------------------------------------------------------------------------
# _prompt_for_genome
# ---------------------------------------------------------------------------

def test_common_genomes_is_nonempty_sequence():
    from schemas import COMMON_GENOMES
    assert len(COMMON_GENOMES) > 0
    assert all(isinstance(g, str) and g for g in COMMON_GENOMES)


def test_prompt_for_genome_returns_selection_by_number(monkeypatch):
    """User entering a valid index returns the corresponding genome string."""
    from preflight import _prompt_for_genome
    from schemas import COMMON_GENOMES
    monkeypatch.setattr("sys.stdin.isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda _prompt="": "1")
    result = _prompt_for_genome()
    assert result == COMMON_GENOMES[0]


def test_prompt_for_genome_returns_selection_by_name(monkeypatch):
    """User entering a genome name directly returns that name."""
    from preflight import _prompt_for_genome
    from schemas import COMMON_GENOMES
    monkeypatch.setattr("sys.stdin.isatty", lambda: True)
    target = COMMON_GENOMES[0]
    monkeypatch.setattr("builtins.input", lambda _prompt="": target)
    result = _prompt_for_genome()
    assert result == target


def test_prompt_for_genome_returns_none_in_noninteractive(monkeypatch):
    """When stdin is not a TTY the prompt returns None rather than blocking."""
    from preflight import _prompt_for_genome
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    result = _prompt_for_genome()
    assert result is None


def test_prompt_for_genome_returns_none_on_empty_input(monkeypatch):
    """Pressing Enter (empty input) returns None so the caller can raise MISSING_REFERENCE."""
    from preflight import _prompt_for_genome
    monkeypatch.setattr("sys.stdin.isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda _prompt="": "")
    result = _prompt_for_genome()
    assert result is None


def test_prompt_for_genome_returns_none_on_invalid_number(monkeypatch):
    """An out-of-range number returns None so the caller can raise MISSING_REFERENCE."""
    from preflight import _prompt_for_genome
    monkeypatch.setattr("sys.stdin.isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda _prompt="": "9999")
    result = _prompt_for_genome()
    assert result is None


def test_cellrangermulti_accepts_custom_protocol_string(tmp_path, monkeypatch):
    """cellrangermulti internally uses the cellranger protocol map (Utils.groovy:12).
    Like cellranger after F2 fix, it must accept any protocol string and delegate
    validation to the pipeline."""
    args = _args(tmp_path)
    args.preset = "cellrangermulti"
    args.protocol = "dropseq"
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


# ── mamba profile ─────────────────────────────────────────────────────────────

def test_mamba_profile_accepted_when_mamba_installed(monkeypatch):
    import shutil
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/mamba" if name == "mamba" else None)
    result = preflight._check_profile("mamba")
    assert result["backend_ready"] is True
    assert result["backend_path"] == "/usr/bin/mamba"


def test_mamba_profile_falls_back_to_conda_binary(monkeypatch):
    import shutil
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/conda" if name == "conda" else None)
    result = preflight._check_profile("mamba")
    assert result["backend_ready"] is True
    assert result["backend_path"] == "/usr/bin/conda"


def test_missing_mamba_uses_missing_conda_error(monkeypatch):
    import shutil
    monkeypatch.setattr(shutil, "which", lambda name: None)
    with pytest.raises(SkillError) as exc:
        preflight._check_profile("mamba")
    assert exc.value.error_code == "MISSING_CONDA"


# ── podman profile ────────────────────────────────────────────────────────────

def test_missing_podman_raises_missing_podman_error(monkeypatch):
    import shutil
    monkeypatch.setattr(shutil, "which", lambda name: None)
    with pytest.raises(SkillError) as exc:
        preflight._check_profile("podman")
    assert exc.value.error_code == "MISSING_PODMAN"


def test_podman_not_running_raises_podman_not_running(monkeypatch):
    import shutil
    import subprocess as sp
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/podman" if name == "podman" else None)
    def _raise(*a, **kw):
        raise sp.TimeoutExpired(cmd=["podman"], timeout=30)
    monkeypatch.setattr(sp, "run", _raise)
    with pytest.raises(SkillError) as exc:
        preflight._check_profile("podman")
    assert exc.value.error_code == "PODMAN_NOT_RUNNING"


def test_podman_profile_accepted_when_running(monkeypatch):
    import shutil
    import subprocess as sp
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/podman" if name == "podman" else None)
    mock_result = sp.CompletedProcess(args=["podman", "info"], returncode=0, stdout="", stderr="")
    monkeypatch.setattr(sp, "run", lambda *a, **kw: mock_result)
    result = preflight._check_profile("podman")
    assert result["backend_ready"] is True
    assert result["backend_path"] == "/usr/bin/podman"


# ── HPC profiles (shifter, charliecloud) ─────────────────────────────────────

def test_missing_shifter_raises_missing_hpc_runtime(monkeypatch):
    import shutil
    monkeypatch.setattr(shutil, "which", lambda name: None)
    with pytest.raises(SkillError) as exc:
        preflight._check_profile("shifter")
    assert exc.value.error_code == "MISSING_HPC_RUNTIME"


def test_shifter_profile_accepted_when_installed(monkeypatch):
    import shutil
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/shifter" if name == "shifter" else None)
    result = preflight._check_profile("shifter")
    assert result["backend_ready"] is True
    assert result["backend_path"] == "/usr/bin/shifter"


def test_missing_charliecloud_raises_missing_hpc_runtime(monkeypatch):
    import shutil
    monkeypatch.setattr(shutil, "which", lambda name: None)
    with pytest.raises(SkillError) as exc:
        preflight._check_profile("charliecloud")
    assert exc.value.error_code == "MISSING_HPC_RUNTIME"


def test_charliecloud_profile_accepted_when_ch_run_installed(monkeypatch):
    import shutil
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/ch-run" if name == "ch-run" else None)
    result = preflight._check_profile("charliecloud")
    assert result["backend_ready"] is True
    assert result["backend_path"] == "/usr/bin/ch-run"


# ── macOS + Docker + /tmp warning ─────────────────────────────────────────────

def test_macos_docker_tmp_emits_warning(monkeypatch, capsys):
    monkeypatch.setattr(preflight.sys, "platform", "darwin")
    preflight._warn_if_macos_docker_tmp("docker", Path("/tmp/some_run"))
    captured = capsys.readouterr()
    assert "WARNING" in captured.err
    assert "/tmp" in captured.err


def test_macos_docker_private_tmp_emits_warning(monkeypatch, capsys):
    monkeypatch.setattr(preflight.sys, "platform", "darwin")
    preflight._warn_if_macos_docker_tmp("docker", Path("/private/tmp/some_run"))
    captured = capsys.readouterr()
    assert "WARNING" in captured.err


def test_macos_docker_home_no_warning(monkeypatch, capsys):
    monkeypatch.setattr(preflight.sys, "platform", "darwin")
    preflight._warn_if_macos_docker_tmp("docker", Path("/Users/someone/my_run"))
    captured = capsys.readouterr()
    assert captured.err == ""


def test_linux_docker_tmp_no_warning(monkeypatch, capsys):
    monkeypatch.setattr(preflight.sys, "platform", "linux")
    preflight._warn_if_macos_docker_tmp("docker", Path("/tmp/some_run"))
    captured = capsys.readouterr()
    assert captured.err == ""


def test_macos_singularity_tmp_no_warning(monkeypatch, capsys):
    monkeypatch.setattr(preflight.sys, "platform", "darwin")
    preflight._warn_if_macos_docker_tmp("singularity", Path("/tmp/some_run"))
    captured = capsys.readouterr()
    assert captured.err == ""
