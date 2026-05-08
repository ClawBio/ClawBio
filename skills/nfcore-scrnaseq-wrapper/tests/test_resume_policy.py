from __future__ import annotations

import json
from argparse import Namespace
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import preflight
from errors import SkillError


_LOCAL_SOURCE = {"source_kind": "local_checkout", "source_ref": "scrnaseq-checkout", "resolved_version": "abc123", "branch": "main", "dirty": False}
_SS = {"sample_count": 1, "unknown_columns": []}


def test_resume_requires_manifest(tmp_path, monkeypatch):
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    gtf.write_text("chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"g1\";\n", encoding="utf-8")
    _mock_preflight_infra(monkeypatch)
    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            _make_resume_args(output_dir, fasta, gtf),
            pipeline_source=_LOCAL_SOURCE,
            samplesheet_summary=_SS,
        )
    assert exc.value.error_code == "INVALID_RESUME_STATE"


def _make_resume_args(output_dir, fasta, gtf, *, preset="star", profile="docker"):
    return Namespace(
        output=str(output_dir),
        resume=True,
        preset=preset,
        protocol="10XV3",
        profile=profile,
        demo=False,
        fasta=str(fasta),
        gtf=str(gtf),
        transcript_fasta=None,
        txp2gene=None,
        simpleaf_index=None,
        kallisto_index=None,
        star_index=None,
        cellranger_index=None,
        barcode_whitelist=None,
        genome=None,
        cellrangerarc_reference=None,
        kb_t1c=None, kb_t2c=None,
        motifs=None, cellrangerarc_config=None,
        cellranger_vdj_index=None,
        gex_frna_probe_set=None, gex_target_panel=None,
        gex_cmo_set=None, fb_reference=None,
        vdj_inner_enrichment_primers=None,
        gex_barcode_sample_assignment=None,
        cellranger_multi_barcodes=None,
    )


def _mock_preflight_infra(monkeypatch):
    monkeypatch.setattr(preflight, "_check_java", lambda: {"path": "/usr/bin/java", "version": "17.0.8"})
    monkeypatch.setattr(preflight, "_check_nextflow", lambda: {"path": "/usr/bin/nextflow", "version": "25.04.0"})
    monkeypatch.setattr(preflight, "_check_profile", lambda profile: {"profile": profile, "backend_path": "/usr/bin/docker", "backend_ready": True})


def _write_manifest(output_dir, *, preset="star", profile="docker", source_kind="local_checkout", version="abc123"):
    (output_dir / "reproducibility").mkdir(parents=True, exist_ok=True)
    (output_dir / "reproducibility" / "manifest.json").write_text(
        json.dumps({"preset": preset, "profile": profile, "pipeline_source": {"source_kind": source_kind, "resolved_version": version}}),
        encoding="utf-8",
    )


def test_resume_allows_manifest(tmp_path, monkeypatch):
    output_dir = tmp_path / "out"
    _write_manifest(output_dir)
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    gtf.write_text("chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"g1\";\n", encoding="utf-8")
    _mock_preflight_infra(monkeypatch)
    result = preflight.run_preflight(
        _make_resume_args(output_dir, fasta, gtf),
        pipeline_source=_LOCAL_SOURCE,
        samplesheet_summary=_SS,
    )
    assert result["ok"] is True


def test_resume_preset_mismatch_raises(tmp_path, monkeypatch):
    """Resume with a different preset must raise INVALID_RESUME_STATE."""
    output_dir = tmp_path / "out"
    _write_manifest(output_dir, preset="kallisto")  # previous run used kallisto
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    gtf.write_text("chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"g1\";\n", encoding="utf-8")
    _mock_preflight_infra(monkeypatch)
    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            _make_resume_args(output_dir, fasta, gtf, preset="star"),  # now requesting star
            pipeline_source=_LOCAL_SOURCE,
            samplesheet_summary=_SS,
        )
    assert exc.value.error_code == "INVALID_RESUME_STATE"


def test_resume_profile_mismatch_raises(tmp_path, monkeypatch):
    """Resume with a different profile must raise INVALID_RESUME_STATE."""
    output_dir = tmp_path / "out"
    _write_manifest(output_dir, profile="conda")  # previous run used conda
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    gtf.write_text("chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"g1\";\n", encoding="utf-8")
    _mock_preflight_infra(monkeypatch)
    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            _make_resume_args(output_dir, fasta, gtf, profile="docker"),  # now requesting docker
            pipeline_source=_LOCAL_SOURCE,
            samplesheet_summary=_SS,
        )
    assert exc.value.error_code == "INVALID_RESUME_STATE"


def test_resume_source_kind_mismatch_raises(tmp_path, monkeypatch):
    """Resume when source_kind changed (local → remote) must raise INVALID_RESUME_STATE."""
    output_dir = tmp_path / "out"
    _write_manifest(output_dir, source_kind="local_checkout", version="abc123")
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    gtf.write_text("chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"g1\";\n", encoding="utf-8")
    _mock_preflight_infra(monkeypatch)
    remote_source = {"source_kind": "remote_repo", "source_ref": "nf-core/scrnaseq", "resolved_version": "4.1.0", "branch": "", "dirty": False}
    with pytest.raises(SkillError) as exc:
        preflight.run_preflight(
            _make_resume_args(output_dir, fasta, gtf),
            pipeline_source=remote_source,  # different source_kind
            samplesheet_summary=_SS,
        )
    assert exc.value.error_code == "INVALID_RESUME_STATE"
