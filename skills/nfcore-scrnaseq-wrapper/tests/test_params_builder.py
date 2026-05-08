from __future__ import annotations

from argparse import Namespace
from pathlib import Path
import sys

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from params_builder import build_params_file


def _base_args(tmp_path: Path, **overrides) -> Namespace:
    defaults = dict(
        demo=False,
        preset="star",
        protocol=None,
        email=None,
        multiqc_title=None,
        skip_cellbender=False,
        skip_fastqc=False,
        skip_emptydrops=False,
        skip_multiqc=False,
        fasta=None,
        gtf=None,
        transcript_fasta=None,
        txp2gene=None,
        simpleaf_index=None,
        kallisto_index=None,
        star_index=None,
        cellranger_index=None,
        barcode_whitelist=None,
        star_feature=None,
        star_ignore_sjdbgtf=False,
        seq_center=None,
        simpleaf_umi_resolution=None,
        kb_workflow=None,
        kb_t1c=None,
        kb_t2c=None,
        skip_cellranger_renaming=False,
        motifs=None,
        cellrangerarc_config=None,
        cellrangerarc_reference=None,
        cellranger_vdj_index=None,
        skip_cellrangermulti_vdjref=False,
        gex_frna_probe_set=None,
        gex_target_panel=None,
        gex_cmo_set=None,
        fb_reference=None,
        vdj_inner_enrichment_primers=None,
        gex_barcode_sample_assignment=None,
        cellranger_multi_barcodes=None,
        genome=None,
        save_reference=False,
        save_align_intermeds=False,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


def test_build_params_file_standard(tmp_path):
    samplesheet = tmp_path / "reproducibility" / "samplesheet.valid.csv"
    samplesheet.parent.mkdir()
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, preset="standard", protocol="10XV2", fasta=str(tmp_path / "g.fa"), gtf=str(tmp_path / "g.gtf"))
    path, payload = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert payload["aligner"] == "simpleaf"
    # Keep --input schema-safe even when the absolute output path contains whitespace.
    assert loaded["input"] == "reproducibility/samplesheet.valid.csv"
    assert loaded["protocol"] == "10XV2"


def test_input_path_is_relative_when_output_dir_contains_spaces(tmp_path):
    output_dir = tmp_path / "output with spaces"
    samplesheet = output_dir / "reproducibility" / "samplesheet.valid.csv"
    samplesheet.parent.mkdir(parents=True)
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, preset="star")

    path, _payload = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=output_dir)

    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["input"] == "reproducibility/samplesheet.valid.csv"
    assert " " not in loaded["input"]


def test_build_params_file_demo_omits_input(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, demo=True, preset="star")
    path, _payload = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert "input" not in loaded


def test_build_params_file_produces_valid_yaml(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    kb_index = tmp_path / "kb_index"
    kb_index.mkdir()
    whitelist = tmp_path / "whitelist.txt"
    whitelist.write_text("ACGT\n", encoding="utf-8")
    args = _base_args(
        tmp_path,
        preset="kallisto",
        protocol="10XV3",
        skip_cellbender=True,
        kallisto_index=str(kb_index),
        barcode_whitelist=str(whitelist),
    )
    path, payload = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    assert path.suffix == ".yaml"
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["aligner"] == "kallisto"
    assert loaded["skip_cellbender"] is True
    # reference paths must be absolute POSIX paths (resolved, forward slashes)
    assert loaded["kallisto_index"] == kb_index.resolve().as_posix()
    assert loaded["barcode_whitelist"] == whitelist.resolve().as_posix()


def test_params_paths_use_forward_slashes(tmp_path):
    """No backslashes in any path written to params.yaml (cross-platform safety)."""
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    fasta.write_text(">c\nACGT\n", encoding="utf-8")
    gtf.write_text("", encoding="utf-8")
    args = _base_args(tmp_path, preset="star", fasta=str(fasta), gtf=str(gtf))
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    text = path.read_text(encoding="utf-8")
    assert "\\" not in text, f"Backslash found in params.yaml:\n{text}"


def test_reference_paths_are_resolved_to_absolute(tmp_path):
    """Relative reference paths must be resolved to absolute before writing to params."""
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">c\nACGT\n", encoding="utf-8")
    args = _base_args(tmp_path, preset="star", fasta=str(fasta))
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert Path(loaded["fasta"]).is_absolute()


def test_outdir_uses_forward_slashes(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, demo=True, preset="star")
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert "\\" not in loaded["outdir"]
    assert "/" in loaded["outdir"]


# ── skip flags ────────────────────────────────────────────────────────────────

def test_skip_fastqc_appears_in_params(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, skip_fastqc=True)
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["skip_fastqc"] is True


def test_skip_emptydrops_maps_to_canonical_skip_cellbender_param(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, skip_emptydrops=True)
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["skip_cellbender"] is True
    assert "skip_emptydrops" not in loaded


def test_skip_multiqc_appears_in_params(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, skip_multiqc=True)
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["skip_multiqc"] is True


def test_skip_flags_absent_when_false(tmp_path):
    """Skip flags must not pollute params.yaml when not requested."""
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, skip_fastqc=False, skip_emptydrops=False, skip_multiqc=False)
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert "skip_fastqc" not in loaded
    assert "skip_emptydrops" not in loaded
    assert "skip_multiqc" not in loaded


# ── aligner tuning ────────────────────────────────────────────────────────────

def test_star_feature_appears_in_params(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, preset="star", star_feature="GeneFull")
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["star_feature"] == "GeneFull"


def test_simpleaf_umi_resolution_appears_in_params(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, preset="standard", simpleaf_umi_resolution="cr-like-em")
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["simpleaf_umi_resolution"] == "cr-like-em"


def test_kb_workflow_appears_in_params(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, preset="kallisto", kb_workflow="lamanno")
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["kb_workflow"] == "lamanno"


def test_aligner_tuning_absent_when_none(tmp_path):
    """Aligner tuning params must not appear in params.yaml when not set."""
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, star_feature=None, simpleaf_umi_resolution=None, kb_workflow=None)
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert "star_feature" not in loaded
    assert "simpleaf_umi_resolution" not in loaded
    assert "kb_workflow" not in loaded


# ── reference management ──────────────────────────────────────────────────────

def test_genome_shortcut_appears_in_params(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, genome="GRCh38")
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["genome"] == "GRCh38"


def test_save_reference_appears_in_params(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, save_reference=True)
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["save_reference"] is True


def test_save_align_intermeds_appears_in_params(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, save_align_intermeds=True)
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["save_align_intermeds"] is True


def test_save_flags_absent_when_false(tmp_path):
    """Neither save_reference nor save_align_intermeds must appear in params.yaml when not
    explicitly requested. Omitting them lets the pipeline's nextflow.config defaults apply
    (save_align_intermeds defaults to true in the upstream pipeline)."""
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, save_reference=False, save_align_intermeds=False)
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert "save_reference" not in loaded
    assert "save_align_intermeds" not in loaded


# ── input/output metadata ─────────────────────────────────────────────────────

def test_email_appears_in_params(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, email="lab@example.com")
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["email"] == "lab@example.com"


def test_multiqc_title_appears_in_params(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, multiqc_title="My Experiment QC")
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["multiqc_title"] == "My Experiment QC"


# ── STARsolo additional options ───────────────────────────────────────────────

def test_star_ignore_sjdbgtf_written_as_string(tmp_path):
    """star_ignore_sjdbgtf schema type is 'string' — must be written as 'true', not boolean True,
    to pass nf-core JSON schema validation."""
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, preset="star", star_ignore_sjdbgtf=True)
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["star_ignore_sjdbgtf"] == "true"
    assert not isinstance(loaded["star_ignore_sjdbgtf"], bool)


def test_seq_center_appears_in_params(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, seq_center="GenomicsCore")
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["seq_center"] == "GenomicsCore"


# ── Kallisto RNA velocity ─────────────────────────────────────────────────────

def test_kb_t1c_appears_in_params(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    t1c = tmp_path / "cdna_t2c.txt"
    t1c.write_text("transcript\tgene\n", encoding="utf-8")
    args = _base_args(tmp_path, preset="kallisto", kb_t1c=str(t1c))
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["kb_t1c"] == t1c.resolve().as_posix()


def test_kb_t2c_appears_in_params(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    t2c = tmp_path / "intron_t2c.txt"
    t2c.write_text("transcript\tgene\n", encoding="utf-8")
    args = _base_args(tmp_path, preset="kallisto", kb_t2c=str(t2c))
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["kb_t2c"] == t2c.resolve().as_posix()


# ── CellRanger ────────────────────────────────────────────────────────────────

def test_skip_cellranger_renaming_appears_in_params(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, preset="cellranger", skip_cellranger_renaming=True)
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["skip_cellranger_renaming"] is True


# ── CellRanger ARC ────────────────────────────────────────────────────────────

def test_cellrangerarc_params_appear_in_params(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    motifs = tmp_path / "jaspar.motifs"
    motifs.write_text("", encoding="utf-8")
    arc_config = tmp_path / "arc.json"
    arc_config.write_text("{}", encoding="utf-8")
    args = _base_args(
        tmp_path,
        preset="cellrangerarc",
        motifs=str(motifs),
        cellrangerarc_config=str(arc_config),
        cellrangerarc_reference="GRCh38-2024-A",
    )
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["motifs"] == motifs.resolve().as_posix()
    assert loaded["cellrangerarc_config"] == arc_config.resolve().as_posix()
    assert loaded["cellrangerarc_reference"] == "GRCh38-2024-A"


# ── CellRanger Multi ──────────────────────────────────────────────────────────

def test_cellranger_multi_path_params_appear_in_params(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    vdj_index = tmp_path / "vdj_idx"
    vdj_index.mkdir()
    barcodes = tmp_path / "barcodes.csv"
    barcodes.write_text("sample,barcode\n", encoding="utf-8")
    cmo_set = tmp_path / "cmo.csv"
    cmo_set.write_text("id,sequence\n", encoding="utf-8")
    args = _base_args(
        tmp_path,
        preset="cellrangermulti",
        cellranger_vdj_index=str(vdj_index),
        cellranger_multi_barcodes=str(barcodes),
        gex_cmo_set=str(cmo_set),
    )
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["cellranger_vdj_index"] == vdj_index.resolve().as_posix()
    assert loaded["cellranger_multi_barcodes"] == barcodes.resolve().as_posix()
    assert loaded["gex_cmo_set"] == cmo_set.resolve().as_posix()


def test_skip_cellrangermulti_vdjref_appears_in_params(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    args = _base_args(tmp_path, preset="cellrangermulti", skip_cellrangermulti_vdjref=True)
    path, _ = build_params_file(args, normalized_samplesheet=samplesheet, output_dir=tmp_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert loaded["skip_cellrangermulti_vdjref"] is True
