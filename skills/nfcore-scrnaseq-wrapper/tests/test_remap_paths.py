from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from remap_paths import cmd_remap, cmd_verify, find_samplesheet, remap_csv, verify_paths

_FASTQ_HEADER = "sample,fastq_1,fastq_2,expected_cells,seq_center\n"


def _write_samplesheet(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample", "fastq_1", "fastq_2", "expected_cells", "seq_center"])
        writer.writeheader()
        writer.writerows(rows)


# ── find_samplesheet ──────────────────────────────────────────────────────────

def test_find_samplesheet_returns_valid_csv_when_present(tmp_path):
    ss = tmp_path / "samplesheet.valid.csv"
    ss.write_text(_FASTQ_HEADER, encoding="utf-8")
    result = find_samplesheet(bundle_dir=tmp_path)
    assert result == ss


def test_find_samplesheet_returns_demo_csv_when_no_valid(tmp_path):
    ss = tmp_path / "samplesheet.demo.csv"
    ss.write_text(_FASTQ_HEADER, encoding="utf-8")
    result = find_samplesheet(bundle_dir=tmp_path)
    assert result == ss


def test_find_samplesheet_prefers_valid_over_demo(tmp_path):
    valid = tmp_path / "samplesheet.valid.csv"
    demo = tmp_path / "samplesheet.demo.csv"
    valid.write_text(_FASTQ_HEADER, encoding="utf-8")
    demo.write_text(_FASTQ_HEADER, encoding="utf-8")
    result = find_samplesheet(bundle_dir=tmp_path)
    assert result == valid


def test_find_samplesheet_returns_none_when_no_csv(tmp_path):
    assert find_samplesheet(bundle_dir=tmp_path) is None


# ── remap_csv ─────────────────────────────────────────────────────────────────

def test_remap_csv_replaces_matching_prefix(tmp_path):
    ss = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(ss, [
        {"sample": "S1", "fastq_1": "/old/data/S1_R1.fastq.gz",
         "fastq_2": "/old/data/S1_R2.fastq.gz", "expected_cells": "", "seq_center": ""},
    ])
    changes = remap_csv(ss, "/old/data", "/new/data", dry_run=False)
    assert len(changes) == 2
    assert changes[0] == ("fastq_1", "/old/data/S1_R1.fastq.gz", "/new/data/S1_R1.fastq.gz")
    assert changes[1] == ("fastq_2", "/old/data/S1_R2.fastq.gz", "/new/data/S1_R2.fastq.gz")
    rows = list(csv.DictReader(ss.read_text(encoding="utf-8").splitlines()))
    assert rows[0]["fastq_1"] == "/new/data/S1_R1.fastq.gz"
    assert rows[0]["fastq_2"] == "/new/data/S1_R2.fastq.gz"


def test_remap_csv_dry_run_does_not_modify_file(tmp_path):
    ss = tmp_path / "samplesheet.valid.csv"
    original = "/old/data/S1_R1.fastq.gz,/old/data/S1_R2.fastq.gz"
    _write_samplesheet(ss, [
        {"sample": "S1", "fastq_1": "/old/data/S1_R1.fastq.gz",
         "fastq_2": "/old/data/S1_R2.fastq.gz", "expected_cells": "", "seq_center": ""},
    ])
    original_text = ss.read_text(encoding="utf-8")
    changes = remap_csv(ss, "/old/data", "/new/data", dry_run=True)
    assert len(changes) == 2
    assert ss.read_text(encoding="utf-8") == original_text
    assert not ss.with_suffix(".bak").exists()


def test_remap_csv_creates_backup_when_modifying(tmp_path):
    ss = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(ss, [
        {"sample": "S1", "fastq_1": "/old/S1_R1.fastq.gz",
         "fastq_2": "/old/S1_R2.fastq.gz", "expected_cells": "", "seq_center": ""},
    ])
    remap_csv(ss, "/old", "/new", dry_run=False)
    assert ss.with_suffix(".bak").exists()


def test_remap_csv_returns_empty_when_no_match(tmp_path):
    ss = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(ss, [
        {"sample": "S1", "fastq_1": "/different/S1_R1.fastq.gz",
         "fastq_2": "/different/S1_R2.fastq.gz", "expected_cells": "", "seq_center": ""},
    ])
    changes = remap_csv(ss, "/old/data", "/new/data", dry_run=False)
    assert changes == []


def test_remap_csv_handles_multiple_samples(tmp_path):
    ss = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(ss, [
        {"sample": "S1", "fastq_1": "/mnt/S1_R1.fastq.gz",
         "fastq_2": "/mnt/S1_R2.fastq.gz", "expected_cells": "", "seq_center": ""},
        {"sample": "S2", "fastq_1": "/mnt/S2_R1.fastq.gz",
         "fastq_2": "/mnt/S2_R2.fastq.gz", "expected_cells": "", "seq_center": ""},
    ])
    changes = remap_csv(ss, "/mnt", "/data", dry_run=False)
    assert len(changes) == 4
    rows = list(csv.DictReader(ss.read_text(encoding="utf-8").splitlines()))
    assert rows[0]["fastq_1"] == "/data/S1_R1.fastq.gz"
    assert rows[1]["fastq_2"] == "/data/S2_R2.fastq.gz"


def test_remap_csv_only_replaces_prefix_not_middle(tmp_path):
    ss = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(ss, [
        {"sample": "S1", "fastq_1": "/data/old/S1_R1.fastq.gz",
         "fastq_2": "/data/old/S1_R2.fastq.gz", "expected_cells": "", "seq_center": ""},
    ])
    changes = remap_csv(ss, "/old", "/new", dry_run=False)
    assert changes == [], "Should not replace /old in the middle of the path"


# ── verify_paths ──────────────────────────────────────────────────────────────

def test_verify_paths_returns_empty_when_all_exist(tmp_path):
    r1 = tmp_path / "S1_R1.fastq.gz"
    r2 = tmp_path / "S1_R2.fastq.gz"
    r1.write_bytes(b"")
    r2.write_bytes(b"")
    ss = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(ss, [
        {"sample": "S1", "fastq_1": str(r1), "fastq_2": str(r2),
         "expected_cells": "", "seq_center": ""},
    ])
    assert verify_paths(ss) == []


def test_verify_paths_returns_missing_paths(tmp_path):
    ss = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(ss, [
        {"sample": "S1", "fastq_1": "/nonexistent/R1.fastq.gz",
         "fastq_2": "/nonexistent/R2.fastq.gz", "expected_cells": "", "seq_center": ""},
    ])
    missing = verify_paths(ss)
    assert len(missing) == 2
    assert "/nonexistent/R1.fastq.gz" in missing
    assert "/nonexistent/R2.fastq.gz" in missing


def test_verify_paths_ignores_empty_fastq_columns(tmp_path):
    ss = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(ss, [
        {"sample": "S1", "fastq_1": "", "fastq_2": "", "expected_cells": "", "seq_center": ""},
    ])
    assert verify_paths(ss) == []


# ── cmd_remap (integration) ───────────────────────────────────────────────────

def test_cmd_remap_succeeds_when_paths_are_remapped_to_existing_files(tmp_path):
    r1 = tmp_path / "fastqs" / "S1_R1.fastq.gz"
    r2 = tmp_path / "fastqs" / "S1_R2.fastq.gz"
    r1.parent.mkdir()
    r1.write_bytes(b"")
    r2.write_bytes(b"")
    ss = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(ss, [
        {"sample": "S1", "fastq_1": "/old/fastqs/S1_R1.fastq.gz",
         "fastq_2": "/old/fastqs/S1_R2.fastq.gz", "expected_cells": "", "seq_center": ""},
    ])
    rc = cmd_remap("/old/fastqs", str(tmp_path / "fastqs"), dry_run=False, bundle_dir=tmp_path)
    assert rc == 0
    rows = list(csv.DictReader(ss.read_text(encoding="utf-8").splitlines()))
    assert rows[0]["fastq_1"] == str(tmp_path / "fastqs" / "S1_R1.fastq.gz")


def test_cmd_remap_returns_nonzero_when_new_paths_missing(tmp_path):
    ss = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(ss, [
        {"sample": "S1", "fastq_1": "/old/S1_R1.fastq.gz",
         "fastq_2": "/old/S1_R2.fastq.gz", "expected_cells": "", "seq_center": ""},
    ])
    rc = cmd_remap("/old", "/nonexistent", dry_run=False, bundle_dir=tmp_path)
    assert rc != 0


def test_cmd_remap_dry_run_returns_zero_without_verifying(tmp_path):
    ss = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(ss, [
        {"sample": "S1", "fastq_1": "/old/S1_R1.fastq.gz",
         "fastq_2": "/old/S1_R2.fastq.gz", "expected_cells": "", "seq_center": ""},
    ])
    rc = cmd_remap("/old", "/nonexistent", dry_run=True, bundle_dir=tmp_path)
    assert rc == 0


def test_cmd_remap_returns_zero_when_no_paths_matched(tmp_path):
    ss = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(ss, [
        {"sample": "S1", "fastq_1": "/other/S1_R1.fastq.gz",
         "fastq_2": "/other/S1_R2.fastq.gz", "expected_cells": "", "seq_center": ""},
    ])
    rc = cmd_remap("/old", "/new", dry_run=False, bundle_dir=tmp_path)
    assert rc == 0


def test_cmd_remap_returns_nonzero_when_no_samplesheet(tmp_path):
    rc = cmd_remap("/old", "/new", dry_run=False, bundle_dir=tmp_path)
    assert rc != 0


# ── cmd_verify (integration) ──────────────────────────────────────────────────

def test_cmd_verify_returns_zero_when_all_paths_exist(tmp_path):
    r1 = tmp_path / "S1_R1.fastq.gz"
    r2 = tmp_path / "S1_R2.fastq.gz"
    r1.write_bytes(b"")
    r2.write_bytes(b"")
    ss = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(ss, [
        {"sample": "S1", "fastq_1": str(r1), "fastq_2": str(r2),
         "expected_cells": "", "seq_center": ""},
    ])
    rc = cmd_verify(bundle_dir=tmp_path)
    assert rc == 0


def test_cmd_verify_returns_nonzero_when_paths_missing(tmp_path):
    ss = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(ss, [
        {"sample": "S1", "fastq_1": "/nonexistent/R1.fastq.gz",
         "fastq_2": "/nonexistent/R2.fastq.gz", "expected_cells": "", "seq_center": ""},
    ])
    rc = cmd_verify(bundle_dir=tmp_path)
    assert rc != 0


def test_cmd_verify_returns_nonzero_when_no_samplesheet(tmp_path):
    rc = cmd_verify(bundle_dir=tmp_path)
    assert rc != 0


# ── find_commands_sh ──────────────────────────────────────────────────────────

def test_find_commands_sh_returns_path_when_exists(tmp_path):
    from remap_paths import find_commands_sh
    cs = tmp_path / "commands.sh"
    cs.write_text("#!/usr/bin/env bash\n", encoding="utf-8")
    assert find_commands_sh(bundle_dir=tmp_path) == cs


def test_find_commands_sh_returns_none_when_absent(tmp_path):
    from remap_paths import find_commands_sh
    assert find_commands_sh(bundle_dir=tmp_path) is None


# ── update_commands_output ────────────────────────────────────────────────────

def test_update_commands_output_replaces_output_flag(tmp_path):
    from remap_paths import update_commands_output
    cs = tmp_path / "commands.sh"
    cs.write_text(
        'python "$SKILL_SCRIPT" \\\n    --output /old/output/dir \\\n    --preset star\n',
        encoding="utf-8",
    )
    update_commands_output(cs, "/new/output/dir")
    content = cs.read_text(encoding="utf-8")
    assert "/new/output/dir" in content
    assert "/old/output/dir" not in content


def test_update_commands_output_creates_backup(tmp_path):
    from remap_paths import update_commands_output
    cs = tmp_path / "commands.sh"
    cs.write_text('    --output /old/path \\\n', encoding="utf-8")
    update_commands_output(cs, "/new/path")
    assert cs.with_suffix(".sh.bak").exists()


def test_update_commands_output_noop_when_no_output_flag(tmp_path):
    from remap_paths import update_commands_output
    cs = tmp_path / "commands.sh"
    original = '#!/usr/bin/env bash\necho hello\n'
    cs.write_text(original, encoding="utf-8")
    update_commands_output(cs, "/new/path")
    assert cs.read_text(encoding="utf-8") == original
    assert not cs.with_suffix(".sh.bak").exists()


# ── CLAWBIO_REPO fallback in commands.sh ─────────────────────────────────────

def test_generated_commands_sh_contains_clawbio_repo_fallback(tmp_path):
    import argparse
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from reporting import write_repro_commands
    args = argparse.Namespace(
        input=str(tmp_path / "samplesheet.csv"),
        output=str(tmp_path / "out"), preset="star", profile="docker",
        pipeline_version="4.1.0", protocol=None, demo=False, check=False,
        resume=False, skip_cellbender=False, skip_fastqc=False,
        skip_emptydrops=False, skip_multiqc=False, star_feature=None,
        star_ignore_sjdbgtf=False, seq_center=None, simpleaf_umi_resolution=None,
        kb_workflow=None, kb_t1c=None, kb_t2c=None, fasta=None, gtf=None,
        transcript_fasta=None, txp2gene=None, simpleaf_index=None,
        kallisto_index=None, star_index=None, cellranger_index=None,
        barcode_whitelist=None, expected_cells=None, genome=None,
        save_reference=False, save_align_intermeds=False,
        skip_cellranger_renaming=False, motifs=None, cellrangerarc_config=None,
        cellrangerarc_reference=None, cellranger_vdj_index=None,
        skip_cellrangermulti_vdjref=False, gex_frna_probe_set=None,
        gex_target_panel=None, gex_cmo_set=None, fb_reference=None,
        vdj_inner_enrichment_primers=None, gex_barcode_sample_assignment=None,
        cellranger_multi_barcodes=None, run_downstream=False,
        email=None, multiqc_title=None,
    )
    write_repro_commands(tmp_path, args=args)
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "CLAWBIO_REPO" in content
