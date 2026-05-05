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
