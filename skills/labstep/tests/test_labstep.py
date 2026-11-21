"""
test_labstep.py — Test suite for the Labstep ELN bridge skill
==============================================================
Run with: pytest skills/labstep/tests/test_labstep.py -v

Uses pre-loaded demo JSON fixtures and mocked labstepPy objects — no
network or API key required.  All assertions are deterministic.
"""

from __future__ import annotations

import json
import sys
from io import StringIO
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from inline_snapshot import snapshot
from clawbio.common import demo_runner

sys.path.insert(0, str(Path(__file__).parent.parent))

from labstep import (
    DEMO_DIR,
    DISCLAIMER,
    _fmt_date,
    _load_demo,
    _write_output,
    format_experiments,
    format_inventory,
    format_protocols,
    run_demo,
)

# ---------------------------------------------------------------------------
# Fixtures — load demo data once
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def experiments():
    return json.loads((DEMO_DIR / "demo_experiments.json").read_text())


@pytest.fixture(scope="module")
def protocols():
    return json.loads((DEMO_DIR / "demo_protocols.json").read_text())


@pytest.fixture(scope="module")
def inventory():
    return json.loads((DEMO_DIR / "demo_inventory.json").read_text())


# ---------------------------------------------------------------------------
# Demo data loading
# ---------------------------------------------------------------------------


def test_demo_experiments_loads():
    """demo_experiments.json loads and has expected structure."""
    data = _load_demo("demo_experiments.json")
    assert isinstance(data, list)
    assert len(data) == 3


def test_demo_protocols_loads():
    """demo_protocols.json loads and has expected structure."""
    data = _load_demo("demo_protocols.json")
    assert isinstance(data, list)
    assert len(data) == 3


def test_demo_inventory_loads():
    """demo_inventory.json loads and has expected structure."""
    data = _load_demo("demo_inventory.json")
    assert isinstance(data, dict)
    assert "resources" in data
    assert "locations" in data


def test_demo_missing_file_exits():
    """_load_demo exits with SystemExit on a missing file."""
    with pytest.raises(SystemExit):
        _load_demo("nonexistent_file.json")


# ---------------------------------------------------------------------------
# _fmt_date
# ---------------------------------------------------------------------------


def test_fmt_date_iso_z():
    assert _fmt_date("2026-01-08T09:14:22Z") == "2026-01-08"


def test_fmt_date_iso_offset():
    assert _fmt_date("2025-11-03T14:07:55+00:00") == "2025-11-03"


def test_fmt_date_short_fallback():
    """Falls back to first 10 chars for bare date strings."""
    assert _fmt_date("2025-09-15") == "2025-09-15"


def test_fmt_date_empty():
    """Empty string returns em-dash."""
    assert _fmt_date("") == "—"


# ---------------------------------------------------------------------------
# format_experiments
# ---------------------------------------------------------------------------


def test_format_experiments_title(experiments):
    md = format_experiments(experiments, title="Test Title")
    assert "# 🔬 Labstep — Test Title" in md


def test_format_experiments_count(experiments):
    md = format_experiments(experiments)
    assert "**3 experiment(s)**" in md


def test_format_experiments_ids(experiments):
    md = format_experiments(experiments)
    assert "## [10241]" in md
    assert "## [10187]" in md
    assert "## [10099]" in md


def test_format_experiments_data_fields_table(experiments):
    md = format_experiments(experiments)
    assert "| Cell Line | MCF7 |" in md
    assert "| MOI | 0.3 |" in md


def test_format_experiments_tags(experiments):
    md = format_experiments(experiments)
    assert "`CRISPR`" in md
    assert "`single-cell`" in md


def test_format_experiments_linked_protocols(experiments):
    md = format_experiments(experiments)
    assert "Lentiviral sgRNA Library Transduction (#3301)" in md


def test_format_experiments_comments(experiments):
    md = format_experiments(experiments)
    assert "D14 replicate 2 had low reads" in md


def test_format_experiments_disclaimer(experiments):
    md = format_experiments(experiments)
    assert DISCLAIMER in md


def test_format_experiments_empty():
    md = format_experiments([])
    assert "**0 experiment(s)**" in md


# ---------------------------------------------------------------------------
# format_protocols
# ---------------------------------------------------------------------------


def test_format_protocols_title(protocols):
    md = format_protocols(protocols, title="My Protocols")
    assert "# 📋 Labstep — My Protocols" in md


def test_format_protocols_count(protocols):
    md = format_protocols(protocols)
    assert "**3 protocol(s)**" in md


def test_format_protocols_version(protocols):
    md = format_protocols(protocols)
    assert "(v3)" in md


def test_format_protocols_steps(protocols):
    md = format_protocols(protocols)
    assert "**Steps (5):**" in md
    assert "**1. Prepare cells**" in md
    assert "**5. Harvest Day 0 sample**" in md


def test_format_protocols_step_body(protocols):
    md = format_protocols(protocols)
    assert "MOI of 0.3" in md
    assert "spinoculation" in md


def test_format_protocols_inventory_fields(protocols):
    md = format_protocols(protocols)
    assert "sgRNA Library (lentiviral) (#5501)" in md
    assert "Puromycin (#5503)" in md


def test_format_protocols_disclaimer(protocols):
    md = format_protocols(protocols)
    assert DISCLAIMER in md


# ---------------------------------------------------------------------------
# format_inventory
# ---------------------------------------------------------------------------


def test_format_inventory_title(inventory):
    md = format_inventory(inventory)
    assert "# 🧪 Labstep — Inventory" in md


def test_format_inventory_resource_count(inventory):
    md = format_inventory(inventory)
    assert "**10 resource(s)**" in md


def test_format_inventory_categories(inventory):
    md = format_inventory(inventory)
    assert "## Viral Library" in md
    assert "## RNA Reagent" in md
    assert "## Solvent" in md


def test_format_inventory_supplier(inventory):
    md = format_inventory(inventory)
    assert "Addgene #73178" in md


def test_format_inventory_hazard(inventory):
    md = format_inventory(inventory)
    assert "GHS06, GHS08" in md


def test_format_inventory_item_locations(inventory):
    md = format_inventory(inventory)
    assert "📍 −80°C Rack B / Box 3" in md
    assert "📍 4°C Fridge / Shelf 3" in md


def test_format_inventory_storage_table(inventory):
    md = format_inventory(inventory)
    assert "## Storage Locations" in md
    assert "| −80°C Rack B / Box 3 |" in md


def test_format_inventory_search_filter(inventory):
    """Search filters resources to only matching names/categories."""
    md = format_inventory(inventory, search="RNA")
    assert "TRIzol Reagent" in md
    assert "RNase-free water" in md
    # Non-RNA items should not appear
    assert "Polybrene" not in md
    assert "Puromycin" not in md


def test_format_inventory_search_title(inventory):
    md = format_inventory(inventory, search="TRIzol")
    assert '"TRIzol"' in md


def test_format_inventory_search_no_storage_table(inventory):
    """Storage location table is omitted when a search filter is active."""
    md = format_inventory(inventory, search="RNA")
    assert "## Storage Locations" not in md


def test_format_inventory_disclaimer(inventory):
    md = format_inventory(inventory)
    assert DISCLAIMER in md


# ---------------------------------------------------------------------------
# run_demo — smoke test (captures stdout)
# ---------------------------------------------------------------------------


def test_run_demo_runs_without_error(capsys):
    """run_demo() completes without raising and prints experiment section."""
    run_demo()
    out = capsys.readouterr().out
    assert "EXPERIMENTS" in out
    assert "PROTOCOL DETAIL" in out
    assert "INVENTORY SNAPSHOT" in out
    assert "INVENTORY SEARCH" in out


def test_run_demo_output_contains_experiments(capsys):
    run_demo()
    out = capsys.readouterr().out
    assert "CRISPR Screen" in out
    assert "scTIP-seq Timecourse" in out
    assert "RNA Extraction Validation" in out


def test_run_demo_output_contains_protocol_steps(capsys):
    run_demo()
    out = capsys.readouterr().out
    assert "spinoculation" in out
    assert "Puromycin selection" in out


def test_run_demo_output_contains_inventory(capsys):
    run_demo()
    out = capsys.readouterr().out
    assert "Brunello sgRNA Library" in out
    assert "TRIzol Reagent" in out


def test_run_demo_disclaimer_present(capsys):
    run_demo()
    out = capsys.readouterr().out
    assert "research and educational tool" in out


# ---------------------------------------------------------------------------
# _write_output — report.md written to output directory
# ---------------------------------------------------------------------------


def test_write_output_creates_report_md(tmp_path):
    """_write_output writes report.md to the given directory."""
    _write_output(tmp_path, "# Test Report\n\nContent here.", label="test", repro_args=["--demo"])
    report = tmp_path / "report.md"
    assert report.exists()
    assert "# Test Report" in report.read_text()


def test_write_output_creates_directory(tmp_path):
    """_write_output creates the output directory if it doesn't exist."""
    out = tmp_path / "new_subdir"
    assert not out.exists()
    _write_output(out, "hello", label="test", repro_args=["--demo"])
    assert out.exists()
    assert (out / "report.md").exists()


def test_run_demo_writes_report_md(tmp_path):
    """run_demo with output_dir writes a combined report.md."""
    run_demo(output_dir=tmp_path)
    report = tmp_path / "report.md"
    assert report.exists()
    content = report.read_text()
    assert "Experiments" in content
    assert "Protocol Detail" in content
    assert "Inventory" in content
    assert "CRISPR Screen" in content
    assert "TRIzol Reagent" in content


def test_run_demo_report_contains_disclaimer(tmp_path):
    run_demo(output_dir=tmp_path)
    content = (tmp_path / "report.md").read_text()
    assert "research and educational tool" in content


# ---------------------------------------------------------------------------
# get_labstep_user — auth error paths
# ---------------------------------------------------------------------------


def test_get_labstep_user_no_key_exits(monkeypatch, tmp_path):
    """Missing API key causes SystemExit."""
    monkeypatch.delenv("LABSTEP_API_KEY", raising=False)
    monkeypatch.chdir(tmp_path)
    import labstep as _ls_module  # noqa: F401 — ensure importable

    with pytest.raises(SystemExit):
        from labstep import get_labstep_user
        with patch.dict("os.environ", {}, clear=True):
            # Patch away .claude/settings.json
            with patch("labstep.Path.exists", return_value=False):
                get_labstep_user()


# ---------------------------------------------------------------------------
# Reproducibility bundle — written via clawbio.common.reproducibility
# ---------------------------------------------------------------------------


def test_run_demo_writes_complete_reproducibility_bundle(tmp_path):
    """The demo bundle carries a portable commands.sh, an env, and checksums."""
    run_demo(output_dir=tmp_path)
    repro = tmp_path / "reproducibility"

    commands = (repro / "commands.sh").read_text(encoding="utf-8")
    assert commands.startswith("#!/usr/bin/env bash")
    assert "CLAWBIO_ROOT" in commands
    assert "$OUTPUT_DIR" in commands
    assert "skills/labstep/labstep.py" in commands
    assert "--demo" in commands
    assert str(tmp_path) not in commands

    environment = (repro / "environment.yml").read_text(encoding="utf-8")
    assert "name: clawbio-labstep" in environment
    assert "labstep" in environment

    checksum_lines = [
        line for line in (repro / "checksums.sha256").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    assert checksum_lines
    labels = set()
    for line in checksum_lines:
        digest, label = line.split("  ", 1)
        assert len(digest) == 64
        assert (tmp_path / label).exists()
        labels.add(label)
    assert {"report.md", "result.json"} <= labels


def test_write_output_records_the_live_mode_flags(tmp_path):
    """A live run replays with its own flags, not the demo ones."""
    _write_output(tmp_path, "# Report", label="experiments", repro_args=["--experiments", "--count", "5"])
    commands = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "--experiments" in commands
    assert "--count" in commands
    assert "\n  5 " in commands
    assert "--demo" not in commands


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/labstep/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("labstep", tmp_path) == snapshot({'report.md': """\
# 🔬 Labstep ELN — Demo Report

## Experiments

# 🔬 Labstep — Recent Experiments (demo)

**3 experiment(s)**

## [10241] CRISPR Screen — BRCA1 sgRNA Library (Batch 3)

- **Created**: 2026-01-08
- **Updated**: 2026-02-17
- **Tags**: `CRISPR`, `functional-screen`, `breast-cancer`

**Data Fields:**

| Field | Value |
|---|---|
| Cell Line | MCF7 |
| Library Size | 18423 |
| MOI | 0.3 |
| Timepoint | D0, D7, D14, D21 |
| Sequencing Platform | Illumina NovaSeq 6000 |
| Analysis Status | Complete — MAGeCK v0.5.9 |

**Linked Protocols**: Lentiviral sgRNA Library Transduction (#3301), Puromycin Selection (48h) (#3302), Genomic DNA Extraction + Amplicon-seq (#3303)

**Comments:**

- *2026-02-10* — D14 replicate 2 had low reads — excluded from final analysis.
- *2026-02-17* — MAGeCK output and volcano plots saved to /data/crispr_screen/batch3/

## [10187] scTIP-seq Timecourse — Doxycycline-inducible MYC

- **Created**: 2025-11-03
- **Updated**: 2026-01-29
- **Tags**: `single-cell`, `chromatin`, `MYC`, `timecourse`

**Data Fields:**

| Field | Value |
|---|---|
| Cell Line | U2OS-iMYC |
| Timepoints | 0h, 6h, 24h, 72h post-dox |
| Cells per Timepoint | 50000 |
| Antibody | anti-H3K27ac (Abcam ab4729) |
| Sequencing Depth | 400000000 |
| Pipeline | Seurat 5.1 / Signac 1.12 |

**Linked Protocols**: scTIP-seq Cell Fixation & Permeabilisation (#3310), Antibody Incubation & Tagmentation (#3311), 10x Chromium Library Prep v3.1 (#3312)

**Comments:**

- *2026-01-22* — 72h timepoint shows clear MYC-driven enhancer activation at TERT locus.

## [10099] RNA Extraction Validation — TRIzol vs RNeasy

- **Created**: 2025-09-15
- **Updated**: 2025-10-02
- **Tags**: `QC`, `RNA`, `method-comparison`

**Data Fields:**

| Field | Value |
|---|---|
| Sample Type | Mouse liver, n=6 |
| RIN (TRIzol) | 9.1 |
| RIN (RNeasy) | 9.4 |
| 260/280 (TRIzol) | 1.92 |
| 260/280 (RNeasy) | 2.05 |
| Yield preference | RNeasy preferred — higher purity |

**Linked Protocols**: TRIzol RNA Extraction (#3280), Qiagen RNeasy Plus Mini Kit (#3281), Agilent Bioanalyzer RNA 6000 Nano QC (#3282)

**Comments:**

- *2025-10-02* — SOP updated: RNeasy is now standard for all downstream sequencing.

---
*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*

## Protocol Detail

# 📋 Labstep — Protocol Detail (demo)

**1 protocol(s)**

## [3301] Lentiviral sgRNA Library Transduction  (v3)

- **Created**: 2025-06-12
- **Updated**: 2026-01-05

**Steps (5):**

**1. Prepare cells**

Seed 2×10⁵ target cells per well in a 6-well plate 24 h before transduction. Use complete DMEM + 10% FBS.

**2. Thaw viral library**

Thaw sgRNA lentiviral library aliquot on ice. Use a MOI of 0.3 to ensure single-integration events.

**3. Transduction**

Add library to cells in fresh medium containing 8 μg/mL polybrene. Centrifuge at 1000 × g, 32°C, 2 h (spinoculation). Replace medium after 6 h.

**4. Puromycin selection**

48 h post-transduction, begin puromycin selection at empirically determined kill curve concentration (typically 2 μg/mL). Maintain for 5 days.

**5. Harvest Day 0 sample**

Pellet ≥500× library-representation cells. Freeze at −80°C for gDNA extraction. This is the T0 reference.

**Inventory Fields**: sgRNA Library (lentiviral) (#5501), Polybrene (#5502), Puromycin (#5503)

---
*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*

## Inventory Snapshot

# 🧪 Labstep — Inventory

**10 resource(s)** across 8 location(s)

## Antibiotic

### [5503] Puromycin dihydrochloride

- **Supplier/Cat#**: Gibco A11138-03
- **Lot**: 2425071
- **Expiry**: 2026-09-15
- **Stock items**: 2
  - Item #9020: 10 mg/mL stock (PBS) | 2 mL | 📍 −20°C Freezer / Box 7
  - Item #9021: 10 mg/mL stock (PBS) | 2 mL | 📍 −20°C Freezer / Box 7

## Beads

### [5512] ConA-coated magnetic beads

- **Supplier/Cat#**: Bangs Labs CP01N
- **Lot**: L0012388
- **Expiry**: 2026-08-01
- **Stock items**: 2
  - Item #9050: Vial 1 | 1 mL | 📍 4°C Fridge / Shelf 3
  - Item #9051: Vial 2 (backup) | 1 mL | 📍 4°C Fridge / Shelf 3

## Buffer

### [5523] RNase-free water

- **Supplier/Cat#**: Thermo AM9932
- **Lot**: 2432101
- **Expiry**: 2028-01-01
- **Stock items**: 2
  - Item #9090: 100 mL bottle | 100 mL | 📍 RT Bench Stock
  - Item #9091: 100 mL bottle | 100 mL | 📍 RT Bench Stock

## Detergent

### [5511] Digitonin (Sigma D141)

- **Supplier/Cat#**: Sigma D141
- **Lot**: SLCK5491
- **Expiry**: 2027-01-01
- **Stock items**: 1
  - Item #9040: 5% stock in DMSO | 200 μL | 📍 −20°C Freezer / Box 12

## Fixative

### [5510] Formaldehyde 37% (RNase-free)

- **Supplier/Cat#**: Thermo 28908
- **Lot**: WJ340291
- **Expiry**: 2026-11-01
- **Hazard**: GHS06, GHS08
- **Stock items**: 1
  - Item #9030: 500 mL bottle | 500 mL | 📍 Fume Hood Cabinet / Shelf 1

## RNA Reagent

### [5520] TRIzol Reagent

- **Supplier/Cat#**: Thermo 15596026
- **Lot**: 2527301
- **Expiry**: 2027-03-01
- **Hazard**: GHS05, GHS06
- **Stock items**: 1
  - Item #9060: 200 mL bottle (in use) | ~120 mL remaining | 📍 Fume Hood Cabinet / Shelf 1

## Solvent

### [5521] Chloroform

- **Supplier/Cat#**: Fisher C606SK
- **Lot**: 199483
- **Expiry**: 2028-01-01
- **Hazard**: GHS06, GHS08
- **Stock items**: 1
  - Item #9070: 1 L bottle | ~800 mL | 📍 Flammable Cabinet / Shelf 2

### [5522] Isopropanol (molecular biology grade)

- **Supplier/Cat#**: Sigma 563935
- **Lot**: BCBV1248
- **Expiry**: 2028-06-01
- **Stock items**: 1
  - Item #9080: 1 L bottle | ~600 mL | 📍 Flammable Cabinet / Shelf 2

## Transfection Reagent

### [5502] Polybrene (hexadimethrine bromide)

- **Supplier/Cat#**: Sigma-Aldrich H9268
- **Lot**: MKBT2912V
- **Expiry**: 2026-12-01
- **Stock items**: 1
  - Item #9010: 10 mg/mL stock | 1 mL | 📍 4°C Fridge / Shelf 2

## Viral Library

### [5501] Brunello sgRNA Library (lentiviral, human)

- **Supplier/Cat#**: Addgene #73178
- **Lot**: LOT-2025-0042
- **Expiry**: 2027-06-01
- **Stock items**: 2
  - Item #9001: Aliquot A — 1e9 TU/mL | 500 μL | 📍 −80°C Rack B / Box 3
  - Item #9002: Aliquot B — 1e9 TU/mL | 500 μL | 📍 −80°C Rack B / Box 3

## Storage Locations

| Location | Items |
|---|---|
| −80°C Rack B / Box 3 | 4 |
| −20°C Freezer / Box 7 | 3 |
| −20°C Freezer / Box 12 | 2 |
| 4°C Fridge / Shelf 2 | 1 |
| 4°C Fridge / Shelf 3 | 2 |
| Fume Hood Cabinet / Shelf 1 | 2 |
| Flammable Cabinet / Shelf 2 | 2 |
| RT Bench Stock | 2 |

---
*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*

## Inventory Search: "RNA"

# 🧪 Labstep — Inventory — "RNA"

**4 resource(s)** across 8 location(s)

## Buffer

### [5523] RNase-free water

- **Supplier/Cat#**: Thermo AM9932
- **Lot**: 2432101
- **Expiry**: 2028-01-01
- **Stock items**: 2
  - Item #9090: 100 mL bottle | 100 mL | 📍 RT Bench Stock
  - Item #9091: 100 mL bottle | 100 mL | 📍 RT Bench Stock

## Fixative

### [5510] Formaldehyde 37% (RNase-free)

- **Supplier/Cat#**: Thermo 28908
- **Lot**: WJ340291
- **Expiry**: 2026-11-01
- **Hazard**: GHS06, GHS08
- **Stock items**: 1
  - Item #9030: 500 mL bottle | 500 mL | 📍 Fume Hood Cabinet / Shelf 1

## RNA Reagent

### [5520] TRIzol Reagent

- **Supplier/Cat#**: Thermo 15596026
- **Lot**: 2527301
- **Expiry**: 2027-03-01
- **Hazard**: GHS05, GHS06
- **Stock items**: 1
  - Item #9060: 200 mL bottle (in use) | ~120 mL remaining | 📍 Fume Hood Cabinet / Shelf 1

## Viral Library

### [5501] Brunello sgRNA Library (lentiviral, human)

- **Supplier/Cat#**: Addgene #73178
- **Lot**: LOT-2025-0042
- **Expiry**: 2027-06-01
- **Stock items**: 2
  - Item #9001: Aliquot A — 1e9 TU/mL | 500 μL | 📍 −80°C Rack B / Box 3
  - Item #9002: Aliquot B — 1e9 TU/mL | 500 μL | 📍 −80°C Rack B / Box 3

---
*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*\
""", 'result.json': """\
{
  "skill": "labstep",
  "version": "0.2.0",
  "completed_at": "<TIMESTAMP>",
  "input_checksum": "",
  "datasets": {},
  "summary": {
    "label": "demo",
    "output": "<OUTPUT_DIR>/report.md"
  },
  "data": {}
}\
"""})
