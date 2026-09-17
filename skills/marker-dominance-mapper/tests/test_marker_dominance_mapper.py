import importlib.util
import json
import subprocess
import sys
from pathlib import Path


SKILL_DIR = Path(__file__).resolve().parents[1]
MODULE_PATH = SKILL_DIR / "marker_dominance_mapper.py"
DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)


def load_module():
    spec = importlib.util.spec_from_file_location("marker_dominance_mapper", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_load_spots_rejects_missing_marker_columns(tmp_path):
    module = load_module()
    bad = tmp_path / "bad.csv"
    bad.write_text("spot_id,x,y,total_counts\nS1,0,0,100\n", encoding="utf-8")
    try:
        module.load_spots(bad)
    except ValueError as exc:
        assert "missing required columns" in str(exc).lower()
    else:
        raise AssertionError("missing marker columns should fail")


def test_empty_input_handled(tmp_path):
    module = load_module()
    empty = tmp_path / "empty.csv"
    empty.write_text("spot_id,x,y,total_counts,EPCAM,PTPRC,COL1A1,MKI67\n", encoding="utf-8")
    try:
        module.load_spots(empty)
    except ValueError as exc:
        assert "no spots" in str(exc).lower()
    else:
        raise AssertionError("empty input should fail")


def test_cli_rejects_malformed_input_without_traceback(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text("spot_id,x,y,total_counts\nS1,0,0,100\n", encoding="utf-8")
    completed = subprocess.run(
        [sys.executable, str(MODULE_PATH), "--input", str(bad), "--output", str(tmp_path / "out")],
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 2
    assert "ERROR:" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_mapping_assigns_regions_and_hotspots():
    module = load_module()
    spots = module.load_spots(SKILL_DIR / "demo_marker_counts.csv")
    result = module.map_spots(spots)
    assert result["summary"]["spot_count"] == 6
    assert result["summary"]["hotspot_count"] == 2
    regions = {spot["spot_id"]: spot["region"] for spot in result["spots"]}
    assert regions["SPOT_A1"] == "immune_edge"
    assert regions["SPOT_B2"] == "tumor_core"
    assert regions["SPOT_C2"] == "stromal_zone"
    assert result["regions"][0]["region"] == "tumor_core"


def test_demo_cli_writes_expected_outputs(tmp_path):
    out = tmp_path / "marker_out"
    completed = subprocess.run(
        [sys.executable, str(MODULE_PATH), "--demo", "--output", str(out)],
        text=True,
        capture_output=True,
        check=True,
    )
    assert "Marker Dominance Mapper" in completed.stdout
    report = (out / "report.md").read_text(encoding="utf-8")
    assert DISCLAIMER in report
    assert "tumor_core" in report
    assert "Synthetic demo data" in report
    result = json.loads((out / "result.json").read_text(encoding="utf-8"))
    assert result["skill"] == "marker-dominance-mapper"
    assert result["summary"]["hotspot_count"] == 2
    assert (out / "tables" / "mapped_spots.csv").exists()
    assert (out / "tables" / "region_summary.csv").exists()
    assert (out / "figures" / "marker_map.svg").exists()
    assert (out / "reproducibility" / "commands.sh").exists()


def test_demo_reproducibility_bundle_is_complete(tmp_path):
    out = tmp_path / "repro_out"
    subprocess.run(
        [sys.executable, str(MODULE_PATH), "--demo", "--output", str(out)],
        text=True,
        capture_output=True,
        check=True,
    )
    repro = out / "reproducibility"

    commands_text = (repro / "commands.sh").read_text(encoding="utf-8")
    assert "CLAWBIO_ROOT" in commands_text
    assert "$OUTPUT_DIR" in commands_text
    assert str(out) not in commands_text

    environment = (repro / "environment.yml").read_text(encoding="utf-8")
    assert "name: clawbio-marker-dominance-mapper" in environment

    checksum_lines = [
        line
        for line in (repro / "checksums.sha256").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    assert checksum_lines
    labels = set()
    for line in checksum_lines:
        digest, label = line.split("  ", 1)
        assert len(digest) == 64
        labels.add(label)
        assert (out / label).exists()
    assert {"report.md", "result.json", "tables/mapped_spots.csv", "tables/region_summary.csv", "figures/marker_map.svg"} <= labels


def test_demo_outputs_are_byte_stable_across_hash_seeds(tmp_path):
    """Checksums are only meaningful if two identical runs produce identical bytes."""
    import os

    digests = []
    for seed in ("0", "12345"):
        out = tmp_path / f"run_{seed}"
        env = dict(os.environ, PYTHONHASHSEED=seed)
        subprocess.run(
            [sys.executable, str(MODULE_PATH), "--demo", "--output", str(out)],
            text=True,
            capture_output=True,
            check=True,
            env=env,
        )
        digests.append((out / "reproducibility" / "checksums.sha256").read_text(encoding="utf-8"))
    assert digests[0] == digests[1]
