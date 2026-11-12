"""generate_demo_bundle.py writes where it is told, never into the checked-in demo."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml

SKILL_DIR = Path(__file__).resolve().parent.parent
DEMO_DIR = SKILL_DIR / "demo"


def _snapshot(root: Path) -> dict[str, tuple[bytes, int]]:
    return {
        p.relative_to(root).as_posix(): (p.read_bytes(), p.stat().st_mtime_ns)
        for p in sorted(root.rglob("*"))
        if p.is_file()
    }


def _generate(target: Path) -> Path:
    subprocess.run(
        [sys.executable, str(SKILL_DIR / "generate_demo_bundle.py"), "--output-dir", str(target)],
        check=True, cwd=target.parent, capture_output=True,
    )
    return target


def test_regenerating_into_a_directory_leaves_the_checked_in_demo_untouched(tmp_path):
    """The test suite regenerated the bundle in place, so every run rewrote the
    checked-in demo CSVs (pandas 3 formats some floats one digit shorter)."""
    before = _snapshot(DEMO_DIR)
    target = _generate(tmp_path / "demo")
    assert _snapshot(DEMO_DIR) == before
    assert set(_snapshot(target)) == set(before)


def test_the_checked_in_demo_holds_the_values_the_generator_writes(tmp_path):
    """--demo and the CLI tests read the committed bundle, so it must still be what
    the generator writes: tables equal to float precision (pandas versions format some
    floats one digit shorter), YAML and JSON equal once parsed."""
    target = _generate(tmp_path / "demo")
    for committed in sorted(p for p in DEMO_DIR.rglob("*") if p.is_file()):
        relative = committed.relative_to(DEMO_DIR)
        fresh = target / relative
        if committed.suffix == ".csv":
            pd.testing.assert_frame_equal(
                pd.read_csv(fresh), pd.read_csv(committed), check_exact=False, rtol=1e-12, atol=0, obj=str(relative)
            )
        elif committed.suffix == ".json":
            assert json.loads(fresh.read_text()) == json.loads(committed.read_text()), relative
        else:
            assert yaml.safe_load(fresh.read_text()) == yaml.safe_load(committed.read_text()), relative
