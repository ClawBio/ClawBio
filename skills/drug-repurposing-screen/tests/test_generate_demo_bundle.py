"""generate_demo_bundle.py writes where it is told, never into the checked-in demo."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent.parent
DEMO_DIR = SKILL_DIR / "demo"


def _snapshot(root: Path) -> dict[str, tuple[bytes, int]]:
    return {
        p.relative_to(root).as_posix(): (p.read_bytes(), p.stat().st_mtime_ns)
        for p in sorted(root.rglob("*"))
        if p.is_file()
    }


def test_regenerating_into_a_directory_leaves_the_checked_in_demo_untouched(tmp_path):
    """The test suite regenerated the bundle in place, so every run rewrote the
    checked-in demo CSVs (pandas 3 formats some floats one digit shorter)."""
    before = _snapshot(DEMO_DIR)
    target = tmp_path / "demo"
    subprocess.run(
        [sys.executable, str(SKILL_DIR / "generate_demo_bundle.py"), "--output-dir", str(target)],
        check=True, cwd=tmp_path, capture_output=True,
    )
    assert _snapshot(DEMO_DIR) == before
    assert set(_snapshot(target)) == set(before)
