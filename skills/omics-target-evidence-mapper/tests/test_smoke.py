from pathlib import Path
import subprocess
import sys


def test_demo_runs(tmp_path: Path) -> None:
    output_dir = tmp_path / "demo_out"
    cmd = [
        sys.executable,
        "skills/omics-target-evidence-mapper/omics_target_evidence_mapper.py",
        "--demo",
        "--output",
        str(output_dir),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert (output_dir / "report.md").exists()
    assert (output_dir / "evidence.json").exists()
    assert (output_dir / "metadata.json").exists()
