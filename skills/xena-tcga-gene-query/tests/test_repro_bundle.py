"""Reproducibility bundle tests for xena-tcga-gene-query.

Every run (demo or live) must write <output_dir>/reproducibility/ with
commands.sh, environment.yml and checksums.sha256 via clawbio.common, and
`cd <output_dir> && sha256sum -c reproducibility/checksums.sha256` must pass.
"""

import hashlib
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "query_tcga_api.py"


def run_demo(out: Path) -> None:
    subprocess.run([sys.executable, str(SCRIPT), "--demo", "--output", str(out)], check=True)


def read_checksums(out: Path) -> dict[str, str]:
    lines = (out / "reproducibility" / "checksums.sha256").read_text().strip().splitlines()
    return {label: digest for digest, label in (line.split("  ", 1) for line in lines)}


def test_demo_writes_bundle(tmp_path):
    out = tmp_path / "out"
    run_demo(out)
    repro = out / "reproducibility"
    for name in ("commands.sh", "environment.yml", "checksums.sha256"):
        assert (repro / name).is_file(), f"missing reproducibility/{name}"


def test_checksums_verify_from_output_dir(tmp_path):
    out = tmp_path / "out"
    run_demo(out)
    entries = read_checksums(out)
    assert {"report.md", "result.json"} <= set(entries)
    for label, digest in entries.items():
        assert not Path(label).is_absolute(), label
        assert digest == hashlib.sha256((out / label).read_bytes()).hexdigest(), label


def test_demo_commands_sh_replays_this_run(tmp_path):
    out = tmp_path / "out"
    run_demo(out)
    content = (out / "reproducibility" / "commands.sh").read_text()
    assert '"$CLAWBIO_ROOT/skills/xena-tcga-gene-query/scripts/query_tcga_api.py"' in content
    assert "--demo" in content
    assert '"$OUTPUT_DIR"' in content
    assert str(out) not in content


def test_live_commands_sh_records_real_args_without_api_key(tmp_path):
    sys.path.insert(0, str(SCRIPT.parent))
    import query_tcga_api

    out = tmp_path / "out"
    fake = {"gene": "EGFR", "cancer": "LUSC", "log2_fold_change": 1.0, "p_value": 0.01}
    with patch.object(query_tcga_api, "_request", return_value=fake):
        query_tcga_api.main([
            "--api-key", "s3cret", "--output", str(out),
            "diff-expr", "--gene", "EGFR", "--cancer", "LUSC",
        ])
    content = (out / "reproducibility" / "commands.sh").read_text()
    assert "diff-expr" in content and "EGFR" in content and "LUSC" in content
    assert "s3cret" not in content
    assert set(read_checksums(out)) >= {"report.md", "result.json"}
