"""Every run writes reproducibility/{commands.sh,environment.yml,checksums.sha256}
and the manifest verifies from the output directory (`sha256sum -c` contract)."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

from clawbio.common.checksums import sha256_file

MODULE_PATH = Path(__file__).parents[1] / "omics_target_evidence_mapper.py"
SPEC = importlib.util.spec_from_file_location("omics_target_evidence_mapper", MODULE_PATH)
mapper = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mapper)


def run(monkeypatch, out: Path, *cli: str) -> Path:
    # Offline: every source reports "unavailable"; the bundle must still be written.
    monkeypatch.setattr(mapper, "safe_request_json", lambda *a, **k: None)
    monkeypatch.setattr(sys, "argv", ["omics_target_evidence_mapper.py", *cli, "--output", str(out)])
    mapper.main()
    return out


def test_demo_writes_bundle_that_verifies(tmp_path, monkeypatch) -> None:
    out = run(monkeypatch, tmp_path / "out", "--demo")
    repro = out / "reproducibility"
    for name in ("commands.sh", "environment.yml", "checksums.sha256"):
        assert (repro / name).is_file(), name

    lines = (repro / "checksums.sha256").read_text().splitlines()
    assert lines
    labels = set()
    for line in lines:
        digest, label = line.split("  ", 1)
        assert not Path(label).is_absolute()
        assert digest == sha256_file(out / label), label
        labels.add(label)
    assert {"reproducibility/commands.sh", "reproducibility/environment.yml", "report.md"} <= labels


def test_commands_sh_records_this_runs_arguments(tmp_path, monkeypatch) -> None:
    out = run(monkeypatch, tmp_path / "out", "--gene", "PCSK9", "--disease", "familial hypercholesterolemia", "--max-papers", "3", "--max-trials", "7")
    # The writer puts one argument per continuation line; join them back.
    content = (out / "reproducibility" / "commands.sh").read_text().replace(" \\\n  ", " ")
    assert '"$CLAWBIO_ROOT/skills/omics-target-evidence-mapper/omics_target_evidence_mapper.py"' in content
    assert "--gene PCSK9" in content
    assert "--disease 'familial hypercholesterolemia'" in content
    assert "--max-papers 3" in content and "--max-trials 7" in content
    assert '--output "$OUTPUT_DIR"' in content
    assert "--demo" not in content
    assert str(tmp_path) not in content
