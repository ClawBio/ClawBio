"""Reproducibility bundle tests for PharmGx Reporter.

Every run writes <output_dir>/reproducibility/ (commands.sh, environment.yml,
checksums.sha256) via clawbio.common.reproducibility, and
`cd <output_dir> && sha256sum -c reproducibility/checksums.sha256` passes.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import pharmgx_reporter
from clawbio.common.checksums import sha256_file


def run_demo(tmp_path, monkeypatch, *extra):
    out = tmp_path / "out"
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setattr(sys, "argv", ["pharmgx_reporter.py", "--demo", "--output", str(out), *extra])
    pharmgx_reporter.main()
    return out


def read_checksums(out):
    lines = (out / "reproducibility" / "checksums.sha256").read_text().strip().splitlines()
    return dict(reversed(line.split("  ", 1)) for line in lines)


def test_demo_writes_bundle(tmp_path, monkeypatch):
    repro = run_demo(tmp_path, monkeypatch) / "reproducibility"
    for name in ("commands.sh", "environment.yml", "checksums.sha256"):
        assert (repro / name).is_file(), name


def test_checksums_verify_from_output_dir(tmp_path, monkeypatch):
    out = run_demo(tmp_path, monkeypatch)
    entries = read_checksums(out)
    assert {"report.md", "report.html", "result.json",
            "reproducibility/commands.sh", "reproducibility/environment.yml"} <= set(entries)
    for label, digest in entries.items():
        assert not Path(label).is_absolute(), label
        assert digest == sha256_file(out / label), label


def test_commands_sh_records_this_run(tmp_path, monkeypatch):
    out = run_demo(tmp_path, monkeypatch)
    content = (out / "reproducibility" / "commands.sh").read_text()
    assert '"$CLAWBIO_ROOT/skills/pharmgx-reporter/pharmgx_reporter.py"' in content
    assert "--demo" in content
    assert '--output \\\n  "$OUTPUT_DIR"' in content
    assert str(out) not in content


def test_real_input_run_records_input_and_flags(tmp_path, monkeypatch):
    out = tmp_path / "out"
    demo = Path(pharmgx_reporter.__file__).resolve().parent / "demo_patient.txt"
    monkeypatch.setattr(sys, "argv", ["pharmgx_reporter.py", "--input", str(demo),
                                      "--output", str(out), "--no-enrich"])
    pharmgx_reporter.main()
    content = (out / "reproducibility" / "commands.sh").read_text()
    assert '"$CLAWBIO_ROOT/skills/pharmgx-reporter/demo_patient.txt"' in content
    assert "--no-enrich" in content
    assert "--demo" not in content
