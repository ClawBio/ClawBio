import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

SKILL_DIR = Path(__file__).resolve().parent.parent
DEMO_FIXTURE = SKILL_DIR.parent / "proteomics-clock" / "data" / "demo_olink_npx.csv.gz"
FORMULA = "predicted_age = intercept + Σ (protein_NPX × coefficient)"

# Runs the script named by the next argument as __main__, with EncodingWarning
# raised for that script and for the clawbio package.
_RUN_WITH_ENCODING_WARNINGS_AS_ERRORS = """
import runpy, sys, warnings
warnings.filterwarnings(
    "error", category=EncodingWarning, module="(__main__|clawbio)([.]|$)"
)
sys.argv = sys.argv[1:]
runpy.run_path(sys.argv[0], run_name="__main__")
"""


def _run_demo_cli(out_dir, *python_options, env=None):
    return subprocess.run(
        [
            sys.executable,
            *python_options,
            str(SKILL_DIR / "organ_aging_studio.py"),
            "--demo",
            "--output",
            str(out_dir),
            "--organs",
            "Heart",
            "--sample-id",
            "DEMO_000",
        ],
        capture_output=True,
        text=True,
        timeout=180,
        env=env,
    )


def test_filter_protein_coefs_top_n():
    from organ_aging_studio import filter_protein_coefs

    coefs = {"A": 3.0, "B": 1.0, "C": 2.0}
    out = filter_protein_coefs(coefs, {"A", "B", "C"}, top_n=2)
    assert set(out) == {"A", "C"}


def test_compute_organ_age_gen1():
    from organ_aging_studio import compute_organ_age

    row = pd.Series({"ProteinA": 2.0, "ProteinB": 1.0, "age": 50.0})
    coefs = {"Intercept": 10.0, "ProteinA": 2.0, "ProteinB": 0.5}
    result = compute_organ_age(row, coefs, generation="gen1", top_n=10)
    # 10 + 2*2 + 0.5*1 = 14.5
    assert result["predicted_age_years"] == 14.5
    assert result["age_delta_years"] == -35.5
    assert len(result["contributions"]) == 2


def test_run_studio_demo(tmp_path):
    from organ_aging_studio import run_studio

    out_dir = tmp_path / "studio"
    result = run_studio(
        input_path=DEMO_FIXTURE,
        output_dir=out_dir,
        organs=["Heart", "Brain"],
        generation="gen1",
        sample_id="DEMO_000",
        top_n=5,
    )
    assert len(result["samples"]) == 1
    assert (out_dir / "report.md").exists()
    assert (out_dir / "result.json").exists()
    assert (out_dir / "tables" / "protein_contributions.csv").exists()

    payload = json.loads((out_dir / "result.json").read_text(encoding="utf-8"))
    assert payload["samples"][0]["organs"]["Heart"]["predicted_age_years"] > 0

    report_text = (out_dir / "report.md").read_text(encoding="utf-8")
    assert "filtered partial-model predictions" in report_text
    assert "Raw delta" in report_text
    assert "Input Sanity" in report_text


def test_cli_demo_smoke(tmp_path):
    out_dir = tmp_path / "cli_out"
    script = SKILL_DIR / "organ_aging_studio.py"
    proc = subprocess.run(
        [
            sys.executable,
            str(script),
            "--demo",
            "--output",
            str(out_dir),
            "--organs",
            "Heart,Brain",
            "--sample-id",
            "DEMO_000",
            "--top-n",
            "3",
        ],
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert proc.returncode == 0, proc.stderr
    assert (out_dir / "report.md").exists()
    assert "Organ Aging Studio complete" in proc.stdout


def test_demo_writes_verifiable_reproducibility_bundle(tmp_path):
    import hashlib

    out_dir = tmp_path / "repro_out"
    proc = subprocess.run(
        [
            sys.executable,
            str(SKILL_DIR / "organ_aging_studio.py"),
            "--demo",
            "--output",
            str(out_dir),
            "--organs",
            "Heart",
            "--sample-id",
            "DEMO_000",
        ],
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert proc.returncode == 0, proc.stderr

    repro = out_dir / "reproducibility"
    for name in ("commands.sh", "environment.yml", "checksums.sha256"):
        assert (repro / name).is_file(), name
    assert not (out_dir / "commands.sh").exists()

    commands = (repro / "commands.sh").read_text(encoding="utf-8")
    assert "python skills/organ-aging-studio/organ_aging_studio.py" in commands
    assert "--demo" in commands
    assert str(SKILL_DIR) not in commands

    lines = (repro / "checksums.sha256").read_text(encoding="utf-8").splitlines()
    labels = {line.split("  ", 1)[1] for line in lines}
    assert {"report.md", "result.json", "reproducibility/environment.yml"} <= labels
    for line in lines:
        digest, label = line.split("  ", 1)
        target = out_dir / label
        assert target.resolve().is_relative_to(out_dir.resolve()), label
        assert hashlib.sha256(target.read_bytes()).hexdigest() == digest, label


def test_cli_demo_does_not_depend_on_locale_encoding(tmp_path):
    # With locale coercion and UTF-8 mode off, the C locale gives Python an
    # ASCII locale encoding on Linux and macOS. Windows ignores LC_ALL and uses
    # its ANSI code page, such as cp1252. Neither can encode Σ or ✓.
    env = {**os.environ, "LC_ALL": "C", "PYTHONCOERCECLOCALE": "0", "PYTHONUTF8": "0"}
    env.pop("PYTHONIOENCODING", None)
    out_dir = tmp_path / "out"
    proc = _run_demo_cli(out_dir, env=env)
    assert proc.returncode == 0, proc.stderr
    assert FORMULA in (out_dir / "report.md").read_text(encoding="utf-8")
    assert "Organ Aging Studio complete" in proc.stdout


def test_cli_demo_names_the_encoding_of_its_text_files(tmp_path):
    # -X warn_default_encoding flags every text open that falls back to the
    # locale encoding. The runner turns the ones in this script and in the
    # clawbio code it calls, such as the reproducibility bundle writers, into
    # errors; other modules only print theirs. -W cannot do this, because its
    # module field must name one module exactly.
    proc = _run_demo_cli(
        tmp_path / "out",
        "-X",
        "warn_default_encoding",
        "-c",
        _RUN_WITH_ENCODING_WARNINGS_AS_ERRORS,
    )
    assert proc.returncode == 0, proc.stderr


def test_cli_requires_input_or_demo(tmp_path):
    script = SKILL_DIR / "organ_aging_studio.py"
    proc = subprocess.run(
        [
            sys.executable,
            str(script),
            "--output",
            str(tmp_path / "out"),
        ],
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert proc.returncode != 0
    assert "Provide --input or use --demo" in proc.stderr
