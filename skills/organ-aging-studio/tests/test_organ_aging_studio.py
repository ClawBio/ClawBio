import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
from inline_snapshot import snapshot

from clawbio.common import demo_runner

sys.path.insert(0, str(Path(__file__).parent.parent))

SKILL_DIR = Path(__file__).resolve().parent.parent
DEMO_FIXTURE = SKILL_DIR.parent / "proteomics-clock" / "data" / "demo_olink_npx.csv.gz"
FORMULA = "predicted_age = intercept + Σ (protein_NPX × coefficient)"


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
    # locale encoding; the filter turns the ones in this script into errors.
    proc = _run_demo_cli(
        tmp_path / "out",
        "-X",
        "warn_default_encoding",
        "-W",
        "error::EncodingWarning:__main__",
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


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/organ-aging-studio/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("organ-aging-studio", tmp_path) == snapshot({'report.md': """\
# Organ Aging Studio Report

**Generated**: <TIMESTAMP> UTC
**Model**: Goeminne et al. (2025) organ-specific proteomic clocks
**Generation**: gen1

Goeminne LJE, Vladimirova A, Eames A, Tyshkovskiy A, Argentieri MA, Ying K, Moqri M, Gladyshev VN (2025). Plasma protein-based organ-specific aging and mortality models unveil diseases as accelerated aging of organismal systems. Cell Metabolism 37(1):205-222.e6. DOI: 10.1016/j.cmet.2024.10.005

## Summary

### Sample `DEMO_000`

| Organ | Predicted age | Chronological | Raw delta | Proteins used |
|-------|---------------|---------------|-----------|---------------|
| Heart | 67.5 | 66.0 | +1.5 | 7/7 |
| Brain | 41.8 | 66.0 | -24.2 | 10/118 |
| Liver | 57.6 | 66.0 | -8.4 | 0/79 |
| Immune | 57.7 | 66.0 | -8.3 | 0/124 |
| Organismal | 57.1 | 66.0 | -8.9 | 0/2347 |

### Sample `DEMO_001`

| Organ | Predicted age | Chronological | Raw delta | Proteins used |
|-------|---------------|---------------|-----------|---------------|
| Heart | 66.3 | 58.0 | +8.3 | 7/7 |
| Brain | 49.3 | 58.0 | -8.7 | 10/118 |
| Liver | 57.6 | 58.0 | -0.4 | 0/79 |
| Immune | 57.7 | 58.0 | -0.3 | 0/124 |
| Organismal | 57.1 | 58.0 | -0.9 | 0/2347 |

### Sample `DEMO_002`

| Organ | Predicted age | Chronological | Raw delta | Proteins used |
|-------|---------------|---------------|-----------|---------------|
| Heart | 62.6 | 68.0 | -5.4 | 7/7 |
| Brain | 39.9 | 68.0 | -28.1 | 10/118 |
| Liver | 57.6 | 68.0 | -10.4 | 0/79 |
| Immune | 57.7 | 68.0 | -10.3 | 0/124 |
| Organismal | 57.1 | 68.0 | -10.9 | 0/2347 |

### Sample `DEMO_003`

| Organ | Predicted age | Chronological | Raw delta | Proteins used |
|-------|---------------|---------------|-----------|---------------|
| Heart | 74.7 | 78.0 | -3.3 | 7/7 |
| Brain | 58.2 | 78.0 | -19.8 | 10/118 |
| Liver | 57.6 | 78.0 | -20.4 | 0/79 |
| Immune | 57.7 | 78.0 | -20.3 | 0/124 |
| Organismal | 57.1 | 78.0 | -20.9 | 0/2347 |

### Sample `DEMO_004`

| Organ | Predicted age | Chronological | Raw delta | Proteins used |
|-------|---------------|---------------|-----------|---------------|
| Heart | 68.4 | 57.0 | +11.4 | 7/7 |
| Brain | 50.5 | 57.0 | -6.5 | 10/118 |
| Liver | 57.6 | 57.0 | +0.6 | 0/79 |
| Immune | 57.7 | 57.0 | +0.7 | 0/124 |
| Organismal | 57.1 | 57.0 | +0.1 | 0/2347 |

### Sample `DEMO_005`

| Organ | Predicted age | Chronological | Raw delta | Proteins used |
|-------|---------------|---------------|-----------|---------------|
| Heart | 70.3 | 57.0 | +13.3 | 7/7 |
| Brain | 55.1 | 57.0 | -1.9 | 10/118 |
| Liver | 57.6 | 57.0 | +0.6 | 0/79 |
| Immune | 57.7 | 57.0 | +0.7 | 0/124 |
| Organismal | 57.1 | 57.0 | +0.1 | 0/2347 |

### Sample `DEMO_006`

| Organ | Predicted age | Chronological | Raw delta | Proteins used |
|-------|---------------|---------------|-----------|---------------|
| Heart | 69.0 | 79.0 | -10.0 | 7/7 |
| Brain | 54.1 | 79.0 | -24.9 | 10/118 |
| Liver | 57.6 | 79.0 | -21.4 | 0/79 |
| Immune | 57.7 | 79.0 | -21.3 | 0/124 |
| Organismal | 57.1 | 79.0 | -21.9 | 0/2347 |

### Sample `DEMO_007`

| Organ | Predicted age | Chronological | Raw delta | Proteins used |
|-------|---------------|---------------|-----------|---------------|
| Heart | 66.2 | 69.0 | -2.8 | 7/7 |
| Brain | 43.0 | 69.0 | -26.0 | 10/118 |
| Liver | 57.6 | 69.0 | -11.4 | 0/79 |
| Immune | 57.7 | 69.0 | -11.3 | 0/124 |
| Organismal | 57.1 | 69.0 | -11.9 | 0/2347 |

### Sample `DEMO_008`

| Organ | Predicted age | Chronological | Raw delta | Proteins used |
|-------|---------------|---------------|-----------|---------------|
| Heart | 69.4 | 54.0 | +15.4 | 7/7 |
| Brain | 54.6 | 54.0 | +0.6 | 10/118 |
| Liver | 57.6 | 54.0 | +3.6 | 0/79 |
| Immune | 57.7 | 54.0 | +3.7 | 0/124 |
| Organismal | 57.1 | 54.0 | +3.1 | 0/2347 |

### Sample `DEMO_009`

| Organ | Predicted age | Chronological | Raw delta | Proteins used |
|-------|---------------|---------------|-----------|---------------|
| Heart | 59.1 | 67.0 | -7.9 | 7/7 |
| Brain | 46.9 | 67.0 | -20.1 | 10/118 |
| Liver | 57.6 | 67.0 | -9.4 | 0/79 |
| Immune | 57.7 | 67.0 | -9.3 | 0/124 |
| Organismal | 57.1 | 67.0 | -9.9 | 0/2347 |

### Sample `DEMO_010`

| Organ | Predicted age | Chronological | Raw delta | Proteins used |
|-------|---------------|---------------|-----------|---------------|
| Heart | 65.1 | 54.0 | +11.1 | 7/7 |
| Brain | 41.6 | 54.0 | -12.4 | 10/118 |
| Liver | 57.6 | 54.0 | +3.6 | 0/79 |
| Immune | 57.7 | 54.0 | +3.7 | 0/124 |
| Organismal | 57.1 | 54.0 | +3.1 | 0/2347 |

### Sample `DEMO_011`

| Organ | Predicted age | Chronological | Raw delta | Proteins used |
|-------|---------------|---------------|-----------|---------------|
| Heart | 69.0 | 54.0 | +15.0 | 7/7 |
| Brain | 47.4 | 54.0 | -6.6 | 10/118 |
| Liver | 57.6 | 54.0 | +3.6 | 0/79 |
| Immune | 57.7 | 54.0 | +3.7 | 0/124 |
| Organismal | 57.1 | 54.0 | +3.1 | 0/2347 |

### Sample `DEMO_012`

| Organ | Predicted age | Chronological | Raw delta | Proteins used |
|-------|---------------|---------------|-----------|---------------|
| Heart | 72.9 | 63.0 | +9.9 | 7/7 |
| Brain | 36.6 | 63.0 | -26.4 | 10/118 |
| Liver | 57.6 | 63.0 | -5.4 | 0/79 |
| Immune | 57.7 | 63.0 | -5.3 | 0/124 |
| Organismal | 57.1 | 63.0 | -5.9 | 0/2347 |

### Sample `DEMO_013`

| Organ | Predicted age | Chronological | Raw delta | Proteins used |
|-------|---------------|---------------|-----------|---------------|
| Heart | 68.3 | 37.0 | +31.3 | 7/7 |
| Brain | 52.2 | 37.0 | +15.2 | 10/118 |
| Liver | 57.6 | 37.0 | +20.6 | 0/79 |
| Immune | 57.7 | 37.0 | +20.7 | 0/124 |
| Organismal | 57.1 | 37.0 | +20.1 | 0/2347 |

### Sample `DEMO_014`

| Organ | Predicted age | Chronological | Raw delta | Proteins used |
|-------|---------------|---------------|-----------|---------------|
| Heart | 72.0 | 39.0 | +33.0 | 7/7 |
| Brain | 36.1 | 39.0 | -2.9 | 10/118 |
| Liver | 57.6 | 39.0 | +18.6 | 0/79 |
| Immune | 57.7 | 39.0 | +18.7 | 0/124 |
| Organismal | 57.1 | 39.0 | +18.1 | 0/2347 |

### Sample `DEMO_015`

| Organ | Predicted age | Chronological | Raw delta | Proteins used |
|-------|---------------|---------------|-----------|---------------|
| Heart | 64.8 | 53.0 | +11.8 | 7/7 |
| Brain | 32.4 | 53.0 | -20.6 | 10/118 |
| Liver | 57.6 | 53.0 | +4.6 | 0/79 |
| Immune | 57.7 | 53.0 | +4.7 | 0/124 |
| Organismal | 57.1 | 53.0 | +4.1 | 0/2347 |

### Sample `DEMO_016`

| Organ | Predicted age | Chronological | Raw delta | Proteins used |
|-------|---------------|---------------|-----------|---------------|
| Heart | 64.0 | 48.0 | +16.1 | 7/7 |
| Brain | 62.6 | 48.0 | +14.6 | 10/118 |
| Liver | 57.6 | 48.0 | +9.6 | 0/79 |
| Immune | 57.7 | 48.0 | +9.7 | 0/124 |
| Organismal | 57.1 | 48.0 | +9.1 | 0/2347 |

### Sample `DEMO_017`

| Organ | Predicted age | Chronological | Raw delta | Proteins used |
|-------|---------------|---------------|-----------|---------------|
| Heart | 68.9 | 64.0 | +4.9 | 7/7 |
| Brain | 49.9 | 64.0 | -14.1 | 10/118 |
| Liver | 57.6 | 64.0 | -6.4 | 0/79 |
| Immune | 57.7 | 64.0 | -6.3 | 0/124 |
| Organismal | 57.1 | 64.0 | -6.9 | 0/2347 |

### Sample `DEMO_018`

| Organ | Predicted age | Chronological | Raw delta | Proteins used |
|-------|---------------|---------------|-----------|---------------|
| Heart | 63.7 | 49.0 | +14.7 | 7/7 |
| Brain | 42.9 | 49.0 | -6.2 | 10/118 |
| Liver | 57.6 | 49.0 | +8.6 | 0/79 |
| Immune | 57.7 | 49.0 | +8.7 | 0/124 |
| Organismal | 57.1 | 49.0 | +8.1 | 0/2347 |

### Sample `DEMO_019`

| Organ | Predicted age | Chronological | Raw delta | Proteins used |
|-------|---------------|---------------|-----------|---------------|
| Heart | 71.3 | 43.0 | +28.3 | 7/7 |
| Brain | 49.1 | 43.0 | +6.1 | 10/118 |
| Liver | 57.6 | 43.0 | +14.6 | 0/79 |
| Immune | 57.7 | 43.0 | +14.7 | 0/124 |
| Organismal | 57.1 | 43.0 | +14.1 | 0/2347 |

## Input Sanity

- No obvious unit issues detected.
## How prediction works

```
predicted_age = intercept + Σ (protein_NPX × coefficient)
```

Coefficients are downloaded from the pinned [organAging](https://github.com/ludgergoeminne/organAging) repository.
Delta values are the raw predicted-minus-chronological gap, not age-residualised acceleration.

ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.\
""", 'tables/protein_contributions.csv': """\
sample_id,organ,protein,npx,coefficient,contribution
DEMO_000,Heart,BMP10,2.099,2.279733,4.7852
DEMO_000,Heart,NPPB,3.386,1.180802,3.9982
DEMO_000,Heart,NMRK2,2.205,-0.665022,-1.4664
DEMO_000,Heart,MYL4,3.542,0.457943,1.622
DEMO_000,Heart,TNNI3,3.783,0.319165,1.2074
DEMO_000,Heart,ITGB1BP2,3.515,-0.090852,-0.3193
DEMO_000,Heart,PXDNL,3.149,0.090366,0.2846
DEMO_000,Brain,BCAN,3.39,-5.156511,-17.4806
DEMO_000,Brain,ADAM22,3.147,1.261917,3.9713
DEMO_000,Brain,ABCA2,1.922,0.605362,1.1635
DEMO_000,Brain,AMBN,3.18,-0.452049,-1.4375
DEMO_000,Brain,ADCYAP1R1,4.426,-0.414915,-1.8364
DEMO_000,Brain,ANGPTL7,4.588,0.315781,1.4488
DEMO_000,Brain,APOD,2.204,0.220026,0.4849
DEMO_000,Brain,ATP1B2,5.531,-0.216323,-1.1965
DEMO_000,Brain,APLP1,2.229,-0.126028,-0.2809
DEMO_000,Brain,ATP6V1G2,3.073,-0.030568,-0.0939
DEMO_001,Heart,BMP10,2.562,2.279733,5.8407
DEMO_001,Heart,NPPB,2.888,1.180802,3.4102
DEMO_001,Heart,NMRK2,3.77,-0.665022,-2.5071
DEMO_001,Heart,MYL4,2.032,0.457943,0.9305
DEMO_001,Heart,TNNI3,3.445,0.319165,1.0995
DEMO_001,Heart,ITGB1BP2,0.355,-0.090852,-0.0323
DEMO_001,Heart,PXDNL,2.245,0.090366,0.2029
DEMO_001,Brain,BCAN,1.644,-5.156511,-8.4773
DEMO_001,Brain,ADAM22,1.84,1.261917,2.3219
DEMO_001,Brain,ABCA2,2.68,0.605362,1.6224
DEMO_001,Brain,AMBN,3.772,-0.452049,-1.7051
DEMO_001,Brain,ADCYAP1R1,2.135,-0.414915,-0.8858
DEMO_001,Brain,ANGPTL7,0.362,0.315781,0.1143
DEMO_001,Brain,APOD,2.136,0.220026,0.47
DEMO_001,Brain,ATP1B2,4.322,-0.216323,-0.9349
DEMO_001,Brain,APLP1,1.411,-0.126028,-0.1778
DEMO_001,Brain,ATP6V1G2,1.754,-0.030568,-0.0536
DEMO_002,Heart,BMP10,2.097,2.279733,4.7806
DEMO_002,Heart,NPPB,0.122,1.180802,0.1441
DEMO_002,Heart,NMRK2,3.146,-0.665022,-2.0922
DEMO_002,Heart,MYL4,3.542,0.457943,1.622
DEMO_002,Heart,TNNI3,3.376,0.319165,1.0775
DEMO_002,Heart,ITGB1BP2,3.486,-0.090852,-0.3167
DEMO_002,Heart,PXDNL,0.674,0.090366,0.0609
DEMO_002,Brain,BCAN,3.958,-5.156511,-20.4095
DEMO_002,Brain,ADAM22,3.037,1.261917,3.8324
DEMO_002,Brain,ABCA2,3.466,0.605362,2.0982
DEMO_002,Brain,AMBN,4.067,-0.452049,-1.8385
DEMO_002,Brain,ADCYAP1R1,1.652,-0.414915,-0.6854
DEMO_002,Brain,ANGPTL7,1.225,0.315781,0.3868
DEMO_002,Brain,APOD,2.587,0.220026,0.5692
DEMO_002,Brain,ATP1B2,2.988,-0.216323,-0.6464
DEMO_002,Brain,APLP1,2.906,-0.126028,-0.3662
DEMO_002,Brain,ATP6V1G2,3.406,-0.030568,-0.1041
DEMO_003,Heart,BMP10,5.778,2.279733,13.1723
DEMO_003,Heart,NPPB,2.96,1.180802,3.4952
DEMO_003,Heart,NMRK2,4.453,-0.665022,-2.9613
DEMO_003,Heart,MYL4,5.307,0.457943,2.4303
DEMO_003,Heart,TNNI3,3.52,0.319165,1.1235
DEMO_003,Heart,ITGB1BP2,2.422,-0.090852,-0.22
DEMO_003,Heart,PXDNL,3.103,0.090366,0.2804
DEMO_003,Brain,BCAN,0.508,-5.156511,-2.6195
DEMO_003,Brain,ADAM22,3.747,1.261917,4.7284
DEMO_003,Brain,ABCA2,5.213,0.605362,3.1558
DEMO_003,Brain,AMBN,1.313,-0.452049,-0.5935
DEMO_003,Brain,ADCYAP1R1,3.738,-0.414915,-1.551
DEMO_003,Brain,ANGPTL7,-0.059,0.315781,-0.0186
DEMO_003,Brain,APOD,-0.453,0.220026,-0.0997
DEMO_003,Brain,ATP1B2,5.22,-0.216323,-1.1292
DEMO_003,Brain,APLP1,4.433,-0.126028,-0.5587
DEMO_003,Brain,ATP6V1G2,2.925,-0.030568,-0.0894
DEMO_004,Heart,BMP10,2.98,2.279733,6.7936
DEMO_004,Heart,NPPB,3.09,1.180802,3.6487
DEMO_004,Heart,NMRK2,1.947,-0.665022,-1.2948
DEMO_004,Heart,MYL4,2.946,0.457943,1.3491
DEMO_004,Heart,TNNI3,1.98,0.319165,0.6319
DEMO_004,Heart,ITGB1BP2,1.985,-0.090852,-0.1803
DEMO_004,Heart,PXDNL,1.407,0.090366,0.1271
DEMO_004,Brain,BCAN,2.901,-5.156511,-14.959
DEMO_004,Brain,ADAM22,5.177,1.261917,6.5329
DEMO_004,Brain,ABCA2,4.286,0.605362,2.5946
DEMO_004,Brain,AMBN,0.699,-0.452049,-0.316
DEMO_004,Brain,ADCYAP1R1,1.02,-0.414915,-0.4232
DEMO_004,Brain,ANGPTL7,2.596,0.315781,0.8198
DEMO_004,Brain,APOD,0.727,0.220026,0.16
DEMO_004,Brain,ATP1B2,3.116,-0.216323,-0.6741
DEMO_004,Brain,APLP1,1.521,-0.126028,-0.1917
DEMO_004,Brain,ATP6V1G2,2.642,-0.030568,-0.0808
DEMO_005,Heart,BMP10,1.413,2.279733,3.2213
DEMO_005,Heart,NPPB,6.695,1.180802,7.9055
DEMO_005,Heart,NMRK2,2.509,-0.665022,-1.6685
DEMO_005,Heart,MYL4,5.347,0.457943,2.4486
DEMO_005,Heart,TNNI3,3.348,0.319165,1.0686
DEMO_005,Heart,ITGB1BP2,3.918,-0.090852,-0.356
DEMO_005,Heart,PXDNL,3.71,0.090366,0.3353
DEMO_005,Brain,BCAN,1.183,-5.156511,-6.1002
DEMO_005,Brain,ADAM22,4.439,1.261917,5.6016
DEMO_005,Brain,ABCA2,2.76,0.605362,1.6708
DEMO_005,Brain,AMBN,4.917,-0.452049,-2.2227
DEMO_005,Brain,ADCYAP1R1,5.747,-0.414915,-2.3845
DEMO_005,Brain,ANGPTL7,4.076,0.315781,1.2871
DEMO_005,Brain,APOD,5.05,0.220026,1.1111
DEMO_005,Brain,ATP1B2,1.708,-0.216323,-0.3695
DEMO_005,Brain,APLP1,3.756,-0.126028,-0.4734
DEMO_005,Brain,ATP6V1G2,1.639,-0.030568,-0.0501
DEMO_006,Heart,BMP10,4.234,2.279733,9.6524
DEMO_006,Heart,NPPB,2.711,1.180802,3.2012
DEMO_006,Heart,NMRK2,2.412,-0.665022,-1.604
DEMO_006,Heart,MYL4,-0.93,0.457943,-0.4259
DEMO_006,Heart,TNNI3,3.44,0.319165,1.0979
DEMO_006,Heart,ITGB1BP2,4.546,-0.090852,-0.413
DEMO_006,Heart,PXDNL,1.621,0.090366,0.1465
DEMO_006,Brain,BCAN,2.022,-5.156511,-10.4265
DEMO_006,Brain,ADAM22,6.23,1.261917,7.8617
DEMO_006,Brain,ABCA2,2.971,0.605362,1.7985
DEMO_006,Brain,AMBN,3.498,-0.452049,-1.5813
DEMO_006,Brain,ADCYAP1R1,4.769,-0.414915,-1.9787
DEMO_006,Brain,ANGPTL7,5.254,0.315781,1.6591
DEMO_006,Brain,APOD,5.467,0.220026,1.2029
DEMO_006,Brain,ATP1B2,5.285,-0.216323,-1.1433
DEMO_006,Brain,APLP1,2.205,-0.126028,-0.2779
DEMO_006,Brain,ATP6V1G2,2.135,-0.030568,-0.0653
DEMO_007,Heart,BMP10,1.169,2.279733,2.665
DEMO_007,Heart,NPPB,3.452,1.180802,4.0761
DEMO_007,Heart,NMRK2,0.805,-0.665022,-0.5353
DEMO_007,Heart,MYL4,4.233,0.457943,1.9385
DEMO_007,Heart,TNNI3,1.928,0.319165,0.6154
DEMO_007,Heart,ITGB1BP2,4.397,-0.090852,-0.3995
DEMO_007,Heart,PXDNL,5.325,0.090366,0.4812
DEMO_007,Brain,BCAN,3.071,-5.156511,-15.8356
DEMO_007,Brain,ADAM22,1.849,1.261917,2.3333
DEMO_007,Brain,ABCA2,1.496,0.605362,0.9056
DEMO_007,Brain,AMBN,1.877,-0.452049,-0.8485
DEMO_007,Brain,ADCYAP1R1,2.296,-0.414915,-0.9526
DEMO_007,Brain,ANGPTL7,3.111,0.315781,0.9824
DEMO_007,Brain,APOD,2.626,0.220026,0.5778
DEMO_007,Brain,ATP1B2,3.808,-0.216323,-0.8238
DEMO_007,Brain,APLP1,1.811,-0.126028,-0.2282
DEMO_007,Brain,ATP6V1G2,4.133,-0.030568,-0.1263
DEMO_008,Heart,BMP10,3.313,2.279733,7.5528
DEMO_008,Heart,NPPB,2.948,1.180802,3.481
DEMO_008,Heart,NMRK2,3.444,-0.665022,-2.2903
DEMO_008,Heart,MYL4,3.131,0.457943,1.4338
DEMO_008,Heart,TNNI3,5.799,0.319165,1.8508
DEMO_008,Heart,ITGB1BP2,1.741,-0.090852,-0.1582
DEMO_008,Heart,PXDNL,1.825,0.090366,0.1649
DEMO_008,Brain,BCAN,1.709,-5.156511,-8.8125
DEMO_008,Brain,ADAM22,4.308,1.261917,5.4363
DEMO_008,Brain,ABCA2,2.972,0.605362,1.7991
DEMO_008,Brain,AMBN,5.327,-0.452049,-2.4081
DEMO_008,Brain,ADCYAP1R1,0.43,-0.414915,-0.1784
DEMO_008,Brain,ANGPTL7,5.443,0.315781,1.7188
DEMO_008,Brain,APOD,3.865,0.220026,0.8504
DEMO_008,Brain,ATP1B2,1.444,-0.216323,-0.3124
DEMO_008,Brain,APLP1,2.839,-0.126028,-0.3578
DEMO_008,Brain,ATP6V1G2,3.751,-0.030568,-0.1147
DEMO_009,Heart,BMP10,0.06,2.279733,0.1368
DEMO_009,Heart,NPPB,1.247,1.180802,1.4725
DEMO_009,Heart,NMRK2,3.392,-0.665022,-2.2558
DEMO_009,Heart,MYL4,2.551,0.457943,1.1682
DEMO_009,Heart,TNNI3,3.711,0.319165,1.1844
DEMO_009,Heart,ITGB1BP2,2.536,-0.090852,-0.2304
DEMO_009,Heart,PXDNL,2.517,0.090366,0.2275
DEMO_009,Brain,BCAN,2.423,-5.156511,-12.4942
DEMO_009,Brain,ADAM22,3.275,1.261917,4.1328
DEMO_009,Brain,ABCA2,2.567,0.605362,1.554
DEMO_009,Brain,AMBN,3.174,-0.452049,-1.4348
DEMO_009,Brain,ADCYAP1R1,5.031,-0.414915,-2.0874
DEMO_009,Brain,ANGPTL7,0.93,0.315781,0.2937
DEMO_009,Brain,APOD,3.467,0.220026,0.7628
DEMO_009,Brain,ATP1B2,2.714,-0.216323,-0.5871
DEMO_009,Brain,APLP1,1.447,-0.126028,-0.1824
DEMO_009,Brain,ATP6V1G2,1.534,-0.030568,-0.0469
DEMO_010,Heart,BMP10,1.008,2.279733,2.298
DEMO_010,Heart,NPPB,4.714,1.180802,5.5663
DEMO_010,Heart,NMRK2,3.008,-0.665022,-2.0004
DEMO_010,Heart,MYL4,3.138,0.457943,1.437
DEMO_010,Heart,TNNI3,1.213,0.319165,0.3871
DEMO_010,Heart,ITGB1BP2,3.497,-0.090852,-0.3177
DEMO_010,Heart,PXDNL,4.22,0.090366,0.3813
DEMO_010,Brain,BCAN,4.509,-5.156511,-23.2507
DEMO_010,Brain,ADAM22,6.285,1.261917,7.9311
DEMO_010,Brain,ABCA2,3.484,0.605362,2.1091
DEMO_010,Brain,AMBN,4.769,-0.452049,-2.1558
DEMO_010,Brain,ADCYAP1R1,2.828,-0.414915,-1.1734
DEMO_010,Brain,ANGPTL7,0.445,0.315781,0.1405
DEMO_010,Brain,APOD,7.618,0.220026,1.6762
DEMO_010,Brain,ATP1B2,1.687,-0.216323,-0.3649
DEMO_010,Brain,APLP1,2.17,-0.126028,-0.2735
DEMO_010,Brain,ATP6V1G2,3.149,-0.030568,-0.0963
DEMO_011,Heart,BMP10,3.295,2.279733,7.5117
DEMO_011,Heart,NPPB,4.128,1.180802,4.8744
DEMO_011,Heart,NMRK2,2.648,-0.665022,-1.761
DEMO_011,Heart,MYL4,0.019,0.457943,0.0087
DEMO_011,Heart,TNNI3,3.985,0.319165,1.2719
DEMO_011,Heart,ITGB1BP2,4.463,-0.090852,-0.4055
DEMO_011,Heart,PXDNL,1.154,0.090366,0.1043
DEMO_011,Brain,BCAN,2.135,-5.156511,-11.0092
DEMO_011,Brain,ADAM22,1.788,1.261917,2.2563
DEMO_011,Brain,ABCA2,1.759,0.605362,1.0648
DEMO_011,Brain,AMBN,3.101,-0.452049,-1.4018
DEMO_011,Brain,ADCYAP1R1,4.857,-0.414915,-2.0152
DEMO_011,Brain,ANGPTL7,2.917,0.315781,0.9211
DEMO_011,Brain,APOD,4.679,0.220026,1.0295
DEMO_011,Brain,ATP1B2,0.926,-0.216323,-0.2003
DEMO_011,Brain,APLP1,1.203,-0.126028,-0.1516
DEMO_011,Brain,ATP6V1G2,4.127,-0.030568,-0.1262
DEMO_012,Heart,BMP10,4.108,2.279733,9.3651
DEMO_012,Heart,NPPB,4.187,1.180802,4.944
DEMO_012,Heart,NMRK2,0.877,-0.665022,-0.5832
DEMO_012,Heart,MYL4,2.67,0.457943,1.2227
DEMO_012,Heart,TNNI3,1.538,0.319165,0.4909
DEMO_012,Heart,ITGB1BP2,2.281,-0.090852,-0.2072
DEMO_012,Heart,PXDNL,3.341,0.090366,0.3019
DEMO_012,Brain,BCAN,4.254,-5.156511,-21.9358
DEMO_012,Brain,ADAM22,1.74,1.261917,2.1957
DEMO_012,Brain,ABCA2,3.779,0.605362,2.2877
DEMO_012,Brain,AMBN,6.091,-0.452049,-2.7534
DEMO_012,Brain,ADCYAP1R1,0.608,-0.414915,-0.2523
DEMO_012,Brain,ANGPTL7,3.576,0.315781,1.1292
DEMO_012,Brain,APOD,2.808,0.220026,0.6178
DEMO_012,Brain,ATP1B2,4.389,-0.216323,-0.9494
DEMO_012,Brain,APLP1,5.947,-0.126028,-0.7495
DEMO_012,Brain,ATP6V1G2,0.496,-0.030568,-0.0152
DEMO_013,Heart,BMP10,3.257,2.279733,7.4251
DEMO_013,Heart,NPPB,1.636,1.180802,1.9318
DEMO_013,Heart,NMRK2,2.369,-0.665022,-1.5754
DEMO_013,Heart,MYL4,3.536,0.457943,1.6193
DEMO_013,Heart,TNNI3,4.181,0.319165,1.3344
DEMO_013,Heart,ITGB1BP2,2.722,-0.090852,-0.2473
DEMO_013,Heart,PXDNL,4.961,0.090366,0.4483
DEMO_013,Brain,BCAN,1.305,-5.156511,-6.7292
DEMO_013,Brain,ADAM22,2.101,1.261917,2.6513
DEMO_013,Brain,ABCA2,5.299,0.605362,3.2078
DEMO_013,Brain,AMBN,5.633,-0.452049,-2.5464
DEMO_013,Brain,ADCYAP1R1,2.101,-0.414915,-0.8717
DEMO_013,Brain,ANGPTL7,2.951,0.315781,0.9319
DEMO_013,Brain,APOD,1.567,0.220026,0.3448
DEMO_013,Brain,ATP1B2,5.864,-0.216323,-1.2685
DEMO_013,Brain,APLP1,3.053,-0.126028,-0.3848
DEMO_013,Brain,ATP6V1G2,3.815,-0.030568,-0.1166
DEMO_014,Heart,BMP10,2.827,2.279733,6.4448
DEMO_014,Heart,NPPB,5.104,1.180802,6.0268
DEMO_014,Heart,NMRK2,2.486,-0.665022,-1.6532
DEMO_014,Heart,MYL4,5.217,0.457943,2.3891
DEMO_014,Heart,TNNI3,4.738,0.319165,1.5122
DEMO_014,Heart,ITGB1BP2,1.34,-0.090852,-0.1217
DEMO_014,Heart,PXDNL,0.589,0.090366,0.0532
DEMO_014,Brain,BCAN,3.795,-5.156511,-19.569
DEMO_014,Brain,ADAM22,-0.186,1.261917,-0.2347
DEMO_014,Brain,ABCA2,2.837,0.605362,1.7174
DEMO_014,Brain,AMBN,2.627,-0.452049,-1.1875
DEMO_014,Brain,ADCYAP1R1,3.008,-0.414915,-1.2481
DEMO_014,Brain,ANGPTL7,-0.101,0.315781,-0.0319
DEMO_014,Brain,APOD,0.59,0.220026,0.1298
DEMO_014,Brain,ATP1B2,0.902,-0.216323,-0.1951
DEMO_014,Brain,APLP1,1.95,-0.126028,-0.2458
DEMO_014,Brain,ATP6V1G2,2.006,-0.030568,-0.0613
DEMO_015,Heart,BMP10,2.548,2.279733,5.8088
DEMO_015,Heart,NPPB,0.897,1.180802,1.0592
DEMO_015,Heart,NMRK2,1.797,-0.665022,-1.195
DEMO_015,Heart,MYL4,2.223,0.457943,1.018
DEMO_015,Heart,TNNI3,1.769,0.319165,0.5646
DEMO_015,Heart,ITGB1BP2,1.206,-0.090852,-0.1096
DEMO_015,Heart,PXDNL,3.277,0.090366,0.2961
DEMO_015,Brain,BCAN,5.162,-5.156511,-26.6179
DEMO_015,Brain,ADAM22,2.211,1.261917,2.7901
DEMO_015,Brain,ABCA2,3.603,0.605362,2.1811
DEMO_015,Brain,AMBN,4.457,-0.452049,-2.0148
DEMO_015,Brain,ADCYAP1R1,3.07,-0.414915,-1.2738
DEMO_015,Brain,ANGPTL7,2.866,0.315781,0.905
DEMO_015,Brain,APOD,3.305,0.220026,0.7272
DEMO_015,Brain,ATP1B2,3.844,-0.216323,-0.8315
DEMO_015,Brain,APLP1,3.321,-0.126028,-0.4185
DEMO_015,Brain,ATP6V1G2,3.856,-0.030568,-0.1179
DEMO_016,Heart,BMP10,0.782,2.279733,1.7828
DEMO_016,Heart,NPPB,3.88,1.180802,4.5815
DEMO_016,Heart,NMRK2,2.758,-0.665022,-1.8341
DEMO_016,Heart,MYL4,1.787,0.457943,0.8183
DEMO_016,Heart,TNNI3,4.445,0.319165,1.4187
DEMO_016,Heart,ITGB1BP2,4.219,-0.090852,-0.3833
DEMO_016,Heart,PXDNL,3.39,0.090366,0.3063
DEMO_016,Brain,BCAN,-0.707,-5.156511,3.6457
DEMO_016,Brain,ADAM22,1.861,1.261917,2.3484
DEMO_016,Brain,ABCA2,4.035,0.605362,2.4426
DEMO_016,Brain,AMBN,3.968,-0.452049,-1.7937
DEMO_016,Brain,ADCYAP1R1,2.325,-0.414915,-0.9647
DEMO_016,Brain,ANGPTL7,1.043,0.315781,0.3294
DEMO_016,Brain,APOD,1.865,0.220026,0.4103
DEMO_016,Brain,ATP1B2,2.024,-0.216323,-0.4378
DEMO_016,Brain,APLP1,2.832,-0.126028,-0.3569
DEMO_016,Brain,ATP6V1G2,1.855,-0.030568,-0.0567
DEMO_017,Heart,BMP10,1.92,2.279733,4.3771
DEMO_017,Heart,NPPB,6.286,1.180802,7.4225
DEMO_017,Heart,NMRK2,3.606,-0.665022,-2.3981
DEMO_017,Heart,MYL4,2.247,0.457943,1.029
DEMO_017,Heart,TNNI3,3.619,0.319165,1.1551
DEMO_017,Heart,ITGB1BP2,5.034,-0.090852,-0.4573
DEMO_017,Heart,PXDNL,4.173,0.090366,0.3771
DEMO_017,Brain,BCAN,1.805,-5.156511,-9.3075
DEMO_017,Brain,ADAM22,3.226,1.261917,4.0709
DEMO_017,Brain,ABCA2,2.398,0.605362,1.4517
DEMO_017,Brain,AMBN,5.053,-0.452049,-2.2842
DEMO_017,Brain,ADCYAP1R1,3.934,-0.414915,-1.6323
DEMO_017,Brain,ANGPTL7,4.005,0.315781,1.2647
DEMO_017,Brain,APOD,0.867,0.220026,0.1908
DEMO_017,Brain,ATP1B2,2.269,-0.216323,-0.4908
DEMO_017,Brain,APLP1,2.669,-0.126028,-0.3364
DEMO_017,Brain,ATP6V1G2,0.293,-0.030568,-0.009
DEMO_018,Heart,BMP10,2.309,2.279733,5.2639
DEMO_018,Heart,NPPB,1.514,1.180802,1.7877
DEMO_018,Heart,NMRK2,5.829,-0.665022,-3.8764
DEMO_018,Heart,MYL4,4.373,0.457943,2.0026
DEMO_018,Heart,TNNI3,4.233,0.319165,1.351
DEMO_018,Heart,ITGB1BP2,2.892,-0.090852,-0.2627
DEMO_018,Heart,PXDNL,1.145,0.090366,0.1035
DEMO_018,Brain,BCAN,3.866,-5.156511,-19.9351
DEMO_018,Brain,ADAM22,3.513,1.261917,4.4331
DEMO_018,Brain,ABCA2,3.336,0.605362,2.0195
DEMO_018,Brain,AMBN,1.553,-0.452049,-0.702
DEMO_018,Brain,ADCYAP1R1,1.399,-0.414915,-0.5805
DEMO_018,Brain,ANGPTL7,3.55,0.315781,1.121
DEMO_018,Brain,APOD,2.03,0.220026,0.4467
DEMO_018,Brain,ATP1B2,2.111,-0.216323,-0.4567
DEMO_018,Brain,APLP1,3.921,-0.126028,-0.4942
DEMO_018,Brain,ATP6V1G2,0.559,-0.030568,-0.0171
DEMO_019,Heart,BMP10,4.586,2.279733,10.4549
DEMO_019,Heart,NPPB,2.151,1.180802,2.5399
DEMO_019,Heart,NMRK2,3.262,-0.665022,-2.1693
DEMO_019,Heart,MYL4,3.493,0.457943,1.5996
DEMO_019,Heart,TNNI3,5.845,0.319165,1.8655
DEMO_019,Heart,ITGB1BP2,4.505,-0.090852,-0.4093
DEMO_019,Heart,PXDNL,1.019,0.090366,0.0921
DEMO_019,Brain,BCAN,2.695,-5.156511,-13.8968
DEMO_019,Brain,ADAM22,5.814,1.261917,7.3368
DEMO_019,Brain,ABCA2,3.019,0.605362,1.8276
DEMO_019,Brain,AMBN,4.029,-0.452049,-1.8213
DEMO_019,Brain,ADCYAP1R1,2.786,-0.414915,-1.156
DEMO_019,Brain,ANGPTL7,1.59,0.315781,0.5021
DEMO_019,Brain,APOD,1.378,0.220026,0.3032
DEMO_019,Brain,ATP1B2,1.704,-0.216323,-0.3686
DEMO_019,Brain,APLP1,4.136,-0.126028,-0.5213
DEMO_019,Brain,ATP6V1G2,3.072,-0.030568,-0.0939
"""})
