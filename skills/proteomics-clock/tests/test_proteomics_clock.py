import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from hypothesis import example, given, settings, strategies as st

sys.path.insert(0, str(Path(__file__).parent.parent))

SKILL_DIR = Path(__file__).resolve().parent.parent
DEMO_FIXTURE = SKILL_DIR / "data" / "demo_olink_npx.csv.gz"


# ---------------------------------------------------------------------------
# Unit tests
# ---------------------------------------------------------------------------


def test_parse_organ_list_defaults():
    from proteomics_clock import DEFAULT_ORGANS, parse_organ_list

    organs = parse_organ_list(None)
    assert organs == DEFAULT_ORGANS
    assert "Bladder" not in organs
    assert "Heart" in organs
    assert "Brain" in organs


def test_parse_organ_list_custom():
    from proteomics_clock import parse_organ_list

    organs = parse_organ_list("Heart,Brain,Kidney")
    assert organs == ["Heart", "Brain", "Kidney"]


def test_parse_organ_list_empty_raises():
    from proteomics_clock import parse_organ_list

    with pytest.raises(ValueError):
        parse_organ_list("")


def test_mortality_to_years():
    from proteomics_clock import mortality_to_years

    # Relative hazard 0 is the cohort average: the UK Biobank mean age per the organAging README.
    assert mortality_to_years(np.array([0.0]))[0] == pytest.approx(57.29426, abs=1e-3)


@given(age=st.floats(min_value=0, max_value=120))
@example(age=0.0)
# Run on its own, the first example also pays for importing the skill (pandas, matplotlib).
@settings(deadline=None)
def test_mortality_to_years_round_trips_gompertz_age(age):
    from proteomics_clock import (
        GOMPERTZ_AVG_HAZARD,
        GOMPERTZ_INTERCEPT,
        GOMPERTZ_SLOPE,
        mortality_to_years,
    )

    relative = GOMPERTZ_INTERCEPT + GOMPERTZ_SLOPE * age - GOMPERTZ_AVG_HAZARD
    assert mortality_to_years(np.array([relative]))[0] == pytest.approx(age, abs=1e-9)


def test_predict_gen1_simple():
    from proteomics_clock import predict_organ_ages

    df = pd.DataFrame({
        "sample_id": ["S1", "S2"],
        "ProteinA": [1.0, 2.0],
        "ProteinB": [3.0, 4.0],
    })
    coefficients = {
        "TestOrgan": {"Intercept": 10.0, "ProteinA": 2.0, "ProteinB": 0.5},
    }
    preds, missing = predict_organ_ages(df, coefficients, ["TestOrgan"], "gen1")
    assert "TestOrgan" in preds.columns
    # S1: 10 + 2*1 + 0.5*3 = 13.5
    assert abs(preds.loc[0, "TestOrgan"] - 13.5) < 1e-6
    # S2: 10 + 2*2 + 0.5*4 = 16.0
    assert abs(preds.loc[1, "TestOrgan"] - 16.0) < 1e-6
    assert missing.empty or missing[missing["organ"] == "TestOrgan"].empty


def test_predict_gen2_no_intercept():
    from proteomics_clock import predict_organ_ages

    df = pd.DataFrame({
        "sample_id": ["S1"],
        "ProteinA": [1.0],
        "ProteinB": [3.0],
    })
    # Gen2 coefficients should NOT have Intercept used
    coefficients = {
        "TestOrgan": {"ProteinA": 2.0, "ProteinB": 0.5},
    }
    preds, missing = predict_organ_ages(df, coefficients, ["TestOrgan"], "gen2")
    # S1: 2*1 + 0.5*3 = 3.5
    assert abs(preds.loc[0, "TestOrgan"] - 3.5) < 1e-6


def test_missing_proteins_reported():
    from proteomics_clock import predict_organ_ages

    df = pd.DataFrame({
        "sample_id": ["S1"],
        "ProteinA": [1.0],
        # ProteinB and ProteinC are missing from input
    })
    coefficients = {
        "TestOrgan": {"Intercept": 5.0, "ProteinA": 1.0, "ProteinB": 2.0, "ProteinC": 3.0},
    }
    preds, missing = predict_organ_ages(df, coefficients, ["TestOrgan"], "gen1")
    assert not missing.empty
    missing_for_organ = missing[missing["organ"] == "TestOrgan"]
    missing_proteins = set(missing_for_organ["protein"].tolist())
    assert "ProteinB" in missing_proteins
    assert "ProteinC" in missing_proteins
    assert "ProteinA" not in missing_proteins


def test_download_coefficients_selects_fold(monkeypatch):
    import proteomics_clock as pc

    text = "Intercept,ProteinA,ProteinB\n1,2,3\n4,5,6\n"
    monkeypatch.setattr(pc, "_download_text", lambda url, cache_name: text)

    coeffs = pc.download_coefficients(["Heart"], "gen1", fold=2)
    assert coeffs["Heart_gen1"]["Intercept"] == 4.0
    assert coeffs["Heart_gen1"]["ProteinA"] == 5.0
    assert coeffs["Heart_gen1"]["ProteinB"] == 6.0

    with pytest.raises(ValueError):
        pc.download_coefficients(["Heart"], "gen1", fold=3)


class _FakeResponse:
    """Serves fixed bytes in place of requests.Response."""

    def __init__(self, content: bytes):
        self.content = content
        self.text = content.decode("utf-8")

    def raise_for_status(self) -> None:
        pass


@st.composite
def coefficient_tables(draw):
    """An Intercept and 1-3 proteins, 1-5 fold rows of finite coefficients, and a 1-based fold."""
    proteins = draw(
        st.lists(st.sampled_from(["NPPB", "BMP10", "UMOD", "GFAP"]), min_size=1, max_size=3, unique=True)
    )
    header = ["Intercept", *proteins]
    coefficient = st.floats(allow_nan=False, allow_infinity=False)
    row = st.lists(coefficient, min_size=len(header), max_size=len(header))
    rows = draw(st.lists(row, min_size=1, max_size=5))
    fold = draw(st.integers(min_value=1, max_value=len(rows)))
    return header, rows, fold


def _coefficient_csv(header, rows, newline, blank_after=()):
    """Serialise a table the way the organAging files are laid out: quoted header, one row per fold.

    Each number in blank_after adds a blank line after that line (1 is the header).
    """
    lines = [",".join(f'"{name}"' for name in header)]
    lines += [",".join(repr(value) for value in row) for row in rows]
    for line_number in sorted(blank_after, reverse=True):
        lines.insert(line_number, "")
    return "".join(line + newline for line in lines).encode("utf-8")


@given(table=coefficient_tables(), newline=st.sampled_from(["\n", "\r\n", "\r\r\n"]))
# Minimal failure found on Windows, where text-mode caching wrote CRLF as CR CR LF.
@example(table=(["Intercept", "NPPB"], [[0.0, 0.0]], 1), newline="\r\n")
# The cache those runs left behind; the old code misreads it on every OS.
@example(table=(["Intercept", "NPPB"], [[0.0, 0.0]], 1), newline="\r\r\n")
@settings(deadline=None)
def test_download_coefficients_selects_fold_from_download_and_cache(tmp_path_factory, table, newline):
    import proteomics_clock as pc

    header, rows, fold = table
    served = _coefficient_csv(header, rows, newline)
    # tmp_path and monkeypatch would be shared by every example, so each example gets a fresh cache.
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("CLAWBIO_CACHE", str(tmp_path_factory.mktemp("cache")))
        mp.setattr(pc.requests, "get", lambda url, **kwargs: _FakeResponse(served))
        downloaded = pc.download_coefficients(["Heart"], "gen1", fold=fold)
        cached = pc.download_coefficients(["Heart"], "gen1", fold=fold)

    expected = {"Heart_gen1": dict(zip(header, rows[fold - 1]))}
    assert downloaded == expected
    assert cached == expected


@pytest.mark.parametrize("fold_row", ["NA,NA,NA", ",,", "57.33,NA,NA"])
def test_download_coefficients_rejects_fold_without_protein_coefficients(monkeypatch, fold_row):
    import proteomics_clock as pc

    text = f"Intercept,NPPB,BMP10\n{fold_row}\n57.36,2.28,-0.09\n"
    monkeypatch.setattr(pc, "_download_text", lambda url, cache_name: text)

    # Unchecked, fold 1 gives the organ no protein weights, so every prediction is a silent NaN.
    with pytest.raises(ValueError, match="no numeric protein coefficients"):
        pc.download_coefficients(["Heart"], "gen1", fold=1)
    assert pc.download_coefficients(["Heart"], "gen1", fold=2) == {
        "Heart_gen1": {"Intercept": 57.36, "NPPB": 2.28, "BMP10": -0.09}
    }


@st.composite
def tables_with_blank_lines(draw):
    """A coefficient table and 1-3 line numbers (1 is the header) to put a blank line after."""
    header, rows, fold = draw(coefficient_tables())
    blank_after = draw(st.lists(st.integers(min_value=1, max_value=len(rows) + 1), min_size=1, max_size=3))
    return header, rows, fold, blank_after


@given(table=tables_with_blank_lines(), newline=st.sampled_from(["\n", "\r\n"]))
# Without the filter, a blank line before fold 2 made fold 3 silently read fold 2's coefficients.
@example(table=(["Intercept", "NPPB"], [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]], 3, [2]), newline="\r\n")
# Minimal failure without the filter: a blank line after the header was read as fold 1.
@example(table=(["Intercept", "NPPB"], [[0.0, 0.0]], 1, [1]), newline="\n")
# Run on its own, the first example also pays for importing the skill (pandas, matplotlib).
@settings(deadline=None)
def test_download_coefficients_skips_blank_lines_between_fold_rows(table, newline):
    import proteomics_clock as pc

    header, rows, fold, blank_after = table
    text = _coefficient_csv(header, rows, newline, blank_after).decode("utf-8")
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(pc, "_download_text", lambda url, cache_name: text)
        assert pc.download_coefficients(["Heart"], "gen1", fold=fold) == {
            "Heart_gen1": dict(zip(header, rows[fold - 1]))
        }
        # A blank line counted as a fold would add a fold that does not exist.
        with pytest.raises(ValueError, match=f"only has {len(rows)} fold rows"):
            pc.download_coefficients(["Heart"], "gen1", fold=len(rows) + 1)


def test_download_text_returns_same_text_from_download_and_cache(tmp_path, monkeypatch):
    import proteomics_clock as pc

    served = b'"Intercept","NPPB"\r\n57.36,2.28\r\n57.33,2.09\r\n'
    gets = []

    def fake_get(url, **kwargs):
        gets.append(url)
        return _FakeResponse(served)

    monkeypatch.setenv("CLAWBIO_CACHE", str(tmp_path))
    monkeypatch.setattr(pc.requests, "get", fake_get)
    url = pc.COEF_URL_TEMPLATES["gen1"].format(organ="Heart")

    downloaded = pc._download_text(url, "Heart_gen1_fold1.csv")
    cached = pc._download_text(url, "Heart_gen1_fold1.csv")

    cache_file = tmp_path / "proteomics-clock" / "Heart_gen1_fold1.csv"
    assert gets == [url]
    assert cache_file.read_bytes() == served
    assert pc._checksum_path(cache_file).read_text().strip() == hashlib.sha256(served).hexdigest()
    assert cached == downloaded == served.decode("utf-8")


def test_load_input_csv(tmp_path):
    from proteomics_clock import load_input

    csv_path = tmp_path / "test.csv"
    df = pd.DataFrame({"sample_id": ["S1"], "ProteinA": [1.5]})
    df.to_csv(csv_path, index=False)
    loaded = load_input(csv_path)
    assert loaded.shape == (1, 2)
    assert "ProteinA" in loaded.columns


# ---------------------------------------------------------------------------
# Integration tests
# ---------------------------------------------------------------------------


def test_run_analysis_with_demo(tmp_path):
    from proteomics_clock import run_analysis

    out_dir = tmp_path / "proteomics_report"
    result = run_analysis(
        output_dir=out_dir,
        organs=["Heart", "Brain", "Kidney"],
        generation="both",
        fold=1,
        convert_mortality_to_years=True,
        age_column="age",
        verbose=False,
        input_path=DEMO_FIXTURE,
    )

    assert result["n_samples"] >= 1
    assert (out_dir / "report.md").exists()
    assert (out_dir / "tables" / "predictions_gen1.csv").exists()
    assert (out_dir / "tables" / "predictions_gen2.csv").exists()
    assert (out_dir / "tables" / "prediction_summary.csv").exists()
    assert (out_dir / "tables" / "missing_proteins.csv").exists()
    assert (out_dir / "tables" / "clock_metadata.json").exists()

    gen1 = pd.read_csv(out_dir / "tables" / "predictions_gen1.csv")
    assert "sample_id" in gen1.columns

    report_text = (out_dir / "report.md").read_text()
    assert "ClawBio is a research and educational tool" in report_text
    assert "Goeminne" in report_text
    assert "Input Sanity" in report_text


def test_cli_smoke_with_demo(tmp_path):
    out_dir = tmp_path / "proteomics_cli"
    script = SKILL_DIR / "proteomics_clock.py"

    cmd = [
        sys.executable,
        str(script),
        "--demo",
        "--output",
        str(out_dir),
        "--organs",
        "Heart,Brain,Kidney",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

    assert proc.returncode == 0, proc.stderr
    data = json.loads(proc.stdout)
    assert data["n_samples"] >= 1
    assert (out_dir / "tables" / "predictions_gen1.csv").exists()
    assert (out_dir / "reproducibility" / "checksums.sha256").exists()
