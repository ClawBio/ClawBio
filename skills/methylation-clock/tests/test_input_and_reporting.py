"""Input and report contracts that run without downloading clock weights."""

import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from hypothesis import given, strategies as st

sys.path.insert(0, str(Path(__file__).parent.parent))

import methylation_clock as clock  # noqa: E402


def test_aggregation_sends_only_probe_columns_to_pyaging(monkeypatch):
    seen = []

    def fake_aggregation(probes, *, verbose):
        seen.extend(probes.columns)
        return probes.rename(columns=lambda col: col.split("_")[0])

    monkeypatch.setattr(
        clock,
        "pya",
        SimpleNamespace(pp=SimpleNamespace(epicv2_probe_aggregation=fake_aggregation)),
    )
    frame = pd.DataFrame({"tissue_type": ["blood"], "cg00000108_BC21": [0.4]})

    result = clock._aggregate_epicv2_probes(frame, verbose=False)

    assert seen == ["cg00000108_BC21"]
    assert result["tissue_type"].tolist() == ["blood"]


def test_geo_download_stays_under_output_directory(tmp_path, monkeypatch):
    def fake_download(geo_id, *, dir, verbose):
        destination = Path(dir) / f"{geo_id}.pkl"
        destination.write_bytes(b"fixture")
        return str(destination)

    monkeypatch.setattr(
        clock,
        "pya",
        SimpleNamespace(data=SimpleNamespace(download_example_data=fake_download)),
    )

    downloaded = clock._resolve_geo_pickle("GSE139307", verbose=False, output_dir=tmp_path)

    assert downloaded == tmp_path / "input" / "GSE139307.pkl"
    assert downloaded.read_bytes() == b"fixture"


def test_csv_sample_ids_become_index(tmp_path):
    source = tmp_path / "samples.csv"
    pd.DataFrame(
        {"sample_id": ["person-A", "person-B"], "cg00000029": [0.2, 0.8]}
    ).to_csv(source, index=False)

    loaded = clock._load_dataframe(source)

    assert loaded.index.tolist() == ["person-A", "person-B"]
    assert "sample_id" not in loaded.columns


def test_duplicate_sample_ids_are_rejected(tmp_path):
    source = tmp_path / "samples.csv"
    pd.DataFrame(
        {"sample_id": ["person-A", "person-A"], "cg00000029": [0.2, 0.8]}
    ).to_csv(source, index=False)

    with pytest.raises(ValueError, match="sample_id.*unique"):
        clock._load_dataframe(source)


@given(st.floats(min_value=1.001, max_value=10, allow_nan=False, allow_infinity=False))
def test_beta_values_above_one_are_rejected(value):
    frame = pd.DataFrame({"cg00000029": [0.2, value]}, index=["s1", "s2"])
    with pytest.raises(ValueError, match="outside.*0.*1"):
        clock._validate_methylation_frame(frame, [])


def test_nonnumeric_probe_is_rejected_before_pyaging():
    frame = pd.DataFrame({"cg00000029": ["0.2", "unknown"]})
    with pytest.raises(ValueError, match="numeric"):
        clock._validate_methylation_frame(frame, [])


def test_no_probe_columns_are_rejected():
    frame = pd.DataFrame({"gender": ["F", "M"]})
    with pytest.raises(ValueError, match="probe"):
        clock._validate_methylation_frame(frame, ["gender"])


def test_input_qc_counts_observed_and_missing_beta_values():
    frame = pd.DataFrame(
        {"cg00000029": [0.1, None], "cg00000108": [0.8, 0.4]},
        index=["s1", "s2"],
    )

    qc = clock._validate_methylation_frame(frame, [])

    assert qc == {
        "n_samples": 2,
        "n_probe_columns": 2,
        "n_observed_values": 3,
        "n_missing_values": 1,
        "beta_min": 0.1,
        "beta_max": 0.8,
    }


def test_gender_column_accepts_full_female_label():
    frame = pd.DataFrame({"gender": ["Female", "M", "unknown"]})
    result = clock._ensure_female_column(frame)
    assert result["female"].iloc[:2].tolist() == [1, 0]
    assert pd.isna(result["female"].iloc[2])


def test_report_shows_clock_units_and_only_existing_figures(tmp_path):
    predictions = pd.DataFrame({"sample_id": ["s1", "s2"], "horvath2013": [26.0, 28.0]})
    metadata = {"Horvath2013": {"unit": ["years"], "predicts": ["chronological age"]}}
    (tmp_path / "figures").mkdir()
    (tmp_path / "figures" / "clock_distributions.png").touch()

    clock.write_report(
        tmp_path,
        "example.csv",
        ["Horvath2013"],
        predictions,
        pd.DataFrame(columns=["clock", "feature"]),
        metadata,
    )

    report = (tmp_path / "report.md").read_text()
    assert "| Clock | Unit |" in report
    assert "| horvath2013 | years |" in report
    assert "clock_distributions.png" in report
    assert "clock_correlation.png" not in report


def test_run_analysis_writes_structured_result_and_qc(tmp_path, monkeypatch):
    source = tmp_path / "samples.csv"
    pd.DataFrame(
        {"sample_id": ["s1", "s2"], "cg00000029": [0.1, 0.8]}
    ).to_csv(source, index=False)

    class FakeAdata:
        def __init__(self, frame):
            self.obs = pd.DataFrame(index=frame.index)
            self.uns = {"horvath2013_metadata": {"unit": ["years"]}}

    def fake_predict(adata, *args, **kwargs):
        adata.obs.loc[:, "horvath2013"] = [26.0, 28.0]

    monkeypatch.setattr(clock, "pya", SimpleNamespace())
    monkeypatch.setattr(clock, "_df_to_adata", lambda df, *args, **kwargs: FakeAdata(df))
    monkeypatch.setattr(clock, "_predict_age", fake_predict)

    output = tmp_path / "report"
    result = clock.run_analysis(
        output,
        ["Horvath2013"],
        [],
        "mean",
        False,
        input_path=source,
        skip_epicv2_aggregation=True,
    )

    assert result["n_samples"] == 2
    assert (output / "result.json").exists()
    assert (output / "tables" / "input_qc.json").exists()
    assert pd.read_csv(output / "tables" / "predictions.csv")["sample_id"].tolist() == ["s1", "s2"]
    assert "Input QC" in (output / "report.md").read_text()


def test_reproducibility_checksums_cover_report_and_result(tmp_path):
    (tmp_path / "tables").mkdir()
    (tmp_path / "figures").mkdir()
    (tmp_path / "report.md").write_text("report")
    (tmp_path / "result.json").write_text("{}")
    (tmp_path / "tables" / "predictions.csv").write_text("sample_id\ns1\n")

    clock.write_reproducibility(
        tmp_path,
        input_path=tmp_path / "source.csv",
        geo_id=None,
        demo=False,
        clocks=["Horvath2013"],
        metadata_cols=[],
        imputer_strategy="mean",
        skip_epicv2_aggregation=True,
    )

    checksums = (tmp_path / "reproducibility" / "checksums.sha256").read_text()
    assert "  report.md\n" in checksums
    assert "  result.json\n" in checksums
    assert "  tables/predictions.csv\n" in checksums
    command = (tmp_path / "reproducibility" / "commands.sh").read_text()
    assert '"$OUTPUT_DIR/replay"' in command
