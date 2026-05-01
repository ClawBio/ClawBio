"""Tests for the data-extractor skill."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Add skill root to sys.path
_SKILL_DIR = Path(__file__).resolve().parent.parent
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))


def test_models_import():
    """Test that core models can be imported."""
    from data_extractor_core.models import PlotType, Confidence, Figure, DataSeries, ExtractedData

    assert PlotType.BAR.value == "bar"
    assert Confidence.HIGH.value == "high"
    assert len(PlotType) == 26


def test_plot_type_enum():
    """Test all expected plot types exist."""
    from data_extractor_core.models import PlotType

    expected = [
        "scatter", "bar", "line", "box", "violin", "histogram",
        "heatmap", "forest", "kaplan_meier", "dot_strip", "stacked_bar",
        "funnel", "roc", "volcano", "waterfall", "bland_altman",
        "paired", "bubble", "area", "dose_response", "manhattan",
        "correlation_matrix", "error_bar", "table", "other", "non_data",
    ]
    for pt in expected:
        assert PlotType(pt), f"PlotType missing: {pt}"


def test_extracted_data_model():
    """Test ExtractedData model creation."""
    from data_extractor_core.models import ExtractedData, PlotType, Confidence, DataSeries

    data = ExtractedData(
        figure_id="test_fig",
        plot_type=PlotType.BAR,
        title="Test bar chart",
        x_label="Category",
        y_label="Value",
        series=[
            DataSeries(
                name="Group A",
                x_values=["Cat1", "Cat2", "Cat3"],
                y_values=[10.5, 20.3, 15.7],
                error_bars_lower=[1.2, 2.1, 1.5],
                error_bars_upper=[1.3, 1.9, 1.6],
            ),
        ],
        confidence=Confidence.HIGH,
    )

    assert data.figure_id == "test_fig"
    assert len(data.series) == 1
    assert data.series[0].name == "Group A"
    assert len(data.series[0].y_values) == 3


def test_cv_calibration_import():
    """Test that cv_calibration module can be imported."""
    from data_extractor_core.cv_calibration import (
        PlotRegion, DetectedMarker, DetectedBar, CalibrationResult,
        calibrate_image, format_calibration_prompt,
    )

    # Test CalibrationResult creation
    result = CalibrationResult()
    assert result.markers == []
    assert result.bars == []


def test_digitizer_prompts():
    """Test that digitizer prompts are defined."""
    from data_extractor_core.digitizer import BASE_PROMPT, PLOT_GUIDANCE, DEFAULT_GUIDANCE

    assert "STEP 1" in BASE_PROMPT
    assert "STEP 4" in BASE_PROMPT
    assert len(PLOT_GUIDANCE) > 10  # Should have many plot types


def test_export_csv(tmp_path):
    """Test CSV export."""
    from data_extractor_core.models import ExtractedData, PlotType, Confidence, DataSeries
    from data_extractor_core.digitizer import export_csv

    results = [
        ExtractedData(
            figure_id="test",
            plot_type=PlotType.SCATTER,
            series=[
                DataSeries(
                    name="Series 1",
                    x_values=[1.0, 2.0, 3.0],
                    y_values=[4.0, 5.0, 6.0],
                ),
            ],
            confidence=Confidence.MEDIUM,
        ),
    ]

    csv_path = export_csv(results, str(tmp_path / "test.csv"))
    assert Path(csv_path).exists()

    content = Path(csv_path).read_text()
    assert "Series 1" in content
    assert "scatter" in content


def test_export_json(tmp_path):
    """Test JSON export."""
    import json
    from data_extractor_core.models import ExtractedData, PlotType, Confidence, DataSeries
    from data_extractor_core.digitizer import export_json

    results = [
        ExtractedData(
            figure_id="test",
            plot_type=PlotType.BAR,
            series=[
                DataSeries(name="A", x_values=["x"], y_values=[1.0]),
            ],
            confidence=Confidence.HIGH,
        ),
    ]

    json_path = export_json(results, str(tmp_path / "test.json"))
    assert Path(json_path).exists()

    data = json.loads(Path(json_path).read_text())
    assert len(data) == 1
    assert data[0]["plot_type"] == "bar"


def test_api_run_no_input():
    """Test api.run() with no input returns error."""
    from data_extractor_api import run

    result = run()
    assert result["success"] is False
    assert "No image" in result["error"]


def test_validate_extraction():
    """Test heuristic validation."""
    from data_extractor_core.models import ExtractedData, PlotType, Confidence, DataSeries
    from data_extractor_core.digitizer import validate_extraction

    # Valid result — should pass
    result = ExtractedData(
        figure_id="test",
        plot_type=PlotType.BAR,
        y_min=0, y_max=100,
        series=[
            DataSeries(name="A", x_values=["x"], y_values=[50.0]),
        ],
        confidence=Confidence.HIGH,
    )
    validated = validate_extraction(result)
    assert validated.confidence == Confidence.HIGH

    # Out-of-range value — should flag
    result2 = ExtractedData(
        figure_id="test",
        plot_type=PlotType.BAR,
        y_min=0, y_max=100,
        series=[
            DataSeries(name="A", x_values=["x"], y_values=[250.0]),
        ],
        confidence=Confidence.HIGH,
    )
    validated2 = validate_extraction(result2)
    assert "outside axis range" in (validated2.notes or "")
