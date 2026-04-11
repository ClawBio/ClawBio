"""Tests for Karpathy-style progress plotter."""
from __future__ import annotations

from pathlib import Path

import pytest

from skills.clawpathy_autoresearch.plotter import plot_progress, ExperimentRecord


@pytest.fixture
def experiment_history() -> list[ExperimentRecord]:
    return [
        ExperimentRecord(experiment=1, score=3.2, kept=True, label="baseline"),
        ExperimentRecord(experiment=2, score=2.8, kept=False, label="worse workflow order"),
        ExperimentRecord(experiment=3, score=4.1, kept=True, label="added gwas-lookup chaining"),
        ExperimentRecord(experiment=4, score=3.5, kept=False, label="removed validation step"),
        ExperimentRecord(experiment=5, score=5.3, kept=True, label="created fine-mapping skill"),
        ExperimentRecord(experiment=6, score=4.9, kept=False, label="aggressive gotcha pruning"),
        ExperimentRecord(experiment=7, score=6.1, kept=True, label="added PheWAS cross-check"),
    ]


def test_plot_creates_png(experiment_history, tmp_path: Path):
    output = tmp_path / "progress.png"
    plot_progress(experiment_history, output)
    assert output.exists()
    assert output.stat().st_size > 1000  # Non-trivial PNG


def test_plot_with_single_experiment(tmp_path: Path):
    history = [ExperimentRecord(experiment=1, score=3.0, kept=True, label="baseline")]
    output = tmp_path / "progress.png"
    plot_progress(history, output)
    assert output.exists()


def test_plot_with_no_kept(tmp_path: Path):
    history = [
        ExperimentRecord(experiment=1, score=3.0, kept=True, label="baseline"),
        ExperimentRecord(experiment=2, score=2.0, kept=False, label="bad change"),
        ExperimentRecord(experiment=3, score=1.5, kept=False, label="worse change"),
    ]
    output = tmp_path / "progress.png"
    plot_progress(history, output)
    assert output.exists()


def test_experiment_record_fields():
    rec = ExperimentRecord(experiment=1, score=5.0, kept=True, label="test")
    assert rec.experiment == 1
    assert rec.score == 5.0
    assert rec.kept is True
    assert rec.label == "test"
