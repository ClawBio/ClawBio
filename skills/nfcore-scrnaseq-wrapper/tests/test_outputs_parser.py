from __future__ import annotations

from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from errors import SkillError
from outputs_parser import build_selection_log, parse_outputs


def test_parse_outputs_prefers_cellbender_combined(tmp_path):
    upstream = tmp_path / "upstream" / "results"
    (upstream / "pipeline_info").mkdir(parents=True)
    (upstream / "multiqc" / "simpleaf").mkdir(parents=True)
    (upstream / "multiqc" / "simpleaf" / "multiqc_report.html").write_text("ok", encoding="utf-8")
    preferred = upstream / "simpleaf" / "mtx_conversions" / "combined_cellbender_filter_matrix.h5ad"
    preferred.parent.mkdir(parents=True)
    preferred.write_text("h5ad", encoding="utf-8")
    raw = upstream / "simpleaf" / "mtx_conversions" / "combined_raw_matrix.h5ad"
    raw.write_text("h5ad", encoding="utf-8")
    result = parse_outputs(tmp_path)
    assert result["preferred_h5ad"] == str(preferred)
    assert result["handoff_available"] is True
    assert result["cellbender_used"] is True


def test_parse_outputs_single_sample_fallback(tmp_path):
    upstream = tmp_path / "upstream" / "results"
    upstream.mkdir(parents=True)
    sample = upstream / "star" / "mtx_conversions" / "sampleA_filtered_matrix.h5ad"
    sample.parent.mkdir(parents=True)
    sample.write_text("h5ad", encoding="utf-8")
    result = parse_outputs(tmp_path)
    assert result["preferred_h5ad"] == str(sample)
    assert result["handoff_available"] is True


def test_parse_outputs_no_canonical_h5ad_when_multiple_samples(tmp_path):
    upstream = tmp_path / "upstream" / "results"
    upstream.mkdir(parents=True)
    for sample_name in ("sampleA", "sampleB"):
        p = upstream / "star" / "mtx_conversions" / f"{sample_name}_filtered_matrix.h5ad"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("h5ad", encoding="utf-8")
    result = parse_outputs(tmp_path)
    assert result["preferred_h5ad"] == ""
    assert result["handoff_available"] is False
    assert len(result["h5ad_candidates"]) == 2


def test_parse_outputs_raises_when_upstream_dir_missing(tmp_path):
    with pytest.raises(SkillError) as exc:
        parse_outputs(tmp_path)
    assert exc.value.error_code == "EXPECTED_OUTPUTS_NOT_FOUND"


def test_parse_outputs_detects_multiqc_and_pipeline_info(tmp_path):
    upstream = tmp_path / "upstream" / "results"
    info = upstream / "pipeline_info"
    info.mkdir(parents=True)
    mqc = upstream / "multiqc" / "star" / "multiqc_report.html"
    mqc.parent.mkdir(parents=True)
    mqc.write_text("report", encoding="utf-8")
    result = parse_outputs(tmp_path)
    assert result["multiqc_report"] == str(mqc)
    assert result["pipeline_info_dir"] == str(info)


def test_parse_outputs_exact_name_match_not_suffix(tmp_path):
    upstream = tmp_path / "upstream" / "results"
    upstream.mkdir(parents=True)
    subdir = upstream / "star"
    subdir.mkdir(parents=True)
    # Both files: the correct preferred name and a file whose path only suffix-matches but name differs
    (subdir / "not_combined_filtered_matrix.h5ad").write_text("h5ad", encoding="utf-8")
    real_preferred = subdir / "combined_filtered_matrix.h5ad"
    real_preferred.write_text("h5ad", encoding="utf-8")
    result = parse_outputs(tmp_path)
    assert result["preferred_h5ad"] == str(real_preferred)


def test_parse_outputs_filtered_preferred_over_raw(tmp_path):
    """combined_filtered_matrix.h5ad must be preferred over combined_raw_matrix.h5ad."""
    upstream = tmp_path / "upstream" / "results"
    upstream.mkdir(parents=True)
    for name in ("combined_filtered_matrix.h5ad", "combined_raw_matrix.h5ad"):
        (upstream / name).write_text("h5ad", encoding="utf-8")
    result = parse_outputs(tmp_path)
    assert Path(result["preferred_h5ad"]).name == "combined_filtered_matrix.h5ad"
    assert result["handoff_available"] is True




def test_parse_outputs_does_not_choose_between_duplicate_canonical_h5ad(tmp_path):
    upstream = tmp_path / "upstream" / "results"
    for aligner in ("star", "simpleaf"):
        candidate = upstream / aligner / "mtx_conversions" / "combined_filtered_matrix.h5ad"
        candidate.parent.mkdir(parents=True, exist_ok=True)
        candidate.write_text("h5ad", encoding="utf-8")
    result = parse_outputs(tmp_path)
    assert result["preferred_h5ad"] == ""
    assert result["handoff_available"] is False
    assert len(result["h5ad_candidates"]) == 2


def test_parse_outputs_sample_names_exclude_h5ad_suffix(tmp_path):
    upstream = tmp_path / "upstream" / "results"
    upstream.mkdir(parents=True)
    # A file that doesn't match any known suffix pattern
    unknown = upstream / "star" / "mysample.h5ad"
    unknown.parent.mkdir(parents=True)
    unknown.write_text("h5ad", encoding="utf-8")
    result = parse_outputs(tmp_path)
    # mysample.h5ad → after replacements still "mysample.h5ad" → must be filtered out
    assert "mysample.h5ad" not in result["samples_detected"]
    # Only one file, so it becomes the preferred fallback
    assert result["preferred_h5ad"] == str(unknown)


def test_selection_log_describes_selected_file():
    candidates = ["/out/star/combined_filtered_matrix.h5ad", "/out/star/combined_raw_matrix.h5ad"]
    log = build_selection_log(candidates, candidates[0])
    assert "combined_filtered_matrix.h5ad" in log
    assert "2 candidate(s)" in log


def test_selection_log_no_h5ad():
    assert build_selection_log([], "") == "No h5ad output files found."


def test_selection_log_ambiguous():
    candidates = ["/out/a/combined_filtered_matrix.h5ad", "/out/b/combined_filtered_matrix.h5ad"]
    log = build_selection_log(candidates, "")
    assert "ambiguous" in log
    assert "2 candidate(s)" in log


def test_parse_outputs_includes_selection_log(tmp_path):
    upstream = tmp_path / "upstream" / "results"
    upstream.mkdir(parents=True)
    h5ad = upstream / "combined_filtered_matrix.h5ad"
    h5ad.write_text("h5ad", encoding="utf-8")
    result = parse_outputs(tmp_path)
    assert "preferred_h5ad_selection_log" in result
    assert "combined_filtered_matrix.h5ad" in result["preferred_h5ad_selection_log"]


def test_preferred_h5ad_order_contains_no_phantom_entries():
    """combined_matrix.h5ad is never produced by the pipeline — concat_h5ad.py always
    writes <id>_<input_type>_matrix.h5ad where input_type is raw/filtered/cellbender_filter."""
    from schemas import PREFERRED_H5AD_ORDER
    assert "combined_matrix.h5ad" not in PREFERRED_H5AD_ORDER
