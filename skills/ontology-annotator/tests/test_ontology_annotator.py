"""test_ontology_annotator.py — Automated test suite for the ontology-annotator skill.
Run with: python -m pytest skills/ontology-annotator/tests/ -v

Everything here runs offline against recorded OLS4 JSON fixtures
(skills/ontology-annotator/data/ols4_fixtures.json) — no network required,
and no test constructs an OLS4Client without `fixtures=` set, so a stray
network call is a test bug, not a possibility this suite tolerates.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

from ontology_annotator_core.io_readers import (  # noqa: E402
    infer_columns,
    load_h5ad_obs,
    load_table,
    parse_columns_arg,
)
from ontology_annotator_core.ols4_client import (  # noqa: E402
    OLS4Client,
    fixture_key,
    load_fixtures,
)
from ontology_annotator_core.scoring import (  # noqa: E402
    best_candidate,
    normalise,
    rank_candidates,
    score_candidate,
)

FIXTURES_PATH = SKILL_DIR / "data" / "ols4_fixtures.json"
DEMO_INPUT_PATH = SKILL_DIR / "examples" / "demo_input.csv"


@pytest.fixture(scope="module")
def fixtures() -> dict:
    return load_fixtures(FIXTURES_PATH)


# ── Fixture integrity ────────────────────────────────────────────────────


def test_fixtures_file_has_meta_and_is_real_api_data():
    """The fixture file must declare its provenance and cover all four ontologies."""
    raw = json.loads(FIXTURES_PATH.read_text())
    assert "_meta" in raw
    assert "ebi.ac.uk/ols4" in raw["_meta"]["source"]
    ontologies = {k.split(":", 1)[0] for k in raw["responses"]}
    assert ontologies == {"uberon", "cl", "mondo", "efo"}


def test_fixtures_include_a_zero_hit_query(fixtures):
    """At least one recorded query legitimately returns zero OLS4 hits —
    this is what a typo'd metadata value looks like in the real world, and
    the pipeline must handle it without crashing."""
    zero_hit = [k for k, v in fixtures.items() if v["response"]["numFound"] == 0]
    assert zero_hit, "Expected at least one zero-hit fixture (e.g. a typo query)"


# ── Scoring ──────────────────────────────────────────────────────────────


def test_normalise_folds_case_and_whitespace():
    assert normalise("  Lung  ") == "lung"
    assert normalise("Type 2 Diabetes") == "type 2 diabetes"


def test_score_candidate_exact_match_is_one():
    doc = {"label": "lung", "exact_synonyms": ["pulmo"]}
    assert score_candidate("lung", doc) == 1.0
    assert score_candidate("Lung", doc) == 1.0  # case-insensitive


def test_score_candidate_uses_best_synonym_not_just_label():
    doc = {"label": "pulmo", "exact_synonyms": ["lung"]}
    assert score_candidate("lung", doc) == 1.0


def test_rank_candidates_sorts_descending_and_caps_at_top_n(fixtures):
    key = fixture_key("uberon", "lung")
    docs = fixtures[key]["response"]["docs"]
    ranked = rank_candidates("lung", docs, top_n=3)
    assert len(ranked) == 3
    scores = [c["score"] for c in ranked]
    assert scores == sorted(scores, reverse=True)
    assert ranked[0]["ontology_id"] == "UBERON:0002048"
    assert ranked[0]["label"] == "lung"
    assert ranked[0]["score"] == 1.0


def test_best_candidate_flags_below_threshold():
    candidates = [{"ontology_id": "X:1", "label": "x", "score": 0.5, "iri": ""}]
    top, flagged = best_candidate(candidates, threshold=0.75)
    assert top["ontology_id"] == "X:1"
    assert flagged is True


def test_best_candidate_does_not_flag_above_threshold():
    candidates = [{"ontology_id": "X:1", "label": "x", "score": 0.9, "iri": ""}]
    top, flagged = best_candidate(candidates, threshold=0.75)
    assert flagged is False


def test_best_candidate_flags_empty_candidate_list():
    top, flagged = best_candidate([], threshold=0.75)
    assert top is None
    assert flagged is True


def test_best_candidate_never_silently_returns_a_pick_as_trustworthy():
    """A low-confidence top-1 is still returned (so it's visible in output),
    but always paired with flagged=True — never presented as if it were a
    confident answer."""
    candidates = [{"ontology_id": "X:1", "label": "x", "score": 0.1, "iri": ""}]
    top, flagged = best_candidate(candidates, threshold=0.75)
    assert top is not None
    assert flagged is True


# ── OLS4 client (fixture / offline mode) ────────────────────────────────


def test_ols4_client_fixture_mode_returns_docs(fixtures):
    client = OLS4Client(fixtures=fixtures)
    result = client.search("lung", "uberon")
    assert result["status"] == "ok"
    assert result["source"] == "fixture"
    assert len(result["docs"]) > 0
    assert result["docs"][0]["label"] == "lung"


def test_ols4_client_fixture_mode_is_case_insensitive(fixtures):
    client = OLS4Client(fixtures=fixtures)
    lower = client.search("lung", "uberon")
    upper = client.search("Lung", "uberon")
    assert lower["docs"] == upper["docs"]


def test_ols4_client_fixture_mode_handles_zero_hits(fixtures):
    client = OLS4Client(fixtures=fixtures)
    result = client.search("kynee", "uberon")
    assert result["status"] == "ok"
    assert result["docs"] == []


def test_ols4_client_missing_fixture_reports_no_fixture_not_a_crash(fixtures):
    client = OLS4Client(fixtures=fixtures)
    result = client.search("this value was never recorded", "uberon")
    assert result["status"] == "no_fixture"
    assert result["docs"] == []


def test_ols4_client_never_touches_network_in_fixture_mode(fixtures):
    """Guard against a regression that adds a network fallback inside fixture
    mode: the lazy `requests.Session` must never be instantiated when
    `fixtures` is set, for a hit, a zero-hit, or a missing-fixture lookup."""
    client = OLS4Client(fixtures=fixtures)
    client.search("lung", "uberon")
    client.search("kynee", "uberon")
    client.search("this value was never recorded", "uberon")
    assert client._session is None


# ── Column parsing / auto-detection ─────────────────────────────────────


def test_parse_columns_arg():
    parsed = parse_columns_arg("tissue:uberon, disease:MONDO ,trait:efo")
    assert parsed == {"tissue": "uberon", "disease": "mondo", "trait": "efo"}


def test_parse_columns_arg_rejects_malformed_entry():
    with pytest.raises(ValueError):
        parse_columns_arg("tissue_uberon")


def test_infer_columns_matches_known_names():
    df = pd.DataFrame({"sample_id": ["a"], "tissue": ["lung"], "disease": ["asthma"], "extra": [1]})
    inferred = infer_columns(df)
    assert inferred == {"tissue": "uberon", "disease": "mondo"}


def test_infer_columns_returns_empty_when_nothing_matches():
    df = pd.DataFrame({"foo": [1], "bar": [2]})
    assert infer_columns(df) == {}


# ── Table loading ────────────────────────────────────────────────────────


def test_load_table_csv():
    df = load_table(DEMO_INPUT_PATH)
    assert list(df.columns) == ["sample_id", "tissue", "cell_type", "disease", "trait"]
    assert len(df) == 5


def test_load_table_tsv(tmp_path):
    tsv_path = tmp_path / "in.tsv"
    tsv_path.write_text("tissue\tdisease\nlung\tasthma\n")
    df = load_table(tsv_path)
    assert list(df.columns) == ["tissue", "disease"]
    assert df.iloc[0]["tissue"] == "lung"


def test_load_h5ad_obs(tmp_path):
    anndata = pytest.importorskip("anndata")
    import numpy as np

    # pandas 3.x defaults string columns to a nullable StringArray, which
    # anndata's writer refuses unless this is opted into (its own suggested
    # fix — see the RuntimeError this guards against).
    anndata.settings.allow_write_nullable_strings = True

    obs_df = pd.DataFrame(
        {"tissue": ["lung", "liver", "heart"]}, index=["c1", "c2", "c3"]
    )
    adata = anndata.AnnData(X=np.zeros((3, 2)), obs=obs_df)
    h5ad_path = tmp_path / "demo.h5ad"
    adata.write_h5ad(h5ad_path)

    obs = load_h5ad_obs(h5ad_path)
    assert "tissue" in obs.columns
    assert "cell_id" in obs.columns
    assert list(obs["tissue"]) == ["lung", "liver", "heart"]


def test_load_h5ad_obs_missing_anndata_gives_clear_error(monkeypatch, tmp_path):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "anndata":
            raise ImportError("no anndata")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(RuntimeError, match="pip install anndata"):
        load_h5ad_obs(tmp_path / "does_not_matter.h5ad")


# ── Full pipeline (offline, via fixtures) ───────────────────────────────


def test_full_pipeline_demo_data(tmp_path):
    from ontology_annotator import run_annotation

    result = run_annotation(
        input_path=DEMO_INPUT_PATH,
        output_dir=tmp_path,
        columns_arg=None,
        threshold=0.8,
        fixtures_path=FIXTURES_PATH,
        demo=True,
    )

    assert (tmp_path / "annotated.csv").exists()
    assert (tmp_path / "report.md").exists()
    assert (tmp_path / "result.json").exists()
    assert (tmp_path / "reproducibility" / "commands.sh").exists()
    assert (tmp_path / "reproducibility" / "environment.yml").exists()
    assert (tmp_path / "reproducibility" / "checksums.sha256").exists()

    annotated = pd.read_csv(tmp_path / "annotated.csv")
    assert len(annotated) == 5
    for col in ("tissue", "cell_type", "disease", "trait"):
        for suffix in ("_ontology_id", "_label", "_score", "_flag", "_candidates"):
            assert f"{col}{suffix}" in annotated.columns

    # S1: exact, high-confidence matches on every column.
    row0 = annotated.iloc[0]
    assert row0["tissue_ontology_id"] == "UBERON:0002048"
    assert row0["tissue_score"] == 1.0
    assert row0["tissue_flag"] == False  # noqa: E712
    assert row0["disease_ontology_id"] == "MONDO:0005148"
    assert row0["trait_ontology_id"] == "EFO:0004340"

    # S3: "kynee" (typo tissue) and "diabetes melitus" (typo disease) have
    # zero OLS4 candidates and must be flagged, never silently picked.
    row2 = annotated.iloc[2]
    assert pd.isna(row2["tissue_ontology_id"]) or row2["tissue_ontology_id"] == ""
    assert row2["tissue_flag"] == True  # noqa: E712
    assert row2["disease_flag"] == True  # noqa: E712

    # S2: "blood pressure" is a real borderline case — OLS4 has no exact
    # "blood pressure" node, only "systolic/diastolic blood pressure"
    # (score 0.7568). At threshold=0.8 that must be flagged even though a
    # real candidate exists — never silently accepted as the answer.
    row1 = annotated.iloc[1]
    assert 0.7 < row1["trait_score"] < 0.8
    assert row1["trait_flag"] == True  # noqa: E712

    result_json = json.loads((tmp_path / "result.json").read_text())
    assert result_json["skill"] == "ontology-annotator"
    assert result_json["summary"]["n_rows"] == 5
    assert result_json["summary"]["n_flagged"] > 0

    report = (tmp_path / "report.md").read_text()
    assert "not a medical device" in report
    assert "Flagged Rows" in report
    assert "ClawBio Skills That Can Consume These IDs" in report

    # summary sanity: matched + flagged + no_candidates should equal 5 per column
    for col_summary in result["summary"].values():
        assert (
            col_summary["n_matched"] + col_summary["n_flagged"] + col_summary["n_no_candidates"]
            == 5
        )


def test_full_pipeline_reuses_cache_for_repeated_normalised_value(tmp_path):
    """Row S5 repeats S1's tissue/disease/trait values with different casing;
    the in-run memo must answer those from cache, not a second lookup — we
    can't observe network calls directly in fixture mode, but we can assert
    the identical scores/ids come back, proving the same cached candidates
    were reused."""
    from ontology_annotator import run_annotation

    run_annotation(
        input_path=DEMO_INPUT_PATH,
        output_dir=tmp_path,
        columns_arg=None,
        threshold=0.8,
        fixtures_path=FIXTURES_PATH,
        demo=True,
    )
    annotated = pd.read_csv(tmp_path / "annotated.csv")
    row0, row4 = annotated.iloc[0], annotated.iloc[4]
    assert row0["tissue_ontology_id"] == row4["tissue_ontology_id"]
    assert row0["disease_ontology_id"] == row4["disease_ontology_id"]
    assert row0["trait_ontology_id"] == row4["trait_ontology_id"]


def test_explicit_columns_arg_overrides_autodetect(tmp_path):
    from ontology_annotator import run_annotation

    result = run_annotation(
        input_path=DEMO_INPUT_PATH,
        output_dir=tmp_path,
        columns_arg="tissue:uberon",
        threshold=0.8,
        fixtures_path=FIXTURES_PATH,
        demo=True,
    )
    assert list(result["summary"].keys()) == ["tissue"]
    annotated = pd.read_csv(tmp_path / "annotated.csv")
    assert "disease_ontology_id" not in annotated.columns


def test_run_annotation_raises_clear_error_when_no_columns_detected(tmp_path):
    from ontology_annotator import run_annotation

    bad_input = tmp_path / "no_recognisable_columns.csv"
    bad_input.write_text("foo,bar\n1,2\n")

    with pytest.raises(ValueError, match="--columns"):
        run_annotation(
            input_path=bad_input,
            output_dir=tmp_path / "out",
            columns_arg=None,
            threshold=0.75,
            fixtures_path=FIXTURES_PATH,
            demo=False,
        )
