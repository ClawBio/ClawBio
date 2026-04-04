"""Tests for the swappable fine-mapping benchmark pipeline."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from finemapping_benchmark import (
    make_benchmark_locus,
    run_method_abf,
    run_method_susie,
    score_method,
    run_benchmark,
    CAUSAL_INDICES,
    METHODS,
)


class TestBenchmarkLocus:
    """Test synthetic locus generation."""

    def test_locus_shape(self):
        df, ld, causal = make_benchmark_locus()
        assert len(df) == 200
        assert ld.shape == (200, 200)
        assert causal == [60, 140]

    def test_locus_has_required_columns(self):
        df, _, _ = make_benchmark_locus()
        for col in ["rsid", "chr", "pos", "beta", "se", "pvalue", "z", "n"]:
            assert col in df.columns, f"Missing column: {col}"

    def test_ld_is_symmetric(self):
        _, ld, _ = make_benchmark_locus()
        np.testing.assert_array_almost_equal(ld, ld.T, decimal=10)

    def test_ld_diagonal_is_one(self):
        _, ld, _ = make_benchmark_locus()
        np.testing.assert_array_almost_equal(np.diag(ld), np.ones(200), decimal=5)

    def test_causal_variants_have_low_pvalues(self):
        df, _, causal = make_benchmark_locus()
        for c in causal:
            assert df.iloc[c]["pvalue"] < 0.01, f"Causal variant {c} p={df.iloc[c]['pvalue']}"

    def test_deterministic_with_same_seed(self):
        df1, ld1, _ = make_benchmark_locus(seed=99)
        df2, ld2, _ = make_benchmark_locus(seed=99)
        np.testing.assert_array_equal(df1["beta"].values, df2["beta"].values)

    def test_different_with_different_seed(self):
        df1, _, _ = make_benchmark_locus(seed=1)
        df2, _, _ = make_benchmark_locus(seed=2)
        assert not np.allclose(df1["beta"].values, df2["beta"].values)


class TestMethodRunners:
    """Test individual fine-mapping methods."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.df, self.ld, self.causal = make_benchmark_locus()

    def test_abf_returns_pips(self):
        result = run_method_abf(self.df, self.ld)
        assert result["method"] == "ABF"
        assert len(result["pips"]) == 200
        assert all(0 <= p <= 1 for p in result["pips"])

    def test_abf_pips_sum_approximately_to_one(self):
        result = run_method_abf(self.df, self.ld)
        pip_sum = sum(result["pips"])
        # ABF PIPs should sum to approximately 1 (or slightly more due to multiple signals)
        assert 0.5 < pip_sum < 5.0, f"PIP sum {pip_sum} outside expected range"

    def test_abf_has_credible_set(self):
        result = run_method_abf(self.df, self.ld)
        assert "credible_set_indices" in result
        assert result["credible_set_size"] > 0
        assert result["coverage"] >= 0.95

    def test_susie_returns_pips(self):
        result = run_method_susie(self.df, self.ld)
        assert result["method"] == "SuSiE"
        assert len(result["pips"]) == 200
        assert all(0 <= p <= 1 for p in result["pips"])

    def test_susie_finds_credible_sets(self):
        result = run_method_susie(self.df, self.ld)
        assert "credible_sets" in result
        # SuSiE should find at least 1 credible set on this well-powered locus
        assert result["n_credible_sets"] >= 1

    def test_both_methods_find_causal_signal(self):
        """Both methods should rank at least one causal variant highly."""
        for method_name, method_fn in METHODS.items():
            result = method_fn(self.df, self.ld)
            pips = np.array(result["pips"])
            # At least one causal variant should be in top 20 by PIP
            top20 = set(np.argsort(-pips)[:20])
            causal_in_top20 = sum(1 for c in self.causal if c in top20)
            assert causal_in_top20 >= 1, f"{method_name} failed to rank any causal in top 20"


class TestScoring:
    """Test the scoring logic."""

    def test_perfect_score(self):
        """A method that puts all PIP on causal variants."""
        pips = [0.0] * 200
        pips[60] = 0.5
        pips[140] = 0.5
        result = {"method": "perfect", "pips": pips, "credible_set_indices": [60, 140],
                  "credible_set_size": 2, "elapsed": 0}
        scored = score_method(result, [60, 140], 200)
        assert scored["recall"] == 1.0
        assert scored["precision"] == 1.0
        assert scored["composite_score"] > 0.8

    def test_zero_score(self):
        """A method that misses all causal variants."""
        pips = [1.0 / 200] * 200
        result = {"method": "random", "pips": pips, "credible_set_indices": [10, 20, 30],
                  "credible_set_size": 3, "elapsed": 0}
        scored = score_method(result, [60, 140], 200)
        assert scored["recall"] == 0.0
        assert scored["composite_score"] < 0.3

    def test_partial_capture(self):
        """A method that captures one of two causal variants."""
        pips = [0.0] * 200
        pips[60] = 0.8
        pips[10] = 0.2
        result = {"method": "partial", "pips": pips, "credible_set_indices": [60, 10],
                  "credible_set_size": 2, "elapsed": 0}
        scored = score_method(result, [60, 140], 200)
        assert scored["recall"] == 0.5
        assert scored["precision"] == 0.5


class TestBenchmarkRunner:
    """Test the full benchmark pipeline."""

    def test_full_benchmark_runs(self):
        result = run_benchmark(methods=["abf", "susie"], seed=42)
        assert "methods" in result
        assert len(result["methods"]) == 2
        assert result["winner"] in ("ABF", "SuSiE")
        assert result["winner_score"] > 0

    def test_benchmark_with_output(self, tmp_path):
        result = run_benchmark(methods=["abf"], seed=42, output_dir=tmp_path)
        assert (tmp_path / "finemapping_benchmark.json").exists()

    def test_both_methods_produce_valid_scores(self):
        result = run_benchmark(methods=["abf", "susie"], seed=42)
        for m in result["methods"]:
            assert "composite_score" in m
            assert m["composite_score"] > 0
            assert m["recall"] >= 0
            assert m["precision"] >= 0


class TestSuSiEInfMethod:
    """Test SuSiE-inf as a swappable fine-mapping method in the benchmark."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.df, self.ld, self.causal = make_benchmark_locus()

    def test_susieinf_registered_in_methods(self):
        """susieinf is registered as a runnable method."""
        assert "susieinf" in METHODS

    def test_susieinf_returns_required_keys(self):
        """run_method_susieinf returns dict with pips, credible_set_indices, method."""
        from finemapping_benchmark import run_method_susieinf
        result = run_method_susieinf(self.df, self.ld)
        for key in ("method", "pips", "credible_set_indices", "credible_set_size", "elapsed"):
            assert key in result, f"Missing key: {key}"

    def test_susieinf_method_label(self):
        """run_method_susieinf sets method to 'SuSiE-inf'."""
        from finemapping_benchmark import run_method_susieinf
        result = run_method_susieinf(self.df, self.ld)
        assert result["method"] == "SuSiE-inf"

    def test_susieinf_pips_length(self):
        """run_method_susieinf returns one PIP per variant."""
        from finemapping_benchmark import run_method_susieinf
        result = run_method_susieinf(self.df, self.ld)
        assert len(result["pips"]) == len(self.df)

    def test_susieinf_pips_in_range(self):
        """All SuSiE-inf PIPs are in [0, 1]."""
        from finemapping_benchmark import run_method_susieinf
        result = run_method_susieinf(self.df, self.ld)
        assert all(0.0 <= p <= 1.0 for p in result["pips"])

    def test_susieinf_finds_causal_signal_in_top20(self):
        """SuSiE-inf ranks at least one causal variant in the top 20 by PIP."""
        from finemapping_benchmark import run_method_susieinf
        result = run_method_susieinf(self.df, self.ld)
        top20 = set(np.argsort(-np.array(result["pips"]))[:20])
        assert any(c in top20 for c in self.causal)

    def test_susieinf_benchmark_produces_valid_score(self):
        """Full benchmark including susieinf produces a valid composite score."""
        result = run_benchmark(methods=["susieinf"], seed=42)
        assert len(result["methods"]) == 1
        scored = result["methods"][0]
        assert "error" not in scored
        assert scored["composite_score"] > 0
        assert scored["recall"] >= 0
        assert scored["precision"] >= 0

    def test_three_way_benchmark_picks_winner(self):
        """Benchmark across abf, susie, susieinf completes and picks a winner."""
        result = run_benchmark(methods=["abf", "susie", "susieinf"], seed=42)
        assert result["winner"] in ("ABF", "SuSiE", "SuSiE-inf")
        assert result["winner_score"] > 0
