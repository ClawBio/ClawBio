"""Tests for the SuSiE fine-mapping skill."""
import importlib.util
import sys
import warnings
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from inline_snapshot import snapshot

from clawbio.common import demo_runner

# Add skill directory to path
SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

import fine_mapping
from fine_mapping_core.abf import compute_abf, _log_abf, DEFAULT_W
from fine_mapping_core.susie import run_susie
from fine_mapping_core.credible_sets import (
    build_credible_sets_susie,
    build_credible_set_abf,
    _greedy_credible_set,
    _purity,
)
from fine_mapping_core.io import load_sumstats, load_ld
from fine_mapping_core.report import generate_markdown

needs_jax = pytest.mark.skipif(importlib.util.find_spec("jax") is None, reason="needs jax (uv sync --extra fine-mapping)")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _small_locus(n=20, seed=0) -> pd.DataFrame:
    """Synthetic 20-variant locus with one strong signal at index 10."""
    rng = np.random.RandomState(seed)
    z = rng.normal(0, 1, n)
    z[10] = 6.0  # injected signal
    return pd.DataFrame({
        "rsid": [f"rs{i}" for i in range(n)],
        "chr": "1",
        "pos": np.arange(1_000_000, 1_000_000 + n * 1000, 1000),
        "z": z,
        "se": np.full(n, 0.01),
        "n": np.full(n, 5000),
        "p": 2 * (1 - 0.5 * (1 + np.sign(z) * (1 - np.exp(-0.7 * z**2)))),
    })


def _identity_ld(n: int) -> np.ndarray:
    """Identity LD matrix (no LD — each variant independent)."""
    return np.eye(n)


# ---------------------------------------------------------------------------
# TestDemoData
# ---------------------------------------------------------------------------


class TestDemoData:
    """Tests for the synthetic demo locus generator."""

    def test_shape(self):
        """make_demo_data returns 200 variants and a 200x200 LD matrix."""
        df, R = fine_mapping.make_demo_data()
        assert len(df) == 200
        assert R.shape == (200, 200)

    def test_columns(self):
        """Demo DataFrame has required columns."""
        df, _ = fine_mapping.make_demo_data()
        for col in ("rsid", "chr", "pos", "z", "se", "n", "p"):
            assert col in df.columns, f"Missing column: {col}"

    def test_reproducible(self):
        """Same seed produces identical z-scores."""
        df1, _ = fine_mapping.make_demo_data(seed=42)
        df2, _ = fine_mapping.make_demo_data(seed=42)
        np.testing.assert_array_equal(df1["z"].values, df2["z"].values)

    def test_different_seeds(self):
        """Different seeds produce different z-scores."""
        df1, _ = fine_mapping.make_demo_data(seed=42)
        df2, _ = fine_mapping.make_demo_data(seed=99)
        assert not np.allclose(df1["z"].values, df2["z"].values)

    def test_ld_is_symmetric(self):
        """LD matrix is symmetric and has unit diagonal."""
        _, R = fine_mapping.make_demo_data()
        np.testing.assert_allclose(R, R.T, atol=1e-10)
        np.testing.assert_allclose(np.diag(R), 1.0, atol=1e-10)

    def test_ld_is_psd(self):
        """LD matrix is positive semi-definite (all eigenvalues >= 0)."""
        _, R = fine_mapping.make_demo_data()
        eigvals = np.linalg.eigvalsh(R)
        assert eigvals.min() >= -1e-6


# ---------------------------------------------------------------------------
# TestABF
# ---------------------------------------------------------------------------


class TestABF:
    """Tests for Approximate Bayes Factor computation."""

    def test_log_abf_shape(self):
        """_log_abf returns an array matching input length."""
        z = np.array([0.0, 1.0, 3.0, 6.0])
        V = np.full(4, 1.0 / 5000)
        result = _log_abf(z, V, DEFAULT_W)
        assert result.shape == (4,)

    def test_log_abf_monotone_in_z(self):
        """Larger |z| produces larger log ABF."""
        z = np.array([1.0, 2.0, 4.0, 6.0])
        V = np.full(4, 1.0 / 5000)
        log_abf = _log_abf(z, V, DEFAULT_W)
        assert np.all(np.diff(log_abf) > 0)

    def test_compute_abf_sums_to_one(self):
        """ABF PIPs sum to 1.0 (they are a proper probability distribution)."""
        df = _small_locus()
        pip = compute_abf(df)
        assert pytest.approx(pip.sum(), abs=1e-9) == 1.0

    def test_compute_abf_signal_has_highest_pip(self):
        """Injected signal at index 10 receives the highest PIP."""
        df = _small_locus()
        pip = compute_abf(df)
        assert pip.argmax() == 10

    def test_compute_abf_uses_se_column(self):
        """When 'se' column is present, it is used for V rather than n."""
        df = _small_locus()
        df_no_n = df.drop(columns=["n"])
        pip = compute_abf(df_no_n)
        assert pytest.approx(pip.sum(), abs=1e-9) == 1.0

    def test_compute_abf_fallback_no_se_no_n(self):
        """Falls back to n_eff=10000 when neither se nor n are present."""
        df = _small_locus()[["rsid", "chr", "pos", "z"]]
        pip = compute_abf(df)
        assert pytest.approx(pip.sum(), abs=1e-9) == 1.0


# ---------------------------------------------------------------------------
# TestSuSiE
# ---------------------------------------------------------------------------


class TestSuSiE:
    """Tests for the sushie-backed SuSiE engine adapter."""

    @needs_jax
    def test_pip_shape(self):
        """run_susie returns PIPs with correct shape."""
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie(z=df["z"].values, R=R, n=5000, L=5)
        assert result["pip"].shape == (20,)

    @needs_jax
    def test_pip_range(self):
        """All PIPs are in [0, 1]."""
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie(z=df["z"].values, R=R, n=5000, L=5)
        assert result["pip"].min() >= 0.0
        assert result["pip"].max() <= 1.0

    @needs_jax
    def test_signal_recovers_high_pip(self):
        """The injected signal at index 10 has the highest PIP."""
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie(z=df["z"].values, R=R, n=5000, L=5)
        assert result["pip"].argmax() == 10

    @needs_jax
    def test_alpha_contains_only_active_signals(self):
        """alpha keeps only signals sushie retained as credible sets.

        The small locus has one injected signal, so with L=3 requested only
        one alpha row survives pruning. Each surviving row is a probability
        distribution over the p variants.
        """
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie(z=df["z"].values, R=R, n=5000, L=3)
        assert result["alpha"].shape == (1, 20)
        np.testing.assert_allclose(result["alpha"].sum(axis=1), 1.0, atol=1e-4)

    @needs_jax
    def test_two_signal_demo_keeps_two_alpha_rows(self):
        """Demo locus (two causal variants) keeps exactly two active signals."""
        df, R = fine_mapping.make_demo_data(seed=42)
        result = run_susie(z=df["z"].values, R=R, n=5000, L=10)
        assert result["alpha"].shape == (2, 200)

    @needs_jax
    def test_engine_is_sushie(self):
        """The engine tag identifies the sushie implementation."""
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie(z=df["z"].values, R=R, n=5000, L=3)
        assert result["engine"] == "sushie"

    @needs_jax
    def test_converges_on_clean_signal(self):
        """The engine converges for a clean single signal and reports n_iter."""
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie(z=df["z"].values, R=R, n=5000, L=5)
        assert result["converged"] is True
        assert result["n_iter"] >= 1

    @needs_jax
    def test_elbo_tracked(self):
        """ELBO history is recorded and contains finite values."""
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie(z=df["z"].values, R=R, n=5000, L=5)
        elbo = result["elbo"]
        assert len(elbo) >= 1
        assert all(np.isfinite(e) for e in elbo)

    @needs_jax
    def test_two_signals_recovered(self):
        """Demo locus with two injected causal variants: BOTH are recovered.

        `60 in top5 or 140 in top5` would pass while the engine missed one of
        the two signals entirely, which is the exact failure a fine-mapper is
        supposed to be judged on. Require both, each as its own credible set.
        """
        df, R = fine_mapping.make_demo_data(seed=42)
        result = run_susie(z=df["z"].values, R=R, n=5000, L=10)
        top5 = set(int(i) for i in np.argsort(-result["pip"])[:5])
        assert {60, 140} <= top5, (
            f"causal indices 60 and 140 not both in the top-5 PIPs: {sorted(top5)}"
        )
        assert result["pip"][60] > 0.8 and result["pip"][140] > 0.8
        leads = sorted(int(np.argmax(row)) for row in result["alpha"])
        assert leads == [60, 140], f"credible-set leads were {leads}"

    @needs_jax
    def test_null_locus_no_phantom_pip(self):
        """Null locus (z=0 everywhere) produces near-zero PIPs."""
        p = 50
        z = np.zeros(p)
        R = np.eye(p)
        result = run_susie(z=z, R=R, n=10000, L=5)
        assert result["pip"].max() < 0.1, (
            f"Null locus phantom PIP: max PIP = {result['pip'].max():.3f}, "
            f"expected < 0.1"
        )
        # No signal survives pruning, so no alpha rows remain
        assert result["alpha"].shape == (0, p)

    @needs_jax
    def test_single_strong_signal_high_pip(self):
        """A lone strong signal receives PIP > 0.8."""
        p = 50
        z = np.zeros(p)
        z[25] = 5.0  # strong signal at index 25
        R = np.eye(p)
        result = run_susie(z=z, R=R, n=10000, L=5)
        assert result["pip"].argmax() == 25
        assert result["pip"][25] > 0.8, (
            f"Signal PIP = {result['pip'][25]:.3f}, expected > 0.8"
        )

    def test_nan_z_raises(self):
        """NaN z-scores are rejected (sushie would silently return garbage)."""
        z = np.zeros(20)
        z[3] = np.nan
        with pytest.raises(ValueError, match="NaN"):
            run_susie(z=z, R=_identity_ld(20), n=5000, L=3)

    def test_nonpositive_n_raises(self):
        """Sample size must be positive."""
        df = _small_locus(n=20)
        with pytest.raises(ValueError, match="positive"):
            run_susie(z=df["z"].values, R=_identity_ld(20), n=0, L=3)

    @needs_jax
    def test_deterministic(self):
        """Two runs on identical input give identical PIPs."""
        df = _small_locus(n=20)
        R = _identity_ld(20)
        r1 = run_susie(z=df["z"].values, R=R, n=5000, L=3)
        r2 = run_susie(z=df["z"].values, R=R, n=5000, L=3)
        np.testing.assert_array_equal(r1["pip"], r2["pip"])


# ---------------------------------------------------------------------------
# TestSuSiEBenchmarkContract
# ---------------------------------------------------------------------------


class TestSuSiEBenchmarkContract:
    """run_susie must stay callable by the external clawbio_bench harness.

    The harness (biostochastics/clawbio_bench, drivers/finemapping_driver.py)
    imports core.susie.run_susie and always calls it with the full keyword
    set below — `w` is passed even when a test case's inputs.json omits it,
    because the driver defaults it to 0.04. It then reads `mu` and `mu2` off
    the result. Dropping any of these is scored as a skill failure, not as an
    interface change, so the adapter keeps them.
    """

    HARNESS_KWARGS = dict(w=0.04, max_iter=100, tol=1e-3, min_purity=0.5)

    @needs_jax
    def test_accepts_harness_keyword_set(self):
        """The exact driver call shape does not raise TypeError."""
        df = _small_locus(n=20)
        result = run_susie(
            z=df["z"].values, R=_identity_ld(20), n=5000, L=5,
            **self.HARNESS_KWARGS,
        )
        assert result["pip"].shape == (20,)

    @needs_jax
    def test_returns_posterior_moments(self):
        """mu (post mean) and mu2 (post second moment) are exposed per signal."""
        df = _small_locus(n=20)
        result = run_susie(z=df["z"].values, R=_identity_ld(20), n=5000, L=5)
        n_signals = result["alpha"].shape[0]
        assert result["mu"].shape == (n_signals, 20)
        assert result["mu2"].shape == (n_signals, 20)

    @needs_jax
    def test_mu2_is_second_moment_not_variance(self):
        """mu2 >= mu**2 elementwise — it is E[b^2], not Var(b).

        clawbio_bench case fm_14_susie_moment_mislabeled scores exactly this
        confusion.
        """
        df = _small_locus(n=20)
        result = run_susie(z=df["z"].values, R=_identity_ld(20), n=5000, L=5)
        assert np.all(result["mu2"] >= result["mu"] ** 2 - 1e-9)

    @needs_jax
    def test_moments_align_with_alpha_rows(self):
        """mu/mu2 are pruned to the same kept signals as alpha."""
        p = 50
        z = np.zeros(p)
        result = run_susie(z=z, R=np.eye(p), n=10000, L=5)
        assert result["alpha"].shape == (0, p)
        assert result["mu"].shape == (0, p)
        assert result["mu2"].shape == (0, p)

    @needs_jax
    def test_min_purity_forwarded_to_engine(self, monkeypatch):
        """min_purity sets sushie's own purity pruning threshold.

        Without forwarding, the CLI's --min-purity could only tighten the
        engine's hardcoded 0.5 downstream, never loosen it.
        """
        from sushie.infer_ss import infer_sushie_ss as real_infer

        seen = {}

        def _spy(*args, **kwargs):
            seen.update(kwargs)
            return real_infer(*args, **kwargs)

        monkeypatch.setattr("sushie.infer_ss.infer_sushie_ss", _spy)
        df = _small_locus(n=20)
        run_susie(z=df["z"].values, R=_identity_ld(20), n=5000, L=3,
                  min_purity=0.8)
        assert seen["purity"] == 0.8

    @needs_jax
    def test_prior_variance_seeds_effect_var(self, monkeypatch):
        """w seeds sushie's prior effect variance rather than being discarded."""
        from sushie.infer_ss import infer_sushie_ss as real_infer

        seen = {}

        def _spy(*args, **kwargs):
            seen.update(kwargs)
            return real_infer(*args, **kwargs)

        monkeypatch.setattr("sushie.infer_ss.infer_sushie_ss", _spy)
        df = _small_locus(n=20)
        run_susie(z=df["z"].values, R=_identity_ld(20), n=5000, L=3, w=0.09)
        assert seen["effect_var"] is not None
        assert float(np.ravel(seen["effect_var"])[0]) == pytest.approx(0.09)

    @needs_jax
    def test_nonconvergence_warns_and_flags(self):
        """Hitting max_iter warns (as susieR does) and reports converged=False.

        PIPs are still returned finite, so the flag is the consumer's cue;
        clawbio_bench fm_12 scores the converged/pip pair and records the
        warning text in its payload.
        """
        df = _small_locus(n=20)
        with pytest.warns(RuntimeWarning, match="did not converge"):
            result = run_susie(z=df["z"].values, R=_identity_ld(20), n=5000,
                               L=5, max_iter=1, tol=1e-3)
        assert result["converged"] is False
        assert result["n_iter"] == 1
        assert np.all(np.isfinite(result["pip"]))

    @needs_jax
    def test_coverage_forwarded_to_engine(self, monkeypatch):
        """coverage sets sushie's credible-set threshold at fit time.

        Otherwise sushie prunes at its hardcoded 0.95 while credible_sets.py
        builds sets at the user's --coverage, and the two disagree.
        """
        from sushie.infer_ss import infer_sushie_ss as real_infer

        seen = {}

        def _spy(*args, **kwargs):
            seen.update(kwargs)
            return real_infer(*args, **kwargs)

        monkeypatch.setattr("sushie.infer_ss.infer_sushie_ss", _spy)
        df = _small_locus(n=20)
        run_susie(z=df["z"].values, R=_identity_ld(20), n=5000, L=3,
                  coverage=0.9)
        assert seen["threshold"] == 0.9

    @needs_jax
    def test_locus_smaller_than_L_still_runs(self):
        """A 9-variant locus under the CLI default L=10 must not error.

        sushie rejects a fit whose min_snps guard is below L, so passing L
        through unchanged turned loci that the old hand-rolled IBSS handled
        into hard failures. run_susie clamps L to p and says so.
        """
        p_small = 9
        z = np.zeros(p_small)
        z[3] = 6.0
        with pytest.warns(RuntimeWarning, match="clamping to L=9"):
            result = run_susie(z=z, R=np.eye(p_small), n=5000, L=10)
        assert result["pip"].shape == (p_small,)
        assert int(result["pip"].argmax()) == 3

    @needs_jax
    def test_L_at_or_below_p_is_not_clamped(self):
        """The clamp is a guard, not a rewrite: L <= p passes through silently."""
        df = _small_locus(n=20)
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            result = run_susie(z=df["z"].values, R=_identity_ld(20), n=5000, L=5)
        assert result["pip"].shape == (20,)

    @needs_jax
    def test_x64_flag_is_restored(self):
        """jax.config is process-global; run_susie must not leave it flipped.

        Any other JAX user in the same interpreter would otherwise silently
        change precision the first time a locus is fine-mapped.
        """
        import jax

        before = jax.config.read("jax_enable_x64")
        try:
            jax.config.update("jax_enable_x64", False)
            df = _small_locus(n=20)
            run_susie(z=df["z"].values, R=_identity_ld(20), n=5000, L=3)
            assert jax.config.read("jax_enable_x64") is False
        finally:
            jax.config.update("jax_enable_x64", before)

    @needs_jax
    def test_max_iter_echoed_back(self):
        """The report needs the budget, not just the iterations used."""
        df = _small_locus(n=20)
        result = run_susie(z=df["z"].values, R=_identity_ld(20), n=5000, L=3,
                           max_iter=123)
        assert result["max_iter"] == 123

    @pytest.mark.parametrize("bad", [0.0, 1.0, 1.5])
    def test_coverage_out_of_open_interval_raises(self, bad):
        """sushie requires 0 < threshold < 1; reject before reaching it."""
        df = _small_locus(n=20)
        with pytest.raises(ValueError, match="coverage"):
            run_susie(z=df["z"].values, R=_identity_ld(20), n=5000, L=3,
                      coverage=bad)

# ---------------------------------------------------------------------------
# TestSuSiEKnownAnswer
# ---------------------------------------------------------------------------


_KNOWN_ANSWER_PATH = (
    Path(__file__).resolve().parent / "fixtures" / "demo_susie_known_answer.json"
)


@needs_jax
class TestSuSiEKnownAnswer:
    """Pin the engine's actual numbers, not just their shapes.

    Every other SuSiE test here asserts a shape, an argmax or a `> 0.8`
    threshold. All of them keep passing if a future sushie release shifts PIP
    magnitudes, reorders alpha rows against `CSIndex`, or rescales mu/mu2 --
    exactly the drift an engine swap makes possible. These assertions compare
    against values recorded from the seeded demo locus at the commit that
    introduced the sushie backend (see the fixture's own `_comment`).
    """

    @pytest.fixture(scope="class")
    def known(self):
        with open(_KNOWN_ANSWER_PATH) as fh:
            return json.load(fh)

    @pytest.fixture(scope="class")
    def fit(self):
        df, R = fine_mapping.make_demo_data(seed=42)
        return df, R, run_susie(z=df["z"].values, R=R, n=5000, L=10)

    def test_engine_version_matches_recording(self, known, fit):
        """A sushie bump must re-record the fixture, not silently drift past it."""
        _, _, result = fit
        assert result["engine_version"] == known["recorded_from"]["engine_version"], (
            "sushie version changed since these numbers were recorded; re-record "
            f"{_KNOWN_ANSWER_PATH.name} in the same commit as the bump"
        )

    def test_convergence_matches_recording(self, known, fit):
        _, _, result = fit
        assert result["converged"] is known["converged"]
        assert result["n_iter"] == known["n_iter"]

    def test_credible_set_count_and_leads(self, known, fit):
        """Exactly the recorded number of sets, led by the recorded variants.

        This is the assertion the phantom-credible-set regression would trip:
        the hand-rolled IBSS reported 10 sets here, 8 of them purity-0 phantoms.
        """
        _, _, result = fit
        assert result["alpha"].shape[0] == known["n_credible_sets"]
        leads = [int(np.argmax(row)) for row in result["alpha"]]
        assert leads == known["alpha_row_leads"]

    def test_pips_match_recorded_values(self, known, fit):
        """PIP magnitudes, not just ordering."""
        _, _, result = fit
        atol = known["tolerance"]["pip_atol"]
        for idx, want in known["pip_by_index"].items():
            got = float(result["pip"][int(idx)])
            assert got == pytest.approx(want, abs=atol), (
                f"PIP at index {idx}: recorded {want}, got {got}"
            )
        assert int((result["pip"] >= 0.1).sum()) == known["n_pip_at_least_0.1"]

    def test_lead_moments_match_recorded_scale(self, known, fit):
        """mu/mu2 stay on the scale this adapter documents.

        clawbio_bench fm_14 reads `mu` expecting susieR's z-unit `r * z`; sushie
        reports standardised effect-size moments instead (mu at index 60 is
        ~0.22 against an injected beta of 0.25). If a sushie release changes
        that scale, the docstring in susie.py and the bench conversation both
        need revisiting -- so fail here rather than let it pass unnoticed.
        """
        _, _, result = fit
        rtol = known["tolerance"]["moment_rtol"]
        for rec in known["lead_moments"]:
            k, i = rec["alpha_row"], rec["variant_index"]
            assert float(result["alpha"][k][i]) == pytest.approx(rec["alpha"], rel=rtol)
            assert float(result["mu"][k][i]) == pytest.approx(rec["mu"], rel=rtol)
            assert float(result["mu2"][k][i]) == pytest.approx(rec["mu2"], rel=rtol)

    def test_credible_set_membership_matches_recording(self, known, fit):
        """CSIndex k -> alpha row k-1 is a mapping no other test pins.

        `build_credible_sets_susie` consumes the rows in order, so an ordering
        change inside sushie would relabel every set without changing a shape.
        """
        df, R, result = fit
        df = df.copy()
        df["pip"] = result["pip"]
        cs_list = build_credible_sets_susie(alpha=result["alpha"], df=df, R=R)
        assert len(cs_list) == known["n_credible_sets"]
        want_leads = known["alpha_row_leads"]
        for cs, want_idx in zip(cs_list, want_leads):
            assert cs["size"] == 1, f"{cs['cs_id']} grew to {cs['size']} variants"
            assert cs["lead_rsid"] == df["rsid"].iloc[want_idx]
            assert cs["purity"] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# TestCredibleSets
# ---------------------------------------------------------------------------


class TestCredibleSets:
    """Tests for credible set construction."""

    def test_greedy_coverage_reached(self):
        """_greedy_credible_set accumulates enough weight to meet coverage."""
        weights = np.array([0.1, 0.3, 0.05, 0.4, 0.15])
        cs = _greedy_credible_set(weights, coverage=0.95)
        assert sum(weights[i] for i in cs) >= 0.95

    def test_greedy_single_dominant(self):
        """When one variant has weight > coverage, CS is size 1."""
        weights = np.zeros(10)
        weights[3] = 0.99
        cs = _greedy_credible_set(weights, coverage=0.95)
        assert cs == [3]

    def test_purity_single_variant(self):
        """Purity of a single-variant CS is 1.0 by definition."""
        R = np.eye(5)
        assert _purity([2], R) == 1.0

    def test_purity_identity_ld(self):
        """Purity of any multi-variant CS under identity LD is 0.0."""
        R = np.eye(5)
        assert _purity([0, 1, 2], R) == pytest.approx(0.0)

    def test_purity_perfect_ld(self):
        """Purity of multi-variant CS under perfect LD is 1.0."""
        R = np.ones((5, 5))
        np.fill_diagonal(R, 1.0)
        assert _purity([0, 1, 2], R) == pytest.approx(1.0)

    def test_abf_credible_set_structure(self):
        """build_credible_set_abf returns a list with one dict with expected keys."""
        df = _small_locus()
        df["pip"] = compute_abf(df)
        cs_list = build_credible_set_abf(pip=df["pip"].values, df=df)
        assert len(cs_list) == 1
        cs = cs_list[0]
        for key in ("cs_id", "size", "coverage", "lead_rsid", "variants"):
            assert key in cs, f"Missing key: {key}"

    def test_abf_credible_set_coverage(self):
        """ABF credible set coverage meets the requested threshold."""
        df = _small_locus()
        df["pip"] = compute_abf(df)
        cs_list = build_credible_set_abf(pip=df["pip"].values, df=df, coverage=0.95)
        assert cs_list[0]["coverage"] >= 0.95

    @needs_jax
    def test_susie_credible_sets_structure(self):
        """build_credible_sets_susie returns list of dicts with expected keys."""
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie(z=df["z"].values, R=R, n=5000, L=3)
        df["pip"] = result["pip"]
        cs_list = build_credible_sets_susie(alpha=result["alpha"], df=df, R=R)
        for cs in cs_list:
            for key in ("cs_id", "signal_index", "size", "coverage", "lead_rsid", "variants"):
                assert key in cs, f"Missing key: {key}"

    @needs_jax
    def test_susie_purity_flag(self):
        """Sets with purity < min_purity are flagged pure=False."""
        df = _small_locus(n=20)
        R = _identity_ld(20)  # all variants independent → purity = 0
        result = run_susie(z=df["z"].values, R=R, n=5000, L=3)
        df["pip"] = result["pip"]
        cs_list = build_credible_sets_susie(
            alpha=result["alpha"], df=df, R=R, min_purity=0.9
        )
        # Under identity LD, all multi-variant CS should be impure
        multi = [cs for cs in cs_list if cs["size"] > 1]
        assert all(not cs["pure"] for cs in multi)


# ---------------------------------------------------------------------------
# TestIO
# ---------------------------------------------------------------------------


class TestIO:
    """Tests for sumstats and LD matrix loading."""

    def test_load_sumstats_z_column(self, tmp_path):
        """load_sumstats accepts a file with a direct 'z' column."""
        f = tmp_path / "stats.tsv"
        f.write_text("rsid\tchr\tpos\tz\tse\tn\n"
                     "rs1\t1\t1000\t3.5\t0.01\t5000\n"
                     "rs2\t1\t2000\t1.2\t0.01\t5000\n")
        df = load_sumstats(f)
        assert list(df["rsid"]) == ["rs1", "rs2"]
        assert pytest.approx(df.loc[0, "z"]) == 3.5

    def test_load_sumstats_beta_se(self, tmp_path):
        """load_sumstats computes z = beta/se when 'z' column is absent."""
        f = tmp_path / "stats.tsv"
        f.write_text("rsid\tchr\tpos\tbeta\tse\n"
                     "rs1\t1\t1000\t0.05\t0.01\n"
                     "rs2\t1\t2000\t0.02\t0.01\n")
        df = load_sumstats(f)
        assert pytest.approx(df.loc[0, "z"]) == 5.0
        assert pytest.approx(df.loc[1, "z"]) == 2.0

    def test_load_sumstats_alias_snp(self, tmp_path):
        """'snp' column is aliased to 'rsid'."""
        f = tmp_path / "stats.tsv"
        f.write_text("snp\tz\n"
                     "rs1\t3.5\n"
                     "rs2\t1.2\n")
        df = load_sumstats(f)
        assert "rsid" in df.columns

    def test_load_sumstats_missing_z_and_beta(self, tmp_path):
        """Raises ValueError when neither 'z' nor 'beta'+'se' are present."""
        f = tmp_path / "bad.tsv"
        f.write_text("rsid\tchr\tpos\tp\n"
                     "rs1\t1\t1000\t0.05\n")
        with pytest.raises(ValueError, match="z.*column|beta.*se"):
            load_sumstats(f)

    def test_load_sumstats_missing_rsid(self, tmp_path):
        """Raises ValueError when no rsid-equivalent column is present."""
        f = tmp_path / "bad.tsv"
        f.write_text("chr\tpos\tz\n"
                     "1\t1000\t3.5\n")
        with pytest.raises(ValueError, match="rsid"):
            load_sumstats(f)

    def test_load_sumstats_locus_filter(self, tmp_path):
        """Locus window filter keeps only variants in the specified range."""
        f = tmp_path / "stats.tsv"
        f.write_text("rsid\tchr\tpos\tz\n"
                     "rs1\t1\t500000\t1.0\n"
                     "rs2\t1\t1000000\t2.0\n"
                     "rs3\t1\t2000000\t3.0\n"
                     "rs4\t2\t1000000\t4.0\n")
        df = load_sumstats(f, chr_="1", start=900_000, end=1_100_000)
        assert list(df["rsid"]) == ["rs2"]

    def test_load_sumstats_csv(self, tmp_path):
        """load_sumstats handles comma-delimited files."""
        f = tmp_path / "stats.csv"
        f.write_text("rsid,chr,pos,z\n"
                     "rs1,1,1000,3.5\n"
                     "rs2,1,2000,1.2\n")
        df = load_sumstats(f)
        assert len(df) == 2

    def test_load_ld_npy(self, tmp_path):
        """load_ld loads a square .npy LD matrix correctly."""
        R = np.eye(5)
        path = tmp_path / "ld.npy"
        np.save(path, R)
        R_loaded = load_ld(path, 5)
        np.testing.assert_allclose(R_loaded, R)

    def test_load_ld_wrong_size(self, tmp_path):
        """load_ld raises ValueError when matrix size != n_variants."""
        R = np.eye(5)
        path = tmp_path / "ld.npy"
        np.save(path, R)
        with pytest.raises(ValueError, match="dimension"):
            load_ld(path, 10)

    def test_load_ld_symmetrised(self, tmp_path):
        """load_ld symmetrises the matrix and sets diagonal to 1."""
        R = np.array([[1.0, 0.6], [0.4, 1.0]])  # slightly asymmetric
        path = tmp_path / "ld.npy"
        np.save(path, R)
        R_loaded = load_ld(path, 2)
        np.testing.assert_allclose(R_loaded, R_loaded.T, atol=1e-10)
        np.testing.assert_allclose(np.diag(R_loaded), 1.0, atol=1e-10)

    def test_load_ld_tsv_with_string_header(self, tmp_path):
        """load_ld skips a string header row (e.g. rsIDs) in a TSV LD matrix."""
        R = np.eye(3)
        path = tmp_path / "ld.tsv"
        # Write with rsID header row
        lines = ["rs1\trs2\trs3\n"] + ["\t".join(f"{v:.1f}" for v in row) + "\n" for row in R]
        path.write_text("".join(lines))
        R_loaded = load_ld(path, 3)
        np.testing.assert_allclose(R_loaded, R, atol=1e-6)

    def test_load_ld_tsv_without_header(self, tmp_path):
        """load_ld loads a numeric-only TSV LD matrix correctly."""
        R = np.array([[1.0, 0.8], [0.8, 1.0]])
        path = tmp_path / "ld.tsv"
        lines = ["\t".join(f"{v:.1f}" for v in row) + "\n" for row in R]
        path.write_text("".join(lines))
        R_loaded = load_ld(path, 2)
        np.testing.assert_allclose(R_loaded, R, atol=1e-6)

    def test_load_sumstats_alias_pval(self, tmp_path):
        """'pval' column is aliased to 'p'."""
        f = tmp_path / "stats.tsv"
        f.write_text("rsid\tchr\tpos\tz\tpval\n"
                     "rs1\t1\t1000\t3.5\t0.001\n"
                     "rs2\t1\t2000\t1.2\t0.23\n")
        df = load_sumstats(f)
        assert "p" in df.columns
        assert pytest.approx(df.loc[0, "p"]) == 0.001

    def test_load_sumstats_alias_p_value(self, tmp_path):
        """'p_value' column is aliased to 'p'."""
        f = tmp_path / "stats.tsv"
        f.write_text("rsid\tz\tp_value\n"
                     "rs1\t3.5\t0.001\n")
        df = load_sumstats(f)
        assert "p" in df.columns

    def test_load_sumstats_alias_pvalue(self, tmp_path):
        """'pvalue' column is aliased to 'p'."""
        f = tmp_path / "stats.tsv"
        f.write_text("rsid\tz\tpvalue\n"
                     "rs1\t3.5\t0.001\n")
        df = load_sumstats(f)
        assert "p" in df.columns


# ---------------------------------------------------------------------------
# TestRunFinemapping
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def demo_run(tmp_path_factory):
    """Run the demo pipeline once and share results across all tests."""
    tmp = tmp_path_factory.mktemp("demo_run")
    results = fine_mapping.run_finemapping(
        sumstats_path=None,
        ld_path=None,
        output_dir=tmp,
        demo=True,
        make_figures=False,
    )
    return results, tmp


class TestRunFinemapping:
    """Integration tests for the top-level pipeline."""

    @needs_jax
    def test_demo_mode_produces_outputs(self, demo_run):
        """Demo mode writes report.md, fine_mapping.json, and tables/pips.tsv."""
        _, tmp = demo_run
        assert (tmp / "report.md").exists()
        assert (tmp / "fine_mapping.json").exists()
        assert (tmp / "tables" / "pips.tsv").exists()

    @needs_jax
    def test_demo_results_structure(self, demo_run):
        """Demo run returns a dict with expected top-level keys."""
        results, _ = demo_run
        for key in ("method", "n_variants", "n_credible_sets", "lead_rsid", "lead_pip"):
            assert key in results, f"Missing key: {key}"

    @needs_jax
    def test_demo_susie_method(self, demo_run):
        """Demo mode includes a synthetic LD matrix, so it uses SuSiE."""
        results, _ = demo_run
        assert results["method"] == "SuSiE"

    @needs_jax
    def test_susie_method_with_ld(self, tmp_path):
        """When LD is provided, the pipeline uses SuSiE."""
        df, R = fine_mapping.make_demo_data(seed=42)
        ss_path = tmp_path / "sumstats.tsv"
        df.to_csv(ss_path, sep="\t", index=False)
        ld_path = tmp_path / "ld.npy"
        np.save(ld_path, R)

        results = fine_mapping.run_finemapping(
            sumstats_path=ss_path,
            ld_path=ld_path,
            output_dir=tmp_path / "out",
            make_figures=False,
        )
        assert results["method"] == "SuSiE"

    def test_abf_method_without_ld(self, tmp_path):
        """Without LD matrix, pipeline uses ABF."""
        df, _ = fine_mapping.make_demo_data(seed=42)
        ss_path = tmp_path / "sumstats.tsv"
        df.to_csv(ss_path, sep="\t", index=False)

        results = fine_mapping.run_finemapping(
            sumstats_path=ss_path,
            ld_path=None,
            output_dir=tmp_path / "out",
            make_figures=False,
        )
        assert results["method"] == "ABF"

    @needs_jax
    def test_lead_pip_in_range(self, demo_run):
        """Lead PIP is a valid probability in [0, 1]."""
        results, _ = demo_run
        assert 0.0 <= results["lead_pip"] <= 1.0

    @needs_jax
    def test_report_contains_disclaimer(self, demo_run):
        """Generated report contains the required safety disclaimer."""
        _, tmp = demo_run
        report = (tmp / "report.md").read_text()
        assert "ClawBio is a research and educational tool" in report

    @needs_jax
    def test_results_json_parseable(self, demo_run):
        """fine_mapping.json is valid JSON with expected keys."""
        _, tmp = demo_run
        data = json.loads((tmp / "fine_mapping.json").read_text())
        assert "credible_sets" in data
        assert "n_variants" in data

    def test_empty_locus_raises(self, tmp_path):
        """Raises ValueError when locus filter returns no variants."""
        df, _ = fine_mapping.make_demo_data()
        ss_path = tmp_path / "sumstats.tsv"
        df.to_csv(ss_path, sep="\t", index=False)

        with pytest.raises(ValueError, match="No variants"):
            fine_mapping.run_finemapping(
                sumstats_path=ss_path,
                ld_path=None,
                output_dir=tmp_path / "out",
                chr_="99",
                start=1,
                end=2,
                make_figures=False,
            )

    def test_locus_window_filter(self, tmp_path):
        """Locus window filter restricts analysis to specified region."""
        df, _ = fine_mapping.make_demo_data()
        ss_path = tmp_path / "sumstats.tsv"
        df.to_csv(ss_path, sep="\t", index=False)

        results = fine_mapping.run_finemapping(
            sumstats_path=ss_path,
            ld_path=None,
            output_dir=tmp_path / "out",
            chr_="1",
            start=109_000_000,
            end=109_050_000,
            make_figures=False,
        )
        # Only ~25% of the 200-variant locus falls in this window
        assert results["n_variants"] < 200

    @needs_jax
    def test_gene_track_param_accepted(self, tmp_path):
        """run_finemapping accepts gene_track=True without error (figures skipped)."""
        results = fine_mapping.run_finemapping(
            sumstats_path=None,
            ld_path=None,
            output_dir=tmp_path,
            demo=True,
            make_figures=False,
            gene_track=True,
        )
        assert "method" in results

    def test_susie_path_rejects_coverage_one_early(self, tmp_path):
        """coverage=1.0 must fail at the CLI, not deep inside run_susie.

        The ABF path still accepts 1.0; only the SuSiE path carries sushie's
        open-interval constraint, so the check is scoped to it. (coverage=0 is
        already rejected for both methods by the shared check above it.)
        """
        with pytest.raises(ValueError, match="SuSiE path"):
            fine_mapping.run_finemapping(
                sumstats_path=None, ld_path=None, output_dir=tmp_path,
                demo=True, make_figures=False, coverage=1.0,
            )

    def test_nonpositive_coverage_rejected_for_both_methods(self, tmp_path):
        with pytest.raises(ValueError, match="Coverage must be in"):
            fine_mapping.run_finemapping(
                sumstats_path=None, ld_path=None, output_dir=tmp_path,
                demo=True, make_figures=False, coverage=0.0,
            )

    def test_susie_path_rejects_zero_min_purity_early(self, tmp_path):
        """Same for min_purity=0, which the old [0, 1] check let through."""
        with pytest.raises(ValueError, match="SuSiE path"):
            fine_mapping.run_finemapping(
                sumstats_path=None, ld_path=None, output_dir=tmp_path,
                demo=True, make_figures=False, min_purity=0.0,
            )

    def test_abf_path_still_accepts_coverage_one(self, tmp_path):
        """The tightened check must not leak onto ABF, which has no such limit."""
        df, _ = fine_mapping.make_demo_data(seed=42)
        ss_path = tmp_path / "sumstats.tsv"
        df.to_csv(ss_path, sep="\t", index=False)
        results = fine_mapping.run_finemapping(
            sumstats_path=ss_path, ld_path=None, output_dir=tmp_path,
            make_figures=False, coverage=1.0,
        )
        assert results["method"] == "ABF"


# ---------------------------------------------------------------------------
# TestReportProvisionalBanner
# ---------------------------------------------------------------------------


class TestReportProvisionalBanner:
    """A non-converged SuSiE fit must say so above its own numbers.

    SKILL.md Gotcha 4 tells the reader not to report PIPs from a run that did
    not converge, but the only signal in the report was one row of a parameter
    table, well below the credible-set tables a reader actually reads.
    """

    def _markdown(self, **params):
        df = pd.DataFrame({
            "rsid": ["rs1", "rs2"], "chr": ["1", "1"],
            "pos": [1000, 2000], "z": [5.0, 0.2],
            "p": [1e-7, 0.8], "pip": [0.9, 0.1],
        })
        base = {"L": 10, "coverage": 0.95, "min_purity": 0.5, "n_iter": 500,
                "max_iter": 500}
        base.update(params)
        return generate_markdown(df=df, credible_sets=[], method="SuSiE",
                                 params=base)

    def test_banner_present_when_not_converged(self):
        md = self._markdown(converged=False)
        assert "Provisional results" in md
        assert md.index("Provisional results") < md.index("## Credible Sets")

    def test_banner_absent_when_converged(self):
        assert "Provisional results" not in self._markdown(converged=True)

    def test_banner_absent_for_abf(self):
        """ABF has no convergence concept, so it must never carry the banner."""
        df = pd.DataFrame({
            "rsid": ["rs1"], "chr": ["1"], "pos": [1000],
            "z": [5.0], "p": [1e-7], "pip": [0.9],
        })
        md = generate_markdown(df=df, credible_sets=[], method="ABF",
                               params={"coverage": 0.95, "w": 0.04})
        assert "Provisional results" not in md

    def test_banner_reports_the_iteration_budget(self):
        md = self._markdown(converged=False, n_iter=1, max_iter=500)
        assert "did not converge in 1 iterations" in md
        assert "max_iter=500" in md


# ---------------------------------------------------------------------------
# TestSuSiEInf
# ---------------------------------------------------------------------------


class TestSuSiEInf:
    """Tests for the SuSiE-inf IBSS algorithm with infinitesimal component."""

    def test_pip_shape(self):
        """run_susie_inf returns a 1-D PIP array of length p."""
        from fine_mapping_core.susie_inf import run_susie_inf
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie_inf(z=df["z"].values, R=R, n=5000, L=5)
        assert result["pip"].shape == (20,)

    def test_pip_range(self):
        """All PIPs returned by run_susie_inf are in [0, 1]."""
        from fine_mapping_core.susie_inf import run_susie_inf
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie_inf(z=df["z"].values, R=R, n=5000, L=5)
        assert result["pip"].min() >= 0.0
        assert result["pip"].max() <= 1.0

    def test_alpha_shape(self):
        """run_susie_inf returns alpha as a (p, L) matrix."""
        from fine_mapping_core.susie_inf import run_susie_inf
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie_inf(z=df["z"].values, R=R, n=5000, L=3)
        assert result["alpha"].shape == (20, 3)

    def test_alpha_columns_sum_to_at_most_one(self):
        """Each column of alpha sums to <= 1 (remainder goes to null component)."""
        from fine_mapping_core.susie_inf import run_susie_inf
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie_inf(z=df["z"].values, R=R, n=5000, L=3)
        col_sums = result["alpha"].sum(axis=0)
        assert np.all(col_sums <= 1.0 + 1e-5), f"alpha columns exceed 1: {col_sums}"
        assert np.all(col_sums > 0.0), f"alpha columns are zero: {col_sums}"

    def test_signal_recovers_high_pip(self):
        """The injected signal at index 10 receives the highest PIP."""
        from fine_mapping_core.susie_inf import run_susie_inf
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie_inf(z=df["z"].values, R=R, n=5000, L=5)
        assert result["pip"].argmax() == 10

    def test_variance_components_returned(self):
        """run_susie_inf returns sigmasq and tausq scalar variance components."""
        from fine_mapping_core.susie_inf import run_susie_inf
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie_inf(z=df["z"].values, R=R, n=5000, L=5)
        assert np.isfinite(result["sigmasq"])
        assert np.isfinite(result["tausq"])
        assert result["sigmasq"] > 0

    def test_ssq_length_matches_L(self):
        """ssq (per-effect prior variances) has length L."""
        from fine_mapping_core.susie_inf import run_susie_inf
        df = _small_locus(n=20)
        R = _identity_ld(20)
        L = 4
        result = run_susie_inf(z=df["z"].values, R=R, n=5000, L=L)
        assert len(result["ssq"]) == L

    def test_two_signals_recovered(self):
        """Demo locus with two causal variants: both appear in top-5 PIPs."""
        from fine_mapping_core.susie_inf import run_susie_inf
        df, R = fine_mapping.make_demo_data(seed=42)
        result = run_susie_inf(z=df["z"].values, R=R, n=5000, L=10)
        top5 = set(np.argsort(-result["pip"])[:5])
        assert 60 in top5 or 140 in top5

    def test_identity_ld_reduces_to_sparse(self):
        """Under identity LD (no LD), tausq should converge near zero."""
        from fine_mapping_core.susie_inf import run_susie_inf
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie_inf(z=df["z"].values, R=R, n=5000, L=5)
        # tausq should be non-negative
        assert result["tausq"] >= 0.0

    def test_cred_inf_returns_list_of_index_lists(self):
        """cred_inf returns a list where each element is a list of SNP indices."""
        from fine_mapping_core.susie_inf import run_susie_inf, cred_inf
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie_inf(z=df["z"].values, R=R, n=5000, L=3)
        csets = cred_inf(result["alpha"], R=R, coverage=0.95, purity=0.0)
        assert isinstance(csets, list)
        for cs in csets:
            assert isinstance(cs, list)
            assert all(isinstance(i, (int, np.integer)) for i in cs)

    def test_cred_inf_captures_injected_signal(self):
        """cred_inf includes index 10 (injected signal) in at least one credible set."""
        from fine_mapping_core.susie_inf import run_susie_inf, cred_inf
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie_inf(z=df["z"].values, R=R, n=5000, L=5)
        csets = cred_inf(result["alpha"], R=R, coverage=0.95, purity=0.0)
        all_indices = {i for cs in csets for i in cs}
        assert 10 in all_indices


# ---------------------------------------------------------------------------
# TestGeneTrack
# ---------------------------------------------------------------------------


class TestGeneTrack:
    """Tests for gene track fetching and rendering."""

    def test_fetch_genes_returns_empty_on_network_error(self):
        """_fetch_genes returns [] when the network request fails."""
        from fine_mapping_core.report import _fetch_genes
        from unittest.mock import patch

        with patch("urllib.request.urlopen", side_effect=OSError("network error")):
            result = _fetch_genes("1", 1_000_000, 2_000_000)
        assert result == []

    def test_fetch_genes_retries_on_429(self):
        """_fetch_genes sleeps for Retry-After seconds and retries on HTTP 429."""
        import json
        import urllib.error
        from fine_mapping_core.report import _fetch_genes
        from unittest.mock import MagicMock, call, patch

        # First two calls raise 429; third succeeds
        gene_data = [{"feature_type": "gene", "gene_id": "ENSG001",
                      "start": 1000, "end": 2000, "strand": 1}]

        headers_429 = {"Retry-After": "0.01"}
        http_429 = urllib.error.HTTPError(
            url="", code=429, msg="Too Many Requests", hdrs=headers_429, fp=None
        )

        success_resp = MagicMock()
        success_resp.read.return_value = json.dumps(gene_data).encode()
        success_resp.__enter__ = lambda s: s
        success_resp.__exit__ = MagicMock(return_value=False)

        side_effects = [http_429, http_429, success_resp]

        with patch("urllib.request.urlopen", side_effect=side_effects), \
             patch("time.sleep") as mock_sleep:
            result = _fetch_genes("1", 1_000_000, 2_000_000)

        assert len(result) == 1
        assert mock_sleep.call_count == 2
        mock_sleep.assert_called_with(pytest.approx(0.01))

    def test_fetch_genes_returns_empty_after_max_retries(self):
        """_fetch_genes returns [] when all 3 attempts receive HTTP 429."""
        import urllib.error
        from fine_mapping_core.report import _fetch_genes
        from unittest.mock import patch

        headers_429 = {"Retry-After": "0.01"}
        http_429 = urllib.error.HTTPError(
            url="", code=429, msg="Too Many Requests", hdrs=headers_429, fp=None
        )

        with patch("urllib.request.urlopen", side_effect=[http_429, http_429, http_429]), \
             patch("time.sleep"):
            result = _fetch_genes("1", 1_000_000, 2_000_000)

        assert result == []

    def test_fetch_genes_returns_empty_on_non_429_http_error(self):
        """_fetch_genes returns [] immediately on non-429 HTTP errors (no retry)."""
        import urllib.error
        from fine_mapping_core.report import _fetch_genes
        from unittest.mock import patch

        http_500 = urllib.error.HTTPError(
            url="", code=500, msg="Internal Server Error", hdrs={}, fp=None
        )

        with patch("urllib.request.urlopen", side_effect=http_500) as mock_urlopen, \
             patch("time.sleep") as mock_sleep:
            result = _fetch_genes("1", 1_000_000, 2_000_000)

        assert result == []
        assert mock_urlopen.call_count == 1  # no retry
        mock_sleep.assert_not_called()

    def test_fetch_genes_filters_to_gene_feature_type(self):
        """_fetch_genes drops non-gene features (e.g. transcripts)."""
        import json
        from fine_mapping_core.report import _fetch_genes
        from unittest.mock import MagicMock, patch

        mock_data = [
            {"feature_type": "gene", "gene_id": "ENSG001",
             "start": 1000, "end": 2000, "strand": 1},
            {"feature_type": "transcript", "gene_id": "ENSG001",
             "start": 1000, "end": 2000, "strand": 1},
        ]
        mock_resp = MagicMock()
        mock_resp.read.return_value = json.dumps(mock_data).encode()
        mock_resp.__enter__ = lambda s: s
        mock_resp.__exit__ = MagicMock(return_value=False)

        with patch("urllib.request.urlopen", return_value=mock_resp):
            result = _fetch_genes("1", 1_000_000, 2_000_000)

        assert len(result) == 1
        assert result[0]["feature_type"] == "gene"

    def test_fetch_genes_strips_chr_prefix(self):
        """_fetch_genes strips 'chr' prefix from chromosome before building URL."""
        import json
        from fine_mapping_core.report import _fetch_genes
        from unittest.mock import MagicMock, call, patch

        mock_resp = MagicMock()
        mock_resp.read.return_value = json.dumps([]).encode()
        mock_resp.__enter__ = lambda s: s
        mock_resp.__exit__ = MagicMock(return_value=False)

        captured = {}

        def capture_urlopen(req, **kwargs):
            captured["url"] = req.full_url
            return mock_resp

        with patch("urllib.request.urlopen", side_effect=capture_urlopen):
            _fetch_genes("chr7", 1_000_000, 2_000_000)

        assert "/7:" in captured["url"], "chr prefix should be stripped from URL"

    def test_plot_gene_track_renders_without_error(self):
        """_plot_gene_track draws strand-aware gene arrows without raising."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from fine_mapping_core.report import _plot_gene_track

        genes = [
            {"external_name": "GENE1", "gene_id": "ENSG001",
             "start": 1_000, "end": 5_000, "strand": 1},
            {"external_name": "GENE2", "gene_id": "ENSG002",
             "start": 3_000, "end": 8_000, "strand": -1},
        ]
        fig, ax = plt.subplots()
        _plot_gene_track(ax, genes, 0, 10_000)
        plt.close(fig)

    def test_plot_gene_track_empty_gene_list(self):
        """_plot_gene_track handles an empty gene list without error."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from fine_mapping_core.report import _plot_gene_track

        fig, ax = plt.subplots()
        _plot_gene_track(ax, [], 0, 10_000)
        plt.close(fig)

    def test_regional_association_gene_track_no_genes(self, tmp_path):
        """gene_track=True with empty gene response produces a single-panel figure."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from fine_mapping_core.report import _plot_regional_association
        from unittest.mock import patch

        df = _small_locus()
        with patch("fine_mapping_core.report._fetch_genes", return_value=[]):
            _plot_regional_association(df, [], tmp_path, plt, gene_track=True)

        assert (tmp_path / "regional_association.png").exists()

    def test_regional_association_gene_track_with_genes(self, tmp_path):
        """gene_track=True with gene data produces a two-panel figure."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from fine_mapping_core.report import _plot_regional_association
        from unittest.mock import patch

        genes = [
            {"feature_type": "gene", "external_name": "GENE1", "gene_id": "ENSG001",
             "start": 1_000_000, "end": 1_005_000, "strand": 1},
            {"feature_type": "gene", "external_name": "GENE2", "gene_id": "ENSG002",
             "start": 1_010_000, "end": 1_015_000, "strand": -1},
        ]
        df = _small_locus()
        with patch("fine_mapping_core.report._fetch_genes", return_value=genes):
            _plot_regional_association(df, [], tmp_path, plt, gene_track=True)

        assert (tmp_path / "regional_association.png").exists()

    def test_regional_association_gene_track_false_by_default(self, tmp_path):
        """_plot_regional_association never calls _fetch_genes when gene_track=False."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from fine_mapping_core.report import _plot_regional_association
        from unittest.mock import patch

        df = _small_locus()
        with patch("fine_mapping_core.report._fetch_genes") as mock_fetch:
            _plot_regional_association(df, [], tmp_path, plt, gene_track=False)

        mock_fetch.assert_not_called()


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/fine-mapping/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@needs_jax
@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("fine-mapping", tmp_path) == snapshot({'fine_mapping.json': """\
{
  "method": "SuSiE",
  "params": {
    "L": 10,
    "coverage": 0.95,
    "min_purity": 0.5,
    "w": 0.04,
    "converged": true,
    "n_iter": 176,
    "max_iter": 500,
    "n_eff": 5000,
    "engine": "sushie",
    "engine_version": "0.20"
  },
  "n_variants": 200,
  "n_credible_sets": 2,
  "lead_rsid": "rs_demo_0061",
  "lead_pip": 0.9999987663493761,
  "credible_sets": [
    {
      "cs_id": "L1",
      "signal_index": 0,
      "size": 1,
      "coverage": 0.9999987663493761,
      "lead_rsid": "rs_demo_0061",
      "lead_alpha": 0.9999987663493761,
      "purity": 1.0,
      "pure": true,
      "variants": [
        {
          "rsid": "rs_demo_0061",
          "chr": "1",
          "pos": 109060301,
          "z": 16.701906539234432,
          "pip": 0.9999987663493761,
          "alpha": 0.9999987663493761,
          "p": 1.345486059166518e-62
        }
      ]
    },
    {
      "cs_id": "L2",
      "signal_index": 1,
      "size": 1,
      "coverage": 0.9919991978808628,
      "lead_rsid": "rs_demo_0141",
      "lead_alpha": 0.9919991978808628,
      "purity": 1.0,
      "pure": true,
      "variants": [
        {
          "rsid": "rs_demo_0141",
          "chr": "1",
          "pos": 109140703,
          "z": 15.440428533490131,
          "pip": 0.9919991978808628,
          "alpha": 0.9919991978808628,
          "p": 9.207641370365281e-54
        }
      ]
    }
  ],
  "pip_above_0_1": 2,
  "pip_above_0_5": 2
}\
""", 'report.md': """\
# SuSiE Fine-Mapping Report

**Date**: <TIMESTAMP> UTC
**Skill**: fine-mapping
**Method**: SuSiE
**Input**: synthetic_demo
**Variants analysed**: 200

---

## Summary

| Parameter | Value |
|-----------|-------|
| Method | SuSiE |
| Variants | 200 |
| Credible sets | 2 |
| Lead variant | rs_demo_0061 (PIP = 1.000) |
| Coverage threshold | 95% |
| Max signals (L) | 10 |
| Min purity | 0.5 |
| Converged | True |
| Iterations | 176 |

---

## Credible Sets

### L1

- **Size**: 1 variants
- **Coverage**: 100.0%
- **Lead variant**: rs_demo_0061 (α = 1.0000)
- **Purity** (min |r|): 1.000

| rsID | Chr | Pos | Z | PIP | α |
|------|-----|-----|---|-----|---|
| rs_demo_0061 | 1 | 109,060,301 | 16.702 | 1.0000 | 1.0000 |

### L2

- **Size**: 1 variants
- **Coverage**: 99.2%
- **Lead variant**: rs_demo_0141 (α = 0.9920)
- **Purity** (min |r|): 1.000

| rsID | Chr | Pos | Z | PIP | α |
|------|-----|-----|---|-----|---|
| rs_demo_0141 | 1 | 109,140,703 | 15.440 | 0.9920 | 0.9920 |

## High-PIP Variants (PIP ≥ 0.1)

| rsID | Chr | Pos | Z | PIP | CS |
|------|-----|-----|---|-----|----|
| rs_demo_0061 | 1 | 109,060,301 | 16.702 | 1.0000 | L1 |
| rs_demo_0141 | 1 | 109,140,703 | 15.440 | 0.9920 | L2 |

## Methodology

Fine-mapping was performed using the **SuSiE** (Sum of Single Effects) algorithm (Wang et al. 2020, JRSS-B). SuSiE models the GWAS signal as a sum of L single effects, each with a sparse prior, and estimates posterior inclusion probabilities (PIPs) via Iterative Bayesian Stepwise Selection (IBSS). An LD matrix was used to account for linkage disequilibrium between variants.

**Citations**:
- Wang G et al. (2020) *A simple new approach to variable selection in regression, with application to genetic fine mapping*. JRSS-B. doi:10.1111/rssb.12388

---

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*\
""", 'result.json': """\
{
  "skill": "fine-mapping",
  "version": "0.1.0",
  "completed_at": "<TIMESTAMP>",
  "input_checksum": "",
  "datasets": {},
  "summary": {
    "method": "SuSiE",
    "n_variants": 200,
    "n_credible_sets": 2,
    "lead_rsid": "rs_demo_0061",
    "lead_pip": 0.9999987663493761
  },
  "data": {
    "method": "SuSiE",
    "params": {
      "L": 10,
      "coverage": 0.95,
      "min_purity": 0.5,
      "w": 0.04,
      "converged": true,
      "n_iter": 176,
      "max_iter": 500,
      "n_eff": 5000,
      "engine": "sushie",
      "engine_version": "0.20"
    },
    "n_variants": 200,
    "n_credible_sets": 2,
    "lead_rsid": "rs_demo_0061",
    "lead_pip": 0.9999987663493761,
    "credible_sets": [
      {
        "cs_id": "L1",
        "signal_index": 0,
        "size": 1,
        "coverage": 0.9999987663493761,
        "lead_rsid": "rs_demo_0061",
        "lead_alpha": 0.9999987663493761,
        "purity": 1.0,
        "pure": true,
        "variants": [
          {
            "rsid": "rs_demo_0061",
            "chr": "1",
            "pos": 109060301,
            "z": 16.701906539234432,
            "pip": 0.9999987663493761,
            "alpha": 0.9999987663493761,
            "p": 1.345486059166518e-62
          }
        ]
      },
      {
        "cs_id": "L2",
        "signal_index": 1,
        "size": 1,
        "coverage": 0.9919991978808628,
        "lead_rsid": "rs_demo_0141",
        "lead_alpha": 0.9919991978808628,
        "purity": 1.0,
        "pure": true,
        "variants": [
          {
            "rsid": "rs_demo_0141",
            "chr": "1",
            "pos": 109140703,
            "z": 15.440428533490131,
            "pip": 0.9919991978808628,
            "alpha": 0.9919991978808628,
            "p": 9.207641370365281e-54
          }
        ]
      }
    ],
    "pip_above_0_1": 2,
    "pip_above_0_5": 2
  }
}\
""", 'tables/credible_sets.tsv': """\
cs_id	cs_size	cs_coverage	lead_rsid	purity	rsid	chr	pos	z	pip	alpha
L1	1	0.999999	rs_demo_0061	1	rs_demo_0061	1	109060301	16.7019	0.999999	0.999999
L2	1	0.991999	rs_demo_0141	1	rs_demo_0141	1	109140703	15.4404	0.991999	0.991999
""", 'tables/pips.tsv': """\
rsid	chr	pos	z	se	p	pip	cs_membership
rs_demo_0001	1	109000000	1.01185	0.0141421	0.311609	2.64129e-49	
rs_demo_0002	1	109001005	-0.332843	0.0141421	0.739253	1.59779e-49	
rs_demo_0003	1	109002010	1.4576	0.0141421	0.144951	4.83815e-49	
rs_demo_0004	1	109003015	1.64772	0.0141421	0.0994092	6.69823e-49	
rs_demo_0005	1	109004020	-1.24386	0.0141421	0.21355	3.4992e-49	
rs_demo_0006	1	109005025	-0.59499	0.0141421	0.55185	1.82278e-49	
rs_demo_0007	1	109006030	0.835051	0.0141421	0.403689	2.21094e-49	
rs_demo_0008	1	109007035	0.390993	0.0141421	0.695802	1.63718e-49	
rs_demo_0009	1	109008040	0.690467	0.0141421	0.4899	1.95527e-49	
rs_demo_0010	1	109009045	-0.145735	0.0141421	0.88413	1.52207e-49	
rs_demo_0011	1	109010050	1.58676	0.0141421	0.112568	5.96768e-49	
rs_demo_0012	1	109011055	1.48521	0.0141421	0.137487	5.0223e-49	
rs_demo_0013	1	109012060	0.773532	0.0141421	0.439208	2.08318e-49	
rs_demo_0014	1	109013065	1.83484	0.0141421	0.0665292	9.46592e-49	
rs_demo_0015	1	109014070	1.17523	0.0141421	0.239902	3.19652e-49	
rs_demo_0016	1	109015075	-0.151416	0.0141421	0.879648	1.52406e-49	
rs_demo_0017	1	109016080	1.07682	0.0141421	0.281561	2.83376e-49	
rs_demo_0018	1	109017085	-0.250573	0.0141421	0.802145	1.55822e-49	
rs_demo_0019	1	109018090	1.58706	0.0141421	0.112499	5.95667e-49	
rs_demo_0020	1	109019095	0.659061	0.0141421	0.509857	1.90441e-49	
rs_demo_0021	1	109020100	0.524594	0.0141421	0.599865	1.74565e-49	
rs_demo_0022	1	109021105	0.0291144	0.0141421	0.976773	1.50482e-49	
rs_demo_0023	1	109022110	1.93332	0.0141421	0.0531967	1.15396e-48	
rs_demo_0024	1	109023115	0.557167	0.0141421	0.577413	1.77643e-49	
rs_demo_0025	1	109024120	0.284516	0.0141421	0.776015	1.56886e-49	
rs_demo_0026	1	109025125	0.362144	0.0141421	0.717244	1.61047e-49	
rs_demo_0027	1	109026130	1.24838	0.0141421	0.211892	3.47954e-49	
rs_demo_0028	1	109027135	1.08376	0.0141421	0.278473	2.81596e-49	
rs_demo_0029	1	109028140	1.13385	0.0141421	0.25686	2.97931e-49	
rs_demo_0030	1	109029145	1.15543	0.0141421	0.247912	3.04645e-49	
rs_demo_0031	1	109030150	2.87064	0.0141421	0.00409657	1.2704e-47	
rs_demo_0032	1	109031155	0.328047	0.0141421	0.742877	1.58005e-49	
rs_demo_0033	1	109032160	0.229058	0.0141421	0.818824	1.53611e-49	
rs_demo_0034	1	109033165	-0.063924	0.0141421	0.949031	1.51346e-49	
rs_demo_0035	1	109034170	1.37651	0.0141421	0.168663	3.94651e-49	
rs_demo_0036	1	109035175	1.90286	0.0141421	0.0570582	9.64314e-49	
rs_demo_0037	1	109036180	0.67664	0.0141421	0.498634	1.8333e-49	
rs_demo_0038	1	109037185	-0.0286992	0.0141421	0.977104	1.51684e-49	
rs_demo_0039	1	109038190	0.0130637	0.0141421	0.989577	1.51332e-49	
rs_demo_0040	1	109039195	1.52075	0.0141421	0.128323	4.23807e-49	
rs_demo_0041	1	109040201	1.65417	0.0141421	0.0980928	4.93073e-49	
rs_demo_0042	1	109041206	1.50502	0.0141421	0.132318	3.68093e-49	
rs_demo_0043	1	109042211	0.36005	0.0141421	0.71881	1.50916e-49	
rs_demo_0044	1	109043216	1.33733	0.0141421	0.181114	2.55738e-49	
rs_demo_0045	1	109044221	1.32137	0.0141421	0.186377	2.30009e-49	
rs_demo_0046	1	109045226	1.54777	0.0141421	0.121677	2.59359e-49	
rs_demo_0047	1	109046231	2.35601	0.0141421	0.0184725	6.91075e-49	
rs_demo_0048	1	109047236	1.90254	0.0141421	0.0571001	2.73443e-49	
rs_demo_0049	1	109048241	2.60082	0.0141421	0.00930018	5.41678e-49	
rs_demo_0050	1	109049246	2.29319	0.0141421	0.0218372	2.48153e-49	
rs_demo_0051	1	109050251	2.89435	0.0141421	0.00379956	3.41738e-49	
rs_demo_0052	1	109051256	3.71105	0.0141421	0.000206459	6.36975e-49	
rs_demo_0053	1	109052261	2.21577	0.0141421	0.026707	1.63432e-49	
rs_demo_0054	1	109053266	4.09471	0.0141421	4.22922e-05	2.20587e-49	
rs_demo_0055	1	109054271	4.87083	0.0141421	1.11285e-06	2.12697e-49	
rs_demo_0056	1	109055276	5.86085	0.0141421	4.62006e-09	2.07193e-49	
rs_demo_0057	1	109056281	7.67274	0.0141421	1.69868e-14	7.25516e-47	
rs_demo_0058	1	109057286	11.253	0.0141421	2.29776e-29	1.66833e-32	
rs_demo_0059	1	109058291	11.155	0.0141421	6.94721e-29	5.67095e-33	
rs_demo_0060	1	109059296	15.8175	0.0141421	2.48414e-56	1.23365e-06	
rs_demo_0061	1	109060301	16.7019	0.0141421	1.34549e-62	0.999999	L1
rs_demo_0062	1	109061306	14.5144	0.0141421	1.02829e-47	7.93516e-15	
rs_demo_0063	1	109062311	12.1836	0.0141421	3.9248e-34	6.81596e-28	
rs_demo_0064	1	109063316	10.3443	0.0141421	4.53759e-25	1.12888e-36	
rs_demo_0065	1	109064321	8.65911	0.0141421	4.81713e-18	1.81145e-43	
rs_demo_0066	1	109065326	7.29307	0.0141421	3.05241e-13	6.5893e-48	
rs_demo_0067	1	109066331	4.99248	0.0141421	5.97035e-07	2.37826e-49	
rs_demo_0068	1	109067336	3.95203	0.0141421	7.75243e-05	1.95699e-49	
rs_demo_0069	1	109068341	4.5309	0.0141421	5.87881e-06	1.13463e-48	
rs_demo_0070	1	109069346	2.88846	0.0141421	0.00387148	2.1457e-49	
rs_demo_0071	1	109070351	1.32955	0.0141421	0.183668	1.59995e-49	
rs_demo_0072	1	109071356	1.09232	0.0141421	0.274693	1.55218e-49	
rs_demo_0073	1	109072361	1.00245	0.0141421	0.316124	1.50758e-49	
rs_demo_0074	1	109073366	2.17611	0.0141421	0.0295471	3.92715e-49	
rs_demo_0075	1	109074371	1.627	0.0141421	0.103736	2.45237e-49	
rs_demo_0076	1	109075376	2.01957	0.0141421	0.0434279	4.94067e-49	
rs_demo_0077	1	109076381	0.777436	0.0141421	0.436901	1.60278e-49	
rs_demo_0078	1	109077386	1.26371	0.0141421	0.206333	2.37789e-49	
rs_demo_0079	1	109078391	1.65115	0.0141421	0.098708	4.21834e-49	
rs_demo_0080	1	109079396	0.652522	0.0141421	0.514064	1.66416e-49	
rs_demo_0081	1	109080402	1.36769	0.0141421	0.171409	3.26934e-49	
rs_demo_0082	1	109081407	0.208229	0.0141421	0.835051	1.50801e-49	
rs_demo_0083	1	109082412	0.474491	0.0141421	0.63515	1.61661e-49	
rs_demo_0084	1	109083417	0.429719	0.0141421	0.6674	1.60327e-49	
rs_demo_0085	1	109084422	-0.405252	0.0141421	0.685292	1.70346e-49	
rs_demo_0086	1	109085427	1.26086	0.0141421	0.207358	3.33423e-49	
rs_demo_0087	1	109086432	0.291132	0.0141421	0.770951	1.55517e-49	
rs_demo_0088	1	109087437	0.762343	0.0141421	0.445855	2.00941e-49	
rs_demo_0089	1	109088442	1.22205	0.0141421	0.22169	3.28585e-49	
rs_demo_0090	1	109089447	1.18099	0.0141421	0.237606	3.14148e-49	
rs_demo_0091	1	109090452	1.39438	0.0141421	0.163204	4.25998e-49	
rs_demo_0092	1	109091457	0.626128	0.0141421	0.531231	1.85061e-49	
rs_demo_0093	1	109092462	0.297814	0.0141421	0.765845	1.57577e-49	
rs_demo_0094	1	109093467	0.638593	0.0141421	0.523088	1.87774e-49	
rs_demo_0095	1	109094472	-0.971264	0.0141421	0.331417	2.5229e-49	
rs_demo_0096	1	109095477	-0.732835	0.0141421	0.463659	2.01535e-49	
rs_demo_0097	1	109096482	-0.609856	0.0141421	0.541957	1.84021e-49	
rs_demo_0098	1	109097487	-0.285551	0.0141421	0.775222	1.57127e-49	
rs_demo_0099	1	109098492	1.11055	0.0141421	0.266761	2.97217e-49	
rs_demo_0100	1	109099497	-0.195435	0.0141421	0.845052	1.5352e-49	
rs_demo_0101	1	109100502	0.507601	0.0141421	0.611733	1.56489e-49	
rs_demo_0102	1	109101507	2.18546	0.0141421	0.0288549	4.46893e-49	
rs_demo_0103	1	109102512	0.530557	0.0141421	0.595726	1.55469e-49	
rs_demo_0104	1	109103517	1.62507	0.0141421	0.104148	2.23624e-49	
rs_demo_0105	1	109104522	-0.0451448	0.0141421	0.963992	2.17363e-49	
rs_demo_0106	1	109105527	0.684183	0.0141421	0.49386	1.51087e-49	
rs_demo_0107	1	109106532	-0.0589677	0.0141421	0.952978	2.19872e-49	
rs_demo_0108	1	109107537	0.553814	0.0141421	0.579706	1.54438e-49	
rs_demo_0109	1	109108542	1.7354	0.0141421	0.0826706	2.50073e-49	
rs_demo_0110	1	109109547	-0.829431	0.0141421	0.40686	6.14906e-49	
rs_demo_0111	1	109110552	1.33581	0.0141421	0.18161	1.79019e-49	
rs_demo_0112	1	109111557	1.1435	0.0141421	0.25283	1.62215e-49	
rs_demo_0113	1	109112562	0.317088	0.0141421	0.751177	1.68572e-49	
rs_demo_0114	1	109113567	-0.527981	0.0141421	0.597513	3.80697e-49	
rs_demo_0115	1	109114572	0.998755	0.0141421	0.317913	1.54674e-49	
rs_demo_0116	1	109115577	0.407818	0.0141421	0.683407	1.61968e-49	
rs_demo_0117	1	109116582	1.18334	0.0141421	0.236673	1.64769e-49	
rs_demo_0118	1	109117587	0.989216	0.0141421	0.322557	1.54191e-49	
rs_demo_0119	1	109118592	2.59001	0.0141421	0.0095973	9.08958e-49	
rs_demo_0120	1	109119597	0.774966	0.0141421	0.43836	1.50434e-49	
rs_demo_0121	1	109120603	0.0234337	0.0141421	0.981304	2.05597e-49	
rs_demo_0122	1	109121608	1.26697	0.0141421	0.205167	1.71458e-49	
rs_demo_0123	1	109122613	1.35864	0.0141421	0.17426	1.80824e-49	
rs_demo_0124	1	109123618	2.56152	0.0141421	0.0104215	8.56027e-49	
rs_demo_0125	1	109124623	2.11706	0.0141421	0.0342546	4.00947e-49	
rs_demo_0126	1	109125628	1.73834	0.0141421	0.0821514	2.49229e-49	
rs_demo_0127	1	109126633	2.96916	0.0141421	0.00298625	2.08822e-48	
rs_demo_0128	1	109127638	0.472593	0.0141421	0.636504	1.58185e-49	
rs_demo_0129	1	109128643	1.21597	0.0141421	0.223996	1.67486e-49	
rs_demo_0130	1	109129648	1.17211	0.0141421	0.241154	1.6421e-49	
rs_demo_0131	1	109130653	2.01257	0.0141421	0.0441593	3.49804e-49	
rs_demo_0132	1	109131658	1.40532	0.0141421	0.159925	1.8749e-49	
rs_demo_0133	1	109132663	3.89169	0.0141421	9.95871e-05	3.0795e-47	
rs_demo_0134	1	109133668	3.62748	0.0141421	0.00028627	1.30113e-47	
rs_demo_0135	1	109134673	3.12035	0.0141421	0.00180648	3.07754e-48	
rs_demo_0136	1	109135678	4.5024	0.0141421	6.72526e-06	2.99141e-46	
rs_demo_0137	1	109136683	6.99002	0.0141421	2.76631e-12	1.92327e-40	
rs_demo_0138	1	109137688	8.96278	0.0141421	3.21185e-19	7.48291e-34	
rs_demo_0139	1	109138693	11.9316	0.0141421	8.34998e-33	9.01059e-21	
rs_demo_0140	1	109139698	15.1115	0.0141421	1.42915e-51	0.0080008	
rs_demo_0141	1	109140703	15.4404	0.0141421	9.20764e-54	0.991999	L2
rs_demo_0142	1	109141708	11.2081	0.0141421	3.82119e-29	2.86205e-24	
rs_demo_0143	1	109142713	7.8028	0.0141421	6.11186e-15	6.43634e-38	
rs_demo_0144	1	109143718	8.39237	0.0141421	4.82126e-17	6.5091e-36	
rs_demo_0145	1	109144723	5.86356	0.0141421	4.54519e-09	2.04989e-43	
rs_demo_0146	1	109145728	5.10262	0.0141421	3.35566e-07	4.19418e-45	
rs_demo_0147	1	109146733	3.97906	0.0141421	6.92188e-05	4.16916e-47	
rs_demo_0148	1	109147738	3.70891	0.0141421	0.00020821	1.68124e-47	
rs_demo_0149	1	109148743	4.32252	0.0141421	1.54371e-05	1.46809e-46	
rs_demo_0150	1	109149748	2.93906	0.0141421	0.00329226	1.96395e-48	
rs_demo_0151	1	109150753	2.24375	0.0141421	0.0248485	4.92697e-49	
rs_demo_0152	1	109151758	1.06303	0.0141421	0.287768	1.57591e-49	
rs_demo_0153	1	109152763	0.181041	0.0141421	0.856336	1.82224e-49	
rs_demo_0154	1	109153768	1.17505	0.0141421	0.239976	1.64354e-49	
rs_demo_0155	1	109154773	1.45208	0.0141421	0.14648	1.93427e-49	
rs_demo_0156	1	109155778	3.1494	0.0141421	0.00163622	3.28135e-48	
rs_demo_0157	1	109156783	1.41222	0.0141421	0.157884	1.87737e-49	
rs_demo_0158	1	109157788	2.18508	0.0141421	0.0288832	4.46394e-49	
rs_demo_0159	1	109158793	0.63935	0.0141421	0.522595	1.5197e-49	
rs_demo_0160	1	109159798	-0.0972506	0.0141421	0.922527	2.28356e-49	
rs_demo_0161	1	109160804	0.0818144	0.0141421	0.934794	1.95837e-49	
rs_demo_0162	1	109161809	0.289096	0.0141421	0.772508	1.71349e-49	
rs_demo_0163	1	109162814	3.1167	0.0141421	0.00182899	3.0091e-48	
rs_demo_0164	1	109163819	0.145977	0.0141421	0.883939	1.87232e-49	
rs_demo_0165	1	109164824	1.78916	0.0141421	0.0735894	2.63665e-49	
rs_demo_0166	1	109165829	0.0343838	0.0141421	0.972571	2.03663e-49	
rs_demo_0167	1	109166834	1.8582	0.0141421	0.0631407	2.85381e-49	
rs_demo_0168	1	109167839	1.30303	0.0141421	0.192565	1.75067e-49	
rs_demo_0169	1	109168844	0.7546	0.0141421	0.450489	1.50473e-49	
rs_demo_0170	1	109169849	0.865005	0.0141421	0.387036	1.51088e-49	
rs_demo_0171	1	109170854	0.246603	0.0141421	0.805216	1.75397e-49	
rs_demo_0172	1	109171859	1.34396	0.0141421	0.178961	1.79624e-49	
rs_demo_0173	1	109172864	0.440105	0.0141421	0.659861	1.59982e-49	
rs_demo_0174	1	109173869	-0.332759	0.0141421	0.739317	2.94985e-49	
rs_demo_0175	1	109174874	-0.386883	0.0141421	0.698843	3.15927e-49	
rs_demo_0176	1	109175879	1.06221	0.0141421	0.288141	1.57264e-49	
rs_demo_0177	1	109176884	2.46756	0.0141421	0.0136036	7.16922e-49	
rs_demo_0178	1	109177889	1.04755	0.0141421	0.294848	1.56448e-49	
rs_demo_0179	1	109178894	0.768182	0.0141421	0.442379	1.50446e-49	
rs_demo_0180	1	109179899	1.17206	0.0141421	0.241173	1.63551e-49	
rs_demo_0181	1	109180904	2.19177	0.0141421	0.0283964	4.46859e-49	
rs_demo_0182	1	109181909	1.10477	0.0141421	0.269259	1.59256e-49	
rs_demo_0183	1	109182914	0.474159	0.0141421	0.635386	1.58492e-49	
rs_demo_0184	1	109183919	1.99113	0.0141421	0.0464659	3.34867e-49	
rs_demo_0185	1	109184924	1.31341	0.0141421	0.189045	1.75536e-49	
rs_demo_0186	1	109185929	2.42026	0.0141421	0.0155096	6.5394e-49	
rs_demo_0187	1	109186934	1.06761	0.0141421	0.285696	1.57344e-49	
rs_demo_0188	1	109187939	-0.340191	0.0141421	0.733713	2.99672e-49	
rs_demo_0189	1	109188944	-0.48396	0.0141421	0.628414	3.61754e-49	
rs_demo_0190	1	109189949	2.53506	0.0141421	0.0112427	8.10953e-49	
rs_demo_0191	1	109190954	2.60775	0.0141421	0.00911399	9.34843e-49	
rs_demo_0192	1	109191959	0.704526	0.0141421	0.481105	1.50912e-49	
rs_demo_0193	1	109192964	0.500825	0.0141421	0.616494	1.56988e-49	
rs_demo_0194	1	109193969	2.34543	0.0141421	0.0190051	5.77372e-49	
rs_demo_0195	1	109194974	-0.22308	0.0141421	0.823474	2.60565e-49	
rs_demo_0196	1	109195979	-0.0107774	0.0141421	0.991401	2.114e-49	
rs_demo_0197	1	109196984	1.52726	0.0141421	0.126696	2.04938e-49	
rs_demo_0198	1	109197989	0.489321	0.0141421	0.624615	1.57374e-49	
rs_demo_0199	1	109198994	0.878795	0.0141421	0.379512	1.51312e-49	
rs_demo_0200	1	109200000	0.720468	0.0141421	0.471237	1.50684e-49	
"""})
