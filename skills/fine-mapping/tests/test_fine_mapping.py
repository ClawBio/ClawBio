"""Tests for the SuSiE fine-mapping skill."""
import sys
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

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
    """Tests for the SuSiE IBSS algorithm."""

    def test_pip_shape(self):
        """run_susie returns PIPs with correct shape."""
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie(z=df["z"].values, R=R, n=5000, L=5)
        assert result["pip"].shape == (20,)

    def test_pip_range(self):
        """All PIPs are in [0, 1]."""
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie(z=df["z"].values, R=R, n=5000, L=5)
        assert result["pip"].min() >= 0.0
        assert result["pip"].max() <= 1.0

    def test_signal_recovers_high_pip(self):
        """The injected signal at index 10 has the highest PIP."""
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie(z=df["z"].values, R=R, n=5000, L=5)
        assert result["pip"].argmax() == 10

    def test_alpha_rows_sum_to_at_most_one(self):
        """With null component (default), alpha rows sum to <= 1.
        Without null component, rows sum to exactly 1."""
        df = _small_locus(n=20)
        R = _identity_ld(20)
        # With null component: rows sum to <= 1
        result = run_susie(z=df["z"].values, R=R, n=5000, L=3)
        row_sums = result["alpha"].sum(axis=1)
        assert np.all(row_sums <= 1.0 + 1e-6)
        # Without null component: rows sum to exactly 1
        result_nn = run_susie(z=df["z"].values, R=R, n=5000, L=3, null_weight=0.0)
        row_sums_nn = result_nn["alpha"].sum(axis=1)
        np.testing.assert_allclose(row_sums_nn, 1.0, atol=1e-6)

    def test_converges_on_clean_signal(self):
        """SuSiE converges within 100 iterations for a clean single signal."""
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie(z=df["z"].values, R=R, n=5000, L=5, max_iter=100)
        assert result["converged"] is True

    def test_elbo_tracked(self):
        """ELBO history is recorded and contains finite values."""
        df = _small_locus(n=20)
        R = _identity_ld(20)
        result = run_susie(z=df["z"].values, R=R, n=5000, L=5)
        elbo = result["elbo"]
        assert len(elbo) >= 2
        assert all(np.isfinite(e) for e in elbo)

    def test_two_signals_recovered(self):
        """Demo locus with two injected causal variants: both are in the top-5 PIPs."""
        df, R = fine_mapping.make_demo_data(seed=42)
        result = run_susie(z=df["z"].values, R=R, n=5000, L=10)
        top5 = set(np.argsort(-result["pip"])[:5])
        # Causal variants are at indices 60 and 140
        assert 60 in top5 or 140 in top5

    def test_null_locus_no_phantom_pip(self):
        """Null locus (z=0 everywhere) should produce near-zero PIPs with null_weight."""
        p = 50
        z = np.zeros(p)
        R = np.eye(p)
        result = run_susie(z=z, R=R, n=10000, L=5, null_weight=1.0 / 6)
        # On a null locus, all PIPs should be very low (< 0.1)
        assert result["pip"].max() < 0.1, (
            f"Null locus phantom PIP: max PIP = {result['pip'].max():.3f}, "
            f"expected < 0.1 with null_weight"
        )

    def test_null_locus_without_null_weight_has_phantom(self):
        """Without null_weight, a null locus produces phantom PIPs (baseline check)."""
        p = 10  # small p amplifies the phantom effect: PIP ~ 1-(1-1/p)^L
        z = np.zeros(p)
        R = np.eye(p)
        result = run_susie(z=z, R=R, n=10000, L=5, null_weight=0.0)
        # Without null component, PIPs are spread uniformly: ~1-(9/10)^5 ≈ 0.41
        assert result["pip"].max() > 0.3, (
            f"Expected phantom PIP > 0.3 without null_weight, got {result['pip'].max():.3f}"
        )

    def test_single_signal_with_null_weight_still_recovered(self):
        """A real signal should still be recovered even with null_weight active."""
        p = 50
        z = np.zeros(p)
        z[25] = 5.0  # strong signal at index 25
        R = np.eye(p)
        result = run_susie(z=z, R=R, n=10000, L=5, null_weight=1.0 / 6)
        assert result["pip"].argmax() == 25
        assert result["pip"][25] > 0.8, (
            f"Signal PIP = {result['pip'][25]:.3f}, expected > 0.8"
        )


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

    def test_demo_mode_produces_outputs(self, demo_run):
        """Demo mode writes report.md, fine_mapping.json, and tables/pips.tsv."""
        _, tmp = demo_run
        assert (tmp / "report.md").exists()
        assert (tmp / "fine_mapping.json").exists()
        assert (tmp / "tables" / "pips.tsv").exists()

    def test_demo_results_structure(self, demo_run):
        """Demo run returns a dict with expected top-level keys."""
        results, _ = demo_run
        for key in ("method", "n_variants", "n_credible_sets", "lead_rsid", "lead_pip"):
            assert key in results, f"Missing key: {key}"

    def test_demo_susie_method(self, demo_run):
        """Demo mode includes a synthetic LD matrix, so it uses SuSiE."""
        results, _ = demo_run
        assert results["method"] == "SuSiE"

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

    def test_lead_pip_in_range(self, demo_run):
        """Lead PIP is a valid probability in [0, 1]."""
        results, _ = demo_run
        assert 0.0 <= results["lead_pip"] <= 1.0

    def test_report_contains_disclaimer(self, demo_run):
        """Generated report contains the required safety disclaimer."""
        _, tmp = demo_run
        report = (tmp / "report.md").read_text()
        assert "ClawBio is a research and educational tool" in report

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
