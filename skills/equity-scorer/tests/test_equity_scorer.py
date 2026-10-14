"""
test_equity_scorer.py — Automated test suite for Equity Scorer

Run with: pytest skills/equity-scorer/tests/test_equity_scorer.py -v

Uses the FIXED demo data (demo_populations.vcf + demo_population_map.csv)
so that all assertions are deterministic and reproducible.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from equity_scorer import (
    GLOBAL_PROPORTIONS,
    DEFAULT_WEIGHTS,
    parse_vcf,
    load_population_map,
    parse_ancestry_csv,
    compute_allele_frequencies,
    compute_heterozygosity,
    compute_pairwise_fst,
    compute_pca,
    compute_representation_index,
    compute_heterozygosity_balance,
    compute_fst_coverage,
    compute_geographic_spread,
    compute_heim_score,
)

PROJ = Path(__file__).resolve().parents[3]
DEMO_VCF = PROJ / "examples" / "demo_populations.vcf"
DEMO_MAP = PROJ / "examples" / "demo_population_map.csv"
DEMO_CSV = PROJ / "examples" / "sample_ancestry.csv"


# ── VCF Parsing ───────────────────────────────────────────────────────────────

def test_parse_vcf_returns_correct_shape():
    samples, variant_ids, geno = parse_vcf(DEMO_VCF)
    assert len(samples) == 50
    assert geno.shape[0] == 50  # n_samples
    assert geno.shape[1] == len(variant_ids)
    assert geno.shape[1] > 0


def test_parse_vcf_genotype_range():
    _, _, geno = parse_vcf(DEMO_VCF)
    unique_vals = set(np.unique(geno))
    assert unique_vals.issubset({-1, 0, 1, 2}), (
        f"Unexpected genotype values: {unique_vals}"
    )


def test_parse_vcf_sample_names():
    samples, _, _ = parse_vcf(DEMO_VCF)
    # Expect samples with population prefixes
    prefixes = {s.split("_")[0] for s in samples}
    assert "AFR" in prefixes
    assert "EUR" in prefixes


# ── Population Map ─────────────────────────────────────────────────────────────

def test_load_population_map_from_csv():
    samples, _, _ = parse_vcf(DEMO_VCF)
    pop_map = load_population_map(DEMO_MAP, samples)
    assert len(pop_map) == len(samples)
    pops = set(pop_map.values())
    assert "AFR" in pops
    assert "EUR" in pops


def test_load_population_map_infers_from_ids():
    samples = ["AFR_001", "EUR_002", "EAS_003"]
    pop_map = load_population_map(None, samples)
    assert pop_map["AFR_001"] == "AFR"
    assert pop_map["EUR_002"] == "EUR"
    assert pop_map["EAS_003"] == "EAS"


def test_parse_ancestry_csv():
    df = parse_ancestry_csv(DEMO_CSV)
    assert "population" in df.columns
    assert "sample_id" in df.columns
    assert len(df) > 0


# ── Allele Frequencies ─────────────────────────────────────────────────────────

def _load_demo_data():
    """Load demo VCF + map for reuse."""
    samples, vids, geno = parse_vcf(DEMO_VCF)
    pop_map = load_population_map(DEMO_MAP, samples)
    pop_indices = {}
    for i, s in enumerate(samples):
        pop = pop_map[s]
        pop_indices.setdefault(pop, []).append(i)
    return samples, vids, geno, pop_indices


def test_allele_frequencies_range():
    _, _, geno, pop_idx = _load_demo_data()
    afs = compute_allele_frequencies(geno, pop_idx)
    for pop, freqs in afs.items():
        assert np.all(freqs >= 0.0), f"{pop} has negative AF"
        assert np.all(freqs <= 1.0), f"{pop} has AF > 1"


# ── Heterozygosity ─────────────────────────────────────────────────────────────

def test_heterozygosity_values():
    _, _, geno, pop_idx = _load_demo_data()
    obs, exp, _, _ = compute_heterozygosity(geno, pop_idx)
    for pop in obs:
        assert 0.0 <= obs[pop] <= 1.0, f"Obs het out of range for {pop}"
        assert 0.0 <= exp[pop] <= 0.5, f"Exp het out of range for {pop}"


def test_heterozygosity_all_pops_present():
    _, _, geno, pop_idx = _load_demo_data()
    obs, exp, _, _ = compute_heterozygosity(geno, pop_idx)
    assert set(obs.keys()) == set(pop_idx.keys())
    assert set(exp.keys()) == set(pop_idx.keys())


# ── FST ────────────────────────────────────────────────────────────────────────

def test_fst_symmetric():
    _, _, geno, pop_idx = _load_demo_data()
    fst_df, _ = compute_pairwise_fst(geno, pop_idx)
    assert np.allclose(fst_df.values, fst_df.values.T), "FST matrix not symmetric"


def test_fst_non_negative():
    _, _, geno, pop_idx = _load_demo_data()
    fst_df, _ = compute_pairwise_fst(geno, pop_idx)
    assert np.all(fst_df.values >= 0.0), "Negative FST values found"


def test_fst_diagonal_zero():
    _, _, geno, pop_idx = _load_demo_data()
    fst_df, _ = compute_pairwise_fst(geno, pop_idx)
    assert np.allclose(np.diag(fst_df.values), 0.0), "Diagonal FST not zero"


# ── PCA ────────────────────────────────────────────────────────────────────────

def test_pca_shape():
    _, _, geno, _ = _load_demo_data()
    coords, var_ratio = compute_pca(geno, n_components=5)
    assert coords.shape[0] == 50  # n_samples
    assert coords.shape[1] == 5
    assert len(var_ratio) == 5


def test_pca_variance_sums_to_at_most_one():
    _, _, geno, _ = _load_demo_data()
    _, var_ratio = compute_pca(geno, n_components=10)
    assert np.sum(var_ratio) <= 1.0 + 1e-6


# ── HEIM Score Components ──────────────────────────────────────────────────────

def test_representation_index_range():
    pop_counts = {"AFR": 8, "AMR": 5, "EAS": 7, "EUR": 22, "SAS": 8}
    result = compute_representation_index(pop_counts)
    ri = result["representation_index"]
    assert 0.0 <= ri <= 1.0


def test_representation_index_perfect():
    """Perfect match to global proportions should give RI close to 1."""
    total = 1000
    pop_counts = {k: int(v * total) for k, v in GLOBAL_PROPORTIONS.items()}
    result = compute_representation_index(pop_counts)
    ri = result["representation_index"]
    assert ri > 0.95


def test_heterozygosity_balance_range():
    het_values = {"AFR": 0.35, "EUR": 0.28, "EAS": 0.25}
    hb = compute_heterozygosity_balance(het_values)
    assert 0.0 <= hb <= 1.0


def test_heterozygosity_balance_empty():
    assert compute_heterozygosity_balance({}) == 0.0


def test_fst_coverage():
    assert compute_fst_coverage(5, 10) == 1.0  # 5 pops → 10 pairs → full
    assert compute_fst_coverage(5, 5) == 0.5
    assert compute_fst_coverage(0, 0) == 0.0


def test_geographic_spread():
    assert compute_geographic_spread({"AFR", "EUR", "EAS", "SAS", "AMR"}) == 5 / 7
    assert compute_geographic_spread(set()) == 0.0
    assert compute_geographic_spread(set(GLOBAL_PROPORTIONS.keys())) == 1.0


# ── HEIM Composite Score ──────────────────────────────────────────────────────

def test_heim_score_structure():
    pop_counts = {"AFR": 8, "AMR": 5, "EAS": 7, "EUR": 22, "SAS": 8}
    het_values = {"AFR": 0.35, "AMR": 0.30, "EAS": 0.25, "EUR": 0.28, "SAS": 0.27}
    result = compute_heim_score(pop_counts, het_values, n_pairwise_fst=10)
    assert "heim_score" in result
    assert "rating" in result
    assert "components" in result
    assert 0 <= result["heim_score"] <= 100
    assert result["rating"] in ("Excellent", "Good", "Fair", "Poor", "Critical")


def test_heim_score_range():
    pop_counts = {"AFR": 8, "AMR": 5, "EAS": 7, "EUR": 22, "SAS": 8}
    het_values = {"AFR": 0.35, "AMR": 0.30, "EAS": 0.25, "EUR": 0.28, "SAS": 0.27}
    result = compute_heim_score(pop_counts, het_values, n_pairwise_fst=10)
    assert 0 <= result["heim_score"] <= 100


def test_heim_score_custom_weights():
    pop_counts = {"AFR": 10, "EUR": 10}
    het_values = {"AFR": 0.3, "EUR": 0.3}
    r1 = compute_heim_score(pop_counts, het_values, 1, weights=(1.0, 0, 0, 0))
    r2 = compute_heim_score(pop_counts, het_values, 1, weights=(0, 1.0, 0, 0))
    # Different weights should give different scores
    assert r1["heim_score"] != r2["heim_score"]


def test_heim_score_unnormalized_weights_bounded():
    """U-2/F-27: Unnormalized weights (1,1,1,1) must still produce score in [0, 100]."""
    pop_counts = {"AFR": 10, "AMR": 10, "EAS": 10, "EUR": 10, "SAS": 10}
    het_values = {"AFR": 0.35, "AMR": 0.30, "EAS": 0.25, "EUR": 0.28, "SAS": 0.27}
    result = compute_heim_score(pop_counts, het_values, 10, weights=(1, 1, 1, 1))
    assert 0 <= result["heim_score"] <= 100, (
        f"HEIM score {result['heim_score']} out of [0, 100] with weights (1,1,1,1)"
    )


def test_heim_score_negative_weight_raises():
    """Negative weights must be rejected."""
    import pytest
    pop_counts = {"AFR": 10, "EUR": 10}
    het_values = {"AFR": 0.3, "EUR": 0.3}
    with pytest.raises(ValueError, match="non-negative"):
        compute_heim_score(pop_counts, het_values, 1, weights=(-1, 0.5, 0.5, 0))


def test_heim_score_zero_weights_raises():
    """All-zero weights must be rejected."""
    import pytest
    pop_counts = {"AFR": 10, "EUR": 10}
    het_values = {"AFR": 0.3, "EUR": 0.3}
    with pytest.raises(ValueError, match="zero"):
        compute_heim_score(pop_counts, het_values, 1, weights=(0, 0, 0, 0))


# ── Edge Case Handling (clawbio_bench eq_01-eq_15) ────────────────────────────

def test_fst_identical_populations():
    """eq_02: Identical populations should produce FST = 0."""
    geno = np.array([[0, 1, 2, 0, 1],
                     [0, 1, 2, 0, 1],
                     [0, 1, 2, 0, 1],
                     [0, 1, 2, 0, 1]])
    pop_idx = {"POP_A": [0, 1], "POP_B": [2, 3]}
    fst_df, fst_dict = compute_pairwise_fst(geno, pop_idx)
    val = fst_dict[("POP_A", "POP_B")]
    assert val == 0.0 or np.isnan(val), f"Identical pops should give FST=0, got {val}"


def test_fst_monomorphic_sites():
    """eq_03: Monomorphic sites (all same genotype) should produce NaN FST, not crash."""
    geno = np.array([[0, 0, 0],
                     [0, 0, 0],
                     [0, 0, 0],
                     [0, 0, 0]])
    pop_idx = {"POP_A": [0, 1], "POP_B": [2, 3]}
    fst_df, fst_dict = compute_pairwise_fst(geno, pop_idx)
    val = fst_dict[("POP_A", "POP_B")]
    assert np.isnan(val) or val == 0.0, f"Monomorphic should give NaN or 0, got {val}"


def test_fst_single_sample_per_pop():
    """eq_04: n=1 per population should not crash."""
    geno = np.array([[0, 1, 2], [2, 1, 0]])
    pop_idx = {"POP_A": [0], "POP_B": [1]}
    fst_df, fst_dict = compute_pairwise_fst(geno, pop_idx)
    val = fst_dict[("POP_A", "POP_B")]
    assert isinstance(val, (float, np.floating)), f"Single sample FST should be float, got {type(val)}"


def test_het_all_heterozygous():
    """eq_05: All-het genotypes should produce bounded heterozygosity."""
    geno = np.array([[1, 1, 1, 1],
                     [1, 1, 1, 1],
                     [1, 1, 1, 1],
                     [1, 1, 1, 1]])
    pop_idx = {"POP_A": [0, 1], "POP_B": [2, 3]}
    obs, exp, _, _ = compute_heterozygosity(geno, pop_idx)
    for pop in obs:
        assert 0.0 <= obs[pop] <= 1.0, f"Obs het out of range for {pop}: {obs[pop]}"


def test_het_all_homozygous():
    """eq_06: All-hom genotypes should produce bounded heterozygosity."""
    geno = np.array([[0, 0, 0, 0],
                     [2, 2, 2, 2],
                     [0, 0, 0, 0],
                     [2, 2, 2, 2]])
    pop_idx = {"POP_A": [0, 1], "POP_B": [2, 3]}
    obs, exp, _, _ = compute_heterozygosity(geno, pop_idx)
    for pop in obs:
        assert obs[pop] == 0.0, f"All-hom should give obs_het=0, got {obs[pop]}"


def test_heim_single_population():
    """eq_08: Single population should produce bounded HEIM, not crash."""
    pop_counts = {"EUR": 50}
    het_values = {"EUR": 0.28}
    result = compute_heim_score(pop_counts, het_values, 0)
    assert 0 <= result["heim_score"] <= 100


def test_heim_balanced_five_pops():
    """eq_07: Balanced 5-superpop dataset should produce high HEIM."""
    pop_counts = {"AFR": 10, "AMR": 10, "EAS": 10, "EUR": 10, "SAS": 10}
    het_values = {"AFR": 0.35, "AMR": 0.30, "EAS": 0.25, "EUR": 0.28, "SAS": 0.27}
    result = compute_heim_score(pop_counts, het_values, 10)
    assert result["heim_score"] > 50, f"Balanced dataset should score >50, got {result['heim_score']}"
    assert 0 <= result["heim_score"] <= 100


def test_fst_three_populations():
    """eq_15: Three populations should compute all 3 pairwise FST values."""
    geno = np.array([[0, 1, 0, 1, 0],
                     [0, 1, 0, 1, 0],
                     [2, 2, 2, 2, 2],
                     [2, 2, 2, 2, 2],
                     [0, 0, 1, 1, 0],
                     [0, 0, 1, 1, 0]])
    pop_idx = {"POP_A": [0, 1], "POP_B": [2, 3], "POP_C": [4, 5]}
    fst_df, fst_dict = compute_pairwise_fst(geno, pop_idx)
    assert len(fst_dict) == 3, f"3 pops should give 3 pairs, got {len(fst_dict)}"
    for key, val in fst_dict.items():
        assert isinstance(val, (float, np.floating)), f"FST for {key} should be float"


# ── FST Label Honesty ─────────────────────────────────────────────────────────

def test_fst_report_label_matches_formula():
    """C-06: The report must label FST as Nei's GST, not Hudson FST."""
    import equity_scorer as eq
    import inspect
    # Check the report template string does not say "Hudson FST"
    source = inspect.getsource(eq)
    assert "Hudson FST" not in source, (
        "Honesty failure: code computes Nei's GST but labels it Hudson FST"
    )
    # The docstring should say Nei's GST
    assert "Nei" in eq.compute_pairwise_fst.__doc__


# ── End-to-End with Demo VCF ──────────────────────────────────────────────────

def test_end_to_end_demo_vcf():
    """Full pipeline on demo VCF produces a valid HEIM score."""
    samples, vids, geno = parse_vcf(DEMO_VCF)
    pop_map = load_population_map(DEMO_MAP, samples)
    pop_indices = {}
    pop_counts = {}
    for i, s in enumerate(samples):
        pop = pop_map[s]
        pop_indices.setdefault(pop, []).append(i)
        pop_counts[pop] = pop_counts.get(pop, 0) + 1

    obs_het, _, _, _ = compute_heterozygosity(geno, pop_indices)
    _, fst_dict = compute_pairwise_fst(geno, pop_indices)

    result = compute_heim_score(pop_counts, obs_het, len(fst_dict))
    assert result["heim_score"] > 0
    assert result["n_samples"] == 50
    assert result["n_populations"] == len(pop_indices)


# ── HEIM inflation regressions ────────────────────────────────────────────────

def test_het_balance_ignores_population_with_no_calls():
    """A NaN het (population with no valid calls) must not read as perfect balance."""
    assert compute_heterozygosity_balance({"AFR": 0.0, "EUR": float("nan")}) == 0.0
    assert compute_heterozygosity_balance({"EUR": float("nan")}) == 0.0


@pytest.mark.parametrize("bad", [float("inf"), float("nan")])
def test_heim_score_rejects_non_finite_weights(bad):
    with pytest.raises(ValueError):
        compute_heim_score({"AFR": 5, "EUR": 5}, {"AFR": 0.3, "EUR": 0.3}, 1, weights=(0, 0, 0, bad))


def test_fst_coverage_excludes_undefined_pairs(tmp_path):
    """A pair whose FST is NaN (no data in one population) was not computed."""
    from equity_scorer import run_vcf_pipeline

    vcf = tmp_path / "in.vcf"
    rows = "".join(f"22\t{p}\trs{p}\tT\tA\t.\tPASS\t.\tGT\t0/0\t1/1\t./.\t./.\n" for p in (100, 200, 300))
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA1\tA2\tB1\tB2\n" + rows
    )
    pop_map = tmp_path / "map.csv"
    pop_map.write_text(
        "sample_id,population,superpopulation,country\n"
        "A1,AFR,African,Kenya\nA2,AFR,African,Kenya\nB1,EUR,European,UK\nB2,EUR,European,UK\n"
    )
    result = run_vcf_pipeline(vcf, pop_map, tmp_path / "out", (0, 0, 1, 0))
    assert result["components"]["fst_coverage"] == 0.0
