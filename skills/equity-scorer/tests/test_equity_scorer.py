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
from inline_snapshot import snapshot

from clawbio.common import demo_runner

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


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/equity-scorer/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("equity-scorer", tmp_path) == snapshot({'report.md': """\
# HEIM Equity Report

**Date**: <TIMESTAMP> UTC
**Input**: `demo_populations.vcf`
**Checksum (SHA-256)**: `f847bb32a384508cd062bda6caac0979297dbef3fe14d7e645c2b21f1f9f7c68`
**Samples**: 50
**Populations**: 5
**Variants analysed**: 500

---

## HEIM Equity Score: 76.2/100 (Good)

![HEIM Gauge](figures/heim_gauge.png)

### Score Breakdown

| Component | Value | Weight | Description |
|-----------|-------|--------|-------------|
| Representation Index | 0.720 | 0.35 | Match to global population proportions |
| Heterozygosity Balance | 0.667 | 0.25 | Genetic diversity relative to theoretical max |
| FST Coverage | 1.000 | 0.2 | Fraction of pairwise comparisons computed |
| Geographic Spread | 0.714 | 0.2 | Continental groups represented (out of 7) |

### Key Findings

- **Most represented**: EUR (44.0%, 2.8x global proportion)
- **Least represented**: AMR (10.0%, 0.8x global proportion)
- **Mean observed heterozygosity**: 0.3335 (highest: AFR at 0.3543)

---

## Population Distribution

| Population | Count | Sample % | Global % | Ratio | Obs Het | Exp Het |
|------------|-------|-----------|-----------|-------|---------|---------|
| AFR | 8 | 16.0% | 17.0% | 0.94x | 0.3543 | 0.3338 |
| AMR | 5 | 10.0% | 13.0% | 0.77x | 0.3472 | 0.3064 |
| EAS | 7 | 14.0% | 22.0% | 0.64x | 0.3014 | 0.2788 |
| EUR | 22 | 44.0% | 16.0% | 2.75x | 0.3139 | 0.3075 |
| SAS | 8 | 16.0% | 26.0% | 0.62x | 0.3508 | 0.3220 |

![Ancestry Distribution](figures/ancestry_bar.png)

## Heterozygosity

![Heterozygosity](figures/heterozygosity.png)

## Pairwise FST

| Comparison | Nei's GST (FST) |
|------------|-----------------|
| AFR vs AMR | 0.0773 |
| AFR vs EAS | 0.1011 |
| AFR vs EUR | 0.0590 |
| AFR vs SAS | 0.0619 |
| AMR vs EAS | 0.0898 |
| AMR vs EUR | 0.0425 |
| AMR vs SAS | 0.0692 |
| EAS vs EUR | 0.0439 |
| EAS vs SAS | 0.0782 |
| EUR vs SAS | 0.0467 |

![FST Heatmap](figures/fst_heatmap.png)


## Principal Component Analysis

![PCA](figures/pca_plot.png)

- PC1 explains 7.6% of variance
- PC2 explains 4.9% of variance
- Top 5 components explain 23.0% of total variance


---

## Methods

- **Tool**: ClawBio Equity Scorer v0.1.0
- **HEIM framework**: Health Equity Index for Minorities (Corpas, 2026)
- **Heterozygosity**: Observed = proportion of heterozygous genotypes per site, averaged across 500 variants. Expected = 2pq from population allele frequencies.
- **FST**: Nei's GST (HT-HS)/HT, ratio of averages across sites. Values floored at 0.
- **PCA**: scikit-learn PCA on mean-imputed genotype matrix (0/1/2 encoding).
- **Global reference**: Approximate continental proportions from the 1000 Genomes Project.

## Reproducibility

```bash
# Re-run this analysis
python equity_scorer.py --input demo_populations.vcf --output demo_output
```

**Input checksum**: `f847bb32a384508cd062bda6caac0979297dbef3fe14d7e645c2b21f1f9f7c68`

## References

- Corpas, M. (2026). ClawBio. https://github.com/ClawBio/ClawBio
- Nei, M. (1973). Analysis of Gene Diversity in Subdivided Populations. PNAS, 70(12), 3321-3323.
- Hudson, R.R., Slatkin, M. & Maddison, W.P. (1992). Estimation of levels of gene flow from DNA sequence data. Genetics, 132(2), 583-589.
- The 1000 Genomes Project Consortium (2015). A global reference for human genetic variation. Nature, 526, 68-74.

---

## Disclaimer

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
""", 'result.json': """\
{
  "skill": "equity-scorer",
  "version": "0.2.0",
  "completed_at": "<TIMESTAMP>",
  "input_checksum": "sha256:f847bb32a384508cd062bda6caac0979297dbef3fe14d7e645c2b21f1f9f7c68",
  "datasets": {},
  "summary": {
    "heim_score": 76.2,
    "rating": "Good",
    "n_samples": 50,
    "n_populations": 5
  },
  "data": {
    "heim_score": 76.2,
    "rating": "Good",
    "components": {
      "representation_index": 0.72,
      "heterozygosity_balance": 0.667,
      "fst_coverage": 1.0,
      "geographic_spread": 0.714
    },
    "weights": {
      "w1": 0.35,
      "w2": 0.25,
      "w3": 0.2,
      "w4": 0.2
    },
    "n_samples": 50,
    "n_populations": 5,
    "population_counts": {
      "AFR": 8,
      "AMR": 5,
      "EAS": 7,
      "EUR": 22,
      "SAS": 8
    },
    "unknown_fraction": 0.0,
    "representation_warning": null,
    "het_source": "computed",
    "fst_estimator": "Nei_GST"
  }
}\
""", 'tables/fst_matrix.csv': """\
,AFR,AMR,EAS,EUR,SAS
AFR,0.0,0.07729026973161958,0.10109176634894292,0.059007353052296485,0.061931475314574345
AMR,0.07729026973161958,0.0,0.08983741986855488,0.042500319633088424,0.06924789557155057
EAS,0.10109176634894292,0.08983741986855488,0.0,0.043865792382887424,0.0781626295766291
EUR,0.059007353052296485,0.042500319633088424,0.043865792382887424,0.0,0.04666476240435645
SAS,0.061931475314574345,0.06924789557155057,0.0781626295766291,0.04666476240435645,0.0
""", 'tables/heim_score.json': """\
{
  "heim_score": 76.2,
  "rating": "Good",
  "components": {
    "representation_index": 0.72,
    "heterozygosity_balance": 0.667,
    "fst_coverage": 1.0,
    "geographic_spread": 0.714
  },
  "weights": {
    "w1": 0.35,
    "w2": 0.25,
    "w3": 0.2,
    "w4": 0.2
  },
  "n_samples": 50,
  "n_populations": 5,
  "population_counts": {
    "AFR": 8,
    "AMR": 5,
    "EAS": 7,
    "EUR": 22,
    "SAS": 8
  },
  "unknown_fraction": 0.0,
  "representation_warning": null,
  "het_source": "computed",
  "fst_estimator": "Nei_GST"
}\
""", 'tables/heterozygosity.csv': """\
population,observed_het,expected_het
AFR,0.35425,0.333796875
AMR,0.34720000000000006,0.30639999999999995
EAS,0.30142857142857143,0.2787959183673469
EUR,0.31390909090909086,0.30745247933884295
SAS,0.35075,0.322015625
""", 'tables/population_summary.csv': """\
population,count,proportion,obs_het,exp_het
AFR,8,0.16,0.35425,0.333796875
AMR,5,0.1,0.34720000000000006,0.30639999999999995
EAS,7,0.14,0.30142857142857143,0.2787959183673469
EUR,22,0.44,0.31390909090909086,0.30745247933884295
SAS,8,0.16,0.35075,0.322015625
"""})
