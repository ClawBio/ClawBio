"""Tests for proteomics-de — differential expression pipeline."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from inline_snapshot import snapshot

from clawbio.common import demo_runner

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from proteomics_de import (
    _sep_for,
    parse_contrast,
    load_metadata,
    align_and_validate,
    log2_transform,
    down_shift_imputation,
    run_pca,
    _bh_fdr,
    run_differential_expression,
)


# ---------------------------------------------------------------------------
# Helpers: synthetic data
# ---------------------------------------------------------------------------

def _make_matrix(n_proteins=10, n_samples=6):
    """Create a synthetic proteomics matrix with some NaN."""
    np.random.seed(42)
    data = np.random.lognormal(mean=20, sigma=1, size=(n_proteins, n_samples))
    # Introduce some missing values
    data[0, 0] = np.nan
    if n_samples > 3:
        data[2, 3] = np.nan
    cols = [f"LFQ intensity S{i+1}" for i in range(n_samples)]
    idx = [f"P{i:05d}" for i in range(n_proteins)]
    return pd.DataFrame(data, index=idx, columns=cols)


def _make_metadata(n_samples=6):
    """Create matching metadata: first half treatment, second half control."""
    half = n_samples // 2
    groups = ["treated"] * half + ["control"] * (n_samples - half)
    df = pd.DataFrame({
        "sample_id": [f"S{i+1}" for i in range(n_samples)],
        "group": groups,
    })
    return df.set_index("sample_id")


# ---------------------------------------------------------------------------
# _sep_for
# ---------------------------------------------------------------------------

class TestSepFor:
    def test_tsv(self):
        assert _sep_for(Path("data.tsv")) == "\t"

    def test_txt(self):
        assert _sep_for(Path("data.txt")) == "\t"

    def test_csv(self):
        assert _sep_for(Path("data.csv")) == ","


# ---------------------------------------------------------------------------
# parse_contrast
# ---------------------------------------------------------------------------

class TestParseContrast:
    def test_valid(self):
        t, c = parse_contrast("treated,control")
        assert t == "treated"
        assert c == "control"

    def test_with_spaces(self):
        t, c = parse_contrast(" treated , control ")
        assert t == "treated"
        assert c == "control"

    def test_invalid_too_few(self):
        with pytest.raises(ValueError, match="Contrast must be"):
            parse_contrast("only_one")

    def test_invalid_too_many(self):
        with pytest.raises(ValueError, match="Contrast must be"):
            parse_contrast("a,b,c")

    def test_invalid_empty_part(self):
        with pytest.raises(ValueError, match="Contrast must be"):
            parse_contrast("treated,")


# ---------------------------------------------------------------------------
# load_metadata
# ---------------------------------------------------------------------------

class TestLoadMetadata:
    def test_valid(self, tmp_path):
        f = tmp_path / "meta.csv"
        f.write_text("sample_id,group\nS1,treated\nS2,control\n")
        meta = load_metadata(f)
        assert "group" in meta.columns
        assert meta.index.name == "sample_id"

    def test_missing_sample_id(self, tmp_path):
        f = tmp_path / "meta.csv"
        f.write_text("name,group\nS1,treated\n")
        with pytest.raises(ValueError, match="sample_id"):
            load_metadata(f)

    def test_missing_group(self, tmp_path):
        f = tmp_path / "meta.csv"
        f.write_text("sample_id,condition\nS1,treated\n")
        with pytest.raises(ValueError, match="group"):
            load_metadata(f)


# ---------------------------------------------------------------------------
# align_and_validate
# ---------------------------------------------------------------------------

class TestAlignAndValidate:
    def test_valid_alignment(self):
        matrix = _make_matrix()
        meta = _make_metadata()
        m, md = align_and_validate(matrix, meta, "treated", "control")
        assert list(m.columns) == list(md.index)

    def test_missing_group_raises(self):
        matrix = _make_matrix()
        meta = _make_metadata()
        with pytest.raises(ValueError, match="not present"):
            align_and_validate(matrix, meta, "treated", "nonexistent")

    def test_lfq_prefix_removed(self):
        matrix = _make_matrix()
        meta = _make_metadata()
        m, _ = align_and_validate(matrix, meta, "treated", "control")
        assert all(not c.startswith("LFQ") for c in m.columns)

    def test_too_few_samples(self):
        matrix = _make_matrix(n_samples=2)
        meta = pd.DataFrame({
            "sample_id": ["S1", "S2"],
            "group": ["treated", "control"],
        }).set_index("sample_id")
        with pytest.raises(ValueError, match="at least 2"):
            align_and_validate(matrix, meta, "treated", "control")


# ---------------------------------------------------------------------------
# log2_transform
# ---------------------------------------------------------------------------

class TestLog2Transform:
    def test_known_values(self):
        df = pd.DataFrame({"a": [1.0, 2.0, 4.0, 8.0]})
        result = log2_transform(df)
        expected = [0.0, 1.0, 2.0, 3.0]
        np.testing.assert_array_almost_equal(result["a"].values, expected)


# ---------------------------------------------------------------------------
# _bh_fdr
# ---------------------------------------------------------------------------

class TestBhFdr:
    def test_known_pvalues(self):
        pvals = np.array([0.01, 0.04, 0.03, 0.20])
        adj = _bh_fdr(pvals)
        assert all(0 <= v <= 1 for v in adj)
        # Adjusted p-values should be >= original
        assert all(a >= p for a, p in zip(adj, pvals))

    def test_monotonicity_of_sorted(self):
        pvals = np.array([0.001, 0.01, 0.05, 0.10, 0.50])
        adj = _bh_fdr(pvals)
        sorted_adj = np.sort(adj)
        # Monotonically non-decreasing when sorted
        assert all(sorted_adj[i] <= sorted_adj[i+1] for i in range(len(sorted_adj)-1))

    def test_single_pvalue(self):
        adj = _bh_fdr(np.array([0.05]))
        assert adj[0] == pytest.approx(0.05)


# ---------------------------------------------------------------------------
# run_differential_expression
# ---------------------------------------------------------------------------

class TestRunDE:
    def test_basic_output(self):
        matrix = _make_matrix()
        meta = _make_metadata()
        matrix, meta = align_and_validate(matrix, meta, "treated", "control")
        log_matrix = log2_transform(matrix)
        imp_matrix = down_shift_imputation(log_matrix)
        result = run_differential_expression(imp_matrix, meta, "treated", "control")
        assert "log2FoldChange" in result.columns
        assert "pvalue" in result.columns
        assert "padj" in result.columns
        assert len(result) == 10

    def test_pvalues_in_range(self):
        matrix = _make_matrix()
        meta = _make_metadata()
        matrix, meta = align_and_validate(matrix, meta, "treated", "control")
        log_matrix = log2_transform(matrix)
        imp_matrix = down_shift_imputation(log_matrix)
        result = run_differential_expression(imp_matrix, meta, "treated", "control")
        assert all(0 <= p <= 1 for p in result["pvalue"].dropna())
        assert all(0 <= p <= 1 for p in result["padj"].dropna())


# ---------------------------------------------------------------------------
# run_pca
# ---------------------------------------------------------------------------

class TestRunPCA:
    def test_output_shape(self):
        matrix = _make_matrix()
        meta = _make_metadata()
        matrix, meta = align_and_validate(matrix, meta, "treated", "control")
        log_matrix = log2_transform(matrix)
        imp_matrix = down_shift_imputation(log_matrix)
        pca_df, var_ratio = run_pca(imp_matrix)
        assert pca_df.shape == (6, 3)  # 6 samples, 3 columns (sample_id, PC1, PC2)
        assert len(var_ratio) == 2
        assert sum(var_ratio) <= 1.0 + 1e-10


# ---------------------------------------------------------------------------
# down_shift_imputation
# ---------------------------------------------------------------------------

class TestDownShiftImputation:
    def test_no_nan_in_output(self):
        matrix = _make_matrix()
        meta = _make_metadata()
        matrix, meta = align_and_validate(matrix, meta, "treated", "control")
        log_matrix = log2_transform(matrix)
        result = down_shift_imputation(log_matrix)
        assert not result.isnull().any().any()

    def test_shape_preserved(self):
        matrix = _make_matrix()
        meta = _make_metadata()
        matrix, meta = align_and_validate(matrix, meta, "treated", "control")
        log_matrix = log2_transform(matrix)
        result = down_shift_imputation(log_matrix)
        assert result.shape == log_matrix.shape


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/proteomics-de/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("proteomics-de", tmp_path) == snapshot({'report.md': """\
# Proteomics Differential Expression Report

## Summary
- Samples: 6
- Proteins pre-filter: 2383
- Proteins post-filter/imputation: 2069
- Contrast: `treated vs control`
- Significant proteins: 138
  - Upregulated: 45
  - Downregulated: 93

## Preprocessing
- Imputation distribution: `figures/imputation_distribution.png`
- PCA (sample clustering): `figures/pca.png`

## Differential Expression Results
- Full results: `tables/de_results.csv`
- Imputed protein matrix: `tables/imputed_proteinGroups.csv`
- Volcano plot: `figures/volcano.png`

### Top 10 Upregulated Proteins (by log2 fold change)

| Protein ID | log2FoldChange | -log10(pvalue) |
|---|---:|---:|
| P16522 | 12.509 | 5.345 |
| Q12379 | 8.234 | 2.729 |
| P09798 | 8.031 | 3.217 |
| P42935 | 6.054 | 2.336 |
| P06633 | 5.861 | 1.452 |
| P53886 | 5.663 | 2.165 |
| Q3E747 | 5.629 | 1.775 |
| P14724 | 5.530 | 1.440 |
| P53912;P25608;P54007 | 5.492 | 2.144 |
| P31111 | 5.474 | 1.993 |

### Top 10 Downregulated Proteins (by log2 fold change)

| Protein ID | log2FoldChange | -log10(pvalue) |
|---|---:|---:|
| Q12341 | -9.621 | 4.294 |
| P38210 | -7.432 | 2.567 |
| P39004;P39003 | -7.085 | 1.839 |
| P38238 | -6.099 | 3.672 |
| P53978 | -5.936 | 3.098 |
| Q03774 | -5.816 | 2.230 |
| P38314 | -5.668 | 2.365 |
| P34163 | -5.450 | 1.919 |
| P15992 | -5.427 | 1.694 |
| Q08723 | -5.425 | 3.365 |

## Reproducibility
- Commands: `reproducibility/commands.sh`
- Environment: `reproducibility/environment.yml`
- Checksums: `reproducibility/checksums.sha256`
- Provenance: `ro-crate-metadata.json`

## Disclaimer
ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.
""", 'ro-crate-metadata.json': """\
{
  "@context": "https://w3id.org/ro/crate/1.1/context",
  "@graph": [
    {
      "@id": "./",
      "@type": "Dataset",
      "datePublished": "<TIMESTAMP>",
      "name": "proteomics-de run",
      "description": "Proteomics differential expression from LFQ or DIA-NN data",
      "version": "0.1.0",
      "license": {
        "@id": "https://spdx.org/licenses/MIT.html"
      },
      "hasPart": [
        {
          "@id": "figures/imputation_distribution.png"
        },
        {
          "@id": "figures/pca.png"
        },
        {
          "@id": "figures/volcano.png"
        },
        {
          "@id": "report.md"
        },
        {
          "@id": "reproducibility/checksums.sha256"
        },
        {
          "@id": "reproducibility/commands.sh"
        },
        {
          "@id": "reproducibility/environment.yml"
        },
        {
          "@id": "tables/de_results.csv"
        },
        {
          "@id": "tables/imputed_proteinGroups.csv"
        }
      ],
      "mentions": {
        "@id": "#run"
      }
    },
    {
      "@id": "ro-crate-metadata.json",
      "@type": "CreativeWork",
      "conformsTo": {
        "@id": "https://w3id.org/ro/crate/1.1"
      },
      "about": {
        "@id": "./"
      }
    },
    {
      "@id": "figures/imputation_distribution.png",
      "@type": "File"
    },
    {
      "@id": "figures/pca.png",
      "@type": "File"
    },
    {
      "@id": "figures/volcano.png",
      "@type": "File"
    },
    {
      "@id": "report.md",
      "@type": "File"
    },
    {
      "@id": "reproducibility/checksums.sha256",
      "@type": "File"
    },
    {
      "@id": "reproducibility/commands.sh",
      "@type": "File"
    },
    {
      "@id": "reproducibility/environment.yml",
      "@type": "File"
    },
    {
      "@id": "tables/de_results.csv",
      "@type": "File"
    },
    {
      "@id": "tables/imputed_proteinGroups.csv",
      "@type": "File"
    },
    {
      "@id": "skills/proteomics-de/proteomics_de.py",
      "@type": "SoftwareSourceCode",
      "name": "proteomics-de",
      "version": "0.1.0",
      "programmingLanguage": {
        "@id": "https://www.python.org/"
      }
    },
    {
      "@id": "#input_type",
      "@type": "PropertyValue",
      "name": "input_type",
      "value": "maxquant"
    },
    {
      "@id": "#contrast",
      "@type": "PropertyValue",
      "name": "contrast",
      "value": "treated,control"
    },
    {
      "@id": "#s0",
      "@type": "PropertyValue",
      "name": "s0",
      "value": "0.1"
    },
    {
      "@id": "#fdr",
      "@type": "PropertyValue",
      "name": "fdr",
      "value": "0.05"
    },
    {
      "@id": "#ttest_df",
      "@type": "PropertyValue",
      "name": "ttest_df",
      "value": "4"
    },
    {
      "@id": "#imputation_shift",
      "@type": "PropertyValue",
      "name": "imputation_shift",
      "value": "1.8"
    },
    {
      "@id": "#imputation_scale",
      "@type": "PropertyValue",
      "name": "imputation_scale",
      "value": "0.3"
    },
    {
      "@id": "#demo",
      "@type": "PropertyValue",
      "name": "demo",
      "value": "True"
    },
    {
      "@id": "#run",
      "@type": "CreateAction",
      "name": "proteomics-de execution",
      "instrument": {
        "@id": "skills/proteomics-de/proteomics_de.py"
      },
      "endTime": "<TIMESTAMP>",
      "result": [
        {
          "@id": "figures/imputation_distribution.png"
        },
        {
          "@id": "figures/pca.png"
        },
        {
          "@id": "figures/volcano.png"
        },
        {
          "@id": "report.md"
        },
        {
          "@id": "reproducibility/checksums.sha256"
        },
        {
          "@id": "reproducibility/commands.sh"
        },
        {
          "@id": "reproducibility/environment.yml"
        },
        {
          "@id": "tables/de_results.csv"
        },
        {
          "@id": "tables/imputed_proteinGroups.csv"
        }
      ],
      "object": [
        {
          "@id": "#input_type"
        },
        {
          "@id": "#contrast"
        },
        {
          "@id": "#s0"
        },
        {
          "@id": "#fdr"
        },
        {
          "@id": "#ttest_df"
        },
        {
          "@id": "#imputation_shift"
        },
        {
          "@id": "#imputation_scale"
        },
        {
          "@id": "#demo"
        }
      ]
    }
  ]
}\
"""})
