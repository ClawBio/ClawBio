import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

import rnaseq_de
from rnaseq_de import (
    _require_gene_column,
    _resolve_shrinkage_coeff,
    _try_de_pydeseq2,
    align_and_validate,
    compute_qc,
    de_simple,
    filter_low_counts,
    load_counts,
    load_metadata,
    parse_contrast,
    parse_formula_terms,
    run_analysis,
)


HERE = Path(__file__).resolve().parent.parent
DEMO_COUNTS = HERE / "examples" / "demo_counts.csv"
DEMO_META = HERE / "examples" / "demo_metadata.csv"
PSEUDO_COUNTS = HERE / "tests" / "fixtures" / "pseudobulk_counts.csv"
PSEUDO_META = HERE / "tests" / "fixtures" / "pseudobulk_metadata.csv"


def test_formula_parsing():
    terms = parse_formula_terms("~ batch + condition")
    assert terms == ["batch", "condition"]


def test_contrast_parsing():
    factor, numerator, denominator = parse_contrast("condition,treated,control")
    assert factor == "condition"
    assert numerator == "treated"
    assert denominator == "control"


def test_loaders_and_alignment():
    counts = load_counts(DEMO_COUNTS)
    metadata = load_metadata(DEMO_META)
    counts, metadata = align_and_validate(
        counts,
        metadata,
        formula_terms=["batch", "condition"],
        factor="condition",
        numerator="treated",
        denominator="control",
    )
    assert counts.shape == (10, 6)
    assert metadata.shape[0] == 6


def test_qc_and_filtering():
    counts = load_counts(DEMO_COUNTS)
    qc = compute_qc(counts)
    assert {"sample_id", "library_size", "detected_genes"}.issubset(set(qc.columns))
    filtered = filter_low_counts(counts, min_count=10, min_samples=2)
    assert filtered.shape[0] <= counts.shape[0]
    assert filtered.shape[0] >= 2


def test_de_simple_detects_direction():
    counts = load_counts(DEMO_COUNTS)
    metadata = load_metadata(DEMO_META)
    results = de_simple(
        counts,
        metadata,
        factor="condition",
        numerator="treated",
        denominator="control",
    )
    by_gene = results.set_index("gene")
    assert by_gene.loc["GeneA", "log2FoldChange"] > 1.0
    assert by_gene.loc["GeneB", "log2FoldChange"] < -1.0


def test_run_analysis_writes_outputs(tmp_path):
    out_dir = tmp_path / "rnaseq_demo"
    result = run_analysis(
        counts_path=DEMO_COUNTS,
        metadata_path=DEMO_META,
        formula="~ batch + condition",
        contrast="condition,treated,control",
        output_dir=out_dir,
        backend="simple",
    )
    assert result["samples"] == 6
    assert (out_dir / "report.md").exists()
    assert (out_dir / "tables" / "de_results.csv").exists()
    assert (out_dir / "figures" / "pca.png").exists()
    assert (out_dir / "figures" / "volcano.png").exists()
    assert (out_dir / "reproducibility" / "checksums.sha256").exists()
    assert (out_dir / "result.json").exists()


def test_run_analysis_pseudobulk_fixture(tmp_path):
    out_dir = tmp_path / "rnaseq_pseudobulk"
    result = run_analysis(
        counts_path=PSEUDO_COUNTS,
        metadata_path=PSEUDO_META,
        formula="~ cell_type + condition",
        contrast="condition,treated,control",
        output_dir=out_dir,
        backend="simple",
    )
    de_df = pd.read_csv(out_dir / "tables" / "de_results.csv")
    assert result["samples"] == 8
    assert result["genes_post"] >= 2
    assert (out_dir / "figures" / "pca.png").exists()
    assert (out_dir / "tables" / "de_results.csv").exists()
    assert {"gene", "log2FoldChange", "padj"}.issubset(set(de_df.columns))
    assert de_df.shape[0] >= 2


def test_result_json_contains_summary(tmp_path):
    out_dir = tmp_path / "rnaseq_result_json"
    run_analysis(
        counts_path=DEMO_COUNTS,
        metadata_path=DEMO_META,
        formula="~ batch + condition",
        contrast="condition,treated,control",
        output_dir=out_dir,
        backend="simple",
    )
    result_data = json.loads((out_dir / "result.json").read_text(encoding="utf-8"))
    assert result_data["skill"] == "rnaseq"
    assert result_data["summary"]["samples"] == 6
    assert result_data["summary"]["contrast"] == "condition,treated,control"


def test_pydeseq2_reports_lfc_shrinkage(tmp_path):
    pytest.importorskip("pydeseq2")
    out_dir = tmp_path / "rnaseq_pydeseq2"
    run_analysis(
        counts_path=DEMO_COUNTS,
        metadata_path=DEMO_META,
        formula="~ batch + condition",
        contrast="condition,treated,control",
        output_dir=out_dir,
        backend="pydeseq2",
    )
    result_data = json.loads((out_dir / "result.json").read_text(encoding="utf-8"))
    assert result_data["summary"]["backend_used"] == "pydeseq2"
    assert result_data["summary"]["lfc_shrinkage_applied"] is True
    assert result_data["summary"]["lfc_shrinkage_coeff"].startswith("condition[")


def test_load_counts_accepts_nfcore_gene_name_column(tmp_path):
    path = tmp_path / "salmon.merged.gene_counts.tsv"
    path.write_text(
        "gene_id\tgene_name\tctrl_1\tctrl_2\ttrt_1\ttrt_2\n"
        "ENSG0001\tGeneA\t48\t52\t310\t295\n"
        "ENSG0002\tGeneB\t250\t240\t42\t38\n",
        encoding="utf-8",
    )
    counts = load_counts(path)
    assert list(counts.index) == ["ENSG0001", "ENSG0002"]
    assert list(counts.columns) == ["ctrl_1", "ctrl_2", "trt_1", "trt_2"]
    assert counts.loc["ENSG0001", "trt_1"] == 310


def test_require_gene_column_keeps_named_index():
    results = pd.DataFrame(
        {"log2FoldChange": [1.2]},
        index=pd.Index(["ENSG0001"], name="gene_id"),
    )
    out = _require_gene_column(results.reset_index())
    assert list(out["gene"]) == ["ENSG0001"]


def test_require_gene_column_fails_when_identifier_missing():
    results = pd.DataFrame({"log2FoldChange": [1.2], "padj": [0.01]})
    with pytest.raises(ValueError, match="missing gene identifiers"):
        _require_gene_column(results)


def test_run_analysis_nfcore_counts_keep_gene_ids(tmp_path):
    counts_path = tmp_path / "counts.tsv"
    meta_path = tmp_path / "metadata.csv"
    counts_path.write_text(
        "gene_id\tgene_name\tctrl_1\tctrl_2\ttrt_1\ttrt_2\n"
        "ENSG0001\tGeneA\t48\t52\t310\t295\n"
        "ENSG0002\tGeneB\t250\t240\t42\t38\n"
        "ENSG0003\tGeneC\t90\t86\t92\t88\n",
        encoding="utf-8",
    )
    meta_path.write_text(
        "sample_id,condition\nctrl_1,control\nctrl_2,control\ntrt_1,treated\ntrt_2,treated\n",
        encoding="utf-8",
    )
    out_dir = tmp_path / "nfcore_handoff"
    run_analysis(
        counts_path=counts_path,
        metadata_path=meta_path,
        formula="~ condition",
        contrast="condition,treated,control",
        output_dir=out_dir,
        backend="simple",
    )
    de_df = pd.read_csv(out_dir / "tables" / "de_results.csv")
    assert set(de_df["gene"]) == {"ENSG0001", "ENSG0002", "ENSG0003"}
    assert de_df["gene"].notna().all()


# ---------------------------------------------------------------------------
# Issue #365, defect 4: the published log2FoldChange must reflect the
# requested contrast. DeseqStats.lfc_shrink(coeff) replaces the published
# LFC with a single coefficient's column, so shrinkage is only safe when
# that coefficient *is* the requested contrast.
# ---------------------------------------------------------------------------


def _make_nb_dataset(
    n_genes: int = 240,
    per_group: int = 4,
    de_fraction: float = 0.25,
    seed: int = 41,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Negative-binomial counts with a known effect direction per gene."""
    rng = np.random.default_rng(seed)
    genes = [f"ENSG{i:06d}" for i in range(n_genes)]
    base = rng.gamma(2.0, 60.0, size=n_genes)
    # even genes up in treatment, odd genes down, the rest near-null via fraction
    de = rng.random(n_genes) < de_fraction
    up = de & (np.arange(n_genes) % 2 == 0)
    down = de & (np.arange(n_genes) % 2 == 1)
    lfc = np.zeros(n_genes)
    lfc[up] = 2.5
    lfc[down] = -2.5
    samples = [f"ctrl_{i}" for i in range(per_group)] + [f"trt_{i}" for i in range(per_group)]
    conds = ["control"] * per_group + ["treated"] * per_group
    counts = np.zeros((n_genes, 2 * per_group))
    for j, cond in enumerate(conds):
        mu = base * (2 ** lfc if cond == "treated" else 1.0)
        counts[:, j] = rng.negative_binomial(10, 10 / (10 + mu))
    cts = pd.DataFrame(counts.astype(int), index=pd.Index(genes, name="gene_id"), columns=samples)
    meta = pd.DataFrame({"condition": conds}, index=pd.Index(samples, name="sample_id"))
    return cts, meta


def test_pydeseq2_lfc_reflects_requested_contrast(tmp_path):
    """End-to-end acceptance test from #365: nf-core native counts in, gene
    column populated out, and published log2FoldChange tracking ratios of the
    normalized group means (correlation near 1, not ~0.33)."""
    pytest.importorskip("pydeseq2")
    cts, meta = _make_nb_dataset()
    counts_path = tmp_path / "salmon.merged.gene_counts.tsv"
    nfcore = cts.reset_index()  # gene_id becomes the first column
    nfcore.insert(1, "gene_name", [f"SYM{i}" for i in range(len(cts))])
    nfcore.to_csv(counts_path, sep="\t", index=False)
    meta_path = tmp_path / "design.csv"
    meta.to_csv(meta_path)

    out_dir = tmp_path / "out"
    result = run_analysis(
        counts_path=counts_path,
        metadata_path=meta_path,
        formula="~ condition",
        contrast="condition,treated,control",
        output_dir=out_dir,
        backend="pydeseq2",
    )
    assert result["backend_used"] == "pydeseq2"

    de_df = pd.read_csv(out_dir / "tables" / "de_results.csv")
    norm = pd.read_csv(out_dir / "tables" / "normalized_counts.csv", index_col=0)
    # defect 1+2 regression: identifiers survive the handoff
    assert de_df["gene"].notna().all()
    assert set(de_df["gene"]) <= set(cts.index)
    assert len(de_df) > 100
    # defect 4: published LFC tracks the normalized group-mean ratios
    meta_groups = meta["condition"].reindex(norm.columns)
    naive = (
        np.log2(norm.loc[:, meta_groups == "treated"].mean(axis=1) + 1)
        - np.log2(norm.loc[:, meta_groups == "control"].mean(axis=1) + 1)
    )
    published = de_df.set_index("gene")["log2FoldChange"].reindex(naive.index)
    assert published.corr(naive) > 0.9
    top5 = de_df.head(5)
    assert (top5["log2FoldChange"].abs() > 0.5).all(), (
        f"top genes by padj carry near-zero effect sizes: {list(top5['log2FoldChange'])}"
    )


def test_pydeseq2_refuses_shrinkage_on_mismatched_coefficient(monkeypatch):
    """A resolved coefficient that is not the requested contrast must never
    replace the published LFC, even if a (buggy or future) resolver hands one
    back that exists in the fitted model."""
    pytest.importorskip("pydeseq2")
    rng = np.random.default_rng(17)
    n_genes = 120
    genes = [f"ENSG{i:06d}" for i in range(n_genes)]
    base = rng.gamma(2.0, 60.0, size=n_genes)
    lfc = np.where(np.arange(n_genes) % 2 == 0, 2.0, -2.0)  # trtA effect only
    samples = (
        [f"ctl_{i}" for i in range(4)]
        + [f"trtA_{i}" for i in range(4)]
        + [f"trtB_{i}" for i in range(4)]
    )
    conds = ["control"] * 4 + ["trtA"] * 4 + ["trtB"] * 4
    counts = np.zeros((n_genes, 12))
    for j, c in enumerate(conds):
        mu = base * (2 ** lfc if c == "trtA" else 1.0)
        counts[:, j] = rng.negative_binomial(10, 10 / (10 + mu))
    cts = pd.DataFrame(counts.astype(int), index=pd.Index(genes, name="gene_id"), columns=samples)
    meta = pd.DataFrame({"condition": conds}, index=pd.Index(samples, name="sample_id"))

    # simulate a resolver that hands back the WRONG level's coefficient
    monkeypatch.setattr(rnaseq_de, "_resolve_shrinkage_coeff", lambda dds, f, n: "condition[T.trtB]")
    with pytest.warns(UserWarning, match="does not reproduce the requested contrast"):
        res, shrinkage = _try_de_pydeseq2(cts, meta, ["condition"], "condition", "trtA", "control")

    assert shrinkage["lfc_shrinkage_applied"] is False
    assert "unshrunk Wald MLE" in shrinkage["lfc_shrinkage_note"]
    # the published LFC is still the requested contrast (trtA vs control):
    # even genes were designed up in trtA, odd genes down
    by_gene = res.set_index("gene")["log2FoldChange"]
    up_genes = [g for i, g in enumerate(genes) if i % 2 == 0]
    down_genes = [g for i, g in enumerate(genes) if i % 2 == 1]
    assert (by_gene.reindex(up_genes) > 0).mean() > 0.9
    assert (by_gene.reindex(down_genes) < 0).mean() > 0.9


def test_resolve_shrinkage_coeff_exact_match_only():
    """The resolver must return '' rather than a substring-matching column
    for a different level (the sign-flip trap from #365)."""
    pytest.importorskip("pydeseq2")

    class FakeLFC:
        pass

    fake_dds = FakeLFC()
    fake_dds.varm = {
        "LFC": pd.DataFrame(
            columns=["Intercept", "condition[T.AB]"],
            index=["g1", "g2"],
        )
    }
    # numerator 'A' is the reference level here (no condition[T.A] column);
    # 'A' is a substring of 'condition[T.AB]' and the old resolver fell for it
    assert _resolve_shrinkage_coeff(fake_dds, "condition", "A") == ""
    assert _resolve_shrinkage_coeff(fake_dds, "condition", "AB") == "condition[T.AB]"
    assert _resolve_shrinkage_coeff(fake_dds, "condition", "missing") == ""


def test_pydeseq2_internal_failure_is_reported_clearly(monkeypatch, tmp_path):
    """A crash inside pydeseq2 must surface as a clear RuntimeError naming the
    backend and the simple-backend escape hatch, not a bare IndexError."""
    pytest.importorskip("pydeseq2")
    from pydeseq2.dds import DeseqDataSet

    def _boom(self):
        raise IndexError("too many indices for array")

    monkeypatch.setattr(DeseqDataSet, "deseq2", _boom)
    cts, meta = _make_nb_dataset(n_genes=60)
    counts_path = tmp_path / "counts.csv"
    meta_path = tmp_path / "meta.csv"
    cts.to_csv(counts_path)
    meta.to_csv(meta_path)
    with pytest.raises(RuntimeError, match="pydeseq2 backend failed.*--backend simple"):
        run_analysis(
            counts_path=counts_path,
            metadata_path=meta_path,
            formula="~ condition",
            contrast="condition,treated,control",
            output_dir=tmp_path / "out",
            backend="pydeseq2",
        )
