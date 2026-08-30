import json
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from rnaseq_de import (
    _require_gene_column,
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
