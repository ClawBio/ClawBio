import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from rnaseq_de import (
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
