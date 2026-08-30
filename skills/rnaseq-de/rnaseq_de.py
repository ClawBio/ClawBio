#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

from scipy.stats import t as _t_dist

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from clawbio.common.report import write_result_json


DISCLAIMER = (
    "ClawBio is a research and educational tool. "
    "It is not a medical device and does not provide clinical diagnoses. "
    "Consult a healthcare professional before making any medical decisions."
)


def _sep_for(path: Path) -> str:
    return "\t" if path.suffix.lower() == ".tsv" else ","


def parse_formula_terms(formula: str) -> list[str]:
    f = formula.strip()
    if not f.startswith("~"):
        raise ValueError("Formula must start with '~', for example: '~ condition'")
    rhs = f[1:].strip()
    if not rhs:
        raise ValueError("Formula must include at least one design factor")
    terms = [term.strip() for term in rhs.split("+") if term.strip()]
    if not terms:
        raise ValueError("Could not parse formula terms")
    return terms


def parse_contrast(contrast: str) -> tuple[str, str, str]:
    parts = [part.strip() for part in contrast.split(",")]
    if len(parts) != 3 or any(not part for part in parts):
        raise ValueError("Contrast must be: factor,numerator,denominator")
    return parts[0], parts[1], parts[2]


_GENE_ID_COLUMNS = ("gene_id", "Geneid", "gene", "GENEID", "GeneID")
_GENE_NAME_COLUMNS = ("gene_name", "gene_symbol", "GeneName", "symbol")


def _select_gene_id_column(columns: pd.Index) -> str:
    for name in _GENE_ID_COLUMNS:
        if name in columns:
            return name
    if len(columns) == 0:
        raise ValueError("Count matrix is missing a gene identifier column")
    return str(columns[0])


def load_counts(path: Path) -> pd.DataFrame:
    sep = _sep_for(path)
    df = pd.read_csv(path, sep=sep)
    if df.shape[1] < 3:
        raise ValueError("Count matrix must have one gene column and at least two samples")
    gene_col = _select_gene_id_column(df.columns)
    drop = [col for col in _GENE_NAME_COLUMNS if col in df.columns and col != gene_col]
    counts = df.drop(columns=drop).set_index(gene_col)
    counts.index.name = gene_col
    counts = counts.apply(pd.to_numeric, errors="coerce")
    if counts.isna().any().any():
        raise ValueError("Count matrix contains non-numeric entries")
    if (counts < 0).any().any():
        raise ValueError("Count matrix contains negative counts")
    counts.index = counts.index.astype(str)
    counts.columns = counts.columns.astype(str)
    return counts


def load_metadata(path: Path) -> pd.DataFrame:
    sep = _sep_for(path)
    metadata = pd.read_csv(path, sep=sep)
    if "sample_id" not in metadata.columns:
        raise ValueError("Metadata must include a 'sample_id' column")
    metadata = metadata.copy()
    metadata["sample_id"] = metadata["sample_id"].astype(str)
    metadata = metadata.set_index("sample_id")
    metadata.index.name = "sample_id"
    return metadata


def align_and_validate(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    formula_terms: list[str],
    factor: str,
    numerator: str,
    denominator: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    missing_samples = [sample for sample in counts.columns if sample not in metadata.index]
    if missing_samples:
        raise ValueError(f"Metadata missing samples: {missing_samples[:5]}")
    metadata = metadata.loc[counts.columns].copy()
    metadata.index = metadata.index.astype(str)
    metadata.index.name = "sample_id"

    for term in formula_terms:
        if term not in metadata.columns:
            raise ValueError(f"Formula term '{term}' not found in metadata")
    if factor not in metadata.columns:
        raise ValueError(f"Contrast factor '{factor}' not found in metadata")

    groups = metadata[factor].astype(str)
    if numerator not in set(groups):
        raise ValueError(f"Contrast numerator '{numerator}' not present in {factor}")
    if denominator not in set(groups):
        raise ValueError(f"Contrast denominator '{denominator}' not present in {factor}")

    n_num = int((groups == numerator).sum())
    n_den = int((groups == denominator).sum())
    if n_num < 2 or n_den < 2:
        raise ValueError("Contrast groups need at least 2 samples each")

    return counts, metadata


def compute_qc(counts: pd.DataFrame) -> pd.DataFrame:
    lib_sizes = counts.sum(axis=0)
    detected_genes = (counts > 0).sum(axis=0)
    qc = pd.DataFrame(
        {
            "sample_id": counts.columns,
            "library_size": lib_sizes.values,
            "detected_genes": detected_genes.values,
        }
    )
    return qc


def filter_low_counts(counts: pd.DataFrame, min_count: int, min_samples: int) -> pd.DataFrame:
    mask = (counts >= min_count).sum(axis=1) >= min_samples
    filtered = counts.loc[mask].copy()
    if filtered.shape[0] < 2:
        raise ValueError("Too few genes after low-count filtering; relax thresholds")
    return filtered


def cpm_log1p(counts: pd.DataFrame) -> pd.DataFrame:
    lib_sizes = counts.sum(axis=0)
    cpm = counts.div(lib_sizes, axis=1) * 1_000_000.0
    return np.log1p(cpm)


def run_pca(norm_log_counts: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    model = PCA(n_components=2)
    matrix = norm_log_counts.T.values
    coords = model.fit_transform(matrix)
    pca_df = pd.DataFrame(
        {
            "sample_id": norm_log_counts.columns,
            "PC1": coords[:, 0],
            "PC2": coords[:, 1],
        }
    )
    return pca_df, model.explained_variance_ratio_


def _bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    n = len(pvalues)
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    adj = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        value = ranked[i] * n / rank
        prev = min(prev, value)
        adj[i] = prev
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adj, 0, 1)
    return out


def de_simple(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    factor: str,
    numerator: str,
    denominator: str,
) -> pd.DataFrame:
    groups = metadata[factor].astype(str)
    num_samples = groups[groups == numerator].index
    den_samples = groups[groups == denominator].index

    cpm = counts.div(counts.sum(axis=0), axis=1) * 1_000_000.0
    log_cpm = np.log2(cpm + 1.0)

    x_num = log_cpm[num_samples]
    x_den = log_cpm[den_samples]

    mean_num = x_num.mean(axis=1)
    mean_den = x_den.mean(axis=1)
    var_num = x_num.var(axis=1, ddof=1)
    var_den = x_den.var(axis=1, ddof=1)

    n_num = float(len(num_samples))
    n_den = float(len(den_samples))
    v_num = var_num / n_num
    v_den = var_den / n_den
    se = np.sqrt(v_num + v_den + 1e-12)
    t_stat = (mean_num - mean_den) / se

    # Welch-Satterthwaite degrees of freedom; clipped to ≥1 to stay on t-distribution
    dof_num = (v_num + v_den) ** 2
    dof_den = (v_num ** 2) / max(n_num - 1, 1) + (v_den ** 2) / max(n_den - 1, 1)
    dof = np.where(dof_den > 0, dof_num / dof_den, 1.0).clip(1.0)
    pvals = 2 * _t_dist.sf(np.abs(t_stat), df=dof)
    padj = _bh_fdr(pvals)

    base_mean = cpm.mean(axis=1)
    log2fc = mean_num - mean_den

    result = pd.DataFrame(
        {
            "gene": counts.index,
            "baseMean": base_mean.values,
            "log2FoldChange": log2fc.values,
            "pvalue": pvals,
            "padj": padj,
        }
    ).sort_values("padj", ascending=True)
    return result


def _try_de_pydeseq2(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    formula_terms: list[str],
    factor: str,
    numerator: str,
    denominator: str,
) -> tuple[pd.DataFrame, dict[str, str | bool]] | None:
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except Exception:
        return None

    md = metadata.copy()
    for term in formula_terms:
        if term == factor:
            observed = pd.Series(md[term].astype(str)).dropna().tolist()
            ordered_levels = [denominator] + [level for level in dict.fromkeys(observed) if level != denominator]
            md[term] = pd.Categorical(md[term].astype(str), categories=ordered_levels, ordered=True)
        else:
            md[term] = md[term].astype("category")

    cts = counts.T.round().astype(int)

    design = "~ " + " + ".join(formula_terms)
    try:
        dds = DeseqDataSet(
            counts=cts,
            metadata=md,
            design=design,
            refit_cooks=True,
        )
    except TypeError:
        dds = DeseqDataSet(
            counts=cts,
            metadata=md,
            design_factors=formula_terms,
            refit_cooks=True,
        )

    dds.deseq2()
    stats = DeseqStats(dds, contrast=[factor, numerator, denominator])
    stats.summary()
    shrinkage = {
        "lfc_shrinkage_applied": False,
        "lfc_shrinkage_coeff": "",
        "lfc_shrinkage_note": "",
    }
    coeff = _resolve_shrinkage_coeff(dds, factor, numerator)
    if coeff:
        try:
            stats.lfc_shrink(coeff=coeff)
            shrinkage["lfc_shrinkage_applied"] = True
            shrinkage["lfc_shrinkage_coeff"] = coeff
        except Exception as exc:
            shrinkage["lfc_shrinkage_note"] = f"Attempted LFC shrinkage with '{coeff}' but failed: {exc}"
    else:
        shrinkage["lfc_shrinkage_note"] = "Could not identify a PyDESeq2 coefficient for LFC shrinkage."
    res = _require_gene_column(stats.results_df.reset_index())
    if "baseMean" not in res.columns:
        base_mean = (counts.div(counts.sum(axis=0), axis=1) * 1_000_000.0).mean(axis=1)
        res = res.merge(base_mean.rename("baseMean"), left_on="gene", right_index=True, how="left")
    cols = ["gene", "baseMean", "log2FoldChange", "pvalue", "padj"]
    keep = [col for col in cols if col in res.columns]
    return res[keep].sort_values("padj", ascending=True), shrinkage


def _require_gene_column(results: pd.DataFrame) -> pd.DataFrame:
    """Return DE results with a non-empty `gene` column, or fail closed."""
    res = results.copy()
    if "gene" not in res.columns:
        for name in (*_GENE_ID_COLUMNS, "index"):
            if name in res.columns:
                res = res.rename(columns={name: "gene"})
                break
    if "gene" not in res.columns:
        raise ValueError("Differential expression results are missing gene identifiers")
    if res["gene"].isna().any():
        raise ValueError("Differential expression results are missing gene identifiers")
    genes = res["gene"].astype(str).str.strip()
    if (genes == "").any() or (genes.str.lower() == "nan").any():
        raise ValueError("Differential expression results are missing gene identifiers")
    return res


def _resolve_shrinkage_coeff(dds: object, factor: str, numerator: str) -> str:
    lfc_table = getattr(getattr(dds, "varm", {}), "get", lambda *_args, **_kwargs: None)("LFC")
    columns = [str(col) for col in getattr(lfc_table, "columns", [])]
    preferred = f"{factor}[T.{numerator}]"
    if preferred in columns:
        return preferred
    suffix = f"[T.{numerator}]"
    for col in columns:
        if col.startswith(f"{factor}[") and col.endswith(suffix):
            return col
    for col in columns:
        if col != "Intercept" and col.startswith(f"{factor}[") and numerator in col:
            return col
    return ""


def run_de(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    formula_terms: list[str],
    factor: str,
    numerator: str,
    denominator: str,
    backend: str,
) -> tuple[pd.DataFrame, str, dict[str, str | bool]]:
    if backend in {"auto", "pydeseq2"}:
        pydeseq_res = _try_de_pydeseq2(counts, metadata, formula_terms, factor, numerator, denominator)
        if pydeseq_res is not None:
            return pydeseq_res[0], "pydeseq2", pydeseq_res[1]
        if backend == "pydeseq2":
            raise RuntimeError("PyDESeq2 backend requested but unavailable in this environment")

    return de_simple(counts, metadata, factor, numerator, denominator), "simple", {
        "lfc_shrinkage_applied": False,
        "lfc_shrinkage_coeff": "",
        "lfc_shrinkage_note": "LFC shrinkage is only available with the PyDESeq2 backend.",
    }


def plot_pca(pca_df: pd.DataFrame, metadata: pd.DataFrame, factor: str, var_ratio: np.ndarray, outpath: Path) -> None:
    plot_df = pca_df.merge(metadata.reset_index(), on="sample_id", how="left")
    plt.figure(figsize=(7, 5))
    for group, group_df in plot_df.groupby(factor):
        plt.scatter(group_df["PC1"], group_df["PC2"], label=str(group), s=50)
    plt.xlabel(f"PC1 ({var_ratio[0] * 100:.1f}%)")
    plt.ylabel(f"PC2 ({var_ratio[1] * 100:.1f}%)")
    plt.title("RNA-seq PCA (pre-DE)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def plot_volcano(de_results: pd.DataFrame, outpath: Path) -> None:
    df = de_results.copy()
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["log2FoldChange", "padj"])
    y = -np.log10(df["padj"].clip(lower=1e-300))
    sig = (df["padj"] < 0.05) & (df["log2FoldChange"].abs() >= 1.0)
    plt.figure(figsize=(7, 5))
    plt.scatter(df.loc[~sig, "log2FoldChange"], y.loc[~sig], s=10, alpha=0.6)
    plt.scatter(df.loc[sig, "log2FoldChange"], y.loc[sig], s=12, alpha=0.8)
    plt.axvline(1.0, linestyle="--", linewidth=1)
    plt.axvline(-1.0, linestyle="--", linewidth=1)
    plt.axhline(-np.log10(0.05), linestyle="--", linewidth=1)
    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10 adjusted p-value")
    plt.title("Volcano Plot")
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def plot_ma(de_results: pd.DataFrame, outpath: Path) -> None:
    df = de_results.copy()
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["baseMean", "log2FoldChange", "padj"])
    x = np.log10(df["baseMean"].clip(lower=1e-6))
    y = df["log2FoldChange"]
    sig = df["padj"] < 0.05
    plt.figure(figsize=(7, 5))
    plt.scatter(x.loc[~sig], y.loc[~sig], s=10, alpha=0.6)
    plt.scatter(x.loc[sig], y.loc[sig], s=12, alpha=0.8)
    plt.axhline(0.0, linestyle="--", linewidth=1)
    plt.xlabel("log10 baseMean")
    plt.ylabel("log2 Fold Change")
    plt.title("MA Plot")
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def write_repro_files(
    output_dir: Path,
    counts_path: Path,
    metadata_path: Path,
    formula: str,
    contrast: str,
    backend: str,
) -> None:
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)

    commands = (
        "python rnaseq_de.py "
        f"--counts {counts_path} "
        f"--metadata {metadata_path} "
        f"--formula \"{formula}\" "
        f"--contrast \"{contrast}\" "
        f"--backend {backend} "
        f"--output {output_dir}\n"
    )
    (repro_dir / "commands.sh").write_text(commands)

    env = """name: clawbio-rnaseq-de
channels:
  - conda-forge
dependencies:
  - python>=3.10
  - pandas
  - numpy
  - matplotlib
  - scikit-learn
  - pydeseq2
"""
    (repro_dir / "environment.yml").write_text(env)

    checksums = []
    for path in [counts_path, metadata_path]:
        checksums.append(f"{_sha256(path)}  {path.name}")
    for path in sorted((output_dir / "tables").glob("*.csv")):
        checksums.append(f"{_sha256(path)}  tables/{path.name}")
    for path in sorted((output_dir / "figures").glob("*.png")):
        checksums.append(f"{_sha256(path)}  figures/{path.name}")
    (repro_dir / "checksums.sha256").write_text("\n".join(checksums) + "\n")


def write_report(
    output_dir: Path,
    n_samples: int,
    n_genes_before: int,
    n_genes_after: int,
    formula: str,
    contrast: str,
    backend_used: str,
    de_results: pd.DataFrame,
    shrinkage_meta: dict[str, str | bool],
) -> None:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    top = de_results.head(10)
    top_rows = "\n".join(
        [
            f"| {row.gene} | {row.log2FoldChange:.3f} | {row.padj:.3e} |"
            for row in top.itertuples(index=False)
        ]
    )
    report = f"""# ClawBio RNA-seq Differential Expression Report

**Date**: {now}
**Samples**: {n_samples}
**Genes (pre-filter)**: {n_genes_before}
**Genes (post-filter)**: {n_genes_after}
**Formula**: `{formula}`
**Contrast**: `{contrast}`
**Backend used**: `{backend_used}`
**LFC shrinkage**: `{"applied" if shrinkage_meta.get("lfc_shrinkage_applied") else "not applied"}`
"""
    if shrinkage_meta.get("lfc_shrinkage_coeff"):
        report += f"**LFC shrinkage coefficient**: `{shrinkage_meta['lfc_shrinkage_coeff']}`\n"
    if shrinkage_meta.get("lfc_shrinkage_note"):
        report += f"**LFC shrinkage note**: {shrinkage_meta['lfc_shrinkage_note']}\n"
    report += f"""

## Pre-DE QC + PCA

- QC summary: `tables/qc_summary.csv`
- PCA figure: `figures/pca.png`

## Differential Expression

- Full results: `tables/de_results.csv`
- Volcano plot: `figures/volcano.png`
- MA plot: `figures/ma_plot.png`

### Top Genes (by adjusted p-value)

| Gene | log2FoldChange | padj |
|---|---:|---:|
{top_rows}

## Reproducibility

- Commands: `reproducibility/commands.sh`
- Environment: `reproducibility/environment.yml`
- Checksums: `reproducibility/checksums.sha256`

## Disclaimer

{DISCLAIMER}
"""
    (output_dir / "report.md").write_text(report)


def run_analysis(
    counts_path: Path,
    metadata_path: Path,
    formula: str,
    contrast: str,
    output_dir: Path,
    backend: str = "auto",
    min_count: int = 10,
    min_samples: int = 2,
) -> dict:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(
            f"Output directory '{output_dir}' is not empty. "
            "Choose a new --output path to avoid overwriting existing reports."
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)

    formula_terms = parse_formula_terms(formula)
    factor, numerator, denominator = parse_contrast(contrast)
    counts = load_counts(counts_path)
    metadata = load_metadata(metadata_path)
    counts, metadata = align_and_validate(counts, metadata, formula_terms, factor, numerator, denominator)

    n_genes_before = counts.shape[0]
    qc = compute_qc(counts)
    filtered = filter_low_counts(counts, min_count=min_count, min_samples=min_samples)
    n_genes_after = filtered.shape[0]
    norm_log = cpm_log1p(filtered)
    pca_df, var_ratio = run_pca(norm_log)
    de_results, backend_used, shrinkage_meta = run_de(
        filtered,
        metadata,
        formula_terms,
        factor,
        numerator,
        denominator,
        backend,
    )

    qc.to_csv(output_dir / "tables" / "qc_summary.csv", index=False)
    norm_log.to_csv(output_dir / "tables" / "normalized_counts.csv")
    de_results.to_csv(output_dir / "tables" / "de_results.csv", index=False)

    plot_pca(pca_df, metadata, factor, var_ratio, output_dir / "figures" / "pca.png")
    plot_volcano(de_results, output_dir / "figures" / "volcano.png")
    plot_ma(de_results, output_dir / "figures" / "ma_plot.png")

    write_repro_files(output_dir, counts_path, metadata_path, formula, contrast, backend)
    write_report(
        output_dir,
        n_samples=counts.shape[1],
        n_genes_before=n_genes_before,
        n_genes_after=n_genes_after,
        formula=formula,
        contrast=contrast,
        backend_used=backend_used,
        de_results=de_results,
        shrinkage_meta=shrinkage_meta,
    )
    write_result_json(
        output_dir=output_dir,
        skill="rnaseq",
        version="0.1.0",
        summary={
            "samples": counts.shape[1],
            "genes_pre": n_genes_before,
            "genes_post": n_genes_after,
            "formula": formula,
            "contrast": contrast,
            "backend_used": backend_used,
            "lfc_shrinkage_applied": bool(shrinkage_meta.get("lfc_shrinkage_applied")),
            "lfc_shrinkage_coeff": str(shrinkage_meta.get("lfc_shrinkage_coeff", "")),
            "lfc_shrinkage_note": str(shrinkage_meta.get("lfc_shrinkage_note", "")),
        },
        data={
            "input": {
                "counts": counts_path.name,
                "metadata": metadata_path.name,
            },
            "tables": ["qc_summary.csv", "normalized_counts.csv", "de_results.csv"],
            "figures": ["pca.png", "volcano.png", "ma_plot.png"],
            "disclaimer": DISCLAIMER,
        },
    )

    return {
        "output_dir": str(output_dir),
        "samples": counts.shape[1],
        "genes_pre": n_genes_before,
        "genes_post": n_genes_after,
        "backend_used": backend_used,
        "lfc_shrinkage_applied": bool(shrinkage_meta.get("lfc_shrinkage_applied")),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="RNA-seq differential expression (bulk + pseudo-bulk)")
    parser.add_argument("--input", help="Compatibility input as 'counts.csv,metadata.csv'")
    parser.add_argument("--counts", help="Path to counts matrix (.csv/.tsv)")
    parser.add_argument("--metadata", help="Path to sample metadata (.csv/.tsv)")
    parser.add_argument("--formula", default="~ condition", help="Design formula, e.g. '~ batch + condition'")
    parser.add_argument("--contrast", default="condition,treated,control", help="Contrast: factor,numerator,denominator")
    parser.add_argument("--backend", choices=["auto", "pydeseq2", "simple"], default="auto")
    parser.add_argument("--min-count", type=int, default=10)
    parser.add_argument("--min-samples", type=int, default=2)
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run with bundled toy dataset")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    here = Path(__file__).resolve().parent
    if args.demo:
        counts_path = here / "examples" / "demo_counts.csv"
        metadata_path = here / "examples" / "demo_metadata.csv"
        formula = "~ batch + condition"
        contrast = "condition,treated,control"
    else:
        if args.input and (not args.counts and not args.metadata):
            parts = [part.strip() for part in args.input.split(",")]
            if len(parts) != 2 or any(not part for part in parts):
                parser.error("--input for rnaseq-de must be 'counts.csv,metadata.csv'")
            counts_path = Path(parts[0])
            metadata_path = Path(parts[1])
        else:
            if not args.counts or not args.metadata:
                parser.error("Either use --demo, provide --counts/--metadata, or use --input 'counts,metadata'")
            counts_path = Path(args.counts)
            metadata_path = Path(args.metadata)
        formula = args.formula
        contrast = args.contrast

    result = run_analysis(
        counts_path=counts_path,
        metadata_path=metadata_path,
        formula=formula,
        contrast=contrast,
        output_dir=Path(args.output),
        backend=args.backend,
        min_count=args.min_count,
        min_samples=args.min_samples,
    )

    print(pd.Series(result).to_json(indent=2))


if __name__ == "__main__":
    main()
