#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA


DISCLAIMER = (
    "ClawBio is a research and educational tool. "
    "It is not a medical device and does not provide clinical diagnoses. "
    "Consult a healthcare professional before making any medical decisions."
)


def _sep_for(path: Path) -> str:
    return "\t" if path.suffix.lower() in (".txt", ".tsv") else ","


def parse_contrast(contrast: str) -> tuple[str, str]:
    parts = [part.strip() for part in contrast.split(",")]
    if len(parts) != 2 or any(not part for part in parts):
        raise ValueError("Contrast must be: treatment,control")
    return parts[0], parts[1]

def load_diann_output(path:Path) -> pd.DataFrame:
    """Load DIA-NN output file and extract LFQ intensity columns"""
    sep = _sep_for(path)
    diann_df = pd.read_csv(path, sep=sep, low_memory=False)
    intensity_cols = [i for i in diann_df.columns if i.endswith('.raw')]
    protein_id_col = ['Protein.Ids']
    diann_df = diann_df[protein_id_col + intensity_cols]
    
    #set index
    diann_df = diann_df.set_index(protein_id_col)
    
    if diann_df.shape[0] < 2:
        raise ValueError("Too few proteins remaining after filtering")
    return diann_df
    
def load_protein_groups(path: Path) -> pd.DataFrame:
    """Load MaxQuant proteinGroups.txt file and apply initial filtering"""
    sep = _sep_for(path)
    proteinGroups_df = pd.read_csv(path, sep=sep, low_memory=False)

    # Auto-detect contaminant column name
    contaminant_col = 'Potential contaminant' if 'Potential contaminant' in proteinGroups_df.columns else 'Contaminant'

    # Select required columns
    need_to_choose = [i for i in proteinGroups_df.columns if i.startswith('LFQ')] +\
        ['Majority protein IDs'] + ['Reverse', 'Only identified by site', contaminant_col]
    proteinGroups_df = proteinGroups_df[need_to_choose]
    proteinGroups_df['Majority protein IDs'] = proteinGroups_df['Majority protein IDs'].fillna('')

    # Contaminant removal
    print(f'Matrix contains {proteinGroups_df.shape[0]} proteins')
    print('Removing reverse, only identified by site, and potential contaminant proteins...')
    proteinGroups_df = proteinGroups_df[~((proteinGroups_df['Reverse'] == '+') |
                        (proteinGroups_df['Only identified by site'] == '+') |
                        (proteinGroups_df[contaminant_col] == '+'))]
    proteinGroups_df = proteinGroups_df.drop(columns=['Reverse', 'Only identified by site', contaminant_col])
    print(f'After contaminant removal, {proteinGroups_df.shape[0]} proteins remain ')

    # Set index
    proteinGroups_df = proteinGroups_df.set_index(['Majority protein IDs'])

    # Convert to numeric, replace 0 with NaN
    proteinGroups_df = proteinGroups_df.apply(pd.to_numeric, errors="coerce")
    proteinGroups_df = proteinGroups_df.replace(0, np.nan)

    # Remove proteins with all missing values
    proteinGroups_df = proteinGroups_df.dropna(how="all")

    if proteinGroups_df.shape[0] < 2:
        raise ValueError("Too few proteins remaining after filtering")

    print('Done!')
    return proteinGroups_df


def load_metadata(path: Path) -> pd.DataFrame:
    sep = _sep_for(path)
    metadata = pd.read_csv(path, sep=sep)
    if "sample_id" not in metadata.columns:
        raise ValueError("Metadata must include a 'sample_id' column")
    if "group" not in metadata.columns:
        raise ValueError("Metadata must include a 'group' column")

    metadata = metadata.copy()
    metadata["sample_id"] = metadata["sample_id"].astype(str)
    metadata = metadata.set_index("sample_id")
    metadata.index.name = "sample_id"
    
    return metadata


def align_and_validate(
    proteomics_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    treatment: str,
    control: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Clean sample names from columns (support both MaxQuant LFQ and DIA-NN raw format)
    sample_names = []
    for col in proteomics_matrix.columns:
        # Handle MaxQuant format: "LFQ intensity SampleName"
        if col.startswith("LFQ intensity "):
            sample_name = col.replace("LFQ intensity ", "")
        # Handle DIA-NN format: "SampleName.raw"
        elif col.endswith(".raw"):
            sample_name = Path(col).stem
        else:
            sample_name = col
        sample_names.append(sample_name)
    proteomics_matrix.columns = sample_names

    # Clean metadata sample ids (remove path and .raw suffix if present)
    metadata.index = metadata.index.map(lambda x: Path(x).stem if x.endswith('.raw') else x)

    # Check samples are present in metadata
    missing_samples = [sample for sample in sample_names if sample not in metadata.index]
    if missing_samples:
        raise ValueError(f"Metadata missing samples: {missing_samples[:5]}")

    # Reorder metadata to match matrix columns
    metadata = metadata.loc[sample_names].copy()
    metadata.index = metadata.index.astype(str)
    metadata.index.name = "sample_id"

    # Check groups exist
    groups = metadata["group"].astype(str)
    if treatment not in set(groups):
        raise ValueError(f"Treatment group '{treatment}' not present in metadata")
    if control not in set(groups):
        raise ValueError(f"Control group '{control}' not present in metadata")

    n_treat = int((groups == treatment).sum())
    n_ctrl = int((groups == control).sum())
    if n_treat < 2 or n_ctrl < 2:
        raise ValueError("Treatment and control groups need at least 2 samples each")

    return proteomics_matrix, metadata


def log2_transform(df: pd.DataFrame) -> pd.DataFrame:
    """Log2 transform LFQ intensities"""
    return np.log2(df)


def down_shift_imputation(log2_transformed_proteinGroups: pd.DataFrame, down_shift: float = 1.8, width: float = 0.3, plot_path: Path = None) -> pd.DataFrame:
    """Down-shift imputation for missing values (revised version)"""
    import seaborn as sns
    import matplotlib.pyplot as plt

    np.random.seed(42)
    df = log2_transformed_proteinGroups.copy()
    impute_df = pd.DataFrame(np.zeros_like(df), columns=df.columns, index=df.index) ## initialize impute_df with zeros
    print('Calculating distribution of total matrix...')
    df = df.replace(-np.inf, np.nan)
    median = np.nanmedian(df.to_numpy().flatten())
    std = np.nanstd(df.to_numpy().flatten())
    mu = median - down_shift * std
    sigma = width * std
    print(f'Median: {median:.2f}, Std: {std:.2f}')
    print(f'In down-shifted distribution, mean: {mu:.2f}, std: {sigma:.2f}')

    print('Imputing missing values...')
    imputed_vals = [] ## store the random num from down-shifted distribution
    norm_data = []
    only_imputed_vals = []  # Store only imputed values for plotting
    for index, row in df.iterrows():
        for num in row.to_numpy():
            if np.isnan(num):
                hit = np.random.normal(size=1, loc=mu, scale=sigma)[0] ## generate random num from down-shifted distribution
                imputed_vals.append(hit)
                norm_data.append(hit)
                only_imputed_vals.append(hit)  # Collect only imputed values
            else:
                imputed_vals.append(num)
                norm_data.append(num)
        impute_df.loc[index, ] = imputed_vals
        imputed_vals = []
    print('Done!')

    # Plot imputation distribution (combined single plot)
    if plot_path:
        data_array = df.to_numpy().flatten()
        original_data = data_array[~np.isnan(data_array)]

        plt.figure(figsize=(9, 6))

        # Plot original data in blue
        sns.histplot(original_data, kde=True, color='#1f77b4', alpha=0.5, label='Original Data')
        # Plot only imputed data in green
        sns.histplot(only_imputed_vals, kde=True, color='#2ca02c', alpha=0.5, label='Imputed Data')

        plt.title('Data Distribution: Original vs Imputed', fontsize=12)
        plt.xlabel('Log2(LFQ Intensity)', fontsize=10)
        plt.ylabel('Count', fontsize=10)
        plt.legend(fontsize=10)

        plt.tight_layout()
        plt.savefig(plot_path, dpi=200)
        plt.close()

    return impute_df


def run_pca(norm_matrix: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    """Run PCA on normalized proteomics data"""
    model = PCA(n_components=2)
    matrix = norm_matrix.T.values
    coords = model.fit_transform(matrix)
    pca_df = pd.DataFrame({
        "sample_id": norm_matrix.columns,
        "PC1": coords[:, 0],
        "PC2": coords[:, 1],
    })
    return pca_df, model.explained_variance_ratio_


def _bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction"""
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


def s0_based_FDR_correction(fc_pvalue_df: pd.DataFrame, degree_of_freedom: int, fdr: float = 0.05, s0: float = 0.1) -> pd.DataFrame:
    '''
    Calculate FDR based on s0 and ta smooth threshold, and annotate regulation status
    '''
    from scipy.stats import t

    def smooth_threshold(x, ta, s0, df):
        """Calculate smoothed -log10(p) threshold for given fold change"""
        import warnings
        warnings.filterwarnings('ignore', category=RuntimeWarning)

        xp = x[x > (ta * s0)]
        xn = x[x < (-ta * s0)]

        dp = xp / ta - s0
        dn = xn / (-ta) - s0
        dp = s0 / dp
        dp = ta * (1 + dp)
        dn = s0 / dn
        dn = ta * (1 + dn)

        fp = t.cdf(dp, df)
        fn = t.cdf(dn, df)
        yp = -np.log10(2 * (1 - fp))
        yn = -np.log10(2 * (1 - fn))

        return {'negative': pd.DataFrame({'log2(FC)': xn, '-log10(p)': yn}, index=xn.index),
                'positive': pd.DataFrame({'log2(FC)': xp, '-log10(p)': yp}, index=xp.index)}

    # Normalize column names
    if 'log2FoldChange' in fc_pvalue_df.columns:
        fc_col = 'log2FoldChange'
    elif 'fc' in fc_pvalue_df.columns:
        fc_col = 'fc'
    else:
        raise ValueError("fc_pvalue_df must contain either 'log2FoldChange' or 'fc' column")

    # Calculate ta from t-distribution based on FDR
    ta = t.ppf(1 - (fdr / 2), df=degree_of_freedom)
    fc_df = fc_pvalue_df[fc_col]

    # Get threshold values
    correction_dict = smooth_threshold(fc_df, ta=ta, s0=s0, df=degree_of_freedom)

    # Merge thresholds
    thresholds = pd.Series(np.nan, index=fc_pvalue_df.index)
    thresholds.loc[correction_dict['positive'].index] = correction_dict['positive']['-log10(p)']
    thresholds.loc[correction_dict['negative'].index] = correction_dict['negative']['-log10(p)']

    # Calculate -log10(pvalue)
    fc_pvalue_df['-log10(pvalue)'] = -np.log10(fc_pvalue_df['pvalue'].clip(lower=1e-300))
    fc_pvalue_df['s0_corrected_-log10(pvalue)'] = thresholds

    # Annotate regulation status
    conditions = [
        (fc_pvalue_df[fc_col] > ta * s0) & (fc_pvalue_df['-log10(pvalue)'] >= fc_pvalue_df['s0_corrected_-log10(pvalue)']),
        (fc_pvalue_df[fc_col] < -ta * s0) & (fc_pvalue_df['-log10(pvalue)'] >= fc_pvalue_df['s0_corrected_-log10(pvalue)']),
        True
    ]
    choices = ['upregulated', 'downregulated', 'non significant']
    fc_pvalue_df['regulation'] = np.select(conditions, choices)

    return fc_pvalue_df


def run_differential_expression(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    treatment: str,
    control: str,
    s0: float = 0.1,
    fdr: float = 0.05,
    ttest_df: int = 4
) -> pd.DataFrame:
    """Two-sample t-test and s0 calibration"""
    groups = metadata["group"].astype(str)
    treat_samples = groups[groups == treatment].index
    ctrl_samples = groups[groups == control].index

    # Get group values
    treat_vals = matrix[treat_samples].values
    ctrl_vals = matrix[ctrl_samples].values

    # Calculate statistics
    mean_treat = np.nanmean(treat_vals, axis=1)
    mean_ctrl = np.nanmean(ctrl_vals, axis=1)
    log2fc = mean_treat - mean_ctrl

    # Two-sample t-test
    t_stat, pvalues = stats.ttest_ind(treat_vals, ctrl_vals, axis=1, equal_var=False, nan_policy="omit")

    # FDR correction
    padj = _bh_fdr(pvalues)

    # Mean intensity across all samples
    mean_intensity = np.nanmean(matrix.values, axis=1)

    # Build results
    result = pd.DataFrame({
        "protein_id": matrix.index,
        "mean_intensity": mean_intensity,
        "mean_treatment": mean_treat,
        "mean_control": mean_ctrl,
        "log2FoldChange": log2fc,
        "pvalue": pvalues,
        "padj": padj
    }).sort_values("padj", ascending=True)

    # Apply FDR correction with s0
    result = s0_based_FDR_correction(result, degree_of_freedom=ttest_df, fdr=fdr, s0=s0)

    # Remove padj column
    result = result.drop(columns=["padj"])

    return result


def plot_pca(pca_df: pd.DataFrame, metadata: pd.DataFrame, var_ratio: np.ndarray, outpath: Path) -> None:
    """Plot PCA figure"""
    plot_df = pca_df.merge(metadata.reset_index(), on="sample_id", how="left")
    plt.figure(figsize=(7, 5))
    for group, group_df in plot_df.groupby("group"):
        plt.scatter(group_df["PC1"], group_df["PC2"], label=str(group), s=50)
    plt.xlabel(f"PC1 ({var_ratio[0] * 100:.1f}%)")
    plt.ylabel(f"PC2 ({var_ratio[1] * 100:.1f}%)")
    plt.title("Proteomics PCA (post-imputation)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def plot_volcano(de_results: pd.DataFrame, outpath: Path, s0: float = 0.1, fdr: float = 0.05, ttest_df: int = 4) -> None:
    """Plot volcano plot with s0 smooth threshold curve"""
    from scipy.stats import t

    df = de_results.copy()
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["log2FoldChange", "pvalue"])

    # Calculate smooth threshold curve for plotting
    ta = t.ppf(1 - (fdr / 2), df=ttest_df)
    x_range = np.linspace(df["log2FoldChange"].min() - 1, df["log2FoldChange"].max() + 1, 1000)

    def get_threshold_curve(x_values, ta, s0, ttest_df):
        import warnings
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        thresholds = []
        for x in x_values:
            if abs(x) <= ta * s0:
                thresholds.append(np.nan)
            else:
                if x > 0:
                    dp = x / ta - s0
                    dp = s0 / dp
                    dp = ta * (1 + dp)
                    fp = t.cdf(dp, ttest_df)
                    y = -np.log10(2 * (1 - fp))
                else:
                    dn = x / (-ta) - s0
                    dn = s0 / dn
                    dn = ta * (1 + dn)
                    fn = t.cdf(dn, ttest_df)
                    y = -np.log10(2 * (1 - fn))
                thresholds.append(y)
        return np.array(thresholds)

    threshold_curve = get_threshold_curve(x_range, ta=ta, s0=s0, ttest_df=ttest_df)

    plt.figure(figsize=(10, 7))

    # Plot points by regulation status
    colors = {
        'upregulated': '#e63946',
        'downregulated': '#1d3557',
        'non significant': '#a8dadc'
    }

    # Rename for backward compatibility
    df['threshold'] = df['s0_corrected_-log10(pvalue)']

    for status, color in colors.items():
        mask = df['regulation'] == status
        plt.scatter(
            df.loc[mask, 'log2FoldChange'],
            df.loc[mask, '-log10(pvalue)'],
            s=15,
            alpha=0.7,
            color=color,
            label=status
        )

    # Plot s0 threshold curve
    plt.plot(x_range, threshold_curve, color='black', linestyle='-', linewidth=1.5, label=f's0={s0}, FDR={fdr} threshold')

    plt.axvline(0, linestyle='-', linewidth=0.8, color='gray')
    plt.axhline(0, linestyle='-', linewidth=0.8, color='gray')

    plt.xlabel("log2 Fold Change (Treatment / Control)", fontsize=12)
    plt.ylabel("-log10 p-value", fontsize=12)
    plt.title("Proteomics Differential Expression Volcano Plot", fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(alpha=0.2, linestyle='--')
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
    input_path: Path,
    input_type: str,
    metadata_path: Path,
    contrast: str,
    s0: float,
    fdr: float,
    ttest_df: int,
    imputation_shift: float,
    imputation_scale: float,
) -> None:
    """Write reproducibility files (commands, environment, checksums)"""
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)

    commands = (
        "python proteomics_de.py "
        f"--input {input_path} "
        f"--input-type {input_type} "
        f"--metadata {metadata_path} "
        f"--contrast \"{contrast}\" "
        f"--s0 {s0} "
        f"--fdr {fdr} "
        f"--ttest-df {ttest_df} "
        f"--imputation-shift {imputation_shift} "
        f"--imputation-scale {imputation_scale} "
        f"--output {output_dir}\n"
    )
    (repro_dir / "commands.sh").write_text(commands)

    env = """name: clawbio-proteomics-de
channels:
  - conda-forge
dependencies:
  - python>=3.10
  - pandas
  - numpy
  - matplotlib
  - scikit-learn
  - scipy
  - seaborn
"""
    (repro_dir / "environment.yml").write_text(env)

    checksums = []
    for path in [input_path, metadata_path]:
        checksums.append(f"{_sha256(path)}  {path.name}")
    for path in sorted((output_dir / "tables").glob("*.csv")):
        checksums.append(f"{_sha256(path)}  tables/{path.name}")
    for path in sorted((output_dir / "figures").glob("*.png")):
        checksums.append(f"{_sha256(path)}  figures/{path.name}")
    (repro_dir / "checksums.sha256").write_text("\n".join(checksums) + "\n")


def run_analysis(
    input_path: Path,
    input_type: str,
    metadata_path: Path,
    contrast: str,
    output_dir: Path,
    s0: float = 0.1,
    fdr: float = 0.05,
    ttest_df: int = 4,
    imputation_shift: float = 1.8,
    imputation_scale: float = 0.3,
) -> dict:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(
            f"Output directory '{output_dir}' is not empty. "
            "Choose a new --output path to avoid overwriting existing reports."
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)

    treatment, control = parse_contrast(contrast)

    # Load data
    raw_input = pd.read_csv(input_path, sep=_sep_for(input_path), low_memory=False)
    n_proteins_before = raw_input.shape[0]

    if input_type == "maxquant":
        protein_matrix = load_protein_groups(input_path)
    elif input_type == "diann":
        protein_matrix = load_diann_output(input_path)
    else:
        raise ValueError(f"Unsupported input type: {input_type}")

    metadata = load_metadata(metadata_path)
    protein_matrix, metadata = align_and_validate(protein_matrix, metadata, treatment, control)

    # Preprocessing
    log2_matrix = log2_transform(protein_matrix)

    # Imputation + plot
    imputed_matrix = down_shift_imputation(
        log2_matrix,
        down_shift=imputation_shift,
        width=imputation_scale,
        plot_path=output_dir / "figures" / "imputation_distribution.png"
    )
    n_proteins_after = imputed_matrix.shape[0]

    # PCA + plot
    pca_df, var_ratio = run_pca(imputed_matrix)
    plot_pca(pca_df, metadata, var_ratio, output_dir / "figures" / "pca.png")

    # Differential expression
    de_results = run_differential_expression(
        imputed_matrix,
        metadata,
        treatment,
        control,
        s0=s0,
        fdr=fdr,
        ttest_df=ttest_df
    )

    # Volcano plot
    plot_volcano(de_results, output_dir / "figures" / "volcano.png", s0=s0, fdr=fdr, ttest_df=ttest_df)

    # Save results
    de_results.to_csv(output_dir / "tables" / "de_results.csv", index=False)
    imputed_matrix.to_csv(output_dir / "tables" / "imputed_proteinGroups.csv")

    # Write reproducibility files
    write_repro_files(
        output_dir,
        input_path=input_path,
        input_type=input_type,
        metadata_path=metadata_path,
        contrast=contrast,
        s0=s0,
        fdr=fdr,
        ttest_df=ttest_df,
        imputation_shift=imputation_shift,
        imputation_scale=imputation_scale,
    )

    # Generate top proteins tables
    upregulated = de_results[de_results['regulation'] == 'upregulated'].sort_values('log2FoldChange', ascending=False).head(10)
    downregulated = de_results[de_results['regulation'] == 'downregulated'].sort_values('log2FoldChange', ascending=True).head(10)

    # Format top upregulated table
    up_rows = []
    for _, row in upregulated.iterrows():
        up_rows.append(f"| {row['protein_id']} | {row['log2FoldChange']:.3f} | {row['-log10(pvalue)']:.3f} |")
    up_table = "\n".join(up_rows) if up_rows else "| No upregulated proteins | | |"

    # Format top downregulated table
    down_rows = []
    for _, row in downregulated.iterrows():
        down_rows.append(f"| {row['protein_id']} | {row['log2FoldChange']:.3f} | {row['-log10(pvalue)']:.3f} |")
    down_table = "\n".join(down_rows) if down_rows else "| No downregulated proteins | | |"

    # Simple report
    report = f"""# Proteomics Differential Expression Report

## Summary
- Samples: {protein_matrix.shape[1]}
- Proteins pre-filter: {n_proteins_before}
- Proteins post-filter/imputation: {n_proteins_after}
- Contrast: `{treatment} vs {control}`
- Significant proteins: {((de_results['regulation'] != 'non significant').sum())}
  - Upregulated: {len(de_results[de_results['regulation'] == 'upregulated'])}
  - Downregulated: {len(de_results[de_results['regulation'] == 'downregulated'])}

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
{up_table}

### Top 10 Downregulated Proteins (by log2 fold change)

| Protein ID | log2FoldChange | -log10(pvalue) |
|---|---:|---:|
{down_table}

## Reproducibility
- Commands: `reproducibility/commands.sh`
- Environment: `reproducibility/environment.yml`
- Checksums: `reproducibility/checksums.sha256`

## Disclaimer
{DISCLAIMER}
"""
    (output_dir / "report.md").write_text(report)

    print(f"Analysis complete! Output saved to: {output_dir}")
    return {
        "output_dir": str(output_dir),
        "samples": protein_matrix.shape[1],
        "proteins_pre": n_proteins_before,
        "proteins_post": n_proteins_after,
        "significant_proteins": int((de_results['regulation'] != 'non significant').sum()),
    }


# ==== Argument Parsing ====
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Minimal proteomics differential expression analysis"
    )

    parser.add_argument("--input", help="Path to input protein quantification file (MaxQuant proteinGroups.txt or DIA-NN output)")
    parser.add_argument("--input-type", choices=["maxquant", "diann"], default="maxquant",
                        help="Input file type: maxquant (default) or diann")
    parser.add_argument("--metadata", help="Path to sample metadata (.csv/.tsv)")
    parser.add_argument("--contrast", default="treated,control", help="Contrast: treatment,control")
    parser.add_argument("--s0", type=float, default=0.1, help="s0 calibration parameter (default: 0.1)")
    parser.add_argument("--fdr", type=float, default=0.05, help="False discovery rate threshold (default: 0.05)")
    parser.add_argument("--ttest-df", type=int, default=4, help="Degree of freedom for t-test (default: 4 for 3+3 replicates)")
    parser.add_argument("--imputation-shift", type=float, default=1.8, help="Down-shift for imputation (default: 1.8)")
    parser.add_argument("--imputation-scale", type=float, default=0.3, help="Scale factor for imputation (default: 0.3)")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run with bundled toy dataset")

    return parser


def validate_args(args, parser):
    if args.demo:
        ## if demo is True, use bundled toy dataset
        here = Path(__file__).resolve().parent
        args.input = str(here / "examples" / "test_proteinGroups.txt")
        args.metadata = str(here / "examples" / "test_metadata.csv")
        args.input_type = "maxquant"
        
    else:
        missing = []
        if args.input is None:
            missing.append("--input")
        if args.metadata is None:
            missing.append("--metadata")
        if args.input_type is None:
            missing.append("--input-type")

        if missing:
            parser.error(f"Missing required arguments: {', '.join(missing)}")

    return args

def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args = validate_args(args, parser)

    result = run_analysis(
        input_path=Path(args.input),
        input_type=args.input_type,
        metadata_path=Path(args.metadata),
        contrast=args.contrast,
        output_dir=Path(args.output),
        s0=args.s0,
        fdr=args.fdr,
        ttest_df=args.ttest_df,
        imputation_shift=args.imputation_shift,
        imputation_scale=args.imputation_scale,
    )

    print(pd.Series(result).to_json(indent=2))


if __name__ == "__main__":
    main()
