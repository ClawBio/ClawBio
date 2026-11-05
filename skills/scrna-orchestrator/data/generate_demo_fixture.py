#!/usr/bin/env python3
"""Generate the bundled demo fixture `demo_pbmc_small.h5ad`.

Not a DAG step — run once, offline, to materialize a small deterministic
synthetic AnnData fixture so `clawbio start skills/scrna_orchestrator` works
without a network download. The cell x gene count matrix uses the same
synthetic-fallback generator (`build_demo_adata()`) the 1.x skill uses when
its real PBMC3k demo download fails (3 gamma/Poisson-distributed clusters
with marker-gene bumps, MT genes for QC filtering) — only the packaging
differs (this port has no `--demo` flag; pipeline.yaml just points
`input_h5ad` at this file).

One deliberate deviation from the 1.x generator: gene *names*. The source
uses placeholder `GeneNNN`/`MT-GENENN` names, which is fine for QC/clustering
but makes the `annotate_celltypist` step vacuous (CellTypist matches on real
human gene symbols; 0 overlap always errors). GENE_NAMES below is 480 real
HGNC symbols — 19 canonical `MT-*` mitochondrial genes (for QC's
`pct_counts_mt` filter) + 461 symbols sampled (seed=0) from CellTypist's
bundled `Immune_All_Low.pkl` feature list — so the bundled fixture exercises
the full `annotate_celltypist` "done" path too, not just clustering/markers.

Usage: python3 generate_demo_fixture.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

# 19 canonical human mitochondrial gene symbols (all `MT-*` entries in
# CellTypist's Immune_All_Low feature list) + 461 real HGNC symbols sampled
# (seed=0) from that same feature list, for CellTypist-overlap fidelity.
GENE_NAMES = [
    "MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ND3",
    "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6", "MT-RNR1", "MT-RNR2", "MT-TE", "MT-TI", "MT-TT", "MT-TW",
    "NUCKS1", "SMIM1", "NSD2", "FXYD3", "SELL", "NPW", "GNA15", "AHCY", "ARID5B", "HSF5", "TCP1",
    "ROMO1", "PELI1", "ASB2", "NSMCE3", "RPL39L", "FHIT", "ZNF106", "IGHM", "FN3K", "HDHD5", "SLC25A6",
    "MAML3", "THBS1", "MYO7A", "TPM3", "NDUFAF3", "TWISTNB", "IL1R1", "RTN3", "MKRN1", "RPS27P23",
    "SNHG7", "IL12RB2", "CARD11", "SLC4A3", "KLHL14", "TP53BP2", "HES1", "SNHG5", "BLVRA", "APOA1",
    "ALDH2", "KLHL23", "DOK3", "ENO3", "RPL32P1", "RIC3", "SLC37A1", "RANBP1", "EIF3H", "NANOS3",
    "PITPNC1", "NOTCH3", "SOX5", "MTMR6", "F2R", "WASHC3", "RFC5", "COMMD7", "INPP5F", "CPNE2",
    "GALNT14", "RIN2", "MX2", "ITGAE", "AMZ2P1", "E2F1", "SLC7A7", "SLIRP", "TAGAP", "PSIP1", "RAB5IF",
    "CLDN5", "SOX6", "NECTIN2", "ARC", "PZP", "SIVA1", "RPL8", "GLMP", "BCL2L14", "IGLC5", "FABP6",
    "CEP128", "ZNF407", "ZFHX3", "MYL12BP1", "MAP2K6", "PTK2B", "RAB38", "ZDBF2", "NETO1", "IL1RAPL1",
    "LEAP2", "GPC4", "LAIR1", "PCLAF", "CISH", "CD80", "P2RY12", "SCART1", "PFDN5", "DDX5", "GLCCI1",
    "CNR1", "SH3BGRL2", "DLEU1", "FCGRT", "CLECL1", "KIAA1522", "PDCD6IPP2", "RXFP1", "NEFH", "GALC",
    "CETP", "GSTM3", "SLC30A5", "MED12L", "CSH1", "SNHG1", "PRAP1", "HMGN1", "CCSER1", "ITPRIP",
    "HPSE", "MPC2", "USP11", "CYCS", "HAS1", "ATP6V0D1", "SLC2A8", "NABP2", "GPRC5C", "YWHAZP2",
    "CTHRC1", "TROAP", "RAD54B", "PTPRO", "MIR29B2CHG", "HIST1H3F", "CHMP1B", "LTB", "KRT72", "PLXNC1",
    "DSN1", "ZHX3", "PCBP3", "CD1D", "PDLIM7", "PMP22", "HDAC2", "ACAA2", "VASH1", "PEG10", "HOTAIRM1",
    "NCAPG", "CXCL11", "GSTM1", "TCIM", "DCAF12", "MVP", "FSTL3", "HOPX", "SCIMP", "CIDEB", "CYP46A1",
    "MCM7", "OSM", "RPS3P6", "WDFY4", "LAT2", "DNAAF1", "BLVRB", "LYL1", "CAST", "CD300LB", "NSG1",
    "CDK17", "BCL2L12", "SRGAP2", "IL16", "RBMS3", "INPP1", "SSBP2", "MTND5P11", "MSMO1", "IDH3A",
    "RPS11", "PPP2R2B", "ZNF891", "PAN3", "GYPB", "ZNF831", "UROD", "SAMD4A", "MIR22HG", "MAPKAPK2",
    "DENND2D", "SOD2", "COL27A1", "EIF4A2", "UBE2E2", "CLDN10", "TIMM9", "TMIGD3", "IGBP1", "DUSP26",
    "HIGD1A", "B4GALT1", "DYNC2H1", "FDXR", "TAL1", "PDCD1LG2", "SERPINB8", "FASN", "TCF7L2", "GAL3ST4",
    "OR2A25", "SCT", "NIPSNAP2", "PCOLCE", "TCEAL8", "TMEM200A", "EPM2A", "ATF3", "DCTPP1", "RDX",
    "PMCH", "ELOVL7", "CDK2AP2", "ARMH1", "CIB2", "PHF14", "IL18BP", "ERG", "CTAG2", "RPL15", "FSCN1",
    "ABRACL", "TCL1A", "IFI44L", "NOXA1", "SH2B3", "FGL1", "BCL6", "LETM2", "MS4A4E", "PARVB", "LPIN2",
    "BRWD3", "RHPN1", "RPS7P3", "GATA1", "TSPAN13", "BEX1", "ARPC4", "BST2", "CASTOR1", "MIR451B",
    "ALAS2", "MSMB", "BAZ2B", "RPS6KA5", "APOC2", "ARG1", "MSC", "PLAT", "ROBO4", "RNASE1", "NME2P1",
    "TERF2IP", "LDOC1", "FOXM1", "TSC22D1", "STARD7", "RPS26", "VTN", "ATP6V1G3", "RPL19", "TUBA1C",
    "IL2RA", "PLEKHG1", "CNOT6L", "HPDL", "PIK3CA", "HMGB1", "EYS", "HS3ST1", "ADAM12", "TRIR",
    "FOXRED2", "MTA3", "FSTL1", "OTOA", "CD55", "ABCA13", "CDC42SE1", "SULF2", "CRLF2", "RPS2P7",
    "COX4I1", "TFF3", "IL3RA", "CUX2", "CCDC34", "CDC25B", "PLD4", "EPB41L2", "BST1", "POGLUT1",
    "HSPA6", "CDKN3", "CHRDL1", "VNN1", "TENT5C", "IER5", "PHB", "CYP2S1", "FGD6", "RXYLT1", "NDUFB4",
    "PDIA6", "THBD", "FBXL22", "ALDOAP1", "DUSP1", "SGO1", "WDR1", "MDS2", "CNTROB", "LPL", "EFHC2",
    "CCDC88A", "FCER2", "TRAV39", "RPL37AP1", "TRA2B", "RPLP1", "ATP5MGP3", "TGFBR2", "TMEM170B",
    "EZR", "TRPV4", "PDE4A", "CHI3L1", "MTLN", "HBZ", "BTBD3", "B3GAT2", "PAPSS2", "IDH2", "SAP18P2",
    "RNGTT", "SRP14P2", "HSPA2", "ARRB2", "GZMB", "MGP", "ADAMTS5", "ADIPOR1", "RPL23AP42", "FAM83A",
    "RHBDD2", "TPST2", "TBATA", "AKAP12", "ZNF66", "PLAGL1", "HBB", "TPSG1", "VEGFB", "SMIM20", "RGS2",
    "CYB5A", "HBA1", "GASK1B", "SOCS3", "RPS12", "TXNDC16", "KDM5B", "FHL2", "SLCO5A1", "HSP90B3P",
    "HSBP1", "IDS", "BFSP2", "PDCD1", "ACSL6", "NPIPB4", "INPP5A", "GABRE", "PKP3", "RRAD", "CYP2U1",
    "STIP1", "UBE2L3", "REEP3", "CRNDE", "XPNPEP2", "RO60", "SYTL4", "BSG", "RPS4Y2", "LYPD2", "ISOC2",
    "EOMES", "ROGDI", "WDR75", "CBLB", "GNPDA1", "MMRN1", "MMP12", "RPL22L1", "SBDS", "ST6GAL2",
    "MYL6P1", "RPS3AP6", "EMC3", "ANKRD34B", "SLC31A2", "KDM6B", "C1QC", "CRTAM", "IL5RA", "INKA1",
    "FBXO17", "DPEP1", "RASSF4", "ANKRD35", "DNAJB9", "MPHOSPH8", "GUCY1B1", "GLRXP1", "RAPGEF6",
    "ADAMTS1", "SCN1B", "HIST1H1B", "TPRG1", "EDN1", "CRELD2", "CPNE5", "NT5C3A", "PDLIM4", "TRAV35",
    "PAM", "GPSM3", "PALD1", "CD93", "ZFAS1", "MS4A6A", "LIMA1", "PHPT1", "BIN1", "IGLL1", "HOXA5",
    "PTGDS", "CTSS", "TXN", "S100A14", "UROS",
]


def build_demo_adata(random_state: int = 0):
    from anndata import AnnData  # type: ignore

    rng = np.random.default_rng(random_state)
    n_cells = 240
    n_genes = len(GENE_NAMES)
    n_clusters = 3
    cells_per_cluster = n_cells // n_clusters

    base_profiles = []
    for i in range(n_clusters):
        base = rng.gamma(shape=2.0, scale=1.2, size=n_genes)
        marker_start = 40 + i * 20
        marker_end = marker_start + 15
        base[marker_start:marker_end] += 6.0
        base_profiles.append(base)

    expr_blocks = []
    labels = []
    for i, base in enumerate(base_profiles):
        lam = np.clip(base, 0.05, None)
        counts = rng.poisson(lam=lam, size=(cells_per_cluster, n_genes))
        libsize_scale = rng.lognormal(mean=0.0, sigma=0.35, size=(cells_per_cluster, 1))
        counts = np.round(counts * libsize_scale).astype(np.int32)
        expr_blocks.append(counts)
        labels.extend([f"cluster_{i}"] * cells_per_cluster)

    x = np.vstack(expr_blocks)
    gene_names = list(GENE_NAMES)

    obs = pd.DataFrame(
        {
            "sample_id": [f"cell_{i:03d}" for i in range(n_cells)],
            "demo_truth": labels,
        },
        index=[f"cell_{i:03d}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=gene_names)
    return AnnData(X=x, obs=obs, var=var)


def main() -> None:
    adata = build_demo_adata(random_state=0)
    out_path = Path(__file__).resolve().parent / "demo_pbmc_small.h5ad"
    adata.write_h5ad(out_path)
    print(f"wrote {out_path} ({adata.n_obs} cells x {adata.n_vars} genes)")


if __name__ == "__main__":
    main()
