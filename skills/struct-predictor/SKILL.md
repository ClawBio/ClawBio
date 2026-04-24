---
name: struct-predictor
description: Protein structure prediction with Boltz-2. Accepts YAML inputs (single protein or multi-chain complex), runs
  boltz predict, extracts per-residue pLDDT and PAE confidence, and writes a markdown report with figures.
license: MIT
metadata:
  version: 0.2.0
  openclaw:
    requires:
      bins:
      - python3
      anyBins:
      - boltz
    always: false
    emoji: 🧱
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    install:
    - kind: uv
      package: boltz
      bins:
      - boltz
      comment: 'GPU: uv pip install ''boltz[cuda]'' -U'
    - kind: uv
      package: numpy
    - kind: uv
      package: matplotlib
    - kind: uv
      package: pyyaml
---

# Struct Predictor

You are the **Struct Predictor**, a specialised agent for protein structure prediction using Boltz-2.

## Core Capabilities

1. **Structure Prediction**: Run Boltz-2 locally on a YAML input
2. **Confidence Extraction**: Per-residue pLDDT (from CIF B-factors) and PAE matrix (from confidence JSON)
3. **Report Generation**: Markdown with pLDDT line plot, PAE heatmap, band breakdown, and reproducibility bundle
4. **Demo Mode**: Trp-cage miniprotein (20 residues, PDB 1L2Y) — runs immediately, no input required

## CLI Reference

```bash
# Single protein or multi-chain complex (YAML)
python skills/struct-predictor/struct_predictor.py \
  --input complex.yaml --output /tmp/struct_out

# Demo (Trp-cage miniprotein, PDB 1L2Y — no input needed)
python skills/struct-predictor/struct_predictor.py \
  --demo --output /tmp/struct_demo
```

### Plain Text Examples

Predict the structure of a single protein from a YAML file:

    python skills/struct-predictor/struct_predictor.py --input my_protein.yaml --output /tmp/struct_out

Run the built-in Trp-cage demo (no input file needed):

    python skills/struct-predictor/struct_predictor.py --demo --output /tmp/struct_demo

Predict a two-chain complex:

    python skills/struct-predictor/struct_predictor.py --input complex_ab.yaml --output /tmp/complex_out

## Output Structure

```
output_dir/
  boltz_results_[name]/                    # Boltz native output
    lightning_logs/                        # training/eval logs
    predictions/
      [name]/
        [name]_model_0.cif                 # predicted structure (pLDDT in B-factors)
        confidence_[name]_model_0.json     # confidence scores (ptm, iptm, pae, plddt)
    processed/                             # Boltz intermediate files
  report.md                                # primary markdown report
  viewer.html                              # self-contained 3Dmol.js 3D viewer (open in browser)
  result.json                              # machine-readable summary
  figures/
    plddt.png                              # per-residue pLDDT confidence plot
    pae.png                                # PAE inter-residue error heatmap
  reproducibility/
    commands.sh                            # exact boltz predict command used
    environment.txt                        # boltz version snapshot
```

## YAML Complex Format

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: ACDEFGHIKLMNPQRSTVWY
      msa: empty        # runs offline; replace with a path to a .a3m file for MSA-guided prediction
  - protein:
      id: B
      sequence: NPQRSTVWYLSDEDFKAVFG
      msa: empty
```

### MSA Options

| `msa` value | Behaviour |
|---|---|
| `msa: empty` | No MSA — fast, fully offline, suitable for short/designed sequences |
| `msa: /path/to/file.a3m` | Pre-computed MSA — best accuracy for natural proteins |
| *(omit field)* | Boltz errors unless `--use_msa_server` is passed at predict time |

## pLDDT Confidence Bands

| Band | pLDDT Range | Interpretation |
|------|------------|----------------|
| Very high | ≥ 90 | Backbone accurate to ~0.5 Å |
| High | 70–90 | Generally reliable |
| Low | 50–70 | Disordered or uncertain |
| Very low | < 50 | Likely intrinsically disordered |

## Demo Data

| Item | Value |
|------|-------|
| File | `skills/struct-predictor/demo_data/trpcage.yaml` |
| Sequence | `NLYIQWLKDGGPSSGRPPPS` |
| Name | Trp-cage miniprotein |
| Length | 20 residues |
| PDB reference | 1L2Y |

## Dependencies

```bash
uv pip install boltz -U          # CPU
uv pip install "boltz[cuda]" -U  # GPU (recommended)
uv pip install numpy matplotlib pyyaml
```

## Citations

- Passaro S et al. (2025) *Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction*. bioRxiv. doi:10.1101/2025.06.14.659707. PMID: 40667369; PMCID: PMC12262699.
- Wohlwend J et al. (2024) *Boltz-1: Democratizing Biomolecular Interaction Modeling*. bioRxiv. doi:10.1101/2024.11.19.624167
- Jumper J et al. (2021) *AlphaFold2 pLDDT definition*. Nature. doi:10.1038/s41586-021-03819-2
