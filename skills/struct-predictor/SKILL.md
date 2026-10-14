---
name: struct-predictor
description: Protein structure prediction with Boltz-2 (default) or OpenFold3. Accepts YAML inputs (single protein or multi-chain complex), runs
  the chosen backend offline, extracts per-residue pLDDT and PAE confidence, and writes a markdown report with figures.
license: MIT
metadata:
  version: 0.3.0
  openclaw:
    requires:
      bins:
      - python3
      anyBins:
      - boltz
      - run_openfold
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
      package: openfold3
      bins:
      - run_openfold
      comment: 'optional, for --backend openfold3: run setup_openfold once to download weights; needs a CUDA GPU'
    - kind: uv
      package: numpy
    - kind: uv
      package: matplotlib
    - kind: uv
      package: pyyaml
---

# Struct Predictor

You are the **Struct Predictor**, a specialised agent for protein structure prediction using Boltz-2 (default) or OpenFold3.

## Core Capabilities

1. **Structure Prediction**: Run Boltz-2 (default) or OpenFold3 (`--backend openfold3`, CUDA GPU required) locally on a YAML input
2. **Confidence Extraction**: Per-residue pLDDT (from CIF B-factors) and PAE matrix (from confidence JSON)
3. **Report Generation**: Markdown with pLDDT line plot, PAE heatmap, band breakdown, and reproducibility bundle
4. **Demo Mode**: Trp-cage miniprotein (20 residues, PDB 1L2Y) — runs immediately, no input required

## CLI Reference

```bash
# Single protein or multi-chain complex (YAML)
python skills/struct-predictor/struct_predictor.py \
  --input complex.yaml --output /tmp/struct_out

# Same input with OpenFold3 instead of the Boltz-2 default
python skills/struct-predictor/struct_predictor.py \
  --input complex.yaml --output /tmp/struct_out --backend openfold3

# Demo (Trp-cage miniprotein, PDB 1L2Y — no input needed)
python skills/struct-predictor/struct_predictor.py \
  --demo --output /tmp/struct_demo
```

Both backends run offline (no MSA server, no templates). The input is the Boltz-style YAML for either backend; for OpenFold3 it is converted to an OpenFold3 query JSON with `use_msas: false`.

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
  predictions/[name]/                      # Boltz native output (default backend)
    [name]_model_0.cif                     # predicted structure (pLDDT in B-factors)
    confidence_[name]_model_0.json         # confidence scores (ptm, iptm, pae, plddt)
  [name]/seed_[n]/                         # OpenFold3 native output (--backend openfold3)
    [name]_seed_[n]_sample_[k]_model.cif   # predicted structure (pLDDT in B-factors); best sample by sample_ranking_score
    [name]_seed_[n]_sample_[k]_confidences.json             # per-atom plddt, pae, pde
    [name]_seed_[n]_sample_[k]_confidences_aggregated.json  # avg_plddt, ptm, iptm, sample_ranking_score
  report.md                                # primary markdown report
  viewer.html                              # self-contained 3Dmol.js 3D viewer (open in browser)
  result.json                              # machine-readable summary
  figures/
    plddt.png                              # per-residue pLDDT confidence plot
    pae.png                                # PAE inter-residue error heatmap
  reproducibility/
    commands.sh                            # exact struct_predictor.py command used
    environment.txt                        # backend package version snapshot
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

The `msa` field is ignored by the OpenFold3 backend, which always runs with `use_msas: false`.

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

## Gotchas

- The model will want to `pip install openfold3` and run it. Do not assume that works on a GPU host: pip can pull a CUDA 13 torch that fails with "NVIDIA driver on your system is too old (found version 12080)" on CUDA 12.8 drivers. Check `nvidia-smi`, then pin a matching build, e.g. `pip install "torch==2.11.0+cu128" --index-url https://download.pytorch.org/whl/cu128 --extra-index-url https://pypi.org/simple`.
- The model will want to call `run_openfold predict` directly. Do not. Its `--use-msa-server` default is not guaranteed offline, so the skill always passes `--use-msa-server=false --use-templates=false`. Dropping those flags can send sequences to the ColabFold server.
- The model will want to say OpenFold3 works anywhere. Do not. It needs a CUDA GPU and a one-time `setup_openfold` weight download; without them stay on the Boltz-2 default.
- The model will want to treat both engines as equally licensed. Do not. Boltz code and weights are MIT. OpenFold3 code is Apache-2.0, but no licence for its model weights is stated in its README or parameters docs; confirm with the OpenFold team before commercial use.

## Dependencies

```bash
uv pip install boltz -U          # CPU
uv pip install "boltz[cuda]" -U  # GPU (recommended)
uv pip install openfold3         # optional, --backend openfold3; CUDA GPU required
setup_openfold                   # once: downloads OpenFold3 weights (~2 GB)
uv pip install numpy matplotlib pyyaml
```

## Citations

- Passaro S et al. (2025) *Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction*. bioRxiv. doi:10.1101/2025.06.14.659707. PMID: 40667369; PMCID: PMC12262699.
- Wohlwend J et al. (2024) *Boltz-1: Democratizing Biomolecular Interaction Modeling*. bioRxiv. doi:10.1101/2024.11.19.624167
- OpenFold Consortium. *OpenFold3*. https://github.com/aqlaboratory/openfold-3
- Abramson J et al. (2024) *Accurate structure prediction of biomolecular interactions with AlphaFold 3*. Nature. doi:10.1038/s41586-024-07487-w
- Jumper J et al. (2021) *AlphaFold2 pLDDT definition*. Nature. doi:10.1038/s41586-021-03819-2
