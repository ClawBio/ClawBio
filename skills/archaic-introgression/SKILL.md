---
name: archaic-introgression
version: 0.1.0
author: Manuel Corpas
description: Detect Neanderthal and Denisovan introgression segments from modern human genomes
domain: genomics
inputs:
  - format: VCF
    description: Modern human genotypes
    required: true
  - format: VCF
    description: Archaic reference genotypes (Neanderthal/Denisovan)
    required: true
outputs:
  - format: JSON
    description: Introgression segments with scores and summary statistics
  - format: BED
    description: Genomic coordinates of introgressed segments
trigger_keywords:
  - archaic introgression
  - Neanderthal DNA
  - Denisovan ancestry
  - IBDmix
  - Sprime
  - hmmix
  - introgressed segments
---

# Archaic Introgression Detector

## Trigger

**Fire when:**
- User asks about Neanderthal or Denisovan DNA in modern humans
- User wants to detect archaic introgression segments
- User mentions IBDmix, Sprime, or hmmix methods
- User has modern + archaic VCF files and wants to compare them

**Do NOT fire when:**
- User asks about general ancestry or population structure (use claw-ancestry-pca)
- User asks about pharmacogenomics or clinical variants
- User wants admixture proportions without segment-level detail

## Why This Exists

Between 1-4% of non-African modern human genomes derive from archaic hominins
(Neanderthals, Denisovans). Identifying these segments matters for understanding
human evolution, disease susceptibility, and immune adaptation. Three complementary
methods exist (IBDmix, Sprime, hmmix), each with different strengths. This skill
wraps all three behind a unified interface with a pure-Python fallback when
external binaries are unavailable.

## Core Capabilities

1. **IBDmix**: Detect segments shared IBD between modern and archaic genomes using LOD score thresholds
2. **Sprime**: Identify introgressed haplotypes without requiring an archaic reference panel
3. **hmmix**: HMM-based detection of archaic ancestry tracts
4. **Pure-Python fallback**: LOD-score-based segment calling when IBDmix binary is not installed
5. **EIGENSTRAT support**: Read .ind/.snp/.geno files (text and binary packed formats)
6. **Summary statistics**: Per-individual introgression burden, segment count, mean length

## Scope

One skill, one task: detect and report archaic introgression segments. Does not
perform downstream functional annotation of introgressed variants (chain with
vcf-annotator for that).

## Input Formats

- **Modern genotypes**: VCF (uncompressed or .vcf.gz with .tbi index)
- **Archaic genotypes**: VCF (Neanderthal/Denisovan reference panel)
- **EIGENSTRAT**: .ind + .snp + .geno files (text or binary packed format)

## Workflow

1. Parse modern VCF to extract sample names and genotype matrix (numpy array)
2. Parse archaic VCF to extract reference genotypes at matching positions
3. Identify shared variant positions between modern and archaic panels
4. Run selected method (IBDmix / Sprime / hmmix) or pure-Python LOD fallback
5. Collect IntrogressionSegment results per sample
6. Compute per-individual summary statistics
7. Write JSON report and optional BED file

## CLI Reference

```bash
# Run with IBDmix on VCF inputs
python archaic_introgression.py \
  --input modern.vcf --archaic archaic.vcf \
  --method ibdmix --output /tmp/introgression

# Run demo with synthetic data
python archaic_introgression.py --demo --output /tmp/introgression_demo

# Filter to specific samples
python archaic_introgression.py \
  --input modern.vcf --archaic archaic.vcf \
  --samples SAMPLE01,SAMPLE02 --output /tmp/introgression

# Adjust LOD threshold
python archaic_introgression.py \
  --input modern.vcf --archaic archaic.vcf \
  --lod 5.0 --output /tmp/introgression
```

## Demo

```bash
python archaic_introgression.py --demo --output /tmp/introgression_demo
```

Runs on bundled `examples/demo_modern.vcf` (3 samples, 10 SNPs on chr22) and
`examples/demo_archaic.vcf` (1 Neanderthal sample, same positions).

## Example Queries

- "How much Neanderthal DNA do I have?"
- "Detect archaic introgression in my VCF"
- "Run IBDmix on my genotype data"
- "Show me introgressed segments from Denisova"
- "Compare my genome against the Vindija Neanderthal"

## Output Structure

```
output_dir/
  introgression_results.json   # Full results with segments and summary
  segments.bed                 # BED file of introgressed regions
```

### JSON structure

```json
{
  "method": "ibdmix",
  "lod_threshold": 3.0,
  "num_samples": 3,
  "segments": [
    {
      "sample": "SAMPLE01",
      "chrom": "chr22",
      "start": 16050075,
      "end": 16051249,
      "archaic_source": "Neanderthal",
      "method": "ibdmix",
      "score": 4.2,
      "num_variants": 6,
      "length": 1174
    }
  ],
  "summary": {
    "SAMPLE01": {
      "total_segments": 1,
      "total_length_bp": 1174,
      "mean_segment_length": 1174.0,
      "mean_score": 4.2
    }
  }
}
```

## Dependencies

- Python 3.11+
- numpy
- Optional: IBDmix binary, Sprime JAR, hmmix binary
- Optional: bcftools (for indexed VCF extraction)

## Gotchas

1. **The model will want to report percentage of genome that is Neanderthal from 10 SNPs.** Do not. Demo data is too sparse for genome-wide estimates. State the segment coordinates and LOD scores only.
2. **The model will want to assume all archaic segments are Neanderthal.** Do not. Check the archaic source label. Denisovan segments have different frequency distributions in different populations.
3. **The model will want to skip the pure-Python fallback and just say IBDmix is required.** Do not. The fallback LOD scoring works well for demonstration and small datasets.
4. **The model will want to interpret introgression as harmful.** Do not. Many introgressed segments are adaptive (e.g., immune genes, altitude adaptation). Present findings neutrally.
5. **The model will want to use text .geno parsing for binary packed files.** Do not. Check for the 'GENO' magic header and switch to binary unpacking (2-bit per genotype).

## Safety

ClawBio is a research and educational tool. It is not a medical device and does
not provide clinical diagnoses. Consult a healthcare professional before making
any medical decisions.

## Agent Boundary

The agent dispatches queries and explains results. The skill executes the
computational pipeline. The agent should not attempt to reimplement IBDmix
LOD scoring outside this module.

## Chaining Partners

- **vcf-annotator**: annotate introgressed variants with ClinVar/gnomAD significance
- **equity-scorer**: assess representation of archaic ancestry detection across populations
- **claw-ancestry-pca**: combine with PCA to contextualise introgression within population structure

## Maintenance

- Review quarterly against new archaic genome releases
- Update when new IBDmix/Sprime/hmmix versions change output formats
- Deprecate if a unified tool supersedes all three methods

## Citations

- Chen L, Wolf AB, Fu W, Li L, Akey JM. Identifying and interpreting apparent Neanderthal ancestry in African individuals. Cell. 2020;180(4):677-687. (IBDmix)
- Browning SR, Browning BL, Zhou Y, et al. Analysis of human sequence data reveals two pulses of archaic Denisovan admixture. Cell. 2018;173(1):53-61. (Sprime)
- Skov L, Hui R, Shchur V, et al. Detecting archaic introgression using an unadmixed outgroup. PLoS Genet. 2018;14(9):e1007641. (hmmix)
