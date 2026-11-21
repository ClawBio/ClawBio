# Skill highlights

Longer write-ups of a few skills. Each skill's own `SKILL.md` is the authoritative reference; the full list is in the [skill directory](https://clawbio.ai/skills.html).

## PharmGx Reporter — *Personal Scale*

Generates a pharmacogenomic report from consumer genetic data (23andMe, AncestryDNA):

- Parses raw genetic data (auto-detects format, including gzip)
- Extracts **31 pharmacogenomic SNPs** across **12 genes** (CYP2C19, CYP2D6, CYP2C9, VKORC1, SLCO1B1, DPYD, TPMT, UGT1A1, CYP3A5, CYP2B6, NUDT15, CYP1A2)
- Calls star alleles and determines metabolizer phenotypes
- Looks up **CPIC drug recommendations** for **51 medications**
- Zero dependencies. Runs in **< 1 second**.

```bash
python pharmgx_reporter.py --input demo_patient.txt --output report
```

**Demo result**: CYP2D6 \*4/\*4 (Poor Metabolizer) → **10 drugs AVOID** (codeine, tramadol, 7 TCAs, tamoxifen), 20 caution, 21 standard.

> ~7% of people are CYP2D6 Poor Metabolizers — codeine gives them zero pain relief. ~0.5% carry DPYD variants where standard 5-FU dose can be lethal. This skill catches both.

## Drug Photo — *Personal Scale*

Snap a photo of any medication in Telegram. ClawBio identifies the drug from the packaging and returns a personalised dosage card against your own genotype.

- Claude vision extracts drug name and visible dose from the photo
- Cross-references your 23andMe genotype against 31 PGx SNPs
- Four-tier classification: **STANDARD DOSING / USE WITH CAUTION / AVOID / INSUFFICIENT DATA**
- Correct VKORC1 complement-strand handling (23andMe reports minus strand for rs9923231)
- Works for warfarin, clopidogrel, codeine, simvastatin, tamoxifen, sertraline, and 20+ others

```bash
python pharmgx_reporter.py --drug warfarin --dose "5mg" --input my_23andme.txt --output report
```

> No command needed in Telegram — send any medication photo and RoboTerri triggers the skill automatically.

## GWAS Lookup — *Population Scale*

Federated variant query across nine genomic databases in a single command:

| Database | What you get |
|----------|-------------|
| GWAS Catalog | Genome-wide significant associations |
| gnomAD | Allele frequencies across 125,748 exomes |
| ClinVar | Clinical significance and condition links |
| Open Targets | Disease-gene evidence scores |
| Ensembl | Functional annotation, regulatory impact |
| GTEx | eQTL data, tissue-specific expression effects |
| LDlink | Linkage disequilibrium across 26 populations |
| UK Biobank PheWAS | Phenome-wide associations across 4,000+ traits |
| LOVD | Variant pathogenicity database |

```bash
python gwas_lookup.py --rsid rs3798220 --output report
python gwas_lookup.py --demo --output /tmp/gwas_lookup_demo
```

## UKB Navigator — *Research Scale*

Semantic search across the UK Biobank schema. Ask in plain English what UK Biobank measures about any phenotype — get field IDs, descriptions, data types, participant counts, and linked publications back instantly.

```bash
python ukb_navigator.py --query "grip strength"   --output report
python ukb_navigator.py --field 21001              --output report   # BMI
python ukb_navigator.py --demo                     --output /tmp/ukb_demo
```

Built on a ChromaDB embedding of the full UKB Data Showcase (22,000+ fields).

## Ancestry PCA — *Population Scale*

Runs principal component analysis on your cohort against the SGDP reference panel (345 samples, 164 global populations):

- Contig normalisation (chr1 vs 1)
- IBD removal (related individuals filtered)
- Common biallelic SNPs only
- Confidence ellipses per population
- Publication-quality **4-panel figure** generated instantly

```bash
python ancestry_pca.py --demo --output ancestry_report
```

**Demo result**: 736 Peruvian samples across 28 indigenous populations. Amazonian groups (Matzes, Awajun, Candoshi) sit in genetic space that no SGDP population occupies — genuinely underrepresented, not just in GWAS, but in the reference panels themselves.

## Semantic Similarity Index — *Systemic Scale*

Computes a Semantic Isolation Index for diseases using 13.1M PubMed abstracts and PubMedBERT embeddings (768-dim):

- **SII** (Semantic Isolation Index): higher = more isolated in literature
- **KTP** (Knowledge Transfer Potential): higher = more cross-disease spillover
- **RCC** (Research Clustering Coefficient): diversity of research approaches
- **Temporal Drift**: how research focus evolves over time
- Publication-quality **4-panel figure**

```bash
python semantic_sim.py --demo --output sem_report
```

**Key finding**: Neglected tropical diseases are **+38% more semantically isolated** (P < 0.0001, Cohen's d = 0.84). 14 of the 25 most isolated diseases are Global South priority conditions. Knowledge silos kill innovation — a malaria immunology breakthrough could help leishmaniasis, but the literatures don't talk to each other.

> Corpas et al. (2026). *HEIM: Health Equity Index for Measuring structural bias in biomedical research.* Under review.

## Galaxy Bridge — *Research Scale*

ClawBio indexes **8,000+ bioinformatics tools** from [usegalaxy.org](https://usegalaxy.org) via the Galaxy Bridge skill. Search by natural language, inspect tool schemas, and execute remotely — all from the CLI.

```bash
# Search Galaxy tools by keyword
python skills/galaxy-bridge/galaxy_bridge.py --search "metagenomics"

# Browse all 86 tool categories
python skills/galaxy-bridge/galaxy_bridge.py --list-categories

# Run a tool on Galaxy (requires GALAXY_API_KEY)
python skills/galaxy-bridge/galaxy_bridge.py --run fastqc --input reads.fq.gz --output results/

# Demo mode (offline, no API key)
python skills/galaxy-bridge/galaxy_bridge.py --demo
```

**Cross-platform chaining**: Galaxy VEP annotates variants → ClawBio PharmGx generates dosage report. Galaxy Kraken2 classifies reads → ClawBio metagenomics profiler. Neither can do this alone.

Built on [BioBlend](https://bioblend.readthedocs.io/) (Galaxy Python SDK). Developed in collaboration with the Galaxy ML SIG.
