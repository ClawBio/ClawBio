---
name: nutrigx-advisor
description: Personalised nutrition report from consumer genetic data (23andMe, AncestryDNA, VCF) — interrogates nutritionally-relevant SNPs and generates actionable dietary guidance, all computed locally.
version: 0.1.0
author: David de Lorenzo
license: MIT
tags: [nutrigenomics, nutrition, diet, genetics, 23andme, ancestrydna, vcf]
metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: "🥦"
    homepage: https://github.com/ClawBio/ClawBio
    os: [darwin, linux]
    install:
      - kind: pip
        package: pandas numpy matplotlib seaborn
        bins: []
      - kind: pip
        package: reportlab
        bins: []
    trigger_keywords:
      - personalised nutrition
      - nutrigenomics
      - diet genetics
      - MTHFR
      - APOE
      - FTO
      - vitamin D
      - caffeine metabolism
      - lactose
      - omega-3
      - folate
      - what should I eat based on my DNA
      - vitamin absorption genetics
      - gluten
---

# NutriGx Advisor — Personalised Nutrition from Genetic Data

**Author**: David de Lorenzo (ClawBio Community)
**Requires**: Python 3.11+, pandas, numpy, matplotlib, seaborn, reportlab (optional)

---

## What This Skill Does

The NutriGx Advisor generates a **personalised nutrition report** from consumer
genetic data (23andMe, AncestryDNA raw files or VCF). It interrogates a curated
set of nutritionally-relevant SNPs drawn from GWAS Catalog, ClinVar, and
peer-reviewed nutrigenomics literature, then translates genotype calls into
actionable dietary and supplementation guidance — all computed locally.

**Key outputs**
- Markdown nutrition report with risk scores and recommendations
- Radar chart of nutrient risk profile
- Gene × nutrient heatmap
- Reproducibility bundle (`commands.sh`, `environment.yml`, SHA-256 checksums)

---

## Trigger Phrases

The Bio Orchestrator should route to this skill when the user says anything like:

- "personalised nutrition", "nutrigenomics", "diet genetics"
- "what should I eat based on my DNA"
- "nutrient metabolism", "vitamin absorption genetics"
- "MTHFR", "APOE", "FTO", "BCMO1", "VDR", "FADS1/2"
- "folate", "omega-3", "vitamin D", "caffeine metabolism", "lactose", "gluten"
- Input files: `.txt` or `.csv` (23andMe), `.csv` (AncestryDNA), `.vcf`

---

## Curated SNP Panel

### Macronutrient Metabolism

| Gene    | SNP        | Nutrient Impact                          | Evidence |
|---------|------------|------------------------------------------|----------|
| FTO     | rs9939609  | Energy balance, fat mass, carb sensitivity | Strong (GWAS) |
| PPARG   | rs1801282  | Fat metabolism, insulin sensitivity      | Moderate |
| APOA5   | rs662799   | Triglyceride response to dietary fat     | Strong |
| TCF7L2  | rs7903146  | Carbohydrate metabolism, T2D risk        | Strong |
| ADRB2   | rs1042713  | Fat oxidation, exercise × diet interaction | Moderate |

### Micronutrient Metabolism

| Gene    | SNP        | Nutrient                | Effect of risk allele            |
|---------|------------|-------------------------|----------------------------------|
| MTHFR   | rs1801133  | Folate / B12            | ↓ 5-MTHF conversion (~70%)       |
| MTHFR   | rs1801131  | Folate / B12            | ↓ enzyme activity (~30%)         |
| MTR     | rs1805087  | B12 / homocysteine      | ↑ homocysteine risk              |
| BCMO1   | rs7501331  | Beta-carotene → Vitamin A | ↓ conversion (~50%)             |
| BCMO1   | rs12934922 | Beta-carotene → Vitamin A | ↓ conversion (compound het)    |
| VDR     | rs2228570  | Vitamin D absorption    | ↓ VDR function                   |
| VDR     | rs731236   | Vitamin D               | ↓ bone mineral density response  |
| GC      | rs4588     | Vitamin D binding       | ↑ deficiency risk                |
| SLC23A1 | rs33972313 | Vitamin C transport     | ↓ renal reabsorption             |
| ALPL    | rs1256335  | Vitamin B6              | ↓ alkaline phosphatase activity  |

### Omega-3 / Fatty Acid Metabolism

| Gene    | SNP        | Nutrient             | Effect                          |
|---------|------------|----------------------|---------------------------------|
| FADS1   | rs174546   | LC-PUFA synthesis    | ↑/↓ EPA/DHA from ALA            |
| FADS2   | rs1535     | LC-PUFA synthesis    | Modulates omega-6:omega-3 ratio |
| ELOVL2  | rs953413   | DHA synthesis        | ↓ elongation of EPA→DHA         |
| APOE    | rs429358   | Saturated fat response | ε4 → ↑ LDL-C on high SFA diet |
| APOE    | rs7412     | Saturated fat response | Combined with rs429358 for ε typing |

### Caffeine & Alcohol

| Gene    | SNP        | Compound    | Effect                         |
|---------|------------|-------------|--------------------------------|
| CYP1A2  | rs762551   | Caffeine    | Slow/Fast metaboliser          |
| AHR     | rs4410790  | Caffeine    | Modulates CYP1A2 induction     |
| ADH1B   | rs1229984  | Alcohol     | Acetaldehyde accumulation risk |
| ALDH2   | rs671       | Alcohol     | Asian flush / toxicity risk    |

### Food Sensitivities

| Gene    | SNP        | Sensitivity          | Effect                          |
|---------|------------|----------------------|---------------------------------|
| MCM6    | rs4988235  | Lactose intolerance  | Non-persistence of lactase      |
| HLA-DQ2 | Proxy SNPs | Coeliac / gluten     | HLA-DQA1/DQB1 risk haplotypes   |

### Antioxidant & Detoxification

| Gene    | SNP        | Pathway              | Effect                          |
|---------|------------|----------------------|---------------------------------|
| SOD2    | rs4880     | Manganese SOD        | ↓ mitochondrial antioxidant     |
| GPX1    | rs1050450  | Selenium / GSH-Px    | ↓ glutathione peroxidase        |
| GSTT1   | Deletion   | Glutathione-S-trans  | Null genotype → ↑ oxidative risk|
| NQO1    | rs1800566  | Coenzyme Q10         | ↓ CoQ10 regeneration            |
| COMT    | rs4680     | Catechol / B vitamins | Met/Val → methylation load     |

---

## Algorithm

### 1. Input Parsing (`parse_input.py`)

Accepts:
- 23andMe `.txt` or `.csv` (tab-separated: rsid, chromosome, position, genotype)
- AncestryDNA `.csv`
- Standard VCF (extracts GT field)

Auto-detects format from header lines. Normalises alleles to forward strand using
a hard-coded reference table (avoids requiring external databases).

### 2. Genotype Extraction (`extract_genotypes.py`)

For each SNP in the panel:
1. Look up rsid in parsed data
2. Return genotype string (e.g. `"AT"`, `"TT"`, `"AA"`)
3. Flag as `"NOT_TESTED"` if absent (common for chip-to-chip variation)

### 3. Risk Scoring (`score_variants.py`)

Each SNP is scored on a **0 / 0.5 / 1.0** scale:
- `0.0` — homozygous reference (lowest risk)
- `0.5` — heterozygous
- `1.0` — homozygous risk allele

Composite **Nutrient Risk Scores** (0–10) are computed per nutrient domain by
summing weighted SNP scores. Weights are derived from reported effect sizes
(beta coefficients or OR) in the primary literature.

Risk categories:
- **0–3**: Low risk — standard dietary advice applies
- **3–6**: Moderate risk — dietary optimisation recommended
- **6–10**: Elevated risk — consider testing and targeted supplementation

> **Important caveat**: These are polygenic risk indicators based on common
> variants. They are not diagnostic. Rare pathogenic variants (e.g. MTHFR
> compound heterozygosity with high homocysteine) require clinical confirmation.

### 4. Report Generation (`generate_report.py`)

Outputs a structured Markdown report with:
- Executive summary (top 3 personalised findings)
- Per-nutrient sections: genotype table → interpretation → recommendation
- Radar chart (matplotlib) of nutrient risk scores
- Gene × nutrient heatmap (seaborn)
- Supplement interactions table
- Disclaimer section
- Reproducibility block

### 5. Reproducibility Bundle (`repro_bundle.py`)

Exports to the output directory (not committed to the repo):
- `commands.sh` — full CLI to reproduce analysis
- `environment.yml` — pinned conda environment
- `checksums.txt` — SHA-256 checksums of input and output files
- `provenance.json` — timestamp and ClawBio version tag

---

## Usage

```bash
# From 23andMe raw data
openclaw "Generate my personalised nutrition report from genome.csv"

# From VCF
openclaw "Run NutriGx analysis on variants.vcf and flag any folate pathway risks"

# Targeted query
openclaw "What does my APOE status mean for my saturated fat intake?"

# Generate a random demo patient and run the report
python examples/generate_patient.py --run
```

---

## File Structure

```
skills/nutrigx-advisor/
├── SKILL.md                      ← this file (agent instructions)
├── nutrigx_advisor.py            ← main entry point
├── parse_input.py                ← multi-format parser
├── extract_genotypes.py          ← SNP lookup engine
├── score_variants.py             ← risk scoring algorithm
├── generate_report.py            ← Markdown + figures
├── repro_bundle.py               ← reproducibility export
├── .gitignore
├── data/
│   └── snp_panel.json            ← curated SNP definitions
├── tests/
│   ├── synthetic_patient.csv     ← fixed 23andMe-format test data (for pytest)
│   └── test_nutrigx.py           ← pytest suite
└── examples/
    ├── generate_patient.py       ← random patient generator (demo use)
    ├── data/                     ← generated patient files land here (gitignored)
    └── output/
        ├── nutrigx_report.md     ← pre-rendered demo report
        ├── nutrigx_radar.png     ← demo radar chart (nutrient risk profile)
        └── nutrigx_heatmap.png   ← demo gene × nutrient heatmap
```

> **Note**: Runtime output directories and randomly generated patient files are
> excluded from version control via `.gitignore`. Only the pre-rendered demo
> report in `examples/output/` is committed.

---

## Privacy

All computation runs **locally**. No genetic data is transmitted. Input files are
read-only; no raw genotype data appears in any output file (reports contain only
gene names, SNP IDs, and risk categories).

---

## Limitations & Disclaimer

1. **Not a medical device.** This skill provides educational, research-oriented
   nutrigenomics analysis. It does not constitute medical advice.
2. **Common variants only.** The panel covers SNPs with MAF > 1% in at least one
   major population. Rare pathogenic variants are out of scope.
3. **Population context.** Effect sizes are predominantly derived from European
   GWAS cohorts. Risk estimates may not generalise equally across all ancestries.
4. **Gene–environment interaction.** Genetic risk scores interact with baseline
   diet, lifestyle, microbiome, and epigenetic state. A "high risk" score does not
   mean a nutrient deficiency is present — it means the individual may benefit from
   monitoring.
5. **Simpson's Paradox note.** Population-level associations used to derive weights
   may not reflect individual trajectories (see Corpas 2025, *Nutrigenomics and
   the Ecological Fallacy*).

---

## Roadmap

- [ ] **v0.2**: Microbiome × genotype interaction module (16S rRNA input)
- [ ] **v0.3**: Longitudinal tracking — compare reports across time
- [ ] **v0.4**: HLA typing for immune-mediated food reactions (coeliac, gluten sensitivity)
- [ ] **v0.5**: Integration with NeoTree neonatal data for maternal nutrition risk scoring
- [ ] **v1.0**: Multi-omics integration (metabolomics + genomics + dietary recall)

---

## References

Key literature underpinning the SNP panel and scoring algorithm:

- Corbin JM & Ruczinski I (2023). Nutrigenomics: current state and future directions. *Annu Rev Nutr*.
- Fenech M et al. (2011). Nutrigenetics and nutrigenomics: viewpoints on the current status. *J Nutrigenet Nutrigenomics*.
- Stover PJ (2006). Influence of human genetic variation on nutritional requirements. *Am J Clin Nutr*.
- Phillips CM (2013). Nutrigenetics and metabolic disease: current status and implications for personalised nutrition. *Nutrients*.
- Minihane AM et al. (2015). APOE genotype, cardiovascular risk and responsiveness to dietary fat manipulation. *Proc Nutr Soc*.
- Frayling TM et al. (2007). A common variant in the FTO gene is associated with body mass index. *Science*.
- Pare G et al. (2010). MTHFR variants and cardiovascular risk. *Hum Genet*.
- Lecerf JM & de Lorgeril M (2011). Dietary cholesterol: from physiology to cardiovascular risk. *Br J Nutr*.
- Tanaka T et al. (2009). Genome-wide association study of plasma polyunsaturated fatty acids in the InCHIANTI Study. *PLoS Genet* (FADS1/2).
- Cornelis MC et al. (2006). Coffee, CYP1A2 genotype, and risk of myocardial infarction. *JAMA*.

---

## Contributing

The SNP panel (`data/snp_panel.json`) is maintained by the skill author.
To suggest additions or corrections, contact David de Lorenzo directly via
GitHub ([@drdaviddelorenzo](https://github.com/drdaviddelorenzo)) or open
an issue tagging him in the main ClawBio repository.
