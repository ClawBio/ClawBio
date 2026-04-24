---
name: wes-clinical-report-es
description: Generates professional clinical PDF reports in Spanish from WES (Whole Exome Sequencing) data with clinical interpretation,
  pharmacogenomic alerts, and follow-up recommendations.
license: PROPRIETARY
metadata:
  version: 1.0.0
  author: Manuel Corpas
  tags:
  - WES
  - exome
  - clinical-report
  - spanish
  - pharmacogenomics
  - PDF
  - ANNOVAR
  inputs:
  - name: WES markdown report
    format: markdown (.md)
    required: true
    description: Structured WES report with sections 1-7 (Exome Summary through Methods)
  - name: Logo left (Predice)
    format: image (JPG/PNG)
    required: false
    description: Predice institutional logo for cover and header
  - name: Logo right (Inbiomedic)
    format: image (JPG/PNG)
    required: false
    description: Inbiomedic institutional logo for cover and header
  outputs:
  - name: Clinical PDF report (Spanish)
    format: PDF (A4)
    description: Professional clinical report in Spanish with interpretation, tables, and disclaimer
  openclaw:
    requires:
      bins:
      - python3
    always: false
    emoji: 🧬
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    min_python: '3.9'
    dependencies:
    - reportlab
    - pandas
    private: true
    trigger_keywords:
    - informe clinico WES
    - clinical report spanish
    - exome PDF report
    - Novogene report
    - Predice
    - Inbiomedic
---

# Informe Clinico WES en Espanol

Skill for generating professional clinical PDF reports in Spanish from
whole exome sequencing (WES) data. Designed for Novogene WES data
(GATK + ANNOVAR pipeline) but adaptable to any WES pipeline with
equivalent annotations. Translates all section headings, table headers,
cell values, and interpretive text from English markdown input to Spanish
PDF output.

## Trigger

**Fire this skill when the user says any of:**
- "informe clinico WES"
- "generar informe exoma en espanol"
- "clinical report in Spanish"
- "Novogene report Spanish"
- "exome PDF report Spanish"
- "Predice/Inbiomedic report"

**Do NOT fire when:**
- User asks for an English report (use `wes-clinical-report-en`)
- User asks for variant annotation only (use `variant-annotation`)
- User asks for ACMG classification only (use `clinical-variant-reporter`)

## Scope

One skill, one task: convert WES markdown reports (English input) into
professional Spanish-language clinical PDFs with full translation and
clinical interpretation.

## Workflow

1. Parse WES markdown report (structured sections 1-7)
2. Extract KPI metrics from Exome Summary
3. Extract pathogenic variants, PGx alerts, rare damaging variants
4. Translate all headings, table headers, cell values, and body text to Spanish
5. Build interpretive summary paragraph in Spanish
6. Render all sections as styled PDF with clinical tables
7. Add ancestry estimation (section 8) if data available
8. Add limitations section (section 9) in Spanish
9. Add Spanish disclaimer and report metadata
10. Output PDF to specified directory

## Capabilities

1. **Clinical interpretation summary (Spanish)**: key findings, high-risk PGx
   alerts, prioritised rare variants, clinical follow-up recommendations.
2. **Full translation layer**: section headings, table headers, metric labels,
   zygosity, classification, consequence, clinical effects, traits, and
   interpretations all translated via lookup dictionaries.
3. **Clinically significant variants**: ClinVar P/LP, ACMG SF v3.2,
   cancer predisposition panel, conflicting variants.
4. **Pharmacogenomics**: CPIC star alleles, clinical effects, affected
   medications with contextualised high-risk alerts.
5. **Fitness and nutrition traits**: genotypes with evidence grades
   (Corpas et al. 2021).
6. **Rare damaging variant prioritisation**: REVEL, CADD, gnomAD AF.
7. **Institutional logos**: Predice (left) and Inbiomedic (right) on cover
   and header.

## Example Output

```
Pagina 1 (portada):
  [Logo Predice]                    [Logo Inbiomedic]
  +---------------------------------------------+
  |  Secuenciacion del Exoma Completo  [SampleN]|
  |  Plataforma / Referencia / Fecha            |
  +---------------------------------------------+
  [KPIs: SNPs Totales | Missense | Stopgain | Raras Patog. | ClinVar P/PP]

  Interpretacion de resultados
  (parrafo interpretativo clinico auto-generado en espanol)

Paginas 2+:
  1. Resumen del Exoma
  2. Variantes de Significancia Clinica
  3. Farmacogenomica
  4. Rasgos de Aptitud Fisica y Nutricion
  5. Variantes Raras Patogenicas Priorizadas
  6. Contexto de Enfermedad y Vias Metabolicas
  7. Metodologia
  8. Estimacion de Ancestria
  9. Limitaciones
  [Aviso legal / Disclaimer]
```

## Usage

```bash
# Generate reports for all samples
python skills/wes-clinical-report-es/wes_clinical_report_es.py \
  --report-dir /path/to/REPORTS/ \
  --output-dir /path/to/PDF-ES/ \
  --logo-left /path/to/logo_predice.jpg \
  --logo-right /path/to/logo_inbiomedic.jpg

# Generate report for a single sample
python skills/wes-clinical-report-es/wes_clinical_report_es.py \
  --report-dir /path/to/REPORTS/ \
  --output-dir /path/to/PDF-ES/ \
  --samples Sample3

# Demo with default Novogene data
python skills/wes-clinical-report-es/wes_clinical_report_es.py --demo
```

## Input format

The skill consumes WES reports in markdown format generated by the
analysis pipeline (scripts 02-12 in `ANALYSIS/SCRIPTS/`). Each markdown
report must follow this structure:

```markdown
# Whole Exome Sequencing Report: SampleN
> **Project** ... | **Platform** ... | ...
## 1. Exome Summary
## 2. Clinically Significant Variants
## 3. Pharmacogenomics
## 4. Fitness and Nutrition Traits
## 5. Prioritised Rare Damaging Variants
## 6. Disease and Pathway Context
## 7. Methods
```

## Gotchas

1. **Input is English, output is Spanish**: the markdown reports are in
   English. The skill translates to Spanish via dictionary lookups. If a
   term is missing from the translation dictionaries, it appears in English
   in the PDF. Check untranslated terms in output.
2. **Logo paths must exist**: if Predice/Inbiomedic logo files are missing,
   the report still generates but without institutional branding. The script
   silently skips missing logos.
3. **Table truncation**: tables with more than 20 rows are truncated in
   the PDF with a Spanish note to consult TSV files. Do not assume all
   data is visible in the PDF.
4. **ClinVar classifications are time-sensitive**: the report reflects
   ClinVar state at annotation time. Classifications are not permanent.
5. **PGx star alleles from SNVs only**: CYP2D6 CNV analysis is not
   included. Do not claim complete metaboliser phenotyping.
6. **Hardcoded project ID**: the footer references project
   X202SC26016276-Z01-F001. Update for other projects.

## Safety

ClawBio is a research and educational tool. It is not a medical device
and does not provide clinical diagnoses. Consult a healthcare professional
before making any medical decisions.

The Spanish disclaimer in the PDF reads: "Este informe se genera
exclusivamente con fines de investigacion y educacion. No es un informe
de diagnostico clinico."

## Agent Boundary

The agent dispatches and explains; the skill executes. The agent should
not modify PDF generation logic inline. All report customisation goes
through CLI flags.

## Chaining Partners

- `variant-annotation`: upstream VCF annotation feeding markdown reports
- `clinical-variant-reporter`: ACMG classification for deeper analysis
- `wes-clinical-report-en`: English language version of the same report

## Maintenance

- Review cadence: quarterly (aligned with ClinVar release cycle)
- Staleness signals: ClinVar version drift, CPIC guideline updates,
  new Spanish medical terminology standards
- Deprecation: if WES is superseded by WGS-only clinical pipelines

## Requirements

- Python 3.9+
- reportlab >= 4.0
- WES markdown reports (see input format above)
- Institutional logos in JPG/PNG (optional)

## Privacy

This skill is **private** and not included in the ClawBio public catalog.
It contains institutional logos and clinical report templates that should
not be distributed publicly.

## References

- Corpas et al. (2021) "Whole Genome Interpretation for a Family of Five"
  *Frontiers in Genetics* 12:535123
- CPIC guidelines for pharmacogenomics
- ClinVar / gnomAD / OMIM / COSMIC / KEGG for variant annotation
