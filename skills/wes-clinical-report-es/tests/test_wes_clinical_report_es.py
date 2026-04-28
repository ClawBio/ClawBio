#!/usr/bin/env python3
"""Tests for the Spanish WES Clinical Report Generator (red/green TDD)."""

import sys
import textwrap
from pathlib import Path

import pytest

# Add skill directory to path
SKILL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL_DIR))


# ── Sample markdown fixture ─────────────────────────────────────────────

SAMPLE_MD = textwrap.dedent("""\
# Whole Exome Sequencing Report: Sample1
> **Project** X202SC26016276-Z01-F001 | **Platform** Illumina NovaSeq PE150 | **Capture** Xplus WES (60.5 Mb) | **Reference** GRCh38/hg38

---

## 1. Exome Summary

| Metric | Count |
|--------|-------|
| Total SNP variants | 25,432 |
| Missense | 11,234 |
| Synonymous | 12,100 |
| Stopgain | 85 |
| Frameshift | 42 |
| Splicing | 120 |
| Loss-of-function (stopgain + frameshift) | 127 |
| Rare coding (gnomAD < 1%) | 3,210 |
| Rare + computationally damaging | 245 |
| ClinVar Pathogenic / Likely Pathogenic | 8 |

The homozygosity to heterozygosity ratio is **2.35**, which is above the expected ~1.5 for outbred populations.

## 2. Clinically Significant Variants

Sample1 carries 8 variant(s) classified as Pathogenic or Likely Pathogenic in ClinVar:

| Gene | Variant | Zygosity | Classification | Consequence | Associated Condition |
|------|---------|----------|----------------|-------------|---------------------|
| BRCA2 | c.5946delT | Het | Pathogenic | frameshift deletion | Hereditary breast/ovarian cancer |
| MUTYH | c.536A>G | Het | Pathogenic | missense SNV | MUTYH-associated polyposis |

An additional 12 variant(s) have conflicting or uncertain classifications.

ACMG SF v3.2 actionable genes: 156 coding variants identified across 73 medically actionable genes, of which 2 have ClinVar P/LP classification.

Cancer predisposition panel: 3 with P/LP classification across 95 cancer predisposition genes.

## 3. Pharmacogenomics

The following pharmacogenomic markers were identified from CPIC-defined star-allele positions. Variants are reported where the genotype differs from the reference allele.

| Gene | Variant | Allele | Zygosity | Clinical Effect | Affected Medications |
|------|---------|--------|----------|-----------------|---------------------|
| CYP2D6 | rs3892097 | *4 | Het | Slow metaboliser | codeine, tramadol, tamoxifen |
| NAT2 | rs1801280 | *5 | Hom | Slow acetylator | isoniazid, hydralazine |

## 4. Fitness and Nutrition Traits

Genotypes at positions associated with fitness and nutrition traits (Corpas et al. 2021, Tables 3-5). Only markers captured by the WES panel and with non-reference genotypes are shown.

Evidence grades: A = strong replication, B = moderate, C = preliminary.

### Fitness

| Gene | Variant | Trait | Interpretation | Ev. |
|------|---------|-------|----------------|-----|
| ACTN3 | rs1815739 | Muscle fibre type (power vs endurance) | XX - endurance phenotype | A |

### Nutrition

| Gene | Variant | Trait | Interpretation | Ev. |
|------|---------|-------|----------------|-----|
| MTHFR | rs1801133 | Folate metabolism (C677T) | CT - 35% reduced | A |

## 5. Prioritised Rare Damaging Variants

245 variants pass all filters: coding, rare (gnomAD AF < 0.01), and computationally predicted damaging (CADD > 20 or REVEL > 0.5). Top 15 ranked by pathogenicity prediction score:

| Gene | Variant | Consequence | Zygosity | REVEL | CADD | gnomAD AF | OMIM Disease |
|------|---------|-------------|----------|-------|------|-----------|-------------|
| ABCA4 | c.5882G>A | missense SNV | Het | 0.92 | 33.0 | 0.0012 | Stargardt disease |
| GJB2 | c.35delG | frameshift deletion | Hom | - | 35.0 | 0.0089 | Deafness, autosomal recessive |

## 6. Disease and Pathway Context

Across the full variant set: 1,245 variants map to OMIM disease entries, 456 overlap GWAS Catalog associations, and 89 have COSMIC somatic mutation records.

KEGG pathways enriched in rare coding variants:
- hsa04010: MAPK signalling pathway (12 variants)
- hsa04151: PI3K-Akt signalling pathway (9 variants)

## 7. Methods

Whole exome sequencing was performed on an Illumina NovaSeq 6000 platform using 150 bp paired-end reads with the Xplus capture kit (60.5 Mb target region). Reads were aligned to the GRCh38/hg38 reference genome using BWA-MEM. Variant calling was performed with GATK HaplotypeCaller v4.3.0 following GATK Best Practices. Functional annotation was performed with ANNOVAR, incorporating ClinVar (2024), gnomAD v3.1.2 (9 population groups), COSMIC, OMIM, SIFT, PolyPhen-2, CADD, REVEL, and 15 additional databases. Pharmacogenomic analysis used CPIC star-allele definitions with evidence enrichment from the ClinPGx API (PharmGKB). Fitness and nutrition trait interpretation followed the evidence framework of Corpas et al. (2021) Frontiers in Genetics 12:535123. Variant prioritisation applied sequential filters: coding consequence, population frequency (gnomAD AF < 0.01), computational pathogenicity (CADD > 20 or REVEL > 0.5).

*Report prepared by ClawBio WES Analysis Pipeline on 2026-04-05*

> **Disclaimer**: This report is generated for research and educational purposes only. It is not a clinical diagnostic report.
""")


@pytest.fixture
def sample_md_path(tmp_path):
    """Write sample markdown to a temp file and return its path."""
    md_file = tmp_path / "Sample1_WES_Report.md"
    md_file.write_text(SAMPLE_MD, encoding="utf-8")
    return md_file


def _require_reportlab():
    pytest.importorskip("reportlab")


# ── Markdown parser tests ───────────────────────────────────────────────

class TestParseMarkdownReport:
    def test_extracts_sample_id(self, sample_md_path):
        from wes_clinical_report_es import parse_markdown_report
        report = parse_markdown_report(sample_md_path)
        assert report["sample_id"] == "Sample1"

    def test_extracts_metadata(self, sample_md_path):
        from wes_clinical_report_es import parse_markdown_report
        report = parse_markdown_report(sample_md_path)
        assert "Project" in report["metadata"]

    def test_extracts_all_sections(self, sample_md_path):
        from wes_clinical_report_es import parse_markdown_report
        report = parse_markdown_report(sample_md_path)
        assert len(report["sections"]) == 7

    def test_section_headings_preserved(self, sample_md_path):
        from wes_clinical_report_es import parse_markdown_report
        report = parse_markdown_report(sample_md_path)
        headings = [s["heading"] for s in report["sections"]]
        assert "1. Exome Summary" in headings
        assert "7. Methods" in headings


class TestParseTable:
    def test_parse_simple_table(self):
        from wes_clinical_report_es import parse_table
        table_text = textwrap.dedent("""\
        | Metric | Count |
        |--------|-------|
        | Total SNP variants | 25,432 |
        | Missense | 11,234 |
        """)
        headers, rows = parse_table(table_text)
        assert headers == ["Metric", "Count"]
        assert len(rows) == 2

    def test_empty_table_returns_empty(self):
        from wes_clinical_report_es import parse_table
        headers, rows = parse_table("")
        assert headers == []
        assert rows == []


class TestExtractTablesAndText:
    def test_separates_text_and_tables(self):
        from wes_clinical_report_es import extract_tables_and_text
        body = textwrap.dedent("""\
        Some intro text.

        | Col1 | Col2 |
        |------|------|
        | a    | b    |

        More text after.
        """)
        parts = extract_tables_and_text(body)
        types = [p[0] for p in parts]
        assert "text" in types
        assert "table" in types


# ── Translation tests ───────────────────────────────────────────────────

class TestTranslation:
    def test_translate_heading(self):
        from wes_clinical_report_es import translate_heading
        result = translate_heading("1. Exome Summary")
        assert result == "1. Resumen del Exoma"

    def test_translate_header(self):
        from wes_clinical_report_es import translate_header
        assert translate_header("Gene") == "Gen"
        assert translate_header("Classification") == "Clasificacion"

    def test_translate_cell_zygosity(self):
        from wes_clinical_report_es import translate_cell
        assert translate_cell("Zygosity", "Het") == "Het"
        assert translate_cell("Zygosity", "Hom") == "Hom"

    def test_translate_cell_classification(self):
        from wes_clinical_report_es import translate_cell
        assert translate_cell("Classification", "Pathogenic") == "Patogenica"

    def test_translate_cell_consequence(self):
        from wes_clinical_report_es import translate_cell
        result = translate_cell("Consequence", "missense SNV")
        assert result == "SNV de cambio de sentido"

    def test_translate_cell_unknown_passes_through(self):
        from wes_clinical_report_es import translate_cell
        result = translate_cell("Gene", "BRCA2")
        assert result == "BRCA2"

    def test_translate_metric(self):
        from wes_clinical_report_es import translate_metric
        result = translate_metric("Total SNP variants")
        assert result == "Variantes SNP totales"


# ── KPI extraction tests ───────────────────────────────────────────────

class TestExtractExomeKPIs:
    def test_extracts_kpi_values(self, sample_md_path):
        from wes_clinical_report_es import parse_markdown_report, extract_exome_kpis_es
        report = parse_markdown_report(sample_md_path)
        kpis = extract_exome_kpis_es(report["sections"][0]["body"])
        assert "SNPs Totales" in kpis
        assert "Missense" in kpis
        assert "ClinVar P/PP" in kpis


# ── Data extraction tests ──────────────────────────────────────────────

class TestExtractPathogenicVariants:
    def test_finds_pathogenic_variants(self, sample_md_path):
        from wes_clinical_report_es import parse_markdown_report, _extract_pathogenic_variants
        report = parse_markdown_report(sample_md_path)
        variants = _extract_pathogenic_variants(report)
        assert len(variants) == 2
        assert variants[0]["gene"] == "BRCA2"

    def test_variant_has_required_fields(self, sample_md_path):
        from wes_clinical_report_es import parse_markdown_report, _extract_pathogenic_variants
        report = parse_markdown_report(sample_md_path)
        variants = _extract_pathogenic_variants(report)
        for v in variants:
            assert "gene" in v
            assert "variant" in v
            assert "zygosity" in v
            assert "condition" in v


class TestExtractPgxAlerts:
    def test_finds_pgx_entries(self, sample_md_path):
        from wes_clinical_report_es import parse_markdown_report, _extract_pgx_alerts
        report = parse_markdown_report(sample_md_path)
        alerts = _extract_pgx_alerts(report)
        assert len(alerts) == 2
        assert alerts[0]["gene"] == "CYP2D6"


class TestExtractMetrics:
    def test_finds_hom_het_ratio(self, sample_md_path):
        from wes_clinical_report_es import parse_markdown_report, _extract_metrics
        report = parse_markdown_report(sample_md_path)
        metrics = _extract_metrics(report)
        assert "hom_het_ratio" in metrics
        assert abs(metrics["hom_het_ratio"] - 2.35) < 0.01


# ── PDF generation tests ───────────────────────────────────────────────

class TestBuildPdf:
    def test_generates_pdf_file(self, sample_md_path, tmp_path):
        _require_reportlab()
        from wes_clinical_report_es import build_sample_pdf_es
        output = tmp_path / "Sample1_Informe_Clinico_WES.pdf"
        result = build_sample_pdf_es(sample_md_path, output)
        assert result.exists()
        assert result.stat().st_size > 1000

    def test_pdf_is_valid(self, sample_md_path, tmp_path):
        _require_reportlab()
        from wes_clinical_report_es import build_sample_pdf_es
        output = tmp_path / "Sample1_Informe_Clinico_WES.pdf"
        build_sample_pdf_es(sample_md_path, output)
        with open(output, "rb") as f:
            header = f.read(5)
        assert header == b"%PDF-"

    def test_pdf_contains_sample_id(self, sample_md_path, tmp_path):
        _require_reportlab()
        from wes_clinical_report_es import build_sample_pdf_es
        output = tmp_path / "Sample1_Informe_Clinico_WES.pdf"
        build_sample_pdf_es(sample_md_path, output)
        content = output.read_bytes()
        assert b"Sample1" in content


# ── Interpretation paragraph tests ──────────────────────────────────────

class TestBuildInterpretation:
    def test_returns_flowables(self, sample_md_path):
        _require_reportlab()
        from wes_clinical_report_es import (
            parse_markdown_report, build_interpretation_paragraph, build_styles
        )
        report = parse_markdown_report(sample_md_path)
        styles = build_styles()
        elements = build_interpretation_paragraph(report, styles)
        assert len(elements) > 0

    def test_interpretation_in_spanish(self, sample_md_path):
        _require_reportlab()
        from wes_clinical_report_es import (
            parse_markdown_report, build_interpretation_paragraph, build_styles
        )
        from reportlab.platypus import Paragraph
        report = parse_markdown_report(sample_md_path)
        styles = build_styles()
        elements = build_interpretation_paragraph(report, styles)
        paragraphs = [e for e in elements if isinstance(e, Paragraph)]
        text = " ".join(p.text for p in paragraphs)
        # Should contain Spanish
        assert "exoma" in text.lower() or "patogenica" in text.lower()


# ── Limitations section tests ───────────────────────────────────────────

class TestLimitationsSection:
    def test_returns_elements(self):
        _require_reportlab()
        from wes_clinical_report_es import _build_limitations_section, build_styles
        styles = build_styles()
        elements = _build_limitations_section(styles)
        assert len(elements) > 5

    def test_limitations_in_spanish(self):
        _require_reportlab()
        from wes_clinical_report_es import _build_limitations_section, build_styles
        from reportlab.platypus import Paragraph
        styles = build_styles()
        elements = _build_limitations_section(styles)
        paragraphs = [e for e in elements if isinstance(e, Paragraph)]
        text = " ".join(p.text for p in paragraphs)
        assert "cobertura" in text.lower() or "limitaciones" in text.lower()


# ── Section rendering tests ─────────────────────────────────────────────

class TestRenderSection:
    def test_renders_section(self, sample_md_path):
        _require_reportlab()
        from wes_clinical_report_es import (
            parse_markdown_report, render_section_es, build_styles
        )
        report = parse_markdown_report(sample_md_path)
        styles = build_styles()
        elements = render_section_es(report["sections"][0], styles)
        assert len(elements) > 0
